"""Minimal pure-Python GDSII reader/writer (no klayout/gdstk available).

Supports BOUNDARY, PATH, BOX, TEXT, SREF, AREF; flattening with
translation, rotation (any angle), x-reflection and magnification.
"""
from __future__ import annotations
import math
import struct


def _real8_dec(b):
    if b == b"\x00" * 8:
        return 0.0
    sign = -1 if b[0] & 0x80 else 1
    exp = (b[0] & 0x7F) - 64
    mant = int.from_bytes(b[1:], "big") / float(1 << 56)
    return sign * mant * (16.0 ** exp)


def _real8_enc(x):
    if x == 0:
        return b"\x00" * 8
    sign = 0x80 if x < 0 else 0
    x = abs(x)
    exp = 0
    while x >= 1:
        x /= 16.0
        exp += 1
    while x < 1 / 16.0:
        x *= 16.0
        exp -= 1
    mant = int(round(x * (1 << 56)))
    if mant >= (1 << 56):
        mant >>= 4
        exp += 1
    return bytes([sign | (exp + 64)]) + mant.to_bytes(7, "big")


class Cell:
    def __init__(self, name):
        self.name = name
        self.polygons = []   # (layer, datatype, [(x,y),...]) in DB units, open ring
        self.paths = []      # (layer, datatype, width, pathtype, pts)
        self.texts = []      # (layer, texttype, (x,y), string, mag)
        self.refs = []       # dict(name, xy, angle, mag, reflect, cols, rows, pitch)


class Library:
    def __init__(self, name="LIB", user_unit=1e-3, db_unit=1e-9):
        self.name = name
        self.user_unit = user_unit
        self.db_unit = db_unit
        self.cells = {}
        self.order = []

    def add(self, cell):
        self.cells[cell.name] = cell
        self.order.append(cell.name)
        return cell

    # ---------------------------------------------------------------- read
    @classmethod
    def read(cls, path):
        data = open(path, "rb").read()
        lib = cls()
        i = 0
        cell = None
        el = None
        while i < len(data):
            ln, rt, dt = struct.unpack(">HBB", data[i:i + 4])
            if ln < 4:
                break
            body = data[i + 4:i + ln]
            i += ln
            tag = (rt << 8) | dt
            if rt == 0x02:  # LIBNAME
                lib.name = body.rstrip(b"\x00").decode()
            elif rt == 0x03:  # UNITS
                lib.user_unit = _real8_dec(body[:8])
                lib.db_unit = _real8_dec(body[8:16])
            elif rt == 0x06:  # STRNAME
                cell = lib.add(Cell(body.rstrip(b"\x00").decode()))
            elif rt in (0x08, 0x09, 0x2D, 0x0C, 0x0A, 0x0B):
                el = {"kind": rt, "mag": 1.0, "angle": 0.0, "reflect": False,
                      "width": 0, "pathtype": 0, "dt": 0}
            elif rt == 0x0D:
                el["layer"] = struct.unpack(">h", body)[0]
            elif rt in (0x0E, 0x16, 0x2E):  # DATATYPE/TEXTTYPE/BOXTYPE
                el["dt"] = struct.unpack(">h", body)[0]
            elif rt == 0x0F:
                el["width"] = struct.unpack(">i", body)[0]
            elif rt == 0x21:
                el["pathtype"] = struct.unpack(">h", body)[0]
            elif rt == 0x10:
                n = len(body) // 4
                v = struct.unpack(">%di" % n, body)
                el["xy"] = list(zip(v[0::2], v[1::2]))
            elif rt == 0x12:
                el["sname"] = body.rstrip(b"\x00").decode()
            elif rt == 0x13:
                el["colrow"] = struct.unpack(">hh", body)
            elif rt == 0x19:
                el["string"] = body.rstrip(b"\x00").decode(errors="replace")
            elif rt == 0x1A:
                el["reflect"] = bool(body[0] & 0x80)
            elif rt == 0x1B:
                el["mag"] = _real8_dec(body)
            elif rt == 0x1C:
                el["angle"] = _real8_dec(body)
            elif rt == 0x11:  # ENDEL
                k = el["kind"]
                if k in (0x08, 0x2D):
                    pts = el["xy"]
                    if len(pts) > 1 and pts[0] == pts[-1]:
                        pts = pts[:-1]
                    cell.polygons.append((el["layer"], el["dt"], pts))
                elif k == 0x09:
                    cell.paths.append((el["layer"], el["dt"], el["width"], el["pathtype"], el["xy"]))
                elif k == 0x0C:
                    cell.texts.append((el["layer"], el["dt"], el["xy"][0], el.get("string", ""), el["mag"]))
                elif k == 0x0A:
                    cell.refs.append(dict(name=el["sname"], xy=el["xy"][0], angle=el["angle"],
                                          mag=el["mag"], reflect=el["reflect"], cols=1, rows=1, pitch=None))
                elif k == 0x0B:
                    c, r = el["colrow"]
                    p0, pc, pr = el["xy"]
                    cell.refs.append(dict(name=el["sname"], xy=p0, angle=el["angle"], mag=el["mag"],
                                          reflect=el["reflect"], cols=c, rows=r,
                                          pitch=((pc[0] - p0[0]) / c, (pc[1] - p0[1]) / c,
                                                 (pr[0] - p0[0]) / r, (pr[1] - p0[1]) / r)))
                el = None
        return lib

    def tops(self):
        used = {r["name"] for c in self.cells.values() for r in c.refs}
        return [n for n in self.order if n not in used]

    # ------------------------------------------------------------- flatten
    def flatten(self, name):
        """Return (polys, texts) in DB units: polys=[(l,d,pts float)], texts=[(l,d,(x,y),s)]."""
        polys, texts = [], []

        def tf_make(xy, ang, mag, refl):
            c, s = math.cos(math.radians(ang)), math.sin(math.radians(ang))

            def f(p):
                x, y = p
                if refl:
                    y = -y
                x, y = x * mag, y * mag
                return (c * x - s * y + xy[0], s * x + c * y + xy[1])
            return f

        def walk(cname, f):
            cell = self.cells[cname]
            for l, d, pts in cell.polygons:
                polys.append((l, d, [f(p) for p in pts]))
            for l, d, w, pt, pts in cell.paths:
                for poly in path_to_polys(pts, w, pt):
                    polys.append((l, d, [f(p) for p in poly]))
            for l, d, xy, s, m in cell.texts:
                texts.append((l, d, f(xy), s))
            for r in cell.refs:
                offs = [(0.0, 0.0)]
                if r["pitch"]:
                    a, b, c2, d2 = r["pitch"]
                    offs = [(i * a + j * c2, i * b + j * d2) for i in range(r["cols"]) for j in range(r["rows"])]
                for ox, oy in offs:
                    g = tf_make((r["xy"][0] + ox, r["xy"][1] + oy), r["angle"], r["mag"], r["reflect"])
                    walk(r["name"], (lambda gg: (lambda p: f(gg(p))))(g))
        walk(name, lambda p: p)
        return polys, texts

    # --------------------------------------------------------------- write
    def write(self, path, max_vertices=4000):
        out = bytearray()

        def rec(rt, dt, body=b""):
            out.extend(struct.pack(">HBB", 4 + len(body), rt, dt))
            out.extend(body)

        def s(txt):
            b = txt.encode()
            return b + (b"\x00" if len(b) % 2 else b"")
        now = (2026, 11, 23, 0, 0, 0)
        rec(0x00, 0x02, struct.pack(">h", 600))
        rec(0x01, 0x02, struct.pack(">12h", *now, *now))
        rec(0x02, 0x06, s(self.name))
        rec(0x03, 0x05, _real8_enc(self.user_unit) + _real8_enc(self.db_unit))
        for n in self.order:
            c = self.cells[n]
            rec(0x05, 0x02, struct.pack(">12h", *now, *now))
            rec(0x06, 0x06, s(c.name))
            for l, d, pts in c.polygons:
                for chunk in split_poly(pts, max_vertices):
                    ip = [(int(round(x)), int(round(y))) for x, y in chunk]
                    ip = dedupe(ip)
                    if len(ip) < 3:
                        continue
                    ip.append(ip[0])
                    rec(0x08, 0x00)
                    rec(0x0D, 0x02, struct.pack(">h", l))
                    rec(0x0E, 0x02, struct.pack(">h", d))
                    rec(0x10, 0x03, struct.pack(">%di" % (2 * len(ip)), *[v for p in ip for v in p]))
                    rec(0x11, 0x00)
            for l, d, xy, txt, mag in c.texts:
                rec(0x0C, 0x00)
                rec(0x0D, 0x02, struct.pack(">h", l))
                rec(0x16, 0x02, struct.pack(">h", d))
                if mag != 1.0:
                    rec(0x1A, 0x01, b"\x00\x00")
                    rec(0x1B, 0x05, _real8_enc(mag))
                rec(0x10, 0x03, struct.pack(">2i", int(round(xy[0])), int(round(xy[1]))))
                rec(0x19, 0x06, s(txt))
                rec(0x11, 0x00)
            for r in c.refs:
                rec(0x0A, 0x00)
                rec(0x12, 0x06, s(r["name"]))
                if r.get("angle") or r.get("reflect") or r.get("mag", 1.0) != 1.0:
                    rec(0x1A, 0x01, b"\x80\x00" if r.get("reflect") else b"\x00\x00")
                    if r.get("mag", 1.0) != 1.0:
                        rec(0x1B, 0x05, _real8_enc(r["mag"]))
                    if r.get("angle"):
                        rec(0x1C, 0x05, _real8_enc(r["angle"]))
                rec(0x10, 0x03, struct.pack(">2i", int(round(r["xy"][0])), int(round(r["xy"][1]))))
                rec(0x11, 0x00)
            rec(0x07, 0x00)
        rec(0x04, 0x00)
        open(path, "wb").write(bytes(out))


def dedupe(pts):
    o = []
    for p in pts:
        if not o or o[-1] != p:
            o.append(p)
    while len(o) > 1 and o[0] == o[-1]:
        o.pop()
    return o


def split_poly(pts, maxv):
    if len(pts) <= maxv:
        return [pts]
    # split by vertical cut at median x (works for simple shapes; used only rarely)
    xs = sorted(p[0] for p in pts)
    xm = xs[len(xs) // 2]
    left = clip_halfplane(pts, xm, True)
    right = clip_halfplane(pts, xm, False)
    return split_poly(left, maxv) + split_poly(right, maxv)


def clip_halfplane(pts, xc, keep_left):
    out = []
    n = len(pts)
    inside = (lambda p: p[0] <= xc) if keep_left else (lambda p: p[0] >= xc)
    for i in range(n):
        a, b = pts[i], pts[(i + 1) % n]
        ia, ib = inside(a), inside(b)
        if ia:
            out.append(a)
        if ia != ib:
            t = (xc - a[0]) / (b[0] - a[0])
            out.append((xc, a[1] + t * (b[1] - a[1])))
    return out


def path_to_polys(pts, w, pathtype):
    hw = w / 2.0
    polys = []
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        dx, dy = x1 - x0, y1 - y0
        L = math.hypot(dx, dy) or 1
        ux, uy = dx / L, dy / L
        ext = hw if pathtype == 2 else 0
        x0e, y0e, x1e, y1e = x0 - ux * ext, y0 - uy * ext, x1 + ux * ext, y1 + uy * ext
        nx, ny = -uy * hw, ux * hw
        polys.append([(x0e + nx, y0e + ny), (x1e + nx, y1e + ny), (x1e - nx, y1e - ny), (x0e - nx, y0e - ny)])
    return polys


def area(pts):
    a = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1]):
        a += x0 * y1 - x1 * y0
    return a / 2.0


def bbox(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return min(xs), min(ys), max(xs), max(ys)
