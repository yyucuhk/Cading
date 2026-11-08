from pathlib import Path
p=Path('2611_design1_L6000_straight_output.py')
s=p.read_text()
a='''    device = m.make_device(STRAIGHT_LENGTH)
    device.flatten()
    require(not kdb.Region(device.kdb_cell.begin_shapes_rec(gf.get_layer((998, 0)))).is_empty(),
            'Original label was not isolated')
    device.remove_layers([(998, 0)], recursive=False)
    original = regions(device.kcl.layout, device.kdb_cell)
    info = device.info
'''
b='''    # The approved GDS is the exact reference, not a new reconstruction.
    # Reading its independently owned polygons avoids export-grid differences
    # in unchanged corners and preserves every existing layer exactly.
    ref_id = base.case_id(width)
    baseline, _ = read_gds(BASE_GDS)
    saved = baseline.cell(ref_id)
    require(saved is not None, 'Reference GDS does not contain the requested width')
    manifest = json.loads(BASE_GDS.with_suffix('.json').read_text())
    require(manifest['gds_sha256'] == hashlib.sha256(BASE_GDS.read_bytes()).hexdigest(),
            'Approved GDS and manifest hashes disagree')
    record = next(v for v in manifest['devices'] if v['case_id'] == ref_id)
    require(record['width_um'] == width and record['straight_length_um'] == STRAIGHT_LENGTH,
            'Approved reference dimensions differ from requested dimensions')
    require(record['straight_pnc_period_um'] == PNC_PITCH,
            'Approved reference has a different straight PNC pitch')
    original = regions(baseline, saved)
    centers0 = record['connection_centers']
    tips0 = record['tip_centers']
    orientations = [math.degrees(math.atan2(c[1]-t[1], c[0]-t[0])) % 360
                    for c, t in zip(centers0, tips0, strict=True)]
    info = dict(connection_centers=centers0, tip_centers=tips0,
                connection_orientations=orientations,
                total_meander_length_um=record['total_meander_um'])
'''
assert a in s
s=s.replace(a,b).replace("f'regenerated_reference_L{layer[0]}_exact'", "f'approved_reference_L{layer[0]}_exact'")
p.write_text(s)
compile(s,str(p),'exec')
