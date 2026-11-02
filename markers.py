"""
Alignment, lithography-caliper and dicing markers.

Each ``*_array`` wrapper stamps one marker into all four quadrants; the
``marker_x``/``marker_y`` pair is the corner *half*-extent (so a marker lands
at every ``(+-marker_x, +-marker_y)``), while `calipermarker_array` takes the
full chip shape instead.
"""
import gdsfactory as gf
from layout.core.layers import Layers
import numpy as np


def _place_at_corners(component, cell, dx, dy):
    """Reference *cell* into *component* at the four (+-dx, +-dy) corners."""
    for x in [-1, 1]:
        for y in [-1, 1]:
            component.add_ref(cell).move((x * dx, y * dy))


@gf.cell
def marker_cross(center = (0,0), arm_length = 150, arm_width = 2, layer = Layers.globalmarker, label = None):
    """Cross marker, with an optional label.

    Parameters
    ----------
    label : str, optional
        Drawn in the bottom right quadrant of the cross.
    """
    C = gf.Component()
    cross = C.add_ref(gf.components.cross(length=arm_length, width=arm_width, layer=layer))

    cross.move((center[0],center[1]))

    if label is not None:
        ref = C.add_ref(gf.components.text(label, size=10, layer=layer))
        ref.center = [center[0] + 20,center[1] -10]

    C.flatten()
    return C

@gf.cell
def marker_cross_array(marker_x = 4500, marker_y = 3300, arm_length = 150, arm_width = 0.5, layer = Layers.globalmarker, label=None):
    """Four `marker_cross` copies at the (+-marker_x, +-marker_y) corners."""
    c = gf.Component()
    _place_at_corners(
        c, marker_cross(arm_length=arm_length, arm_width=arm_width, layer=layer, label=label),
        marker_x, marker_y)
    return c

@gf.cell
def caliperarm(angle=0, num_notches=10, notch_width=1, layer1=Layers.label, layer2=Layers.globalmarker):
    """Lithography caliper arm: two notch rows that walk out of registration.

    Parameters
    ----------
    angle : float
        Rotation, in degrees.
    num_notches : int
        Notch pairs on each side of the middle pair.
    layer1, layer2 : Layer
        Inner and outer notch rows.
    """
    c = gf.Component()

    c << gf.components.litho_calipers(
        notch_size = [notch_width,notch_width*5],
        notch_spacing = 5,
        num_notches =  num_notches,
        offset_per_notch = 0.1,
        row_spacing = notch_width*.5,
        layer1 = layer1,
        layer2 = layer2
    )

    c.rotate(angle)
    return c

@gf.cell
def calipermarker(
    arm_length=100,
    arm_width=2,
    num_notches=None,
    layer1=Layers.label,
    layer2=Layers.globalmarker,
    center_square_size=3,
):
    """Caliper marker: four `caliperarm` rows around a cross with center squares.
    It has lithography calipers and a few other markers.

    Parameters
    ----------
    arm_length, arm_width : float
        Cross arm length and width, in microns.
    num_notches : int, optional
        Notch pairs on each side of the middle pair; derived from the arm
        dimensions when left as ``None``.
    layer1, layer2 : Layer
        Inner and outer notch rows; ``layer2 = None`` drops the arms.

    Returns
    -------
    gf.Component
        Ports ``center_down`` and ``center_up`` at the marker center.
    """
    c = gf.Component()

    if num_notches is None:
        num_notches = int(arm_length / arm_width / 4)-1

    if layer2 is not None:
        for i in range(4):
            ang = -90*i
            arm = c.add_ref(caliperarm(num_notches=num_notches, notch_width=arm_width, layer1=layer1, layer2=layer2))
            arm.move((-arm.center[0], 0))
            arm.rotate(ang)
            arm.move((-arm_length*np.sin(i*np.pi/2), -arm_length*np.cos(i*np.pi/2)))

    # The cross, with two small squares meeting at its center.
    r = c << gf.components.rectangle(size=(arm_length - 5, arm_width), centered=True, layer=layer1)
    r.move([5+(arm_length - 5)/2, 0])
    r = c << gf.components.rectangle(size=(arm_length - 5, arm_width), centered=True, layer=layer1)
    r.move([-5-(arm_length - 5)/2, 0])
    r = c << gf.components.rectangle(size=(arm_width, arm_length - 5), centered=True, layer=layer1)
    r.move([0, 5+(arm_length - 5)/2])
    r = c << gf.components.rectangle(size=(arm_width, arm_length - 5), centered=True, layer=layer1)
    r.move([0, -5-(arm_length - 5)/2])

    # Preserve the original square centers at (+/-1.5, +/-1.5). Increasing
    # center_square_size expands each square about its own fixed center rather
    # than moving the two squares farther away from the marker origin.
    square_offset = 1.5
    r = c << gf.components.rectangle(
        size=(center_square_size, center_square_size),
        centered=True,
        layer=layer1,
    )
    r.move([square_offset, square_offset])

    r = c << gf.components.rectangle(
        size=(center_square_size, center_square_size),
        centered=True,
        layer=layer1,
    )
    r.move([-square_offset, -square_offset])

    c.add_port(name = 'center_down', center = [0.0, 0.0], width = 1, orientation = -90, layer=layer1)
    c.add_port(name = 'center_up', center = [0.0, 0.0], width = 1, orientation = 90, layer=layer1)
    return c

@gf.cell
def calipermarker_array(
    chip_shape=(4500, 3300),
    arm_length=150,
    arm_width=2,
    num_notches=None,
    layer1=Layers.label,
    layer2=Layers.globalmarker,
    center_square_size=3,
):
    """Four `calipermarker` copies at the chip corners (+-chip_shape/2).

    Returns
    -------
    gf.Component
        Ports ``center_down`` and ``center_up`` at the marker center.
    """

    c = gf.Component()

    _place_at_corners(
        c,
        calipermarker(
            arm_length=arm_length,
            arm_width=arm_width,
            num_notches=num_notches,
            layer1=layer1,
            layer2=layer2,
            center_square_size=center_square_size,
        ),
        chip_shape[0]/2, chip_shape[1]/2)

    c.add_port(name = 'center_down', center = [0.0, 0.0], width = 1, orientation = -90, layer=layer1)
    c.add_port(name = 'center_up', center = [0.0, 0.0], width = 1, orientation = 90, layer=layer1)

    return c

@gf.cell
def recframe(arm_width = 2,square_size=50, layer1=Layers.globalmarker):
    """Hollow square frame of outer size ``square_size`` and border ``arm_width``."""
    D1 = gf.Component()
    D1 << gf.components.rectangle(size = (square_size, square_size), layer=layer1)

    # Center the cutout inside D1 so the frame border is uniform on all sides.
    D2 = gf.Component()
    (D2 << gf.components.rectangle(
        size = (square_size-2*arm_width, square_size-2*arm_width), layer=layer1
    )).move((arm_width, arm_width))

    return gf.boolean(A = D1, B = D2, operation = 'not', layer=layer1)

@gf.cell
def DLalignmarker(center = (0,0), arm_length = 150, layer1=Layers.label, layer2=Layers.globalmarker):
    """Box-in-box vernier marker ("DL").

    A cross plus 2x2 grids of `recframe` outlines on ``layer1``, and 2x2 grids of
    solid squares on ``layer2`` at an inner and an outer pitch, so the two
    lithography steps can be read against each other by eye.
    """
    # TODO: simplify.
    arm_width = arm_length/30
    square_size=arm_length/3

    D = gf.Component()
    D << gf.components.cross(length = arm_length*2, width = arm_width, layer = layer1)

    device_list1 = []
    for _ in [1, 2, 3, 4]:
        D2 = gf.components.rectangle(size = (square_size, square_size), layer = layer2)
        device_list1.append(D2)
    G = gf.grid(device_list1,
                spacing = (arm_length/2,arm_length/2),
                shape = (2,2),
                align_x = 'origin',
                align_y = 'origin')

    D << G

    device_list2 = []
    for _ in [1, 2, 3, 4]:
        D3 = recframe(arm_width = arm_width,square_size=square_size+4*arm_width, layer1= layer1)
        device_list2.append(D3)
    G2 = gf.grid(device_list2,
                spacing = ((arm_length)/2-arm_width*4,(arm_length)/2-arm_width*4),
                shape = (2,2),
                align_x = 'center',
                align_y = 'center')
    D<<G2

    device_list1 = []
    a1=0.75
    for _ in [1, 2, 3, 4]:
        D4 = gf.components.shapes.rectangle(size = (square_size-a1*arm_width, square_size-a1*arm_width), layer = layer2)
        device_list1.append(D4)
    G3 = gf.grid(device_list1,
                spacing = ((arm_length-(square_size-a1*arm_width))*2,(arm_length-(square_size-a1*arm_width))*2),
                shape = (2,2),
                align_x = 'center',
                align_y = 'center')
    D<<G3
    D.center = (center[0],center[1])
    return D

@gf.cell
def DLalignmarker_array(marker_x = 4500, marker_y = 3300, arm_length = 150, layer1=Layers.label, layer2=Layers.globalmarker):
    """Four `DLalignmarker` copies at the (+-marker_x, +-marker_y) corners.

    Returns
    -------
    gf.Component
        Ports ``center_down`` and ``center_up`` at the marker center.
    """
    D = gf.Component()
    _place_at_corners(
        D, DLalignmarker(arm_length=arm_length, layer1=layer1, layer2=layer2),
        marker_x, marker_y)

    D.add_port(name = 'center_down', center = [0.0, 0.0], width = 1, orientation = -90, layer=layer1)
    D.add_port(name = 'center_up', center = [0.0, 0.0], width = 1, orientation = 90, layer=layer1)

    return D

@gf.cell
def marker_square_array(marker_x, marker_y, size, layer):
    """Four centered squares at the (+-marker_x, +-marker_y) corners."""
    c = gf.Component()
    _place_at_corners(
        c, gf.components.rectangle(size=(size, size), centered=True, layer=layer),
        marker_x, marker_y)
    return c

@gf.cell
def marker_ebl(chip_shape=(40000, 8000), facet_length1=1000):
    """Dicing and outline marker crosses at the four chip corners.

    The default parameters are based on YFAv7_250413_AsymAMPMPMDC_WG30cm_v2.

    Parameters
    ----------
    chip_shape : (float, float)
        (x-length, y-length) of the chip, in microns.
    facet_length1 : float
        Facet length, in microns; the outline markers sit half of it past the chip
        edge in y, and a fixed 500 um past it in x.
    """
    c = gf.Component()

    marker_outline_x = chip_shape[0]/2+500
    marker_outline_y = chip_shape[1]/2+facet_length1/2

    c.add_ref(marker_cross_array(marker_x = chip_shape[0]/2 + 25, marker_y = chip_shape[1]/2, layer = Layers.wg, label="dice"))
    c.add_ref(marker_cross_array(marker_x = marker_outline_x, marker_y = marker_outline_y, layer = Layers.wg, label="out"))

    c.flatten()

    return c
