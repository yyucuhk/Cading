# Design 2 - offset ridge converter review layout

This is a separate layout. Design 1 source and GDS are unchanged.

## Geometry

One device with an 800 um middle meander straight; the existing left-bend alignment is retained. Both coupler ends use the same complete all-layer factory, related by a 180-degree rotation. The device label reads `DESIGN 2_L800_...`.

All dimensions below are in micrometres. Local x increases from the meander interface toward the outer taper.

| Feature | Local x range | Width | Centerline y |
|---|---|---|---|
| Layer 15 across PNC | 0 to 45 | 1.0 constant | 0 |
| Layer 15 first transition | 45 to 85 | 1.0 to 0.6, inherited quintic profile | 0 |
| Layer 15 new converter | 85 to 145 | 0.6 to 0.2, linear profile | 0 to 0.840, linear shift |
| Layer 10 post-PNC section | 45 to 85 | 3.0 to 1.0, unchanged profile | 0 |
| Layer 10 final taper | 85 to 185 | 1.0 to 0.2, unchanged quintic profile stretched to 100 um | 0 |

The first ridge transition retains the previous 40 um length because no new length was specified for the 1.0-to-0.6 um stage. The 100 um setting applies to the final Layer-10 taper, not the combined 40+100 um Layer-10 transition.

## Separation convention requiring review

`RIDGE_CONVERTER_Y_SHIFT = 0.840` is interpreted as the final CENTERLINE-TO-CENTERLINE offset. At the converter end, Layer 10 is approximately 0.454 um wide, so the measured vertical edge-to-edge gap in the exported geometry is 0.513 um, not 0.840 um.

An intended 0.840 um same-x edge gap would instead require approximately 1.167 um centerline shift for these profiles. Confirm the convention before using this draft for fabrication.

## Files and execution

Run from the existing Cading repository environment:

```bash
python 2611_tapered_yue_design2.py
```

Output: `gds/2611_tapered_yue_design2.gds`

The existing Design 1 output is not overwritten. The pinned generation/test engine was gdsfactory 9.45.0. Shared local helpers such as `tapered_local.py`, `layers.py`, and the tether library remain dependencies.

The complete per-device geometry is constructed by `make_device(meander_straight_length)`. The new ridge polygon is in `design2_offset_ridge_converter()`; user-editable parameters remain near the beginning of the file.

The Layer-16 protection follows the revised ridge. The existing Layer-2/Layer-14 exposure-window generation, PNC cells, supports, pads, and moved Layer-15 alignment ticks are retained. Their preservation is checked against the inherited design with the requested 100 um final Layer-10 taper.

## Checks

`gds/2611_design2_validation.json` records 36 passing checks: actual width/offset measurements, continuous optical ridge, both end assemblies, all required layers, export fidelity, coverage by Layers 14 and 16, retained left-bend alignment, and unchanged original Design 1 source/GDS.

The checks validate layout geometry, not optical conversion efficiency or fabrication process performance.
