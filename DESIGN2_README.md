# Design 2: three-stage slab / two-stage offset ridge

Run `python 2611_tapered_yue_design2.py` from the Cading root.
Output: `gds/2611_tapered_yue_design2.gds`.
Run `python validate_2611_design2.py` for the current polygon-level checks.

One complete device has an 800 um middle meander straight and a DESIGN 2 label.
The existing 45 um PNC section is unchanged; the following distances are
measured outward from its 3 um Layer-10 endpoint.

| Layer | Stage | Start width (um) | End width (um) | Length (um) |
|---|---|---:|---:|---:|
| 10 | 1 | 3.000 | 1.650 | 20 |
| 10 | 2 | 1.650 | 0.980 | 100 |
| 10 | 3 | 0.980 | 0.200 | 10 |
| 15 | 1 | 1.000 | 0.500 | 20 |
| 15 | 2 | 0.500 | 0.200 | 100 |

Layer 10 totals 130 um. Layer 15 totals 120 um and has NO third stage or
terminal straight extension. The second Layer-15 stage shifts linearly
from zero to 0.840 um in local Y. This is the retained CENTERLINE offset,
not an edge gap: the endpoint same-x Y edge gap is 0.250 um because Layer
10 is 0.980 um wide there and Layer 15 is 0.200 um wide.

Layer-10 stages and the first Layer-15 stage retain the quintic profile.
The second Layer-15 stage retains linear width and centerline interpolation.
Both end assemblies use one local-coordinate function and 180-degree rotation.
Stage ports are connected exactly, without clipping or patching devices.
Release/protection masks are regenerated from the changed geometry. The
alignment-tip landmark tracks the end of stage 1; ticks remain on Layer 15.
Tethers/supports and their masks track the new Layer-10 terminal port.
Left-bend alignment remains x=-900 um; the optical meander is unchanged.
Design 1 source and GDS are not modified.

Editable dimensions are grouped in sections 3 and 4 of the Python source:
`TAPER_STAGE1_LENGTH`, `TAPER_STAGE2_LENGTH`, `TAPER_STAGE3_LENGTH`,
`SLAB_STAGE1_END_WIDTH`, `SLAB_STAGE2_END_WIDTH`,
`RIDGE_EXTENSION_START_WIDTH`, `RIDGE_CONVERTER_END_WIDTH`, and
`RIDGE_CONVERTER_Y_SHIFT`.
Geometry validation is not an optical-performance simulation or fabrication DRC.
