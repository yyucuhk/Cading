# Editing Design 1 and Design 2

Both generators use the same section numbers. The latest meander update uses
three equal straight sections joined by two unchanged 180-degree Euler bends.

| Design | Script | Default output |
|---|---|---|
| Design 1: aligned ridge | `2611_tapered_yue_dual_resist.py` | `gds/2611_tapered_yue_dual_resist_L10_L12.gds` |
| Design 2: offset ridge | `2611_tapered_yue_design2.py` | `gds/2611_tapered_yue_design2.gds` |

Keep these scripts in the Cading repository root with `layers.py`, `markers.py`,
`phononics.py`, `tapered_local.py`, `tether_variants.py`,
`extract_layer2_gdsfactory.py`, `text.py`, and `tether_variants.gds`.
The scripts do not import one another. Changes to one script's settings do not
change the other design. Shared support-library changes can affect both.

## 1. Start in Section 01: USER PARAMETERS

Dimensions are in micrometres (um), angles in degrees. Thus 980 nm is `0.980`,
500 nm is `0.500`, and 200 nm is `0.200`.

### 1A: Number of devices and placement

```python
MEANDER_STRAIGHT_LENGTHS = (1400.0,)  # one device; keep the trailing comma
DEVICE_Y_PITCH = 600.0
BOTTOM_DEVICE_Y = -700.0
LEFT_BEND_ALIGNMENT_X = -900.0
```

One tuple entry creates one complete device, ordered bottom to top. For three:

```python
MEANDER_STRAIGHT_LENGTHS = (200.0, 800.0, 1400.0)
```

The tuple is the ONLY device-count control. Both scripts currently retain the
one-device 1400 um setting. Each entry controls ALL THREE straight sections,
not the total optical path length. With L=800, the straights are 800/800/800 um.
The horizontal port coordinates are -L/2 and +L/2; complete end assemblies
follow those ports without any change to their local geometry.
Other Section-01 settings are shared by
all devices in that script; they are not per-device dictionaries.

The device origin is placed at `BOTTOM_DEVICE_Y + index * DEVICE_Y_PITCH`.
This is not the y coordinate of a taper tip. Complete devices are translated
in x so the leftmost Layer-11 PNC bend edge equals `LEFT_BEND_ALIGNMENT_X`.
Every layer, label and support moves together. Nothing is clipped or spliced.

### 1B: Optical taper dimensions

The following stages start AFTER the 15 -> 3 um PNC transition. The PNC region
has its own fixed length (45 um), independent of its PNC pitch.

#### Design 1

| Parameter | Current value | Meaning |
|---|---:|---|
| `RIDGE_TAPER_LENGTH` | 40 | First post-PNC stage on both layers |
| `SLAB_TAPER_LENGTH` | 40 | Final Layer-10 stage only |
| `PNC_TAPER_SUPPORT_CUTOFF_WIDTH` | 3 | Layer-10 first-stage input width |
| `SLAB_EXTENSION_START_WIDTH` | 1 | Layer-10 intermediate width |
| `REVIEW_TETHER_CONNECTION_WIDTH` | 0.2 | Layer-10 final width and butterfly neck |
| `RIDGE_MEANDER_WIDTH` | 1 | Layer-15 meander and taper input width |
| `RIDGE_EXTENSION_START_WIDTH` | 0.1 | Layer-15 taper END width |

The inherited variable names are retained for compatibility. In particular,
`RIDGE_EXTENSION_START_WIDTH` is the narrow post-PNC taper end in the outward
coordinate convention used by this assembly.

Default sequence: Layer 10 is 3 -> 1 -> 0.2 um over 40/40 um. Layer 15 is
1 -> 0.1 um over the first 40 um, then ends. There is no offset ridge converter.

#### Design 2

```python
TAPER_STAGE1_LENGTH = 20.0
TAPER_STAGE2_LENGTH = 100.0
TAPER_STAGE3_LENGTH = 10.0  # Layer 10 only
SLAB_STAGE1_END_WIDTH = 1.65
SLAB_STAGE2_END_WIDTH = 0.980
RIDGE_MEANDER_WIDTH = 1.0
RIDGE_EXTENSION_START_WIDTH = 0.500
RIDGE_CONVERTER_END_WIDTH = 0.200
RIDGE_CONVERTER_Y_SHIFT = 0.840
REVIEW_TETHER_CONNECTION_WIDTH = 0.200
```

| Stage | Length | Layer 10 | Layer 15 |
|---|---:|---|---|
| 1 | 20 um | 3 -> 1.65 um | 1 -> 0.500 um |
| 2 | 100 um | 1.65 -> 0.980 um | 0.500 -> 0.200 um, laterally shifted |
| 3 | 10 um | 0.980 -> 0.200 um | Absent |

Layer 10 totals 130 um; Layer 15 ends at 120 um. Layer-10 stages and the first
Layer-15 stage use a quintic profile. The second Layer-15 stage uses linear
width and centerline interpolation. `RIDGE_CONVERTER_Y_SHIFT` is the final
CENTERLINE shift, not an edge gap. With the current widths the same-x vertical
edge gap at the ridge endpoint is 0.840 - 0.980/2 - 0.200/2 = 0.250 um.

### 1C and 1D: Meander and PNC supports

`S_BEND_COMPACT_180_RADIUS` and `S_BEND_EULER_P` control each Euler bend.
`RIDGE_MEANDER_WIDTH` in 1B is the optical width. `S_BEND_WIDTH` is the wider
construction/support envelope and is NOT that optical width.

`MEANDER_STRAIGHT_PNC_PERIOD` controls the exact center-to-center PNC pitch
on the three meander straights (currently 5 um). `REVIEW_SUPPORT_BEND_PERIOD`
independently controls the nominal fitted-to-arc bend pitch.
`REVIEW_SUPPORT_ANCHOR_PERIOD` is an alias of `MEANDER_STRAIGHT_PNC_PERIOD`.
The taper follows that pitch; `PNC_TAPER_TRANSITION_LENGTH = 45.0` independently fixes its length.
Changing the straight-meander pitch does not change that taper's length. `PHONONIC_SUPPORT_SEGMENTS` and
`PHONONIC_SUPPORT_GAPS` define the library cell's masses and gaps.
`S_BEND_TAPER_LENGTH` belongs to the inherited noncompact route; it is not the
equal-straight sweep control.

### Equal-straight length calculation

In both scripts, `make_device(L)` now draws:

```text
straight L -> 180-degree Euler bend -> straight L
           -> opposite 180-degree Euler bend -> straight L
```

The centerline length between the two meander ports is
`L_total = 3 * L + 2 * L_bend180`. The code measures `L_bend180` with
`gf.path.euler(radius=S_BEND_COMPACT_180_RADIUS, angle=180,
p=S_BEND_EULER_P, use_eff=True).length()`, rather than using the circular
approximation `pi*R`. The device ML label and metadata use this same formula.
PNC/optical tapers outside the meander ports and tethers are NOT included.

For the unchanged radius and Euler fraction, the bend contribution is constant;
increasing L by 100 um increases the total meander length by 300 um.
The left bend alignment and Y pitch are retained. End assemblies move along X
as needed; they are not clipped, stretched or copied between devices.
`S_BEND_HALF_LENGTH` is no longer an editable constant: Section 09 calculates
its per-device replacement as L/2. L must be positive. Keep `TAPER_ANGLE=0`;
rotate an entire device instance when a rotated device is needed.

### 1E: Tethers and support islands

Edit `TETHER_VARIANT`, `TETHER_PAD_SCALE`, `TETHER_CONNECTION_BLEND_LENGTH`,
`SUPPORT_RECTANGLE_SIZE`, `EXTRA_RECTANGLE_GAP` and `SUPPORT_SURROUND_SLOT_WIDTH`.
`REVIEW_TETHER_CONNECTION_WIDTH` couples the slab's final width to its neck.
`EXTRA_RECTANGLE_GAP` is the edge-to-edge distance from each outer taper
endpoint to the nearest edge of its support rectangle, not the distance
between the two devices or between opposite support islands.

The optional `SUPPORT_TIP_*` values describe extra lateral supports disabled
by `add_support_tips=False` inside `make_device()`. They do not control the
small fiber-alignment ticks. The latter remain on Layer 15.

### 1F: Release/protection masks

`POSITIVE_WINDOW_MARGIN` expands the Layer-2 taper window;
`POSITIVE_WINDOW_END_EXTENSION` must not be smaller than that margin.
`RIDGE_PROTECTION_MARGIN` sets the Layer-16 ridge offset. Layer 14 includes the
cleaned Layer-2 mask and solid meander/taper fill; bounded taper holes are filled.

The long-channel angles, start widths and lengths are independent on Layers 2
and 14. `LAYER14_CHANNEL_START_OFFSET` is measured from the first taper-stage
endpoint, not the final ridge/slab endpoint. Companion-window dimensions
`ALIGNMENT_TIP_WINDOW_*` describe areas removed near the ticks, not tick sizes.

Actual tick dimensions are fixed in `tapered_local.tapered_fiber_coupler` and
matched in `alignment_tip_region()` (Section 08). To change the 0.5 x 10 um
ticks or their transverse offsets, change BOTH definitions consistently.
Do not just change the support-tip or companion-window parameters.

### 1G: Labels, chip outline, global markers, output

Label font size/offset, marker sizes, chip dimensions, top-cell naming and
output filename are grouped here. The chip height retains the approved
minimum-height/row-pitch rule. Inspect device fit inside the outline when
adding many rows; a count change does not guarantee fabrication-field fit.
`TOP_CELL_NAME` is only a label; changing it does not alter device length.
There is no exposed dummy `EBL_SUBFIELD_SIZE`: it was unused and drew nothing.
No lithography subfields or write-field borders are created by this cleanup.

## 2. Where to modify actual drawing code

The `# %% [NN]` headers are navigation/cell markers. They do not make later
sections independent of imports, parameters or derived settings above them.

| Section | Contents | Typical edit |
|---|---|---|
| 00 | Imports and repository-relative paths | Runtime/dependencies |
| 01 | All frequently edited parameters, 1A-1G | Dimensions and device count |
| 02 | Layer aliases and derived lengths | Keep formulas linked to Section 01 |
| 03 | Quintic taper and Design-2 offset converter | Taper interpolation law |
| 04 | Meander/taper PNC shields | PNC placement and continuous side bands |
| 05 | Complete two-layer taper | Layer-10/15 stage construction and ports |
| 06 | Coupler packaging | Long slots, glue spots and inherited fiber guide |
| 07 | Layer-2/14/16 masks, steps 7A-7I | Exposure/protection geometry |
| 08 | Alignment-tip matching and layer transfer | Tick geometry/layer assignment |
| 09 | `make_device(L)`, steps 9A-9G | Assemble all layers at both ends |
| 10 | `build_layout()` and `add_chip_outline_and_markers()` | Placement/global marks |
| 11 | Input checks and `export_gds()` | File export and basic input validation |

The connected meander path is implemented in `tapered_local.py` and receives
its dimensions through `common_device_kwargs` in Section 09. Complete stages
are joined through ports; do not replace portions by global-x clipping.

## 3. Run after editing

Use a fresh process after changing Section 01. Function defaults and cached
GDSFactory cells are evaluated at import; rerunning only a later editor cell
can use old values. Run these commands separately in the Cading root:

```bash
python 2611_tapered_yue_dual_resist.py
python 2611_tapered_yue_design2.py
```

The saved path is printed. The existing GDS at the selected path is overwritten.
The cleanup uses the existing GDSFactory 9.45.0 environment and support files;
it does not add a new runtime module dependency.

## 4. Check changes

Run `python validate_2611_equal_straights.py` for the current equal-straight
geometry, complete-end preservation, connected meander/masks, exact length
labels, export and placement checks in both designs. This regenerates both
production GDS files. The original pre-change comparison is available when
the matching baseline is retained in `validation_2611_equal_straights/baseline/`.


```bash
python validate_2611_cleanup.py
```

This count-aware check covers both scripts, all device layers, ridge/protection
connectivity, tip assignment, left-bend placement, repeatability, plus one- and
five-device generation. It also exercises 200/800/1400 um individual devices.
Outputs go to `validation_2611_cleanup/`, leaving the production GDS untouched.

For a before/after refactor comparison, run `python validate_2611_cleanup.py capture`
BEFORE editing, then run the ordinary check AFTER editing. A baseline is not
automatically approved: changing dimensions intentionally can and should fail
an unchanged-geometry comparison. Remove/move obsolete baseline files when
validating a newly approved geometry. Without a baseline, no claim of
before/after preservation is made.

Older design-specific validators may contain assertions for fixed historical
sweeps or stage dimensions. Update their expected values when intentionally
changing those dimensions. `validate_2611_design2.py` remains useful for the
current approved Design-2 stages. Geometry checks do not replace optical
simulation, complete fabrication DRC, or lithography-field planning.

## 5. Preserved invariants

One `make_device(L)` call owns a complete flattened all-layer device. Exposure
mask construction happens before the final move of alignment ticks from
Layer 10 to Layer 15. Do not remove the independent-copy/flatten steps or use
recursive layer deletion on cached child cells. Shared-cell mutation was the
cause of earlier missing structures; these ownership rules are deliberate.

## Independent straight-meander PNC pitch (current)

`STRAIGHT_PNC_PERIOD = 5.0` in Section 01 of `2611_build_combined_sweep.py`
controls the straight PNC spacing for all 18 devices. The two standalone
generators use `MEANDER_STRAIGHT_PNC_PERIOD = 5.0` in Section 01 / 1D.
This is CENTER-TO-CENTER pitch, not edge-to-edge clearance. All three 1400 um
straights use an exact 5 um grid, with 279 interior PNC pairs per straight.
One shared PNC pair at each straight/bend junction is retained.

`REVIEW_SUPPORT_ANCHOR_PERIOD = MEANDER_STRAIGHT_PNC_PERIOD` makes the
PNC taper follow the straight-section pitch. Its physical length is fixed
separately by `PNC_TAPER_TRANSITION_LENGTH = 45.0`; optical tapers do not move.
The bend sweep remains nominal 2/3/5 um with the existing fitted-to-arc
sampler. Labels still report design, post-PNC taper length and bend period
only; the common straight period is recorded separately in CSV/JSON.

Run `python validate_2611_straight_pitch.py` after regenerating the combined
sweep to measure every straight PNC population directly in the exported GDS.

## Taper-junction pitch correction

PNCs now continue on the straight-section grid through each 45 um PNC taper,
including a pair at the connection port. The 1400 um designs use 5 um
pitch on both straights and tapers; the separate 6000 um width sweep uses
10 um on both. Bend samplers and all optical/release/protection geometry
are unchanged. No extra cell is forced onto a non-grid taper endpoint.
Run `python validate_2611_taper_pitch.py` after regenerating the four GDSs
to check both taper ends and the junction-to-straight spacings. Historical
baseline reports describe the previous geometry; the new report is
`gds/2611_taper_pitch_validation.json`.

### Reproducing the approved mask grid

The short standalone designs and 18-device sweep were rebuilt with
GDSFactory 9.45.0 and KLayout 0.30.12. The separate 6000 um sweep uses
GDSFactory 9.45.0 and KLayout 0.30.13, matching its original export.
Pin the matching versions for exact polygon comparisons; changing the
geometry-engine version can change nanometre-grid rounding in offset masks.
The taper-pitch report records both environments.


## Local input/output PNC reinforcement (current)

At both input/output straight-to-PNC-taper connections, only the first 5 um
on the taper side now uses five 1 um longitudinal center-to-center periods.
This means six paired PNC stations at local distances 0, 1, 2, 3, 4, 5 um.
The remaining PNCs retain their original positions, dimensions and grid phase.
Short-layout stations continue at 10, 15, ..., 45 um (regular pitch 5 um).
The 6000 um width-sweep stations continue at 10, 20, 30, 40 um (regular pitch
10 um); its 5-to-10 um transition gap is deliberately retained to avoid moving
any original station. Straight and bend PNC patterns are unchanged.

Edit Section 01 / 1D in either underlying design generator:

```python
PNC_CONNECTION_PITCH = 1.0
PNC_CONNECTION_PERIODS = 5  # intervals; six stations; 0 restores the regular grid
```

Both combined/long builders inherit these settings from the applicable source.
The dense segment is mirrored by the normal whole-end construction; no local
patch is applied to just one end. The 45 um envelope, optical tapers, windows,
other layers, labels and placements stay unchanged. The JSON manifests record
the separate connection pitch/count/span. Run `validate_2611_taper_pitch.py`
after regenerating all four outputs; capture an approved baseline before the
edit to enable exact polygon-difference checks. Increased support density is
verified geometrically; mechanical robustness has not been simulated.
