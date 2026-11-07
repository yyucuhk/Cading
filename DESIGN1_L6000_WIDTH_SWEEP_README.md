> CURRENT REVISION: PNC terminates at 5 um width, at 31.378626347 um on the reference profile. See `PNC5UM_README.txt` and `validate_2611_pnc5um.py`. Earlier 45 um active-length / 3 um cutoff descriptions below are historical.

# Separate Design 1: 6000 um straight sections, 10 um PNCs, three ridge widths

Run from the existing Cading repository:

```bash
python 2611_design1_L6000_P10_width_sweep.py
```

The new generator does not overwrite any existing Design-1, Design-2 or
18-device-sweep source or output. Its files are:

- `gds/2611_design1_L6000_P10_width_sweep.gds`
- The same stem `.json`: measured dimensions, checks, placement and hash.
- The same stem `.csv`: the three device parameters.
- `validation_2611_L6000_width_sweep/`: actual exported previews and build logs.

## Geometry and placement

| Row, bottom to top | Optical ridge width, Layer 15 | Each straight | Straight PNC pitch | Bend PNC pitch |
|---|---:|---:|---:|---:|
| 1 | 1.62 um | 6000 um | 10 um | 10 um |
| 2 | 1.64 um | 6000 um | 10 um | 10 um |
| 3 | 1.68 um | 6000 um | 10 um | 10 um |

Each device has THREE equal straight sections and TWO 180-degree Euler bends.
The inherited effective bend radius is 100 um, Euler fraction 0.2. Meander
centerline length is `3*6000 + 2*actual_Euler180_length`, approximately
18741.822 um (18.741822 mm), excluding the end PNC transitions, optical
tapers and tethers outside the two meander ports. The JSON records measured
lengths, not just the nominal settings.

Devices are translated together with every layer, at 600 um Y pitch with
origins -700, -100 and 500 um. Leftmost Layer-11 bend edges align at X=-900 um.
There is one Layer-2 polygon label per device, e.g. `D1_W1.62_T80_P10`.
The old meander-length label is suppressed before flattening so it cannot
cut an exposure mask. The new wider reference outline encloses fabrication
geometry and inherited global markers; it is not the old 5 mm-wide chip.
Layer 51 remains a nonfabrication fiber guide and is excluded from chip-fit
and fabrication-mask overlap tests.

## PNC pitch convention

10 um is CENTER-TO-CENTER pitch, not edge-to-edge clearance. On straights,
there are 599 interior pairs per 6000 um section, with the inherited shared
pair at each straight/bend junction. On bends, the pitch is measured along
the meander CENTERLINE arc, not the inner or outer offset support curve.
Interior pairs are at exact 10 um intervals. The approximately 370.911 um
bend arc is not an integer multiple of 10 um, so the remaining length is
split symmetrically at the two bend ends: approximately 10.456 um between
each junction pair and the nearest interior pair. There are 36 interior
pairs per bend. The routine does not stretch the interior pitch to fit.

The 45 um PNC TAPER transition now uses the same 10 um pitch as the straights.
Stations are at 0, 10, 20, 30, 40 um from its port. The final 5 um
is an unpopulated remainder, not a compressed last interval.
The request changes straight and bending meander regions, not that taper.

## Retained Design-1 taper topology

Layer 10 remains 3 -> 1 -> 0.2 um over 40 + 40 um after the 45 um PNC
transition. Layer 15 retains its aligned single-stage topology: the swept
ridge width continues through the meander and PNC transition, then tapers
to the existing 0.1 um endpoint over 40 um. No offset or third stage is added.
The optical taper INPUT follows the swept width to avoid a width step.
Layer 16 follows the changed ridge. Support pads, tethers, end windows,
alignment ticks on Layer 15 and all other Design-1 layers are retained.

## Manual editing

Change Section 01 of `2611_design1_L6000_P10_width_sweep.py`:

```python
STRAIGHT_LENGTH = 6000.0
WG_WIDTHS = (1.62, 1.64, 1.68)
STRAIGHT_PNC_PERIOD = 10.0
BEND_PNC_PERIOD = 10.0
ROW_PITCH = 600.0
FIRST_ROW_Y = -700.0
LEFT_BEND_X = -900.0
```

The builder loads `2611_tapered_yue_dual_resist.py` as a parameterized
in-memory copy in a fresh process for each width. This safely binds all
function defaults/derived values before any cell is cached. Only the bend
PNC sampling expression and final old-label layer are adapted in memory;
shared source files are never edited. Do not run isolated later sections
without rerunning the whole script. Run it inside Cading with the existing
support modules and `tether_variants.gds` available.

The build fails rather than exporting success when widths, pitches,
connections, full end symmetry, mask coverage, inter-device overlaps or
all-layer export comparisons fail. These are geometry checks, not a full
fabrication-rule check or an optical-performance simulation.

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

The reinforcement is at the NARROW, 3 um-wide END of the PNC taper, directly
before the optical taper, not at the wide meander-to-PNC transition.
For the unchanged 45 um PNC taper, five 1 um periods occupy the final 5 um:
six paired stations at local x = 40, 41, 42, 43, 44, 45 um.
Coordinates are measured outward from each meander connection port.
The normal whole-end rotation applies the same geometry at both ends.

The erroneously added wide-end stations at x=1,2,3,4 (and x=5 in the long
layout) are removed. The wide end and all regular stations before x=40 keep
their original pitch and phase. Short-layout stations are 0,5,10,...,35,
then 40,41,42,43,44,45 um. The long-layout stations are 0,10,20,30,
then 40,41,42,43,44,45 um. Straight and bend PNCs are unchanged.

Edit Section 01 / 1D in either underlying design generator:

```python
PNC_CONNECTION_PITCH = 1.0
PNC_CONNECTION_PERIODS = 5  # five intervals / six stations; 0 disables
```

The dense group is anchored to the narrow endpoint (cutoff_x); its start is
cutoff_x minus pitch times periods. Both derived builders inherit these
controls and record the region and start/end distances in their manifests.
The 45 um PNC envelope, optical taper dimensions, other layers, labels and
placements remain unchanged. Run validate_2611_taper_pitch.py after
regenerating all four GDSs. Capture a baseline before editing to compare all
non-PNC layers and confine changes to the old/new reinforcement regions.
Mechanical robustness has not been simulated; checks are geometric only.
