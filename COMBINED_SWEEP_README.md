# Design 1 + Design 2: 18-device parameter sweep

Run from your existing Cading folder:

```bash
python 2611_build_combined_sweep.py
python validate_2611_combined_sweep.py
```

Use the existing helper modules and tether_variants.gds. Validated environment:
Python 3.11, gdsfactory 9.45.0, matplotlib. The first command regenerates and
checks the sweep; the second independently checks the exported GDS.

## Settings to edit

Section 01 of `2611_build_combined_sweep.py`:

```python
STRAIGHT_LENGTH = 1400.0
BEND_PNC_PERIODS = (2.0, 3.0, 5.0)
POST_PNC_TAPER_LENGTHS = (80.0, 100.0, 120.0)
DESIGN_X_OFFSET = 5000.0
ROW_PITCH = 600.0
FIRST_ROW_Y = -700.0
LEFT_BEND_X = -900.0
STAGE_RATIOS = {1: (1.0, 1.0), 2: (2.0, 10.0, 1.0)}
```

Each of the THREE meander straights is 1400 um. The two Euler bends retain
100 um effective radius and 0.2 Euler fraction. The centerline length is
approximately 4941.823 um, excluding end tapers and tethers. The straight-region
PNC nominal period remains 15 um. Only bend PNC periods are swept.

## Definition of total taper length

The swept total is the entire Layer-10 taper AFTER the unchanged 45 um PNC
transition. Widths are unchanged. Length is distributed proportionally among
the previous stages. Layer 15 retains its topology and ends before Layer 10;
no extra ridge stage or terminal straight is added.

| Total L10 (um) | Design 1 stages (um) | Design 2 stages (um) |
|---:|---|---|
| 80 | 40 / 40 | 12.307692 / 61.538462 / 6.153846 |
| 100 | 50 / 50 | 15.384615 / 76.923077 / 7.692308 |
| 120 | 60 / 60 | 18.461538 / 92.307692 / 9.230769 |

Design 1 L10 widths: 3 -> 1 -> 0.2 um; L15: 1 -> 0.1 um during stage 1.
Design 2 L10 widths: 3 -> 1.65 -> 0.980 -> 0.200 um; L15: 1 -> 0.500 ->
0.200 um during stages 1 and 2 only. Design 2 retains the 0.840 um terminal
CENTERLINE shift and nominal 0.250 um endpoint vertical edge gap.
GDS coordinates are quantized to the existing 1 nm database grid.

## Rows, bottom to top (same order in both designs)

| Row | Nominal bend PNC period (um) | Post-PNC L10 total (um) | Origin Y (um) |
|---:|---:|---:|---:|
| 1 | 2 | 80 | -700 |
| 2 | 2 | 100 | -100 |
| 3 | 2 | 120 | 500 |
| 4 | 3 | 80 | 1100 |
| 5 | 3 | 100 | 1700 |
| 6 | 3 | 120 | 2300 |
| 7 | 5 | 80 | 2900 |
| 8 | 5 | 100 | 3500 |
| 9 | 5 | 120 | 4100 |

Design 1 is the left column; Design 2 is translated +5000 um in X.
Leftmost PNC bend edges align at X=-900 um and X=4100 um respectively.
All layers move together. Additional labels identify design, row, bend period,
taper total and common straight length, e.g. `D2_05_B3_T100_L1400`.

The unchanged PNC sampler fits integer populations into each Euler arc.
Consequently 2/3/5 um are NOMINAL requested periods. Measured interior spacings
are approximately 2.004980, 3.015796, and 5.012657 um, respectively, with
184, 122, and 73 interior PNC pairs per bend. The shared boundary PNCs remain.
Per-case JSON data record these values for BOTH bends.

## Output files

- `gds/2611_design1_design2_L1400_18device_sweep.gds`: all 18 complete devices.
- Same stem `.csv`: sweep and placement table.
- Same stem `.json`: measured geometry and placement manifest.
- `gds/2611_combined_sweep_validation.json`: independent exported-GDS checks.
- `validation_2611_combined_sweep/combined_sweep.png`: actual exported preview.

The original generator defaults are also changed to one 1400 um device,
but the new combined builder creates the 9 cases per design. Their standalone
GDS files remain separate from the combined sweep output.

## Geometry isolation and verification

Each case loads only a parameter-adjusted in-memory copy of the approved
source in a fresh Python process. Geometry function bodies are unchanged.
Derived values and default arguments are evaluated before any cell is cached.
No loaded globals are changed between devices; no device is clipped/spliced.

The builder checks the actual three straights, two bend populations, taper
junction widths/endpoints, every layer, ridge connectivity, protection and
release coverage, and retained L15 alignment ticks. It re-reads the combined
GDS and compares every layer with the separately generated devices. Non-PNC
layers must not change when only bend period changes. All fabrication-layer
pairwise inter-device overlaps are checked.

The independent validator measures both ends of every exported device and
compares full end-assembly symmetry, excluding only exact label polygons.
It also checks actual instance counts, X/Y placement, all layers and connectivity.

The combined frame encompasses fabrication geometry: X=-3500 to +8400 um,
Y=-5000 to +5000 um (11.9 x 10 mm). This is one encompassing reference outline,
not a claim of fitting inside a previous standalone 5 mm-wide chip frame.
Inherited global markers are added once around this frame. Long Layer-51
fiber guides are nonfabrication geometry; they are retained and may cross
columns. They do not determine the frame or exposure-mask collision checks.

These checks are geometric, not optical simulations or foundry DRC.
Historical single-device reports do not validate this new sweep. Update the
independent validator expectations when changing the requested matrix or widths.
