DESIGN 1: 6000 um MEANDER WITH A 2000 um STRAIGHT OUTPUT
=====================================================

Run from the Cading root:
    python 2611_design1_L6000_straight_output.py

Output (this file is overwritten; other saved layouts are unchanged):
    gds/2611_design1_L6000_P10_straight_output_comparison.gds

The default file contains six devices in one vertical column:
    Bottom three: approved original two-taper devices, W=1.62/1.64/1.68 um.
    Top three: matching W=1.62/1.64/1.68 um straight-output variants.
    Y origins: -700, -100, 500, 1100, 1700, 2300 um.
    Y pitch: 600 um. All left PNC bend edges align at X=-900 um.

The new variants replace the LOWER-RIGHT output assembly only. The opposite
UPPER-LEFT input retains the current 45 um 15-to-5 um PNC transition, five
1 um dense PNC periods at its narrow end, original optical taper, butterfly
tether, support island, alignment tips, release windows. Its Layer-51 fiber guide is omitted from this comparison.

Each original meander straight remains 6000 um. A separate 2000 um straight
is added outward from its original output port (local X=3000, Y=-200 um).
The new optical endpoint is local X=5000, Y=-200 um. This is NOT an extension
measured from the old taper tip, and does not change either Euler bend.

The extension is the same constant-width Layer-15 guide as the meander,
without Layer-16 protection in the exported file. To retain the existing suspended-waveguide support
construction, it continues the 15 um support envelope with 10 um PNC pitch,
Layer-2 side release bands and the complete Layer-14 opening. Layer 10 does
not carry a taper in this extension, just as it does not carry the original
straight meander's temporary broad construction shape. There are no output
taper stages, PNC width transition, dense terminal PNC group, butterfly,
support island, alignment tips, or output fiber guide. The end face is flat.

PNC pairs are at extension coordinates 0,10,...,1990 um. The pair at the
old port joins the existing regular grid; no doubled or half PNC is allowed.
The final support pair is 10 um before the flat optical end. No new
lithography subfield borders or chip-edge coupling performance is implied.

Main parameters are in Section 01. INCLUDE_REFERENCE_DEVICES=False produces
only the three new variants, still using 600 um spacing. The three approved
widths are verified against the saved source GDS; intentionally changing the
reference design requires reviewing/updating the baseline comparison.

Layer construction is in make_straight_output(). It replaces only owned,
flattened geometry beyond the terminal port and BELOW the middle arm.
Both complete Euler bends and all three original straight sections are
outside the replacement region. The original cells are not mutated.

Checks include exact original-device polygon comparison, exact preservation
outside the output assembly on every layer, measured output widths, exact
2000 um length, output PNC populations, connected ridge and release-mask
coverage, retained input ticks, absence of output taper/tether/pads/guides,
correct device/label counts, alignment, spacing and overlaps. The generated
JSON records these checks. Geometry checks are not optical-loss simulation,
mechanical qualification, or complete fabrication DRC.

Reproduction environment:
    Python 3.11
    gdsfactory==9.45.0
    klayout==0.30.13
    numpy==2.4.6
Keep the existing Design-1 and long-width builder, local support modules,
and tether_variants.gds in this directory. The downloadable release also
includes these dependencies; it never changes their contents.

OUTPUT LAYER CLEANUP
====================
Layers 16 and 51 are omitted everywhere in the six-device comparison GDS.
Global alignment markers formerly on Layer 2 are moved (not copied) to
Layer 15 at exactly the same positions, with the same shapes. Layer-10
markers are retained. Device Layer-2 release windows and all six device
labels remain unchanged on Layer 2. Alignment tips already on Layer 15
remain unchanged. No waveguide, PNC, taper, placement, or frame changes.

The policy is explicit in Section 01 and finalize_output_layers(). It runs
after all construction/placement, on an independent KLayout layout, before
writing the deliverable. Intermediate worker masks may include Layer 16
for baseline checks; they are not the final export.

This change applies ONLY to the comparison builder/output. The source
three-device width sweep, standalone designs and 18-device sweep are not
modified. The original baseline GDS must retain its manifest hash.


TWO-COUPLER LAYER POLICY (current)
Only the three two-coupler devices in this comparison omit Layer 11.
Layer 14 is retained only at the two existing terminal assemblies, starting
at the actual Layer-10 taper entrance (local X=-3045 and +3045 um).
Both outward channels and their existing pads/openings are preserved.
The meander and preceding 15-to-5 um PNC-transition openings are removed
from Layer 14; no other device layer is edited. All three straight-output
variants remain identical to the preceding comparison.
The global removal of Layers 16/51 and L2-to-L15 marker move remain active.
Device labels are intentionally unchanged; P10 in the two-coupler labels
is a legacy source identifier, NOT a remaining Layer-11 pattern.
The original three-device source GDS remains unchanged because it is the
shared source for constructing the straight-output variants.
Control: TWO_COUPLER_END_WINDOWS_ONLY in Section 01 of the builder.
Current exported-mask audit: gds/2611_two_coupler_windows_validation.json.
The old layer-cleanup report is superseded by this current report.
