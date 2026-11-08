DESIGN 1: 6000 um MEANDER WITH A 2000 um STRAIGHT OUTPUT
=====================================================

Run from the Cading root:
    python 2611_design1_L6000_straight_output.py

New output (existing layouts are not overwritten):
    gds/2611_design1_L6000_P10_straight_output_comparison.gds

The default file contains six devices in one vertical column:
    Bottom three: approved original two-taper devices, W=1.62/1.64/1.68 um.
    Top three: matching W=1.62/1.64/1.68 um straight-output variants.
    Y origins: -700, -100, 500, 1100, 1700, 2300 um.
    Y pitch: 600 um. All left PNC bend edges align at X=-900 um.

The new variants replace the LOWER-RIGHT output assembly only. The opposite
UPPER-LEFT input retains the current 45 um 15-to-5 um PNC transition, five
1 um dense PNC periods at its narrow end, original optical taper, butterfly
tether, support island, alignment tips, release windows and fiber guide.

Each original meander straight remains 6000 um. A separate 2000 um straight
is added outward from its original output port (local X=3000, Y=-200 um).
The new optical endpoint is local X=5000, Y=-200 um. This is NOT an extension
measured from the old taper tip, and does not change either Euler bend.

The extension is the same constant-width Layer-15 guide as the meander,
with Layer-16 protection. To retain the existing suspended-waveguide support
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
2000 um length, output PNC populations, connected ridge/protection, mask
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
