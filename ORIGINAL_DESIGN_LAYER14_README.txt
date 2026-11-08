ORIGINAL DESIGN 1 AND DESIGN 2 - LAYER 14 ONLY

Run the standalone scripts or 2611_build_combined_sweep.py normally.
Keep original_design_layer14.py beside these scripts.
Section 01: LAYER14_POST_PNC_MARGIN = 1.0; LAYER14_SIMPLE_TERMINAL_WINDOW = True.

The full 45 um PNC transition, its endpoint opening, meander and bend windows
remain unchanged. At the PNC endpoint the optical-taper window starts at the
existing broad width and smoothly narrows across the first optical stage.
The remaining optical-taper opening has 1 um normal expansion per side.
At the final optical-taper tip, a linear channel starts at the current narrow
opening width and reaches the existing outer endpoint and two 300 um circles.
The old rectangular frame is absent on Layer 14; Layer 2 is untouched.

Only standalone exports and the combined 18-device sweep call this final
Layer-14 policy. Raw make_device() construction is unchanged so the L6000
family retains its separate masks and all saved L6000 layouts remain untouched.
The current exported-geometry report is gds/2611_original_design_windows_validation.json.
Prior reports describe earlier GDS revisions. These are geometric checks,
not optical/mechanical simulation or a complete fabrication DRC.

The saved release reuses exact approved non-L14 records to avoid CAD-engine
rounding differences when rebuilding old geometry. The source construction
functions for all other layers are unchanged.
