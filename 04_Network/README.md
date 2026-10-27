# 04 · Network

Finds proteins that move together, whether or not a curated set groups them, and asks whether
those groups respond to training, differ between the arms or track muscle growth.

[`01_modules`](01_modules/README.md) builds a signed WGCNA network on abundance with participant
removed, tests the module eigengenes with the stage 02 designs, enriches each module against the
project set list and checks preservation from pre- to post-training samples. Participant explains
0.61 to 0.69 of the leading components of the raw data. Participant adjustment and module settings
are from Matheus's PR #24.
