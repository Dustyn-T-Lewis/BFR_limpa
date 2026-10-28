# 04 · Network

Finds proteins that move together, whether or not a curated set groups them, and asks whether
those groups respond to training, differ between the arms or track muscle size and echo intensity.

[`01_modules`](01_modules/README.md) builds four signed WGCNA networks: on all samples with
participant removed, on each leg's training change, and on the pre- and post-training samples
alone. It names each module by enrichment and hub protein, tests the main-network eigengenes with
the stage 02 contrasts and the change-network modules for training, and checks whether
pre-training modules hold after training. Participant adjustment and the signed bicor network are
from Matheus's PR #24.

[`02_module_traits`](02_module_traits/README.md) draws the module-trait heatmap of the YvO Figure 5
for each network: contrasts and phenotype associations on one r-equivalent scale, BH and nominal
hits marked.
