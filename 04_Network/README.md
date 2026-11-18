# 04 · Network

Finds proteins that move together, whether or not a curated set groups them, and asks whether
those groups respond to training, differ between the arms or track muscle size and echo intensity.

- [`01_modules`](01_modules/README.md): four signed WGCNA networks, module names, eigengene and
  change tests, preservation. Participant adjustment and the signed bicor network are from
  Matheus's PR #24.
- [`02_module_traits`](02_module_traits/README.md): the module-trait heatmap of the YvO Figure 5
  for each network.

Each network holds 7 or 8 modules (power 8, `deepSplit = 2`, at least 20 proteins), with 39% to
50% of proteins unassigned. Training raises the oxidative phosphorylation module and lowers the
calcium transport module in both arms; no module differs between the arms at BH. Three
mitochondrial module levels track RF size or echo intensity at BH, with nominal hits at chance.
