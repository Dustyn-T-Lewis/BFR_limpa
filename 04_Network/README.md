# 04 · Network

Finds proteins that move together, whether or not a curated set groups them, and asks whether
those groups respond to training, differ between the arms or track muscle size and echo intensity.

- [`01_modules`](01_modules/README.md): which proteins form co-expression modules in four signed
  WGCNA networks (main, change, pre, post), whether the modules move with training or differ
  between the arms, and whether pre-training modules hold after training. Participant adjustment
  and the signed bicor network are from Matheus's PR #24.
- [`02_module_traits`](02_module_traits/README.md): whether any module's level or change tracks
  muscle size or echo intensity, drawn as the module-trait heatmap of the YvO Figure 5.

Each step README gives its inputs, outputs and results.
