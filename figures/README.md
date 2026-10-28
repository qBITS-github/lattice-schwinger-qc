# Frozen reproduction plots and new driven figures

**FROZEN REPRODUCTION — 19 October 2026, source `5497cd4a763b00f92269f7e8e650e2664d5c72f1`.**

The 16 PNGs in the eight `baseline/fig1_*` through `baseline/fig8_*` directories are the
closed reproduction plots. `reproduction.png` is the repository calculation;
`overlay.png` compares it with the reference paper. Their filenames, image bytes
and original captions are preserved. This directory label and the
[hash manifest](../docs/baseline/reproduction_manifest.json) identify their frozen status.
The PNGs are locally read-only; Git does not preserve that permission setting.

The eight selected driven result figures belong in `driven/results/`.
Generate them from saved results, label them **DRIVEN STUDY**, and retain the
input/output hashes and selection rules in their plot manifest. The
[figure index](../docs/driven/figures.md) states what each retained figure
establishes. No new plot belongs in the frozen directories. Other future
static reruns also need a separate output location.
