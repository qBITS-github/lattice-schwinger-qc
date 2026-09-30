# Reproducing the datasets and figures

The reproduction of Chen, Cheng & Guo, [arXiv:2607.02894](https://arxiv.org/abs/2607.02894), is closed for the claim in the README Research Status: Table I, Figs. 1–8 at $N = 8$, and the Fig. 2 critical field at $N = 8, 10, 12, 14, 16, 18$.

The datasets in `data/` and the figures in `figures/` are that record. Their SHA-256 hashes are in [paper_deliverables.md](paper_deliverables.md). Those files are enough to check the quoted numbers. The scripts below are what produced them, and running a script overwrites the file it writes. Appendix A (Figs. 9–11) is outside this record.

Install [requirements.txt](../requirements.txt) on Python 3.14.7. The scripts find their own paths, so the working directory can be anywhere.

## Datasets

Exact diagonalization first. It writes the exact results for Figs. 1–8, including the critical-field scan at $N = 8, 10, 12, 14, 16, 18$. $N = 18$ makes this the long run.

```bash
python scripts/run_exact_diagonalization.py
```

VQE writes Table I at $\varepsilon = 0$, then the Fig. 1 ground-state scan (`RealAmplitudes` `reps=4`, SLSQP `ftol=1e-12`, `maxiter=5000`, set in `scripts/run_vqe.py`). Each later field is one batch: the previous parameters plus 5 random restarts. VQD reads that scan and writes the first-excited states (`reps=10` in `src/vqd.py`), again one batch per field. Trotter does not need either; it writes the quench trajectories (`order=2`, fourteen products inside each stored $\Delta t = 0.1$).

The comparison against the paper is [paper_deliverables.md](paper_deliverables.md). Run VQD after VQE: it uses the ground states that script just wrote.

```bash
python scripts/run_vqe.py
python scripts/run_vqd.py
python scripts/run_trotter.py
```

## Figures

Each paper figure is a directory under `figures/`. `reproduction.png` is that figure from this repository: dashed curves are exact diagonalization, markers are VQE, VQD, or Trotter. `overlay.png` places that result on the published figure: a solid line is this exact diagonalization, markers are this repository's quantum data, a dashed line is the paper, and gray dots are the paper's quantum markers.

The current `figures/` files were drawn from the current `data/` files. Fig. 1 was written after `data/VQD/e1_scan_n8.npz`, Fig. 2 after `data/ED/critical_field_scaling.npz`, and Figs. 3–8 after `data/Trotter/quench_n8.npz` and the exact quench files.

```bash
python visualization/ED.py
```

That writes both files for Figs. 1–8, using the quantum datasets when they exist. The other scripts write the same files for the figures that use that algorithm:

```bash
python visualization/VQE.py      # Fig. 1
python visualization/VQD.py      # Fig. 1
python visualization/Trotter.py  # Figs. 3–8
python visualization/overlay.py  # every overlay
```

Figs. 9–11 are the $N = 12$ appendix and are not produced.
