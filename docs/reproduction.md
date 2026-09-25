# Reproducing the datasets and figures

`data/` and `figures/` are not in the repository. Generate them from the repo root. The scripts find their own paths, so the working directory can be anywhere.

The runs need NumPy, SciPy, Matplotlib, Qiskit, and `qiskit-algorithms`.

## Datasets

Exact diagonalization first. It writes `data/ED/` for Figs. 1–8, including the critical-field scan through $N = 16$. This is the long run.

```bash
python scripts/run_exact_diagonalization.py
```

VQE writes `data/VQE/` (Table I at $\varepsilon = 0$, then the Fig. 1 ground-state scan). VQD reads that file and writes `data/VQD/`. Trotter does not need either; it writes `data/Trotter/`.

```bash
python scripts/run_vqe.py
python scripts/run_vqd.py
python scripts/run_trotter.py
```

## Figures

Each paper figure is a directory under `figures/`. `reproduction.png` is that figure from this repository: dashed curves are exact diagonalization, markers are VQE, VQD, or Trotter. `overlay.png` places that result on the published figure: a solid line is this exact diagonalization, markers are this repository's quantum data, a dashed line is the paper, and gray dots are the paper's quantum markers.

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

Figs. 9–11 are the $N = 12$ appendix and are not produced. Where the markers leave the exact curves, and what to change, is in [matching_paper_figures.md](matching_paper_figures.md).
