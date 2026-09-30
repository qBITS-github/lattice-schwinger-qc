# Reproducing the datasets and figures

Numerical results and `figures/` are generated locally and are not in the repository. Generate them from the repo root. The scripts find their own paths, so the working directory can be anywhere.

The runs need NumPy, SciPy, Matplotlib, Qiskit, and `qiskit-algorithms`.

## Datasets

Exact diagonalization first. It writes the exact results for Figs. 1–8, including the critical-field scan at $N = 8, 10, 12, 14, 16, 18$. $N = 18$ makes this the long run.

```bash
python scripts/run_exact_diagonalization.py
```

VQE writes Table I at $\varepsilon = 0$, then the Fig. 1 ground-state scan (`RealAmplitudes` `reps=4`, SLSQP `ftol=1e-12`, `maxiter=5000`, set in `scripts/run_vqe.py`). VQD reads that scan and writes the first-excited states (`reps=10` in `src/vqd.py`). Trotter does not need either; it writes the quench trajectories (`order=2`, fourteen products inside each stored $\Delta t = 0.1$).

The comparison against the paper is [paper_deliverables.md](paper_deliverables.md). Run VQD after VQE: it uses the ground states that script just wrote.

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

Figs. 9–11 are the $N = 12$ appendix and are not produced.
