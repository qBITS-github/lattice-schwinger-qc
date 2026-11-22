# Contributing

## Branches

Push to `dev` first. Always.

`dev` is where work lands. Commit and push there before anything else. Do not push directly to `main`.

Open a pull request only when a merge into `main` is actually required. That pull request goes from `dev` into `main`. `main` stays the stable line: reviewed exact-diagonalization results and the docs that describe them. Day-to-day changes stay on `dev` until that review is warranted.

## General

- Match the code that is already here. `src/baseline/model.py` builds the Pauli Hamiltonian and observables, and `scripts/baseline/run_exact_diagonalization.py` writes the `.npz` datasets. `visualization/baseline/ed.py`, `VQE.py`, `VQD.py`, and `Trotter.py` write `reproduction.png` for the paper figures that use that algorithm. `visualization/baseline/overlay.py` writes `overlay.png`.
- Lattice units stay $a = m = g = 1$ unless a change is explicitly about those parameters. The quench starts from the zero-field ground state.
- Keep the stagger $(-1)^{n+1}$ on the paper's site label $n = 1 \ldots N$. That is $(-1)$ raised to the Python index $n - 1$, which starts at 0. The extra $+1$ is only that zero-based index. It does not change the physics. Use the same factor in the mass, the Gauss law, and the charge observables. Mixing it with $(-1)^n$ on the same 1-based label is an indexing mismatch.
- A physics or numerics change that affects a published comparison needs the matching note in `docs/` and, when the README states a number, an update there too. State only results that have been recomputed.
- Appendix A, tensor networks, and hardware runs stay out of result claims until they exist. VQE, VQD, Trotter, and the Fig. 2 critical field through $N = 18$ are part of the results. When those numbers change, update the README and `docs/baseline/paper_deliverables.md`. The README's Research Status section is the boundary.
- The static datasets in `data/baseline/ed/`, `data/baseline/vqe/`, `data/baseline/vqd/`, `data/baseline/trotter/` and plots in `figures/baseline/fig*/` are the frozen reproduction record. Preserve their hashes. New driven data belong in `data/driven/<run_id>/`; retain the raw points, controls and source provenance listed in `data/driven/run_inventory.json`. Keep only the selected result figures in `figures/driven/results/`, with the questions and limitations in `docs/driven/figures.md`. Caches and bytecode are disposable.
- Keep commits focused. One change of physics, one dataset, or one doc correction is easier to review than a mixed push.
