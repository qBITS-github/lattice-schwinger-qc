# Contributing

## Branches

Push to `dev` first. Always.

`dev` is where work lands. Commit and push there before anything else. Do not push directly to `main`.

Open a pull request only when a merge into `main` is actually required. That pull request goes from `dev` into `main`. `main` stays the stable line: reviewed exact-diagonalization results and the docs that describe them. Day-to-day changes stay on `dev` until that review is warranted.

## General

- Match the code that is already here. `src/schwinger_model.py` builds the Pauli Hamiltonian and observables, and `scripts/run_exact_diagonalization.py` writes the `.npz` datasets. `visualization/ED.py`, `VQE.py`, `VQD.py`, and `Trotter.py` write `reproduction.png` for the paper figures that use that algorithm. `visualization/overlay.py` writes `overlay.png`.
- Lattice units stay $a = m = g = 1$ unless a change is explicitly about those parameters. The quench starts from the zero-field ground state.
- Keep the stagger $(-1)^{n+1}$. Site index $n = 1 \ldots N$ in the paper is Python index $n - 1$. A sign change here moves the exact-diagonalization curves off the paper.
- A physics or numerics change that affects a published comparison needs the matching note in `docs/` and, when the README states a number, an update there too. State only results that have been recomputed.
- Leave Appendix A, VQE, VQD, Trotter, tensor networks, and hardware runs out of result claims until they exist in the repository. The README's Research Status section is the boundary.
- Do not commit `__pycache__/`, `figures/archived/`, or `data/archived/`. Those paths are gitignored.
- Keep commits focused. One change of physics, one dataset, or one doc correction is easier to review than a mixed push.
