# Lattice Schwinger model

Reproducing and extending Chen, Cheng & Guo,
[arXiv:2607.02894](https://arxiv.org/abs/2607.02894), on nonequilibrium
dynamics in the 1+1D lattice Schwinger model.

Lattice units are $a = m = g = 1$, and the quench starts from the zero-field ground state. `data/` and `figures/` are generated locally; [docs/reproduction.md](docs/reproduction.md) is the run order.

## Current Results

`src/schwinger_model.py` builds the Pauli Hamiltonian and the paper's observables. The stagger is $(-1)^{n+1}$, and with that sign the exact-diagonalization curves sit on the paper's.

For $N = 8$, the spectrum, charge dynamics, electric-field energy, vacuum fidelity, and early-time decay rate match Figs. 1 and 3–8. The zero-field ground-state energy is $-4.63805774$, in agreement with Table I. The critical field for $N = 8, 10, 12, 14, 16$ matches the Fig. 2 markers. A fit through those five sizes intercepts near $0.465$. The paper's quoted intercept, $0.469$, includes $N = 18$.

VQE, VQD, and second-order Trotter are implemented in `src/` and run from `scripts/`. Their markers do not yet sit on the paper's figures. The ground-state scan leaves the curve at the two vacuum jumps, the excited-state scan leaves it over most of $\varepsilon \in [0.45, 2.48]$, and the Trotter quench drifts once $\varepsilon \gtrsim 1.5$. What to change is in [docs/matching_paper_figures.md](docs/matching_paper_figures.md).

Documentation:

- [docs/reproduction.md](docs/reproduction.md) — generate `data/` and `figures/`.
- [docs/matching_paper_figures.md](docs/matching_paper_figures.md) — why the quantum markers drift, and what to change.
- [docs/building_hamiltonian.md](docs/building_hamiltonian.md) — the Pauli Hamiltonian, Eqs. (8)–(11).
- [docs/building_observables.md](docs/building_observables.md) — the observables and the exact time evolution.

## Next Steps

- Random restarts on the VQE and VQD field scans, and a finer Trotter step, so the markers sit on Figs. 1 and 3–8. See [docs/matching_paper_figures.md](docs/matching_paper_figures.md).
- Appendix A at $N = 12$ (Figs. 9–11), and the $N = 18$ ground state so the Fig. 2 extrapolation uses the same sizes as the paper.
- A time-dependent external field.
- Larger lattices, with a tensor-network baseline, then ansatz benchmarks and hardware runs.

## Research Status

The numbers above are the exact-diagonalization baseline. VQE, VQD, and Trotter run, and their disagreement with the published markers is described in [docs/matching_paper_figures.md](docs/matching_paper_figures.md). That disagreement is not yet a result claimed here. The other items under Next Steps are not part of the results yet.
