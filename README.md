# Lattice Schwinger model

Reproducing and extending Chen, Cheng & Guo,
[arXiv:2607.02894](https://arxiv.org/abs/2607.02894), on nonequilibrium
dynamics in the 1+1D lattice Schwinger model.

Lattice units are $a = m = g = 1$, and the quench starts from the zero-field ground state. `data/` and `figures/` are generated locally; [docs/reproduction.md](docs/reproduction.md) is the run order.

## Current Results

`src/schwinger_model.py` builds the Pauli Hamiltonian and the paper's observables. The stagger is $(-1)^{n+1}$, and with that sign the exact-diagonalization curves sit on the paper's.

For $N = 8$, the spectrum, charge dynamics, electric-field energy, vacuum fidelity, and early-time decay rate match Figs. 1 and 3–8. The zero-field ground-state energy is $-4.63805774$, in agreement with Table I. The critical field for $N = 8, 10, 12, 14, 16$ matches the Fig. 2 markers. A fit through those five sizes intercepts near $0.465$. The paper's quoted intercept, $0.469$, includes $N = 18$.

VQE, VQD, and second-order Trotter are implemented in `src/` and run from `scripts/`. At $N = 8$ their markers sit on the exact curves for Figs. 1 and 3–8. On the Fig. 1 grid the ground-state fidelity is at least $0.9969$ and the first-excited fidelity is at least $0.9931$. The quench repeats the second-order product fourteen times inside each stored $\Delta t = 0.1$, and the plotted observables stay within a line width of the exact evolution. The account is in [docs/investigations/investigation.md](docs/investigations/investigation.md).

Documentation:

- [docs/reproduction.md](docs/reproduction.md) — generate `data/` and `figures/`.
- [docs/investigations/investigation.md](docs/investigations/investigation.md) — why the markers sit on Figs. 1 and 3–8. The checks are in [docs/investigations/procedure.md](docs/investigations/procedure.md).
- [docs/building_hamiltonian.md](docs/building_hamiltonian.md) — the Pauli Hamiltonian, Eqs. (8)–(11).
- [docs/building_observables.md](docs/building_observables.md) — the observables and the exact time evolution.

## Next Steps

- Appendix A at $N = 12$ (Figs. 9–11), and the $N = 18$ ground state so the Fig. 2 extrapolation uses the same sizes as the paper.
- A time-dependent external field.
- Larger lattices, with a tensor-network baseline, then ansatz benchmarks and hardware runs.

## Research Status

The exact-diagonalization baseline and the $N = 8$ VQE, VQD, and Trotter markers are the results claimed here. The other items under Next Steps are not part of the results yet.
