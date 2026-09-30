# Lattice Schwinger model

Reproducing and extending Chen, Cheng & Guo,
[arXiv:2607.02894](https://arxiv.org/abs/2607.02894), on nonequilibrium
dynamics in the 1+1D lattice Schwinger model.

Lattice units are $a = m = g = 1$, and the quench starts from the zero-field ground state. The datasets in `data/` and the figures in `figures/` are the reproduction record. [docs/reproduction.md](docs/reproduction.md) is how they were produced.

## Current Results

`src/schwinger_model.py` builds the Pauli Hamiltonian and the paper's observables. The stagger on the paper's site label $n = 1 \ldots N$ is written $(-1)^{n+1}$. The extra $+1$ is only because Python indexes from zero. It does not change the physics. With that factor used in every staggered term, the exact-diagonalization curves sit on the paper's.

For $N = 8$, the spectrum, charge dynamics, electric-field energy, vacuum fidelity, and early-time decay rate match Figs. 1 and 3–8. The maximum deviation from each digitized paper stroke is in [docs/paper_deliverables.md](docs/paper_deliverables.md). The exact zero-field ground-state energy is $-4.63805774$, in agreement with Table I. The VQE energy at that point is $-4.63798396$, fidelity $99.9988\%$ (paper $99.9931\%$). On the Fig. 1 grid the VQE ground-state fidelity is at least $0.999599$ and the VQD first-excited fidelity is at least $0.991098$. The critical field for $N = 8, 10, 12, 14, 16, 18$ matches the Fig. 2 markers. A linear fit through those six sizes has slope $1.7665$ and intercept $0.4688$ (paper $\approx 1.767$ and $\approx 0.469$).

VQE, VQD, and second-order Trotter are implemented in `src/` and run from `scripts/`. At $N = 8$ their markers sit on the exact curves for Figs. 1 and 3–8. The quench repeats the second-order product fourteen times inside each stored $\Delta t = 0.1$, and the plotted observables stay within a line width of the exact evolution. The numbers are in [docs/paper_deliverables.md](docs/paper_deliverables.md).

Documentation:

- [docs/paper_deliverables.md](docs/paper_deliverables.md) — configuration, paper numbers, and this project's numbers.
- [docs/reproduction.md](docs/reproduction.md) — generate the figures.
- [docs/investigations/investigation.md](docs/investigations/investigation.md) — the earlier repair that put the markers on Figs. 1 and 3–8. The checks are in [docs/investigations/procedure.md](docs/investigations/procedure.md).
- [docs/physics/building_hamiltonian.md](docs/physics/building_hamiltonian.md) — the Pauli Hamiltonian, Eqs. (8)–(11).
- [docs/physics/building_observables.md](docs/physics/building_observables.md) — the observables and the exact time evolution.

## Next Steps

- Appendix A at $N = 12$ (Figs. 9–11).
- A time-dependent external field.
- Larger lattices, with a tensor-network baseline, then ansatz benchmarks and hardware runs.

## Research Status

The reproduction phase is closed. The record is the exact-diagonalization baseline, including the Fig. 2 critical field through $N = 18$, and the $N = 8$ VQE, VQD, and Trotter markers, stored in `data/` and `figures/` and compared in [docs/paper_deliverables.md](docs/paper_deliverables.md). The items under Next Steps are later work.
