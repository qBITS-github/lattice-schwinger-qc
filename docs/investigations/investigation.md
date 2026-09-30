# Why the markers left the exact curves

This note is the repair that put the markers back on the curves. The published scans are VQE `reps=4` and VQD `reps=10`, stored in `data/` and compared in [../paper_deliverables.md](../paper_deliverables.md). The fidelities in the body are from that repair. `field_scan` and `excited_scan` take one fixed batch, and the datasets in `data/` are that batch.

Lattice units are $a = m = g = 1$, $N = 8$. Fig. 1 uses $\varepsilon$ on `np.linspace(0, 3, 81)`. Figs. 3–8 use the quench fields already stored by `scripts/run_trotter.py`, with $\Delta t = 0.1$ out to $t = 12$. The Pauli Hamiltonian and the observables are unchanged.

## Two vacuum jumps

The exact ground states on either side of $\varepsilon = 0.675 \to 0.7125$ are orthogonal, and so are the ground states on either side of $\varepsilon = 1.7625 \to 1.800$. Each jump is a first-order change of vacuum. `field_scan` used to carry one parameter vector forward and run one SLSQP minimization of Eq. (18). SLSQP follows the local gradient. The parameters that prepare the old vacuum sit in a separate basin of the new Hamiltonian, so that minimization stays there. The energy it reports is the continuation of the old vacuum, now above the true ground state.

At the first jump the old scan had fidelity $\approx 0$ and an energy error of $0.117$. After the second jump it stayed in the old basin through $\varepsilon = 1.9875$, where the fidelity was $5.2 \times 10^{-9}$ and the energy error was $0.754$. The ends of the same scan were already faithful: fidelity $0.9996$ at $\varepsilon = 0$ and $0.9999$ at $\varepsilon = 3$.

The repair is the one Section III.B already asks for: at each later field, minimize from the previous parameters and from random restarts, and keep the lowest energy. A restart replaces the current vector only when its energy is lower. Most restarts still stop in a higher minimum; at the first jump about one draw in eight reaches the new vacuum. Further draws are tried while the best state is still short of the ground state. The kept vector is the lowest energy. Fidelity only decides whether to draw again.

After that scan the minimum fidelity on the 81-point grid is $0.9969$, at $\varepsilon = 0.7125$, and the maximum $|E_0^{\mathrm{VQE}} - E_0^{\mathrm{ED}}|$ is $1.143 \times 10^{-2}$, at the same field. Both jumps pass. Before the change the minimum fidelity was $5.2 \times 10^{-9}$ and the maximum energy error was $0.754$.

## First excited state

The old excited-state scan had the same warm-start failure on a wider interval. The minimum fidelity was $1.1 \times 10^{-14}$, at $\varepsilon = 1.2375$, the maximum $|E_1^{\mathrm{VQD}} - E_1^{\mathrm{ED}}|$ was $0.894$, at $\varepsilon = 2.475$, and 56 of the 81 fields were below fidelity $0.99$. The penalty in Eq. (19) was already doing its job: the overlap with the VQE reference stayed below $10^{-8}$. One minimization was landing in a different orthogonal state. On the interval where that reference was itself the wrong vacuum, Eq. (19) could not select the first excited state. The reference is now the VQE state at the same field, and only after that state is the ground state on the whole grid.

Restarts plus that reference put every field on the exact curve except the rearrangement of the first excited state inside the $Q = -1$ sector, near $\varepsilon = 0.5$. There the exact excited states at neighboring fields overlap by $0.65$, and a shallow `RealAmplitudes` circuit cannot make the minimizer of Eq. (19) the exact first excited state. Through nine repetitions the lowest penalized cost near $\varepsilon = 0.4875$ and $0.525$ still has fidelity under $0.99$. At ten repetitions that lowest cost is `evecs[:, 1]`. The scan keeps the lowest penalized cost, and draws again while the fidelity is still short of $0.99$. SLSQP is the same optimizer as Eq. (18), and it receives $\partial C_1 / \partial \theta$.

After the scan the minimum fidelity is $0.9931$, at $\varepsilon = 0.4875$, and the maximum energy error is $2.43 \times 10^{-3}$, at that same field. Every stored $\varepsilon$ passes. The overlap with the VQE reference stays below $2 \times 10^{-9}$.

## Why $\Delta t = 0.1$ fails at large $\varepsilon$

The quench already uses the paper's initial state, the $\varepsilon = 0$ ground state, and one symmetric second-order product over the Pauli terms of $H(\varepsilon)$. There are 44 of those terms. The spectral norm of $H$ is $31.9$ at $\varepsilon = 0.5$, $54.9$ at $\varepsilon = 1.5$, and $102.5$ at $\varepsilon = 3$, so $\lVert H \rVert \Delta t$ at the stored step $\Delta t = 0.1$ runs from $3.2$ to $10.3$. A single product of 44 noncommuting factors, each carrying a piece of that step, leaves the exact trajectory once $\varepsilon \gtrsim 1.5$.

With one product the final-time fidelity was $0.9999$, $0.994$, $0.948$, $0.621$, $0.096$, and $0.038$ at $\varepsilon = 0.5, 1, 1.5, 2, 2.5, 3$. At $\varepsilon = 3$ the fidelity along the quench fell to $0.010$ and the energy drift missed the exact curve by about $0.51$.

Grouping the same Hamiltonian into commuting sets and taking one symmetric product at $\Delta t = 0.1$ still missed the bar. The repair that passes is the shorter effective step: Eq. (22) is repeated 14 times inside each stored interval of $0.1$, so the product is taken at $0.1/14$. The stored times are unchanged.

After that change the final-time fidelities are $1$, $1$, $0.999998$, $0.999983$, $0.999930$, and $0.999880$ on the same six fields. The largest plotted observable deviation, in units of the line width of `reproduction.png`, is $0.029$, $0.067$, $0.23$, $0.43$, $0.64$, and $0.81$. All six are under one line width, and every final fidelity is above $0.99$.

## Files and commands

Changed:

- `src/vqe.py` — warm start and random restarts in `field_scan`, lowest Eq. (18) energy kept
- `src/vqd.py` — same pattern for Eq. (19), ten `RealAmplitudes` repetitions, derivative passed to SLSQP
- `src/trotter.py` — fourteen second-order products inside each stored $\Delta t = 0.1$
- `scripts/run_vqe.py`, `scripts/run_vqd.py`, `scripts/run_trotter.py` — the pass/fail table

Regenerate the datasets, then the figures:

```bash
python scripts/run_vqe.py
python scripts/run_vqd.py
python scripts/run_trotter.py
python visualization/VQE.py
python visualization/VQD.py
python visualization/Trotter.py
```

`run_vqd.py` stops if the VQE fidelity on the grid is below $0.99$. Each script prints $\varepsilon$, fidelity, the energy or observable error, and pass or fail, and raises if any row fails.
