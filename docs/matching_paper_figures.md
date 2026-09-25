# Matching the paper's figures

The exact-diagonalization curves already sit on the published figures. The markers that leave those curves are the VQE and VQD field scan, and the second-order Trotter quench at the paper's stated step $\Delta t = 0.1$.

Lattice units are $a = m = g = 1$, $N = 8$. The numbers below are from the current runs in `data/VQE/`, `data/VQD/`, and `data/Trotter/`.

## Fig. 1

`RealAmplitudes` with three repetitions already represents both ends of the scan. The ground-state fidelity is $0.9996$ at $\varepsilon = 0$ and $0.9999$ at $\varepsilon = 3$. The first excited state is essentially exact for $\varepsilon \ge 2.55$.

The ground-state markers leave the curve only at the two first-order jumps. Fidelity is $0$ at $\varepsilon \approx 0.71$, and again from $\varepsilon = 1.80$ to $1.99$, where the energy sits as much as $0.75$ above the exact ground state. `field_scan` carries one parameter vector forward and runs a single SLSQP minimization. Across those jumps the ground state changes to an orthogonal vacuum, and that local minimization stays in the old basin.

The excited-state markers leave over a wider interval, $\varepsilon \approx 0.45$ to $2.48$. The overlap with the exact first excited state drops to $0$ and the energy is as much as $0.89$ too high, while the overlap with the VQE reference stays below $10^{-8}$. The penalty in Eq. (19) is already enforcing orthogonality. One warm-started minimization is landing in a different orthogonal state. On the two intervals where the VQE reference itself is the wrong vacuum, Eq. (19) cannot select the first excited state.

To put the Fig. 1 markers on the curves:

- At each field, keep several random restarts together with the warm start.
- For VQE, keep the lowest energy. For VQD, keep the lowest penalized cost.
- Build the VQD reference from the VQE state at that same field, after the VQE state is the ground state.

More repetitions, a larger penalty, or a different Hamiltonian will not close these gaps. Section III.B already asks for random restarts, and it also says to start each field from the previous parameters. The published markers require both, not the warm start alone.

## Figs. 3–8

The quench already follows Eqs. (21)–(22). The initial state is the $\varepsilon = 0$ ground state. Each step is one symmetric product of $\Delta t = 0.1$, repeated out to $t/a = 12$, over the 44 Pauli terms of $H(\varepsilon)$. Total charge stays at zero. The state fidelity does not:

| $\varepsilon$ | minimum fidelity |
|---:|---:|
| 0.5 | 0.9998 |
| 1.0 | 0.993 |
| 1.5 | 0.946 |
| 2.0 | 0.621 |
| 2.5 | 0.096 |
| 3.0 | 0.010 |

$\lVert H \rVert$ grows from about $32$ at $\varepsilon = 0.5$ to about $103$ at $\varepsilon = 3$, so $\lVert H \rVert \Delta t$ grows from about $3$ to about $10$. A step of that size, split across 44 noncommuting factors, is why the markers leave the exact curves for $\varepsilon \gtrsim 1.5$, and why the vacuum fidelity and the decay rate fall apart for $\varepsilon \ge 2.5$.

To put those markers on the curves, shrink the effective step: take $\Delta t$ below $0.1$, or use more than one symmetric repetition inside each step of $0.1$. Grouping $H$ into a few commuting pieces, instead of one exponential per Pauli string, cuts the same error at $\Delta t = 0.1$. The initial state and the order of the formula are already the ones in Section III.C.
