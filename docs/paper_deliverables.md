# Paper deliverables (quantitative)

**Paper:** Chen, Cheng, Guo, [arXiv:2607.02894](https://arxiv.org/abs/2607.02894) (2026).

## Configuration

| | Paper | This project |
| --- | --- | --- |
| Lattice | $N=8$ (main), open BC; $a=m=g=1$ | same |
| Ground / excited | VQE; VQD with overlap penalty. Paper does not state `reps`. | `RealAmplitudes`, `reverse_linear`, statevector. VQE Table I and Fig. 1 scan `reps=4`. VQD `reps=10`, orthogonalized to those `reps=4` ground states. |
| Optimizer | SLSQP. Paper does not state `ftol`, `maxiter`, or restarts. | SLSQP. VQE: `ftol=1e-12`, `maxiter=5000`. Table I: 10 restarts, `seed=42`, `draw_seed=4042`, winning restart index 7 (`restart_energies` in `data/VQE/table1_n8.npz`). Fig. 1 VQE and VQD: one batch per point, warm start plus 5 random restarts (`seed=42`). The files in `data/` are that batch. Fidelity is recorded after the scan and does not open another batch. |
| ED | exact diagonalization | full $2^N$ after Gauss-law elimination; `eigsh` `which="SA"` (dense `eigh` for quench evolution). Fig. 2: $N=8,10,12,14,16,18$, coarse `linspace(0,1.5,30)`, bisection `tol=1e-4`. |
| Quench | $\varepsilon=0$ vacuum → field; 2nd-order Trotter, $\Delta t=0.1$, $t/a\in[0,12]$ | same 121-point grid. Trotter `order=2`, `dt=0.1`, internal `reps=14` (paper does not state `reps`). Initial state is the ED vacuum. |

---

## Table I — VQE vs ED ($N=8$, $\varepsilon=0$)

| | ED energy | VQE energy | Rel. error | Fidelity |
| --- | ---: | ---: | ---: | ---: |
| **Paper** | $-4.63805774$ | $-4.63766032$ | $8.569\times10^{-5}$ | $99.9931\%$ |
| **Project** | $-4.638057740873512$ | $-4.637983956989514$ | $1.590835822243348\times10^{-5}$ | $0.9999877469280656$ ($99.9988\%$) |

Project ansatz/optimizer: `reps=4`, SLSQP `ftol=1e-12`, `maxiter=5000`, 10 restarts.

---

## Spectrum scan fidelities (Fig. 1 grid)

Config: $N=8$, $a=m=g=1$, $\varepsilon=\mathrm{linspace}(0,3,81)$. VQE Fig. 1 scan: `reps=4`, SLSQP `ftol=1e-12`, `maxiter=5000`, warm start + 5 restarts, `seed=42`. VQD: `reps=10`, SLSQP `maxiter=1000`, class-default `ftol`, warm start + 5 restarts, `seed=42`, orthogonality reference = the `reps=4` ground states. Paper quotes neither fidelity.

| | Paper | Project |
| --- | --- | --- |
| $\min F$ (VQE ground) | not quoted | $0.9995993561009304$ at $\varepsilon=0.7125$ |
| $\max F$ (VQE ground) | not quoted | $0.9999999735077645$ at $\varepsilon=3$ |
| $\min F$ (VQD excited) | not quoted | $0.991098295470183$ at $\varepsilon=0.525$ |
| Max $\|E_0^{\mathrm{VQE}}-E_0^{\mathrm{ED}}\|$ | not quoted | $4.866953001871899\times10^{-4}$ at $\varepsilon=0.7125$ |
| Max $\|E_1^{\mathrm{VQD}}-E_1^{\mathrm{ED}}\|$ | not quoted | $0.014760982524809396$ at $\varepsilon=0.525$ |

Selected energies at scan points (project):

| $\varepsilon$ | ED $\Gamma$ | ED $E_0$ | ED $E_1$ | VQE $E_0$ | VQE $F$ | VQD $E_1$ | VQD $F$ |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| $0$ | $-0.446415$ | $-4.638058$ | $-3.555144$ | $-4.637984$ | $0.999988$ | $-3.554800$ | $0.999932$ |
| $0.525$ | $-0.444074$ | $-3.677328$ | $-2.793564$ | $-3.677246$ | $0.999987$ | $-2.778803$ | $0.991098$ |
| $0.7125$ | $-0.329241$ | $-2.993447$ | $-2.878571$ | $-2.992961$ | $0.999599$ | $-2.878221$ | $0.999912$ |
| $1.8$ | $-0.221985$ | $0.039575$ | $0.071245$ | $0.039760$ | $0.999925$ | $0.071530$ | $0.999958$ |
| $3$ | $-0.223799$ | $8.222978$ | $8.394409$ | $8.222979$ | $1.000000$ | $8.394414$ | $0.999999$ |

Condensate-jump midpoints (ED): $\varepsilon\approx 0.69375$, $1.78125$.

---

## Critical field $\varepsilon_c$ vs $1/N$ (Fig. 2)

Config: $a=m=g=1$, open BC, full $2^N$ sparse `eigsh` (no sector cut). $\varepsilon_c$ = bisection of the condensate jump, `tol=1e-4`, after `linspace(0,1.5,30)`. `data/ED/critical_field_scaling.npz` stores the condensate and $\varepsilon_c$. Paper does not state the grid or the tolerance.

| $N$ | $1/N$ | Paper $\varepsilon_c$ | Project $\varepsilon_c$ |
| ---: | ---: | ---: | ---: |
| $8$ | $0.125000$ | $0.692384$ | $0.692391$ |
| $10$ | $0.100000$ | $0.643062$ | $0.643041$ |
| $12$ | $0.083333$ | $0.613544$ | $0.613542$ |
| $14$ | $0.071429$ | $0.593925$ | $0.593944$ |
| $16$ | $0.062500$ | $0.579960$ | $0.579952$ |
| $18$ | $0.055556$ | $0.569529$ | $0.569546$ |

| Fit | Paper (6 pts, digitized) | Project ($N=8\ldots18$) |
| --- | ---: | ---: |
| Intercept $\varepsilon_c(\infty)$ | $\approx 0.469$ (fit $0.4688$) | $0.468835$ |
| Slope $\mathrm{d}\varepsilon_c/\mathrm{d}(1/N)$ | $\approx 1.767$ | $1.766521$ |

---

## Charge / energy conservation (quench, $N=8$)

Config: $\varepsilon\in\{0.5,1.0,1.5,2.0\}$, $t/a\in[0,12]$.

| $\varepsilon$ | Paper $Q_N$ | Project $\max\|Q_N^{\mathrm{ED}}\|$ | Paper $\Delta E$ | Project $\max\|\Delta E^{\mathrm{ED}}\|$ |
| ---: | --- | ---: | ---: | ---: |
| $0.5$ | remains $0$ | $2.78\times10^{-16}$ | $\approx 0$ | $3.11\times10^{-15}$ |
| $1.0$ | remains $0$ | $3.33\times10^{-16}$ | $\approx 0$ | $8.88\times10^{-16}$ |
| $1.5$ | remains $0$ | $2.78\times10^{-16}$ | $\approx 0$ | $3.11\times10^{-15}$ |
| $2.0$ | remains $0$ | $2.78\times10^{-16}$ | $\approx 0$ | $1.24\times10^{-14}$ |

Trotter $\max\|Q_N\|\sim 3.3\times10^{-16}$.

---

## Vacuum fidelity $P_{\mathrm{vac}}$ / Trotter–ED (quench)

| $\varepsilon$ | Paper | Project min Trotter–ED fidelity on $[0,12]$ |
| ---: | --- | ---: |
| $0.5$ | near unity / slower decay | $1.00000000$ |
| $1.0$ | near unity / slower decay | $0.99999981$ |
| $1.5$ | faster departure | $0.99999816$ |
| $2.0$ | faster departure | $0.99998307$ |
| $2.5$ | strong-field | $0.99992999$ |
| $3.0$ | strongest | $0.99988039$; ED $\min P_{\mathrm{vac}} = 3.028\times10^{-7}$ |

These fidelities are `min` over time in `data/Trotter/quench_n8.npz` (`reps=14`). The exact $P_{\mathrm{vac}}$ minimum is from `data/ED/vacuum_fidelity_n8.npz`.

One symmetric product at the paper's $\Delta t = 0.1$, measured while choosing `reps` and recorded in [investigations/investigation.md](investigations/investigation.md), had final-time fidelity $0.9999$, $0.994$, $0.948$, $0.621$, $0.096$, and $0.038$ at $\varepsilon = 0.5, 1, 1.5, 2, 2.5, 3$. At $\varepsilon = 3$ the fidelity along that quench fell to $0.010$. That single-product trajectory is not in `data/Trotter/`; the stored quench is the 14-fold product.

---

## Effective decay rate $\gamma_{\mathrm{eff}}$ (fit $0\le t/a\le 1$)

| $\varepsilon$ | Paper | Project $\gamma_{\mathrm{eff}}$ |
| ---: | ---: | ---: |
| $0.5$ | not tabulated | $0.031198$ |
| $1.0$ | not tabulated | $0.119958$ |
| $1.5$ | not tabulated | $0.259015$ |
| $2.0$ | not tabulated | $0.442253$ |
| $2.5$ | not tabulated | $0.664476$ |
| $3.0$ | not tabulated | $0.917994$ |

Paper: $\gamma_{\mathrm{eff}}$ increases monotonically with $\varepsilon$. From `data/ED/decay_rate_n8.npz`, the slope of $\gamma_{\mathrm{eff}}$ versus $\varepsilon$ is $0.357187$.

---

## Appendix / missing

| Item | Paper | Project |
| --- | --- | --- |
| Figs. 4–5, 8 ($Q_i$, $H_E$, $q_n$; $N=8$) | curves only | ED series stored; no separate numeric table |
| Figs. 9–11 ($N=12$ dynamics) | shown | not computed |

---

## Overlay residuals

Maximum $|y_{\mathrm{ED}} - y_{\mathrm{paper}}|$ after reading the published SVG strokes with `visualization/overlay.py` and interpolating the current `data/ED/` curves onto those strokes. Fig. 2 compares $\varepsilon_c$ with the digitized markers. Fig. 1's raw maximum sits on the vertical stroke at a condensate jump. Away from $0.08$ of either jump ($\varepsilon \approx 0.69$ and $\varepsilon \approx 1.78$), the condensate residual is $1\times10^{-5}$ and the $E_0$ residual is $1.1\times10^{-3}$. The $E_1$ residual away from those jumps is $2.2\times10^{-2}$, at $\varepsilon = 0.80$ on the digitized stroke.

| Figure | Quantity | max $\|{\mathrm{ED}} - {\mathrm{paper}}\|$ |
| --- | --- | ---: |
| 1 | condensate | $8.0\times10^{-2}$ (at the jump) |
| 1 | $E_0$ | $4.5\times10^{-2}$ (at the jump) |
| 1 | $E_1$ | $4.7\times10^{-2}$ (at the jump) |
| 2 | $\varepsilon_c$ marker | $2.1\times10^{-5}$ |
| 3 | $Q_N$ and $\Delta E$ | $2.7\times10^{-6}$ |
| 4 | $Q_i$ | $8.6\times10^{-5}$ |
| 5 | $H_E$ | $2.2\times10^{-4}$ |
| 6 | $P_{\mathrm{vac}}$ | $1.1\times10^{-4}$ |
| 7 | $\gamma_{\mathrm{eff}}$ | $3.0\times10^{-5}$ |
| 8 | $q_n$ | $2.6\times10^{-4}$ |

The Fig. 2 paper column above matches these SVG markers. The largest marker residual is $N = 10$.

