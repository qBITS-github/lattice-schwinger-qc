# Procedure for putting the markers on the curves

This is the sequence that was actually run on branch `match-paper-markers`. The physics of the two failures is in `investigation.md`. Nothing here was committed.

The bars were fixed in advance: fidelity at least 0.99 to the exact state at every stored field, and, for the quench, a plotted observable error smaller than the line width of `reproduction.png`. The Hamiltonian, the observables, the ansatz family, and the optimizer stayed as they were. SLSQP and `RealAmplitudes` with reverse-linear entanglement were not replaced.

## 1. Check the diagnosis against the code

Read `README.md`, `docs/investigations/matching_paper_figures.md`, `docs/reproduction.md`, `src/schwinger_model.py`, `src/vqe.py`, `src/vqd.py`, `src/trotter.py`, and the three scripts in `scripts/`. The note and the code agreed.

- Table I already called `best_restart`: several random SLSQP starts, lowest Eq. (18) energy kept.
- `field_scan` did not. It carried one parameter vector from the previous field and ran one minimization.
- `excited_scan` drew random starts only at the first field, then warm-started.
- The Trotter step was `SuzukiTrotter(order=2, reps=1)` on all 44 Pauli strings, at the stored $\Delta t = 0.1$.

The existing VQE, VQD, and Trotter results were read before they were overwritten. Those are the “before” numbers in `investigation.md`.

## 2. Locate the fields that fail

Exact diagonalization of the same `build_hamiltonian` on `np.linspace(0, 3, 81)`:

- Ground states on either side of $\varepsilon = 0.675 \to 0.7125$ have overlap 0, and the same is true for $\varepsilon = 1.7625 \to 1.800$.
- The first excited state changes character at several fields. The hard window for a shallow circuit is $\varepsilon \approx 0.49$ to $0.53$, where the neighbor overlap of `evecs[:, 1]` drops to about 0.65, the state is in the $Q = -1$ sector, and the half-chain entropy is about 0.6. A second entropy peak sits near $\varepsilon = 1.54$.

Line widths for Figs. 3–8 were read off the matplotlib layouts in `visualization/deps/paperfigs.py` and `visualization/deps/style.py`, converted to data coordinates. Those limits are the constants in `scripts/run_trotter.py`.

## 3. Ground state: restarts, then more draws

`field_scan` was changed so each later field minimizes Eq. (18) from the previous parameters and from five random vectors (`seed=42`), and keeps the lowest energy. $\varepsilon = 0$ stays the Table I vector.

A spot check at the first jump with only those five draws missed: best fidelity 0.950. A separate count of 24 random SLSQP starts at $\varepsilon = 0.7125$ found the new vacuum about one time in eight. The scan was then changed to draw further batches of five, up to 48 random starts, while the best fidelity is still below 0.99. The kept state is still the lowest energy. Fidelity only opens another batch.

Re-check at $\varepsilon = 0.7125$, the second jump, and a large $\varepsilon$, then:

```bash
python scripts/run_vqe.py
```

Minimum fidelity 0.9969, maximum energy error $1.143 \times 10^{-2}$, both at $\varepsilon = 0.7125$. All 81 rows passed. `scripts/run_vqe.py` prints that table and raises if a row fails. That floor is the `reps=3` scan. The production scan is `reps=4`; see [../paper_deliverables.md](../paper_deliverables.md).

## 4. Trotter: one grouped product, then a shorter step

Two repairs were compared on the quench fields, including $\varepsilon = 3$.

Commuting groups from `SparsePauliOp.group_commuting` gave three groups. One symmetric product of those groups at $\Delta t = 0.1$ did not pass: at $\varepsilon = 3$, final fidelity 0.980 and maximum energy-drift error 0.129.

The change that was kept is inside the existing step. `one_step_unitary` still receives the stored $\Delta t = 0.1$ and still uses `SuzukiTrotter(order=2)`. `reps` was raised until the line-width bar passed. Twelve repetitions were short of the Fig. 6 vacuum-fidelity width at $\varepsilon = 3$. Thirteen sat just under that width. Fourteen are what `SUZUKI_REPS` is set to, so each stored step is the Eq. (22) product at $0.1/14$.

```bash
python scripts/run_trotter.py
```

Final fidelities run from 1 at $\varepsilon = 0.5$ to 0.999880 at $\varepsilon = 3$. The largest observable error is 0.81 line widths, at $\varepsilon = 3$.

## 5. Excited state: restarts were not enough

`excited_scan` was given the same pattern as the ground-state scan: warm start plus five random starts, lowest Eq. (19) cost kept, further batches while the overlap with `evecs[:, 1]` is below 0.99, cap 48. `scripts/run_vqd.py` refuses to start unless every VQE fidelity on the grid is at least 0.99, and it builds the penalty reference from the VQE state at that same field.

Three `RealAmplitudes` repetitions cannot reach fidelity 0.99 on the first excited state near $\varepsilon = 0.6$ even when that overlap is the function being maximized (best about 0.988 over eight seeds). Four repetitions can, on the fields that were spot-checked first ($\varepsilon = 0.34, 0.6, 1.05, 1.8, 2.1$), and the lowest penalized cost was the exact first excited state there. The circuit was set to four repetitions and the full scan was started.

That scan passed the early fields and then missed after using the full budget of draws:

- $\varepsilon = 0.4875$: fidelity 0.953
- $\varepsilon = 0.525$: fidelity 0.912

The process was killed. The saved VQD file was still the old one, because the script writes only at the end.

## 6. What was tested before adding repetitions

These checks used the VQE state as the Eq. (19) reference. They were not left in the repository.

- Maximizing $|\langle \mathrm{evecs}[:,1] | \psi(\theta)\rangle|^2$ directly: four repetitions top out near 0.982 and 0.970 at the two failing fields. Six repetitions can reach about 0.995, but those parameters have a higher energy than the cost minima, which sit near fidelity 0.95.
- Plain cost minimization from random starts at 6, 8, and 10 repetitions: the lowest cost was still below fidelity 0.99 (about 0.95, 0.96, and 0.98).
- Starting cost minimization from a high-fidelity parameter vector: SLSQP walked off it, down to a lower cost and a lower fidelity.
- Walking the field in steps of 0.0125 from a good state at $\varepsilon = 0.4125$: the fidelity fell to 0.959 at 0.4875 and 0.638 at 0.525. A finer warm start follows the variational minimum, and that minimum was not yet `evecs[:, 1]`.

So the selection rule was left alone. The circuit had to be deep enough that the lowest cost was the first excited state.

## 7. Ten repetitions, with the same SLSQP

`trial_ansatz` returns `RealAmplitudes(num_qubits, reps=10)`. Nine repetitions, six random cost minima, and the VQE reference: zero of the six were above fidelity 0.99 at either failing field. Ten repetitions: one of six was, and that start had the lowest cost at both $\varepsilon = 0.4875$ and $0.525$ (fidelity 0.999, energy error about $10^{-3}$).

A ten-repetition SLSQP with finite-difference gradients is too slow for 81 fields and up to 48 starts. The trial state is real, because it is RY and CX applied to $|0\rangle$. `_trial_amplitudes` evaluates that circuit, and `_penalized_cost_and_gradient` returns Eq. (19) and $\partial C_1/\partial\theta$. Both were checked against Qiskit’s `Statevector` and against a finite-difference derivative of `penalized_cost` before they were used. `minimize_excited` passes that derivative to the same SLSQP (`maxiter=1000`, default `ftol`). One start then takes about a second.

Spot check of `excited_scan` itself, not a one-off script, on $\varepsilon = 0, 0.4875, 0.525, 0.7125, 1.05, 1.5375, 1.8, 3$. All eight fidelities were at least 0.993. Then:

```bash
python scripts/run_vqd.py
```

Minimum fidelity 0.9931 at $\varepsilon = 0.4875$, maximum energy error $2.43 \times 10^{-3}$ at that field, overlap with the VQE reference below $2 \times 10^{-9}$. All 81 rows passed. That scan was orthogonalized to the `reps=3` ground states. The production VQD scan uses the `reps=4` ground states; see [../paper_deliverables.md](../paper_deliverables.md).

## 8. Figures

```bash
python visualization/VQE.py
python visualization/VQD.py
python visualization/Trotter.py
```

Fig. 1 is written by the first two, Figs. 3–8 by the third. Each also writes `overlay.png`.

## What was not kept

- A new optimizer, a new ansatz family, a tensor network, or a different Hamiltonian.
- A larger penalty. The old excited-state failures already had overlap with the reference below $10^{-8}$.
- The commuting-group Trotter product. It missed the fidelity bar at $\Delta t = 0.1$.
- More than ten repetitions. Ten is the smallest depth at which the lowest Eq. (19) cost on the failing fields was `evecs[:, 1]`.
