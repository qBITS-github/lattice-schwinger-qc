# VQE fidelity gap

Lattice units are $a = m = g = 1$, $N = 8$, statevector, SLSQP, `RealAmplitudes` with reverse-linear entanglement. The exact energy at $\varepsilon = 0$ matches Table I: $-4.63805774$.

Two gaps were in hand. Table I ($\varepsilon = 0$) was fidelity $0.999604$ ($99.9604\%$) at `reps=3`. The paper's Table I fidelity is $99.9931\%$. The Fig. 1 scan minimum was $0.996927$ at $\varepsilon = 0.7125$, also at `reps=3`. The paper does not state `reps`, `ftol`, `maxiter`, or the restart count.

## Hypothesis

The Table I gap is the optimizer stopping early, or too few random starts. If neither moves the fidelity to about $99.99\%$, the ansatz is too shallow.

## Steps

1. Same `reps=3` circuit and the same initial points. SLSQP `ftol` from $10^{-6}$ to $10^{-12}$, `maxiter` from $1000$ to $5000$. Every exit was a success, and none hit `maxiter`. Table I fidelity went from $0.999604$ to $0.999681$.
2. Twenty random restarts at `reps=3`, with that tighter tolerance. Best fidelity $0.999720$.
3. Raise `reps`. At `reps=4`, ten restarts, `ftol` $10^{-12}$, `maxiter` $5000$, seed $42$, draw seed $4042$, restart index $7$: energy $-4.637983956989514$, relative error $1.591\times10^{-5}$, fidelity $0.9999877469280656$ ($99.9988\%$). That is past the paper's $99.9931\%$, so the depth ladder stopped there.
4. The Fig. 1 scan was then rerun at `reps=4` with a warm start plus five random restarts per field. No further batches were required.

## Result

The gap was ansatz depth. `reps=3` stayed near $99.97\%$ at $\varepsilon = 0$. `reps=4` is the production setting in `scripts/run_vqe.py`.

Table I is fidelity $99.9988\%$. On the 81-point Fig. 1 grid the minimum fidelity is $0.999599$ at $\varepsilon = 0.7125$. Both numbers are in [../paper_deliverables.md](../paper_deliverables.md).
