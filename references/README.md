# References and relationship to this project

The published baseline figure SVGs used by the overlay code are in
[`baseline/figures/`](baseline/figures/).

The frozen baseline reproduces selected results from Chen–Cheng–Guo. The
driven extension adapts protocols from the literature; it does **not**
currently reproduce another paper's numerical results under matching conditions.

## Primary references

- **Chen, Cheng and Guo — [Digital Quantum Simulation of Nonequilibrium Dynamics in the Schwinger Model under a Strong External Electric Field](https://arxiv.org/abs/2607.02894).**
  Source of the baseline finite open-boundary Hamiltonian, static scans and
  constant-field quench experiments. Selected results are reproduced; see the
  [quantitative comparison and limitations](../docs/baseline/paper_deliverables.md).
  The time-dependent pulse experiments extend this setup.

- **Ghim and Honda — [Digital Quantum Simulation for Spectroscopy of Schwinger Model](https://arxiv.org/abs/2404.14788).**
  Source of the spectroscopy idea: scan an oscillating perturbation and measure
  loss of vacuum survival. Our weak-sine experiment adapts this method, with a
  different drive operator, waveform and setup. Their published spectra are
  not reproduced. Read Sections 3–4.

- **Aleksandrov, Sevostyanov and Shabaev — [Schwinger particle production: Rapid switch off of the external field versus dynamical assistance](https://arxiv.org/abs/2210.15626).**
  Source of the two-Sauter-pulse protocol and assistance/switching comparisons.
  Their observable is pair density in prescribed-field 3+1-dimensional QED;
  ours is vacuum-excitation probability in a finite interacting 1+1-dimensional
  lattice. Our additive-response difference is not their enhancement factor,
  which also includes a locally constant field approximation comparison.
  Read Section II, Section V.1 and Eq. (20). This is an adaptation, not a
  numerical replication.

## Supporting references — not reproduced here

- **Schützhold, Gies and Dunne — [Dynamically assisted Schwinger mechanism](https://arxiv.org/abs/0807.0754).**
  Original motivation for combining a strong slow field with a weak fast field.
- **Huang and Taya — [Spin-dependent dynamically assisted Schwinger mechanism](https://arxiv.org/abs/1904.08200).**
  Precedent for treating the strong field nonperturbatively and the weak field
  perturbatively. Our response expansion does not reproduce their spin-resolved
  calculation.
- **Buyens et al. — [Real-time simulation of the Schwinger effect with Matrix Product States](https://arxiv.org/abs/1612.00739).**
  Reference for interacting Schwinger-model quenches, local dynamics and
  interpretation using excitation spectra.
- **D’Anna, Krstić Marinković and Pinto Barros — [Adiabatic state preparation for digital quantum simulations of QED in 1+1D](https://arxiv.org/abs/2411.01079).**
  Reference for charge-sector crossings and limitations of adiabatic preparation.

## What validates the driven results?

The [experiment specification](../docs/driven/experiment.md) defines our actual
model, pulses and observables. Validation uses weak-field perturbation theory,
Hamiltonian comparisons, numerical refinement and independent propagation.
The remaining controls, state populations, preparation and size studies are
calculations for this project, not reproductions of published figures.

There is currently **no matched published dataset validating the complete
driven extension**. In particular, the gap near `2.658 g` is computed from our
own neutral-sector Hamiltonian at `N=8, a=m=g=1`; it is not taken from the
Ghim–Honda spectrum. These references establish provenance and context,
not novelty or independent confirmation of our new results.
