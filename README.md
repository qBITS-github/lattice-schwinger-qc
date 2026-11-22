# Lattice Schwinger model

Reproduce Chen, Cheng & Guo,
[arXiv:2607.02894](https://arxiv.org/abs/2607.02894), then study how a strong
electric pulse changes the neutral vacuum's response to a weak, short pulse.

## Research Status

**Frozen reproduction:** Table I, Figs. 1–8 at N=8, and the Fig. 2 critical
field through N=18, including the N=8 VQE, VQD and Trotter comparisons.
The 12 datasets and 16 plots are locally read-only, with baseline plots under
`figures/baseline/`, and pinned by the [hash manifest](docs/baseline/reproduction_manifest.json).
See [reproduction instructions](docs/baseline/reproduction.md) and
[measured comparisons](docs/baseline/paper_deliverables.md).

**Driven study:** neutral-sector spectra, compact pulses, time-dependent
evolution, response predictions, matched controls, ramp preparation and
selected size checks. The active code produces the eight retained figures.
These are finite-lattice results; tunnelling,
continuum behaviour and hardware readiness remain unestablished.

- [Documentation index](docs/README.md) — baseline records, driven study and background notes.
- [Intuitive walkthrough](docs/driven/results.md) — each experiment, plot and conclusion.
- [Figures](docs/driven/figures.md) — eight results, with formal captions and PDF exports.
- [Running guide](docs/driven/guide.md) — commands, code map and output format.
- [Experiment specification](docs/driven/experiment.md) — conventions and observables.
- [Hamiltonian](docs/physics/building_hamiltonian.md),
  [observables](docs/physics/building_observables.md) and
  [study notes](docs/study/schwinger_model_intuitive.md).

New numerical outputs use `data/driven/<run_id>/`. The eight selected result
figures are in `figures/driven/results/`; regenerate them with
`python visualization/driven/plot.py`.
The [run inventory](data/driven/run_inventory.json) defines the nine retained experiments.
Use `python scripts/driven/run.py --experiment all --suffix _rerun` to recompute them,
then `python visualization/driven/plot.py --run-suffix _rerun` to plot that set.

```bash
python scripts/baseline/check_reproduction.py
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
```

## Repository layout

The baseline reproduction and driven extension use the same directory structure:

| Purpose | Baseline reproduction | Driven extension |
| --- | --- | --- |
| Physics and storage code | [src/baseline](src/baseline/) | [src/driven](src/driven/) |
| Simulation and verification commands | [scripts/baseline](scripts/baseline/) | [scripts/driven](scripts/driven/) |
| Plotting code | [visualization/baseline](visualization/baseline/) | [visualization/driven](visualization/driven/) |
| Saved data | [data/baseline](data/baseline/) | [data/driven](data/driven/) |
| Generated figures | [figures/baseline](figures/baseline/) | [figures/driven](figures/driven/) |
| Documentation | [docs/baseline](docs/baseline/) | [docs/driven](docs/driven/) |

[configs/driven](configs/driven/) holds experiment settings;
[tests](tests/) holds verification checks.
[References](references/README.md) records the papers and their relationship to
this project. Published baseline SVGs are in `references/baseline/figures/`;
historical investigation outputs are kept separately in `artifacts/baseline/`.

## Next Steps

- Explain the driven response mechanistically and assess its relation to prior work.
- Reproduce Appendix A at N=12 (Figs. 9–11).
- Extend to larger lattices with a tensor-network baseline, then ansatz and hardware benchmarks.
