# Running the eight-figure driven study

The active driven workflow contains one simulation entry point, seven small
physics/storage modules, nine experiment configurations, and one renderer.
The original baseline scripts remain unchanged. All commands below run from
the repository root using the existing `requirements.txt` environment.

## Render the saved results

```bash
python visualization/driven/plot.py
```

This reads saved data only. It writes exactly eight PNG/PDF pairs under
`figures/driven/results/`, plus the plotting-source snapshot and hash manifest.
The [walkthrough](results.md) explains the figures; their
[formal captions](../../figures/driven/results/captions.md) define every panel.
PDFs use embedded TrueType fonts and vector lines/text, with dense spectral
scatter rasterized at 300 dpi. Layouts target a 7.2-inch two-column width.
The baseline plotting style is not modified.

## Run the retained experiments

Try one actual main-scan point:

```bash
python scripts/driven/run.py --experiment case_a_n8_final --suffix _trial --limit 1
```

Recompute the complete selected experiment set into fresh directories, then
render those results:

```bash
python scripts/driven/run.py --experiment all --suffix _rerun
python visualization/driven/plot.py --run-suffix _rerun
```

The full set includes long scans and N=16 sparse evolution; it is substantially
more work than a trial point. `--limit 1` limits each selected experiment to
one row. Incomplete trial sets cannot generate the full eight figures.

Resume an interrupted experiment only when its code and configuration match:

```bash
python scripts/driven/run.py --experiment case_a_n8_final --suffix _trial --resume
```

The default suffix is `_rerun`. Existing runs are rejected unless `--resume`
is explicit and the complete configuration and numerical code match. Use a
new suffix after editing physics code. Historical source is preserved inside
each saved study; its original uncommitted implementation can be recovered
from `source.tar.gz` without replacing the current checkout.

The [inventory](../../data/driven/run_inventory.json) is the single recipe list:

| Experiment name | Saved rows | Purpose |
| --- | ---: | --- |
| `spectrum_n8_refined_final` | 79 | Accessible neutral states and charge comparison |
| `weak_checks_n8_v2` | 12 | Small-amplitude and inaccessible-gap checks |
| `calibration_n8_v2` | 15 | Selection of strong backgrounds 0.6 and 1 |
| `case_a_n8_final` | 150 | Main two-Sauter-pulse scan |
| `controls_n8_final` | 12 | Matched controls, response and final populations |
| `independent_case_a_final` | 2 | Independent midpoint steps |
| `case_b_zero_n8_final` | 438 | Zero-background frequency response |
| `plateau_n8_final` | 36 | Ramp preparation and background probe |
| `size_case_a_final` | 9 | Fixed-protocol size/spatial comparison |

There are 753 saved scan rows. Calibration supports the selection of the
main protocols; it does not add an extra figure. A null `indices` field means
the full configuration grid. The control recipe lists its exact 12 historical
indices, rather than silently launching its entire possible grid. The scan
is a Cartesian product of JSON lists, with the first coordinate changing
slowest. Figure-specific subsets are recorded in the plot manifest.

## How the components connect

```text
scripts/driven/run.py + run_inventory.json + configs/driven/*.json
    -> src/driven/study.py
       -> model.py: sector Hamiltonian, vacuum, excitation spectrum
       -> pulses.py: actual compact fields, derivatives and boundaries
       -> evolution.py: chronological adaptive/midpoint evolution
       -> observables.py: excitation, energy, occupation and checks
       -> response.py: weak-field and strong-pulse response predictions
       -> io.py: numeric points, summaries, checksums and provenance
    -> data/driven/<experiment><suffix>/
    -> visualization/driven/plot.py -> eight figures + manifest
```

`study.py` runs the same initial vacuum with no field, strong only,
weak only and both fields. It also executes only the configured controls and
refinements needed by the retained studies. The independent midpoint method
is kept because it directly supports Figure 8. Driven circuit benchmarks,
stationary-probe protocols and unrelated parameter scans are outside this
active workflow. The baseline VQE/VQD/Trotter code remains intact.

Each saved study contains `summary.json`, `metadata.json`, `points/*.npz`
with JSON sidecars, `source.tar.gz`, and `invocations.jsonl`. Load NPZ with
`allow_pickle=False`. Trajectories distinguish output `times` from actual
`integration_times`; large sectors stream observables and retain endpoint
states. `data/driven/_cache/` is disposable: studies contain their own points.
Simulations create no empty figure directories.

## Verify the implementation and frozen baseline

```bash
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
python scripts/baseline/check_reproduction.py
```

The tests compare the sector Hamiltonian with the independent baseline
builder, check spectral anchors and Lanczos moments, compare propagation with
constant-field evolution, test short-pulse resolution and convergence, check
work and response calculations, and exercise provenance/resume/corruption
handling. The baseline checker verifies the original 12 data and 16 image
hashes. It does not regenerate or overwrite them.
