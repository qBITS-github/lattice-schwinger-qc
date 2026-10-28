# The eight driven figures

These figures address how a strong pulse modifies neutral-vacuum response to
a weak pulse. They are separate from the frozen baseline reproduction.

Read the [intuitive walkthrough](results.md) for the experiments and
conclusions, or the [formal captions](../../figures/driven/results/captions.md)
for definitions and plotting conventions.

| Figure | Scientific purpose |
| --- | --- |
| [1. Accessible states](../../figures/driven/results/01_accessible_states.pdf) | Identify charge-allowed gaps and coupling strengths. |
| [2. Weak probe](../../figures/driven/results/02_weak_probe.pdf) | Validate the accessible resonance and small-amplitude response; expose approximation limits. |
| [3. Two-pulse result](../../figures/driven/results/03_two_pulse_result.pdf) | Establish absolute excitation and positive nonadditivity on the main grid. |
| [4. Matched controls](../../figures/driven/results/04_matched_controls.pdf) | Test peak height, input budget, timing and polarity. |
| [5. Response and states](../../figures/driven/results/05_response_and_states.pdf) | Compare the strong-pulse response expansion and resolve the complete signed final-state excess. |
| [6. Background preparation](../../figures/driven/results/06_background_preparation.pdf) | Test ramp fidelity and show frequency-dependent enhancement and suppression. |
| [7. Size dependence](../../figures/driven/results/07_size_dependence.pdf) | Test the selected protocol through N=16 and compare central/edge local response. |
| [8. Numerical checks](../../figures/driven/results/08_numerical_checks.pdf) | Compare the representative effect with independent-method and tail changes. |

Each figure has a PNG preview and PDF export of the same data. Figures 6–8
can serve as supporting figures rather than being crowded into a manuscript's
main narrative. The [plot manifest](../../figures/driven/results/plot_manifest.json)
records selections and exact input/output hashes; the
[inventory](../../data/driven/run_inventory.json) identifies the nine retained
studies, including the calibration that selected the main backgrounds.

```bash
python visualization/driven/plot.py
```

The [running guide](guide.md) covers complete recomputation with the
same selected experiments. New figures should answer a distinct scientific
question supported by a completed experiment.
