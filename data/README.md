# Frozen reproduction data and new driven outputs

**FROZEN REPRODUCTION — 19 October 2026, source `5497cd4a763b00f92269f7e8e650e2664d5c72f1`.**

The 12 existing NPZ files in `baseline/ed/`, `baseline/vqe/`, `baseline/vqd/`
and `baseline/trotter/` are the
closed reproduction record. Keep their contents and names unchanged. The exact
list and SHA-256 hashes are in [the manifest](../docs/baseline/reproduction_manifest.json).
They are locally read-only; Git does not preserve that permission setting.

New driven results belong only in `driven/<run_id>/`: numeric NPZ/JSON points
under `points/`, a readable `summary.json`, and `metadata.json` for provenance.
See the [running guide](../docs/driven/guide.md) and the
[experiment specification](../docs/driven/experiment.md). Other future static
reruns also need a separate output location; do not unlock this baseline.

The [run inventory](driven/run_inventory.json) lists the retained driven
studies and their exact configurations. Each retains its numeric points,
checksums and source archive. `driven/_cache/` is disposable; saved studies
do not depend on it. Superseded development scans are not retained.

Historical metadata and source archives retain the paths used when those runs
were made.
