# Data for the eight driven figures

The inventory retains nine studies and 753 scan rows: the eight studies used
by the figures plus the strong-pulse calibration that selected their main
backgrounds. Their complete original numeric points, metadata and actual
source archives remain unchanged. Off-scope studies are outside the active
repository; they are not inputs to the current conclusions.

Each study has `summary.json`, `metadata.json`, numeric NPZ/JSON points,
`source.tar.gz`, and `invocations.jsonl`. The archived source identifies how
historical data were produced. Current code is a smaller implementation of
the selected workflow; use a fresh suffix when rerunning it. Resume requires
identical configuration and numerical code. `_cache/` is disposable.

See the [running guide](../../docs/driven/guide.md),
[definitions](../../docs/driven/experiment.md), and
[intuitive results walkthrough](../../docs/driven/results.md).
