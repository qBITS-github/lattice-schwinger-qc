"""One entry point for the nine experiments behind the eight retained figures.

The inventory supplies configurations and the actual historical row selection.
A suffix puts reruns in fresh data/driven directories; it never overwrites the
saved evidence or any baseline reproduction output. See docs/driven/guide.md.
"""
import argparse
import json
import os
from pathlib import Path
import sys

# Small dense eigensystems are faster without many BLAS worker threads.
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(name, "1")

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"src"))
from driven.study import run


def main():
    inventory = json.loads((ROOT/"data/driven/run_inventory.json").read_text())
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment", choices=[*inventory["runs"], "all"], default="case_a_n8_final")
    parser.add_argument("--suffix", default="_rerun", help="fresh output suffix (default: _rerun)")
    parser.add_argument("--limit", type=int, help="maximum rows per experiment for a trial run")
    parser.add_argument("--resume", action="store_true", help="resume identical code and settings")
    args = parser.parse_args()
    if args.limit is not None and args.limit < 1: parser.error("--limit must be positive")
    if not args.suffix or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in args.suffix):
        parser.error("--suffix must be nonempty and contain only letters, digits, '_' or '-'")
    names = inventory["runs"] if args.experiment == "all" else [args.experiment]
    for name in names:
        recipe = inventory["recipes"][name]
        path = run(ROOT/inventory["configurations"][name], recipe["kind"], name+args.suffix,
                   recipe["indices"], args.limit, args.resume)
        print(f"Saved {path}", flush=True)


if __name__ == "__main__":
    main()
