"""Check the frozen reproduction hashes; optionally protect those files locally.

Run from any directory. No dataset is regenerated, and this script never updates
the manifest. Git does not preserve read-only permission bits: use --protect
after a fresh checkout to restore the local protection against overwriting.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protect", action="store_true", help="remove write permissions after all hashes pass")
    args = parser.parse_args()
    manifest = json.loads((ROOT / "docs" / "baseline" / "reproduction_manifest.json").read_text())
    failures = []
    paths = [ROOT / name for name in manifest["files"]]
    for name, expected in manifest["files"].items():
        path = ROOT / name
        if not path.is_file() or sha256(path.read_bytes()).hexdigest() != expected:
            failures.append(name)
    if failures:
        raise SystemExit("Frozen reproduction changed or missing:\n" + "\n".join(failures))
    if args.protect:
        for path in paths:
            path.chmod(path.stat().st_mode & ~0o222)
    protected = sum((path.stat().st_mode & 0o222) == 0 for path in paths)
    print(f"PASS: {len(paths)} frozen reproduction files match {manifest['source_commit']}.")
    print(f"Local read-only protection: {protected}/{len(paths)} files.")


if __name__ == "__main__":
    main()
