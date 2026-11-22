"""Save driven scans, their numeric points and the source that produced them.

The shared runner passes plain JSON configurations to StudyStore. This module
handles output isolation, checksums, resume and caching; it runs no physics.
The frozen reproduction uses its original, separate readers and writers.
"""

from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def run_path(run_id):
    """Return a driven-data path, rejecting traversal and symlinks.

    A run ID is a single directory name, e.g. ``case_b_n8_20261019_001``.
    The fixed root separates every new study from the frozen reproduction.
    """
    if not isinstance(run_id, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", run_id):
        raise ValueError("run_id must contain only letters, digits, '_' and '-'")
    path = ROOT / "data" / "driven" / run_id
    if any(parent.is_symlink() for parent in (path, *path.parents)):
        raise ValueError(f"Output paths must not contain symlinks: {path}")
    return path


def _source_metadata():
    """Identify both HEAD and the actual code, including uncommitted additions.

    HEAD alone does not identify an uncommitted experiment. File hashes cover
    the Python code, driven configurations and dependency list used at save time.
    StudyStore archives these exact files to recover uncommitted code later.
    """
    def git(*args):
        return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()

    paths = [ROOT / "requirements.txt"]
    for directory, pattern in (("src", "*.py"), ("scripts", "*.py"),
                               ("visualization", "*.py"), ("configs/driven", "*.json")):
        paths.extend((ROOT / directory).rglob(pattern))
    return {
        "commit": git("rev-parse", "HEAD"),
        "dirty": bool(git("status", "--porcelain", "--untracked-files=all")),
        "file_sha256": {str(p.relative_to(ROOT)): sha256(p.read_bytes()).hexdigest()
                        for p in sorted(paths)},
    }


def json_text(value):
    """Canonical JSON used for exact specification comparisons and cache keys."""
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False,
                      default=lambda x: x.tolist() if isinstance(x, np.ndarray) else x.item())+"\n"


def specification_hash(value):
    return sha256(json_text(value).encode()).hexdigest()


class StudyStore:
    """One scan directory, small NPZ/JSON points, exact-match resume and cache.

    Completed points have a JSON sidecar with their NPZ hash. A crash before
    that sidecar leaves an incomplete point that can be recomputed on resume.
    Neither a changed configuration nor changed numerical code can be mixed
    into an existing run. Selection is recorded for every invocation.
    """
    def __init__(self, run_id, config, kind, resume=False, selection=None):
        self.data = run_path(run_id)
        self.source = _source_metadata()
        # Figure/config additions need not invalidate an existing physics scan.
        # Each scan's actual config is compared in full below.
        self.code = {p: h for p, h in self.source["file_sha256"].items()
                     if p.startswith(("src/", "scripts/")) or p == "requirements.txt"}
        signature = {"kind": kind, "config": config, "code": self.code}
        path = self.data/"metadata.json"
        if self.data.exists():
            if not resume or not path.is_file():
                raise FileExistsError("Run exists; use --resume with the identical specification")
            previous = json.loads(path.read_text())
            if previous["signature"] != json.loads(json_text(signature)):
                raise ValueError("Refusing to mix incompatible configuration or source code")
        else:
            self.data.mkdir(parents=True)
            packages = {}
            for package in ("numpy", "scipy", "matplotlib", "qiskit", "qiskit-algorithms"):
                try: packages[package] = version(package)
                except PackageNotFoundError: packages[package] = None
            path.write_text(json_text({"schema_version": 2, "study": "driven", "run_id": run_id,
                "signature": signature, "source": self.source, "python": sys.version,
                "packages": packages, "created_utc": datetime.now(timezone.utc).isoformat(),
                "units": "Repository numerical units; hbar=1. Config times/frequencies are explicit numerical values."}))
            # HEAD does not contain this uncommitted implementation. Archive it.
            import tarfile
            with tarfile.open(self.data/"source.tar.gz", "w:gz") as archive:
                for name in self.source["file_sha256"]:
                    archive.add(ROOT/name, arcname=name)
        (self.data/"points").mkdir(exist_ok=True)
        with (self.data/"invocations.jsonl").open("a") as stream:
            stream.write(json.dumps({"utc": datetime.now(timezone.utc).isoformat(),
                                     "selection": selection, "resume": resume})+"\n")

    def point(self, specification, compute):
        """Return exact-matching saved data or compute/save it once.

        specification includes preparation, pulse, physical output times,
        solver settings and all extra analysis options. Cache sharing is
        conservative: identical numerical code is also required.
        """
        key = {"specification": specification, "code": self.code}
        identifier = specification_hash(key)
        cache = ROOT/"data"/"driven"/"_cache"
        cache.mkdir(exist_ok=True)
        if cache.is_symlink():
            raise ValueError("Cache must not be a symlink")
        local = self.data/"points"/identifier
        cached = cache/identifier
        for prefix in (local, cached):
            sidecar = prefix.with_suffix(".json")
            numeric = prefix.with_suffix(".npz")
            if sidecar.is_file() and numeric.is_file():
                metadata = json.loads(sidecar.read_text())
                if metadata["key"] != json.loads(json_text(key)):
                    raise ValueError("Cache specification mismatch")
                if sha256(numeric.read_bytes()).hexdigest() != metadata["npz_sha256"]:
                    raise ValueError(f"Corrupt saved point: {numeric}")
                with np.load(numeric, allow_pickle=False) as saved:
                    arrays = {k: saved[k] for k in saved.files}
                if prefix == cached:
                    import shutil
                    shutil.copyfile(numeric, local.with_suffix(".npz"))
                    shutil.copyfile(sidecar, local.with_suffix(".json"))
                return arrays, metadata["result"], identifier
        arrays, result = compute()
        if any(np.asarray(a).dtype.kind not in "biufc" for a in arrays.values()):
            raise ValueError("Only numeric arrays may enter saved points")
        np.savez_compressed(local.with_suffix(".npz"), **arrays)
        metadata = {"key": key, "result": result,
                    "npz_sha256": sha256(local.with_suffix(".npz").read_bytes()).hexdigest()}
        local.with_suffix(".json").write_text(json_text(metadata))
        import shutil
        shutil.copyfile(local.with_suffix(".npz"), cached.with_suffix(".npz"))
        shutil.copyfile(local.with_suffix(".json"), cached.with_suffix(".json"))
        return arrays, result, identifier

    def table(self, rows):
        """Save all selected rows, retaining previously completed scan points."""
        path = self.data/"summary.json"
        previous = json.loads(path.read_text()) if path.exists() else []
        merged = {str(row["index"]): row for row in previous}
        merged.update({str(row["index"]): row for row in rows})
        path.write_text(json_text(sorted(merged.values(), key=lambda r: r["index"])))
