"""VQD markers for paper Fig. 1(b).

Reads ``data/VQD/e1_scan_n8.npz``. The markers are the saved ⟨H⟩ of the
penalized minimization, Eq. (20).

Run from anywhere: ``python visualization/VQD.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "deps"))

import numpy as np

from style import ROOT, param_line, skip_missing


def load_excited_scan(path=None):
    """Fig. 1(b) VQD series: first excited energy, in units of 1/a."""
    path = path or ROOT / "data" / "VQD" / "e1_scan_n8.npz"
    if skip_missing(path, "VQD Fig. 1(b)"):
        return None
    with np.load(path) as data:
        a = float(data["a"])
        return {
            "eps": np.array(data["eps"], dtype=float),
            "E1": np.array(data["E1_vqd"], dtype=float) * a,
            "E1_ed": np.array(data["E1_ed"], dtype=float) * a,
            "E0_ed": np.array(data["E0_ed"], dtype=float) * a,
            "caption": param_line(int(data["N"]), a, float(data["m"]), float(data["g"])),
        }


def main():
    """Write Fig. 1 reproduction.png and overlay.png."""
    from paperfigs import write
    import overlay

    write([1])
    overlay.draw([1])


if __name__ == "__main__":
    main()
