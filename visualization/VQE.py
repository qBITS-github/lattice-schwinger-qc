"""VQE markers for paper Fig. 1.

Reads ``data/VQE/ground_states_n8.npz``. The condensate is Eq. (17) evaluated
on the saved state vectors with ``schwinger_model.chiral_condensate``. The
energies are the saved ⟨H⟩ values.

Run from anywhere: ``python visualization/VQE.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "deps"))

import numpy as np

from style import ROOT, param_line, skip_missing
from schwinger_model import chiral_condensate


def load_ground_states(path=None):
    """Fig. 1 VQE series: condensate and ground-state energy, in units of 1/a."""
    path = path or ROOT / "data" / "VQE" / "ground_states_n8.npz"
    if skip_missing(path, "VQE Fig. 1"):
        return None
    with np.load(path) as data:
        N = int(data["N"])
        a = float(data["a"])
        psis = np.array(data["psi_vqe"])
        gamma = np.array([chiral_condensate(psi, N, a=a) * a for psi in psis])
        return {
            "eps": np.array(data["eps"], dtype=float),
            "gamma": gamma,
            "E0": np.array(data["E_vqe"], dtype=float) * a,
            "E0_ed": np.array(data["E0_ed"], dtype=float) * a,
            "caption": param_line(N, a, float(data["m"]), float(data["g"])),
        }


def main():
    """Write Fig. 1 reproduction.png and overlay.png."""
    from paperfigs import write
    import overlay

    write([1])
    overlay.draw([1])


if __name__ == "__main__":
    main()
