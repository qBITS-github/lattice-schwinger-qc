"""Exact-diagonalization datasets for the paper figures.

The loaders below read ``data/ED/``. Running this file writes
``reproduction.png`` and ``overlay.png`` for Figs. 1–8. Quantum markers are
included when ``data/VQE``, ``data/VQD``, and ``data/Trotter`` exist.

Run from anywhere: ``python visualization/ED.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "deps"))

import numpy as np

from style import ROOT, param_line, skip_missing

ED_DIR = ROOT / "data" / "ED"


def _open(path, what):
    path = Path(path)
    if skip_missing(path, what):
        return None
    return np.load(path)


def load_condensate_spectrum(path=None):
    """Fig. 1 curves: condensate and the lowest two energies, in units of 1/a."""
    data = _open(path or ED_DIR / "condensate_spectrum_n8.npz", "ED Fig. 1")
    if data is None:
        return None
    with data:
        a = float(data["a"])
        return {
            "eps": np.array(data["eps"]),
            "gamma": np.array(data["gamma"]) * a,
            "E0": np.array(data["E0"]) * a,
            "E1": np.array(data["E1"]) * a,
            "caption": param_line(int(data["N"]), a, float(data["m"]), float(data["g"])),
        }


def load_critical_scaling(path=None):
    """Fig. 2: critical field versus 1/N."""
    data = _open(path or ED_DIR / "critical_field_scaling.npz", "ED Fig. 2")
    if data is None:
        return None
    with data:
        return {
            "N": np.array(data["N"], dtype=float),
            "eps_c": np.array(data["eps_c"], dtype=float),
            "a": float(data["a"]),
            "m": float(data["m"]),
            "g": float(data["g"]),
        }


def load_conserved(path=None):
    """Fig. 3: total charge and energy drift of the exact quench."""
    data = _open(path or ED_DIR / "quench_conserved_n8.npz", "ED Fig. 3")
    if data is None:
        return None
    with data:
        a = float(data["a"])
        return {
            "t_over_a": np.array(data["times"]) / a,
            "eps": np.array(data["eps"]),
            "Q_N": np.array(data["Q_N"]),
            "delta_E": np.array(data["delta_E"]) * a,
            "caption": param_line(int(data["N"]), a, float(data["m"]), float(data["g"])),
        }


def load_charges(path=None):
    """Figs. 4 and 8: spatial-point charge and site charge."""
    data = _open(path or ED_DIR / "charge_dynamics_n8.npz", "ED Figs. 4 and 8")
    if data is None:
        return None
    with data:
        a = float(data["a"])
        return {
            "t_over_a": np.array(data["times"]) / a,
            "eps": np.array(data["eps"]),
            "Q_i": np.array(data["Q_i"]),
            "q_n": np.array(data["q_n"]),
            "caption": param_line(int(data["N"]), a, float(data["m"]), float(data["g"])),
        }


def load_field_energy(path=None):
    """Fig. 5: electric-field energy, in units of 1/a."""
    data = _open(path or ED_DIR / "field_energy_n8.npz", "ED Fig. 5")
    if data is None:
        return None
    with data:
        a = float(data["a"])
        return {
            "t_over_a": np.array(data["times"]) / a,
            "eps": np.array(data["eps"]),
            "H_E": np.array(data["H_E"]) * a,
            "caption": param_line(int(data["N"]), a, float(data["m"]), float(data["g"])),
        }


def load_vacuum_fidelity(path=None):
    """Fig. 6: vacuum fidelity of the exact quench."""
    data = _open(path or ED_DIR / "vacuum_fidelity_n8.npz", "ED Fig. 6")
    if data is None:
        return None
    with data:
        a = float(data["a"])
        return {
            "t_over_a": np.array(data["times"]) / a,
            "eps": np.array(data["eps"]),
            "P_vac": np.array(data["P_vac"]),
            "a": a,
            "caption": param_line(int(data["N"]), a, float(data["m"]), float(data["g"])),
        }


def load_decay_rate(path=None):
    """Fig. 7: early-time effective decay rate, in units of 1/a."""
    data = _open(path or ED_DIR / "decay_rate_n8.npz", "ED Fig. 7")
    if data is None:
        return None
    with data:
        a = float(data["a"])
        return {
            "eps": np.array(data["eps"]),
            "gamma_eff": np.array(data["gamma_eff"]) * a,
            "N": int(data["N"]),
            "a": a,
            "t_min": float(data["t_over_a_min"]),
            "t_max": float(data["t_over_a_max"]),
        }


def main():
    """Write reproduction.png and overlay.png for Figs. 1–8."""
    from paperfigs import write
    import overlay

    write()
    overlay.draw()


if __name__ == "__main__":
    main()

