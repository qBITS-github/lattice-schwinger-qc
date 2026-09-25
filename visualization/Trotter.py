"""Second-order Trotter markers for paper Figs. 3–8.

Reads the state trajectories in ``data/Trotter/quench_n8.npz`` and evaluates
the same observables as the exact-diagonalization datasets: total charge and
energy drift, spatial and site charge, electric-field energy, vacuum
fidelity, and the early-time decay rate.

Run from anywhere: ``python visualization/Trotter.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "deps"))

import numpy as np

from style import ROOT, param_line, skip_missing
from schwinger_model import (
    build_hamiltonian,
    field_energy,
    local_charge,
    total_charge,
    vacuum_fidelity,
)


def load_quench(path=None):
    """State trajectories saved by ``scripts/run_trotter.py``."""
    path = path or ROOT / "data" / "Trotter" / "quench_n8.npz"
    if skip_missing(path, "Trotter Figs. 3–8"):
        return None
    with np.load(path) as data:
        return {
            "N": int(data["N"]),
            "a": float(data["a"]),
            "m": float(data["m"]),
            "g": float(data["g"]),
            "eps": np.array(data["eps"], dtype=float),
            "times": np.array(data["times"], dtype=float),
            "psi0": np.array(data["psi0"]),
            "psi_trotter": np.array(data["psi_trotter"]),
            "caption": param_line(int(data["N"]), float(data["a"]), float(data["m"]), float(data["g"])),
        }


def measure(quench):
    """Observables of the Trotter trajectory, in the same units as the ED plots.

    Charge, field energy, and fidelity call the functions in
    ``schwinger_model``. The energy drift is ⟨H(ε)⟩(t) − ⟨H(ε)⟩(0), the same
    difference ``scripts/run_exact_diagonalization.py`` stores as ``delta_E``.
    """
    N = quench["N"]
    a = quench["a"]
    m = quench["m"]
    g = quench["g"]
    eps = quench["eps"]
    psis = quench["psi_trotter"]
    psi0 = quench["psi0"]
    n_eps, n_times, _dim = psis.shape
    Q_N = np.empty((n_eps, n_times))
    delta_E = np.empty((n_eps, n_times))
    H_E = np.empty((n_eps, n_times))
    P_vac = np.empty((n_eps, n_times))
    q_n = np.empty((n_eps, n_times, N))
    Q_i = np.empty((n_eps, n_times, N // 2))
    for k, eps_k in enumerate(eps):
        print(f"    observables at eps={eps_k:g}")
        H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps_k))
        energies = np.empty(n_times)
        P_vac[k] = vacuum_fidelity(psis[k], psi0)
        for i, psi in enumerate(psis[k]):
            Q_N[k, i] = total_charge(psi, N)
            energies[i] = float(np.real(np.vdot(psi, H @ psi)))
            H_E[k, i] = field_energy(psi, N, a=a, g=g, eps=float(eps_k))
            q = local_charge(psi, N)
            q_n[k, i] = q
            Q_i[k, i] = q[0::2] + q[1::2]
        delta_E[k] = energies - energies[0]
    return {
        "t_over_a": quench["times"] / a,
        "times": quench["times"],
        "eps": eps,
        "Q_N": Q_N,
        "delta_E": delta_E * a,
        "H_E": H_E * a,
        "P_vac": P_vac,
        "q_n": q_n,
        "Q_i": Q_i,
        "a": a,
        "caption": quench["caption"],
    }


def main():
    """Write reproduction.png and overlay.png for Figs. 3–8."""
    from paperfigs import write
    import overlay

    write([3, 4, 5, 6, 7, 8])
    overlay.draw([3, 4, 5, 6, 7, 8])


if __name__ == "__main__":
    main()

