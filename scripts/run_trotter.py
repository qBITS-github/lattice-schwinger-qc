"""Write the second-order Trotter quench dataset.

The evolution lives in ``src/trotter.py``. This file chooses the fields and
times used for Figs. 3–6 and writes ``data/Trotter/``.

Run from anywhere: ``python scripts/run_trotter.py``.
"""

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from schwinger_model import (
    build_hamiltonian,
    field_energy,
    local_charge,
    total_charge,
    vacuum_fidelity,
)
from trotter import SUZUKI_REPS, quench_compare

# Line width of reproduction.png, in the data coordinates of that panel.
_FIG3_LIMIT = 7.512019e-4
_FIG4_LIMIT = 1.386808e-2
_FIG8_LIMIT = 1.437154e-2
_FIG6_LIMIT = 8.141762e-3
_FIG7_LIMIT = 5.933533e-3
_FIG5_LIMIT = {
    0.5: 7.301044e-4,
    1.0: 3.285470e-3,
    1.5: 1.237231e-2,
    2.0: 2.878979e-2,
}
_QUENCH_EPS = (0.5, 1.0, 1.5, 2.0)


def _linewidth_ratio(eps, psi_trotter, psi_ed, psi0, times, N):
    """Largest plotted |O − O_exact| in units of that panel's line width.

    Figs. 3–5 and 8 use the quench fields 0.5, 1, 1.5, and 2. Figs. 6 and 7
    use every stored ε. run_trotter.py passes a field when this ratio is
    below 1 and the final-time state fidelity is at least 0.99.
    """
    ratios = []
    p_t = vacuum_fidelity(psi_trotter, psi0)
    p_e = vacuum_fidelity(psi_ed, psi0)
    ratios.append(float(np.max(np.abs(p_t - p_e)) / _FIG6_LIMIT))
    mask = (times >= 0.0) & (times <= 1.0)
    gamma_t = float(np.polyfit(times[mask], -np.log(p_t[mask]), 1)[0])
    gamma_e = float(np.polyfit(times[mask], -np.log(p_e[mask]), 1)[0])
    ratios.append(abs(gamma_t - gamma_e) / _FIG7_LIMIT)
    if any(np.isclose(eps, value) for value in _QUENCH_EPS):
        q_t = np.array([total_charge(psi, N) for psi in psi_trotter])
        q_e = np.array([total_charge(psi, N) for psi in psi_ed])
        ratios.append(float(np.max(np.abs(q_t - q_e)) / _FIG3_LIMIT))
        H = build_hamiltonian(N, eps=float(eps))
        e_t = np.array([float(np.real(np.vdot(psi, H @ psi))) for psi in psi_trotter])
        e_e = np.array([float(np.real(np.vdot(psi, H @ psi))) for psi in psi_ed])
        ratios.append(float(np.max(np.abs((e_t - e_t[0]) - (e_e - e_e[0]))) / _FIG3_LIMIT))
        he_t = np.array([field_energy(psi, N, eps=float(eps)) for psi in psi_trotter])
        he_e = np.array([field_energy(psi, N, eps=float(eps)) for psi in psi_ed])
        ratios.append(float(np.max(np.abs(he_t - he_e)) / _FIG5_LIMIT[float(eps)]))
        dq = 0.0
        d_spatial = 0.0
        for psi_t, psi_exact in zip(psi_trotter, psi_ed):
            qn_t = local_charge(psi_t, N)
            qn_e = local_charge(psi_exact, N)
            dq = max(dq, float(np.max(np.abs(qn_t - qn_e))))
            spatial_t = qn_t[0::2] + qn_t[1::2]
            spatial_e = qn_e[0::2] + qn_e[1::2]
            d_spatial = max(d_spatial, float(np.max(np.abs(spatial_t - spatial_e))))
        ratios.append(dq / _FIG8_LIMIT)
        ratios.append(d_spatial / _FIG4_LIMIT)
    return max(ratios)


def _print_result_table(eps, fidelity, error):
    """Print ε, final fidelity, observable error, and pass/fail.

    observable_error is the largest plotted deviation in units of the line
    width. Pass means the final fidelity is at least 0.99 and that deviation
    is below one line width.
    """
    print(f"{'eps':>8}  {'fidelity':>10}  {'observable_error':>18}  result")
    failed = False
    for e, f, err in zip(eps, fidelity, error):
        ok = float(f) >= 0.99 and float(err) < 1.0
        failed = failed or not ok
        print(f"{float(e):8.4f}  {float(f):10.6f}  {float(err):18.3e}  {'pass' if ok else 'fail'}")
    return failed


def main():
    """Quench the ε = 0 ED ground state and write both trajectories."""
    N = 8
    a = 1.0
    m = 1.0
    g = 1.0
    dt = 0.1
    # scripts/run_exact_diagonalization.py passes these fields to evolve() for Figs. 3–6.
    # 0.5, 1.0, 1.5, 2.0 are the quench set; 2.5 and 3.0 are added for the
    # vacuum-fidelity figure. times match that script: linspace(0, 12, 121).
    eps_values = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    times = np.linspace(0.0, 12.0, 121)

    print(
        f"Trotter: N={N}, dt={dt}, order=2, reps={SUZUKI_REPS}, "
        f"times 0..12 ({len(times)} points)"
    )
    psi0, psi_trotter, psi_ed, fidelity = quench_compare(
        N, a, m, g, eps_values, times, dt
    )

    out_dir = ROOT / "data" / "Trotter"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "quench_n8.npz"
    np.savez(
        out_path,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
        dt=np.float64(dt),
        order=np.int64(2),
        reps=np.int64(SUZUKI_REPS),
        eps=eps_values,
        times=times,
        psi0=psi0,
        psi_trotter=psi_trotter,
        psi_ed=psi_ed,
        fidelity=fidelity,
    )
    print(f"    saved {out_path}")
    print("    observable errors in units of the reproduction line width")
    ratios = np.empty(len(eps_values))
    for i, eps in enumerate(eps_values):
        print(f"    observables at eps={eps:g}", flush=True)
        ratios[i] = _linewidth_ratio(
            float(eps), psi_trotter[i], psi_ed[i], psi0, times, N
        )
    final_fidelity = fidelity[:, -1]
    failed = _print_result_table(eps_values, final_fidelity, ratios)
    if failed:
        raise RuntimeError("Trotter markers left the exact curves.")


if __name__ == "__main__":
    main()
