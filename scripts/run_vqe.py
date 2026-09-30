"""Write the Table I and Fig. 1 VQE datasets.

The minimization lives in ``src/vqe.py``. This file only chooses the run
and writes ``data/VQE/``.

Run from anywhere: ``python scripts/run_vqe.py``.
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

from vqe import (
    best_restart,
    ed_ground_pair,
    field_scan,
    random_parameters,
    relative_energy_error,
    trial_ansatz,
)


def _print_result_table(eps, fidelity, error, error_name):
    """Print ε, fidelity, energy error, and pass/fail against fidelity 0.99.

    run_vqe.py calls this on the Fig. 1 ground-state scan. Pass means the
    stored state overlaps evecs[:, 0] at that ε by at least 0.99.
    """
    print(f"{'eps':>8}  {'fidelity':>10}  {error_name:>14}  result")
    failed = False
    for e, f, err in zip(eps, fidelity, error):
        ok = float(f) >= 0.99
        failed = failed or not ok
        print(f"{float(e):8.4f}  {float(f):10.6f}  {float(err):14.3e}  {'pass' if ok else 'fail'}")
    return failed


def main():
    """Table I, then the Fig. 1 ground-state scan."""
    N = 8
    a = 1.0
    m = 1.0
    g = 1.0
    n_restarts = 5
    seed = 42
    # Fig. 1 spans ε ∈ [0, 3]. The paper does not state the sample count;
    # this is the grid scripts/run_exact_diagonalization.py stores in data/ED/condensate_spectrum_n8.npz.
    eps_values = np.linspace(0.0, 3.0, 81)
    out_dir = ROOT / "data" / "VQE"

    ansatz = trial_ansatz(N)
    print(
        f"Table I: N={N}, a=m=g=1, eps=0, RealAmplitudes reps={ansatz.reps}, "
        f"{n_restarts} restarts"
    )
    _H0, energy_ed, psi_ed = ed_ground_pair(N, a, m, g, 0.0)
    print(f"    ED E0 = {energy_ed:.8f}")
    initial_points = random_parameters(ansatz, n_restarts, seed)
    best, restart_energies = best_restart(ansatz, _H0, psi_ed, initial_points)
    rel = relative_energy_error(best["energy"], energy_ed)
    print(f"    best E_VQE = {best['energy']:.8f}")
    print(f"    relative error = {rel:.6e}")
    print(f"    fidelity = {best['fidelity']:.6f}")
    if best["energy"] < energy_ed - 1e-7:
        raise RuntimeError("VQE energy fell below diagonalize() evals[0].")

    out_dir.mkdir(parents=True, exist_ok=True)
    table_path = out_dir / "table1_n8.npz"
    np.savez(
        table_path,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
        eps=np.float64(0.0),
        reps=np.int64(ansatz.reps),
        seed=np.int64(seed),
        E_vqe=np.float64(best["energy"]),
        E0_ed=np.float64(energy_ed),
        relative_error=np.float64(rel),
        fidelity=np.float64(best["fidelity"]),
        theta=best["theta"],
        psi_vqe=best["psi"],
        psi_ed=psi_ed,
        restart_energies=restart_energies,
    )
    print(f"    saved {table_path}")

    print("Fig. 1 ground states, warm start plus random restarts...")
    energies, energies_ed, fidelities, thetas, psis = field_scan(
        ansatz, best["theta"], eps_values, N, a, m, g, n_restarts, seed
    )
    scan_path = out_dir / "ground_states_n8.npz"
    np.savez(
        scan_path,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
        reps=np.int64(ansatz.reps),
        eps=eps_values,
        E_vqe=energies,
        E0_ed=energies_ed,
        fidelity=fidelities,
        theta=thetas,
        psi_vqe=psis,
    )
    print(f"    saved {scan_path}")
    energy_error = np.abs(energies - energies_ed)
    print(f"    max |E_VQE - E0_ED| = {np.max(energy_error):.3e}")
    print(f"    min fidelity = {float(fidelities.min()):.4f}")
    failed = _print_result_table(eps_values, fidelities, energy_error, "energy_error")
    if failed:
        raise RuntimeError("Fig. 1 ground-state fidelity fell below 0.99.")


if __name__ == "__main__":
    main()
