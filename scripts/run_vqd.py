"""Write the Fig. 1(b) VQD excited-state dataset.

The minimization lives in ``src/vqd.py``. This file loads the VQE ground
states, runs the scan, and writes ``data/VQD/``. It does not run VQE.

Run from anywhere: ``python scripts/run_vqd.py``. Requires ``python scripts/run_vqe.py``
first.
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

from vqd import excited_scan, trial_ansatz


def _print_result_table(eps, fidelity, error, error_name):
    """Print ε, fidelity, energy error, and pass/fail against fidelity 0.99.

    run_vqd.py calls this on the Fig. 1 excited-state scan. Pass means the
    stored state overlaps evecs[:, 1] at that ε by at least 0.99.
    """
    print(f"{'eps':>8}  {'fidelity':>10}  {error_name:>14}  result")
    failed = False
    for e, f, err in zip(eps, fidelity, error):
        ok = float(f) >= 0.99
        failed = failed or not ok
        print(f"{float(e):8.4f}  {float(f):10.6f}  {float(err):14.3e}  {'pass' if ok else 'fail'}")
    return failed


def main():
    """Load ψ0^VQE, compute E1(ε), write data/VQD/e1_scan_n8.npz."""
    n_restarts = 5
    seed = 42
    ground_path = ROOT / "data" / "VQE" / "ground_states_n8.npz"
    if not ground_path.is_file():
        raise FileNotFoundError(
            f"{ground_path} is missing. Run python scripts/run_vqe.py before scripts/run_vqd.py."
        )
    with np.load(ground_path) as data:
        N = int(data["N"])
        a = float(data["a"])
        m = float(data["m"])
        g = float(data["g"])
        eps = np.array(data["eps"], dtype=float)
        psi_vqe = np.array(data["psi_vqe"])
        ground_fidelity = np.array(data["fidelity"], dtype=float)

    if float(ground_fidelity.min()) < 0.99:
        raise RuntimeError(
            "VQE state is not the ground state on the whole field grid. "
            "Run python scripts/run_vqe.py before scripts/run_vqd.py."
        )
    reps = trial_ansatz(N).reps
    print(
        f"VQD E1 scan: N={N}, RealAmplitudes reps={reps}, "
        f"{len(eps)} fields, warm start plus {n_restarts} restarts"
    )
    energies, energies_ed, ground_ed, fidelities, overlaps, betas, thetas, psis = excited_scan(
        N, a, m, g, eps, psi_vqe, n_restarts, seed
    )
    out_dir = ROOT / "data" / "VQD"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "e1_scan_n8.npz"
    np.savez(
        out_path,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
        reps=np.int64(reps),
        eps=eps,
        E1_vqd=energies,
        E1_ed=energies_ed,
        E0_ed=ground_ed,
        fidelity=fidelities,
        overlap_with_vqe_ground=overlaps,
        beta=betas,
        theta=thetas,
        psi_vqd=psis,
    )
    print(f"    saved {out_path}")
    energy_error = np.abs(energies - energies_ed)
    print(f"    max |E1_VQD - E1_ED| = {np.max(energy_error):.3e}")
    print(f"    min fidelity vs evecs[:, 1] = {float(fidelities.min()):.4f}")
    print(f"    max |⟨ψ0^VQE|ψ1^VQD⟩|^2 = {float(np.max(overlaps)):.3e}")
    failed = _print_result_table(eps, fidelities, energy_error, "energy_error")
    if failed:
        raise RuntimeError("Fig. 1 excited-state fidelity fell below 0.99.")


if __name__ == "__main__":
    main()
