"""Write the Fig. 1(b) VQD excited-state dataset.

The minimization lives in ``src/vqd.py``. This file loads the VQE ground
states, runs the scan, and writes ``data/VQD/``. It does not run VQE.

Run from anywhere: ``python scripts/run_vqd.py``. Requires ``python scripts/run_vqe.py``
first.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from vqd import excited_scan


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
        reps = int(data["reps"])
        eps = np.array(data["eps"], dtype=float)
        psi_vqe = np.array(data["psi_vqe"])

    print(
        f"VQD E1 scan: N={N}, RealAmplitudes reps={reps}, "
        f"{len(eps)} fields, {n_restarts} restarts at eps={eps[0]:g}"
    )
    energies, energies_ed, ground_ed, fidelities, overlaps, betas, thetas, psis = excited_scan(
        N, a, m, g, reps, eps, psi_vqe, n_restarts, seed
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
    print(f"    max |E1_VQD - E1_ED| = {np.max(np.abs(energies - energies_ed)):.3e}")
    print(f"    min fidelity vs evecs[:, 1] = {float(fidelities.min()):.4f}")
    print(f"    max |⟨ψ0^VQE|ψ1^VQD⟩|^2 = {float(np.max(overlaps)):.3e}")


if __name__ == "__main__":
    main()
