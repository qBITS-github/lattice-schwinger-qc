"""Write the second-order Trotter quench dataset.

The evolution lives in ``src/trotter.py``. This file chooses the fields and
times used for Figs. 3–6 and writes ``data/Trotter/``.

Run from anywhere: ``python scripts/run_trotter.py``.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from trotter import quench_compare


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

    print(f"Trotter: N={N}, dt={dt}, order=2, reps=1, times 0..12 ({len(times)} points)")
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
        reps=np.int64(1),
        eps=eps_values,
        times=times,
        psi0=psi0,
        psi_trotter=psi_trotter,
        psi_ed=psi_ed,
        fidelity=fidelity,
    )
    print(f"    saved {out_path}")


if __name__ == "__main__":
    main()
