# NOTE: Table I and all "QC" / VQE / VQD / Trotter markers from
# Chen, Cheng & Guo (arXiv:2607.02894) are intentionally NOT reproduced.
# This pipeline is exact-diagonalization only.

"""Exact-diagonalization datasets for the lattice Schwinger model.

Each section writes one ``.npz`` dataset under ``data/ED/``. ``visualization/ED.py``
only loads those files. No figure is produced here.

Run from anywhere: ``python scripts/run_exact_diagonalization.py``. Paths are relative to the
repository root, not the working directory.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np

from schwinger_model import (
    build_hamiltonian,
    chiral_condensate,
    diagonalize,
    evolve,
    field_energy,
    local_charge,
    spatial_point_charge,
    total_charge,
    vacuum_fidelity,
)


def _write_npz(path, **arrays):
    """Write one dataset consumed by the exact-diagonalization plots in ``visualization/ED.py``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **arrays)


def vacuum_state(N, a, m, g):
    """Eps = 0 ground state of Eqs. (9)-(11), the reference vacuum for Eq. (16) and Figs. 3-6."""
    H0 = build_hamiltonian(N, a=a, m=m, g=g, eps=0.0)
    _evals, evecs = diagonalize(H0, 1)
    return np.asarray(evecs[:, 0]).reshape(-1)


def condensate_scan(N, a, m, g, eps_values):
    """Eqs. (9)-(11) and (17): ground-state Gamma, E0, and E1 versus eps (Figs. 1 and 2)."""
    eps_values = np.asarray(eps_values, dtype=float)
    gamma = np.empty(eps_values.size)
    E0 = np.empty(eps_values.size)
    E1 = np.empty(eps_values.size)
    for i, eps in enumerate(eps_values):
        H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
        evals, evecs = diagonalize(H, 2)
        gamma[i] = chiral_condensate(evecs[:, 0], N, a=a)
        E0[i] = evals[0]
        E1[i] = evals[1]
    return gamma, E0, E1


# ---- Section 1: Chiral condensate & spectrum vs external field (paper Fig. 1) ----

def condensate_and_spectrum(N, a, m, g, eps_values, out_path):
    """Paper Fig. 1 / Eq. (17): ground-state chiral condensate and the lowest two energies versus eps."""
    eps_values = np.asarray(eps_values, dtype=float)
    gamma, E0, E1 = condensate_scan(N, a, m, g, eps_values)
    _write_npz(
        out_path,
        eps=eps_values,
        gamma=gamma,
        E0=E0,
        E1=E1,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
    )
    i0 = int(np.argmin(np.abs(eps_values)))
    print(f"    saved {out_path}")
    print(f"    ground state energy       E0(eps={eps_values[i0]:g}) = {E0[i0]:.8f}")
    print(f"    first excited state energy E1(eps={eps_values[i0]:g}) = {E1[i0]:.8f}")
    print(f"    Gamma(0)={gamma[0]:.4f}  Gamma(eps_max)={gamma[-1]:.4f}")
    return eps_values, gamma, E0, E1


# ---- Section 2: Finite-size critical field scaling (paper Fig. 2) ----

def _refine_condensate_jump(N, a, m, g, lo, hi, gamma_lo, gamma_hi, tol=1e-4):
    """Bisect the Eq. (17) jump inside the coarse interval [lo, hi].

    A 30-point grid from 0 to 1.5 only names the interval. N=12, 14, and 16
    share one interval, so the interval midpoint assigns them the same eps_c.
    The crossing is the field where the ground-state condensate leaves gamma_lo
    and joins gamma_hi. That location is what Fig. 2 plots.
    """
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(mid))
        _evals, evecs = diagonalize(H, k=2)
        gamma_mid = chiral_condensate(evecs[:, 0], N, a=a)
        if abs(gamma_mid - gamma_lo) <= abs(gamma_mid - gamma_hi):
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def critical_field_scaling(N_values, a, m, g, eps_values, out_path):
    """Paper Fig. 2: eps_c(N) at the Eq. (17) condensate jump, refined inside the coarse interval."""
    eps_values = np.asarray(eps_values, dtype=float)
    N_values = np.asarray(N_values, dtype=int)
    eps_c = np.empty(N_values.size)
    gamma_scan = np.empty((N_values.size, eps_values.size))
    for i, N in enumerate(N_values):
        print(f"    N={int(N)}: scanning condensate")
        gamma, _E0, _E1 = condensate_scan(int(N), a, m, g, eps_values)
        gamma_scan[i] = gamma
        jump_idx = int(np.argmax(np.abs(np.diff(gamma))))
        lo = float(eps_values[jump_idx])
        hi = float(eps_values[jump_idx + 1])
        eps_c[i] = _refine_condensate_jump(
            int(N), a, m, g, lo, hi, float(gamma[jump_idx]), float(gamma[jump_idx + 1])
        )
        print(
            f"    N={int(N)}: eps_c={eps_c[i]:.4f}  "
            f"bracket=[{lo:.4f}, {hi:.4f}]  "
            f"|Delta Gamma|={abs(gamma[jump_idx + 1] - gamma[jump_idx]):.4f}"
        )
    _write_npz(
        out_path,
        N=N_values,
        eps_c=eps_c,
        eps_grid=eps_values,
        gamma=gamma_scan,
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
    )
    print(f"    saved {out_path}")
    return N_values, eps_c


def quench_evolved_states(N, a, m, g, eps_values, times, states_path):
    """Figs. 3-5: exact e^{-i H(eps) t}|psi_0> from the eps=0 vacuum. Reuses ``states_path`` when it matches."""
    eps_values = np.asarray(eps_values, dtype=float)
    times = np.asarray(times, dtype=float)
    states_path = Path(states_path)
    if states_path.is_file():
        with np.load(states_path) as data:
            match = (
                int(data["N"]) == int(N)
                and np.isclose(float(data["a"]), a)
                and np.isclose(float(data["m"]), m)
                and np.isclose(float(data["g"]), g)
                and data["eps"].shape == eps_values.shape
                and np.allclose(data["eps"], eps_values)
                and data["times"].shape == times.shape
                and np.allclose(data["times"], times)
            )
            if match:
                print(f"    loaded evolved states from {states_path} (evolve not rerun)")
                return np.array(data["psi0"]), np.array(data["psis"])

    psi0 = vacuum_state(N, a, m, g)
    psis = np.empty((eps_values.size, times.size, psi0.size), dtype=complex)
    print(f"    evolving {eps_values.size} quench fields x {times.size} times (N={N})")
    for k, eps in enumerate(eps_values):
        H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
        psis[k] = evolve(psi0, H, times)
    _write_npz(
        states_path,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
        eps=eps_values,
        times=times,
        psi0=psi0,
        psis=psis,
    )
    return psi0, psis


# ---- Section 3: Quench sanity check: conserved quantities (paper Fig. 3, repurposed) ----

def quench_conserved(N, a, m, g, eps_values, times, states_path, out_path):
    """Paper Fig. 3 sanity check: Q_N(t) from Eq. (12) and Delta E(t) = <H_eps>(t) - E(0)."""
    eps_values = np.asarray(eps_values, dtype=float)
    times = np.asarray(times, dtype=float)
    _psi0, psis = quench_evolved_states(N, a, m, g, eps_values, times, states_path)
    Q_N = np.empty((eps_values.size, times.size))
    delta_E = np.empty((eps_values.size, times.size))
    for k, eps in enumerate(eps_values):
        H_eps = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
        energies = np.empty(times.size)
        for i, psi in enumerate(psis[k]):
            Q_N[k, i] = total_charge(psi, N)
            energies[i] = np.real(psi.conj() @ (H_eps @ psi))
        delta_E[k] = energies - energies[0]
    q_max = float(np.max(np.abs(Q_N)))
    e_max = float(np.max(np.abs(delta_E)))
    print(f"    max |Q_N|={q_max:.3e}  max |Delta E|={e_max:.3e}")
    if q_max > 1e-8 or e_max > 1e-8:
        raise RuntimeError(
            "Conserved quantities drifted above 1e-8; exact evolution is not conserving Q_N or energy."
        )
    _write_npz(
        out_path,
        times=times,
        eps=eps_values,
        Q_N=Q_N,
        delta_E=delta_E,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
    )
    print(f"    saved {out_path}")
    return times, Q_N, delta_E


# ---- Section 4: Spatial and site-resolved charge dynamics (paper Figs. 4, 8) ----

def charge_dynamics(N, a, m, g, eps_values, times, states_path, out_path):
    """Paper Figs. 4 and 8: Q_i(t) from Eq. (14) and q_n(t) from Eq. (13) on the Fig. 3 trajectories."""
    eps_values = np.asarray(eps_values, dtype=float)
    times = np.asarray(times, dtype=float)
    _psi0, psis = quench_evolved_states(N, a, m, g, eps_values, times, states_path)
    n_eps, n_times, _dim = psis.shape
    Q_i = np.empty((n_eps, n_times, N // 2))
    q_n = np.empty((n_eps, n_times, N))
    for k, eps in enumerate(eps_values):
        print(f"    site and spatial charges for eps={eps:g}")
        for i, psi in enumerate(psis[k]):
            q_n[k, i] = local_charge(psi, N)
            Q_i[k, i] = spatial_point_charge(psi, N)
    _write_npz(
        out_path,
        times=times,
        eps=eps_values,
        Q_i=Q_i,
        q_n=q_n,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
    )
    print(f"    saved {out_path}  Q_i {Q_i.shape}  q_n {q_n.shape}")
    return Q_i, q_n


# ---- Section 5: Electric field energy dynamics (paper Fig. 5) ----

def field_energy_dynamics(N, a, m, g, eps_values, times, states_path, out_path):
    """Paper Fig. 5 / Eq. (15): electric-field energy H_E(t) on the Fig. 3 trajectories."""
    eps_values = np.asarray(eps_values, dtype=float)
    times = np.asarray(times, dtype=float)
    _psi0, psis = quench_evolved_states(N, a, m, g, eps_values, times, states_path)
    H_E = np.empty((eps_values.size, times.size))
    for k, eps in enumerate(eps_values):
        print(f"    H_E(t) for eps={eps:g}")
        for i, psi in enumerate(psis[k]):
            H_E[k, i] = field_energy(psi, N, a=a, g=g, eps=float(eps))
    _write_npz(
        out_path,
        times=times,
        eps=eps_values,
        H_E=H_E,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
    )
    print(f"    saved {out_path}")
    return H_E


# ---- Section 6: Vacuum-state fidelity dynamics (paper Fig. 6) ----

def vacuum_fidelity_dynamics(N, a, m, g, eps_values, times, out_path):
    """Paper Fig. 6 / Eq. (16): P_vac(t) = |<psi_0|psi(t)>|^2 for each quench field."""
    eps_values = np.asarray(eps_values, dtype=float)
    times = np.asarray(times, dtype=float)
    psi0 = vacuum_state(N, a, m, g)
    P_vac = np.empty((eps_values.size, times.size))
    for k, eps in enumerate(eps_values):
        H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
        psis = evolve(psi0, H, times)
        P_vac[k] = vacuum_fidelity(psis, psi0)
        print(f"    eps={eps:g}: P_vac(0)={P_vac[k, 0]:.6f}  min P_vac={P_vac[k].min():.4f}")
    _write_npz(
        out_path,
        times=times,
        eps=eps_values,
        P_vac=P_vac,
        N=np.int64(N),
        a=np.float64(a),
        m=np.float64(m),
        g=np.float64(g),
    )
    print(f"    saved {out_path}")
    return P_vac


# ---- Section 7: Effective decay rate vs field strength (paper Fig. 7) ----

def effective_decay_rate(fid_path, a, out_path, t_over_a_min=0.0, t_over_a_max=1.0):
    """Paper Fig. 7: gamma_eff from -ln P_vac(t) = gamma_eff * t + c on the Eq. (16) early-time window."""
    with np.load(fid_path) as data:
        times = np.array(data["times"], dtype=float)
        eps_values = np.array(data["eps"], dtype=float)
        P_vac = np.array(data["P_vac"], dtype=float)
        N = int(data["N"])
        stored_a = float(data["a"])
    if not np.isclose(stored_a, a):
        raise ValueError(f"Requested a={a} does not match a={stored_a} stored in {fid_path}")
    mask = (times / a >= t_over_a_min) & (times / a <= t_over_a_max)
    t_fit = times[mask]
    if t_fit.size < 2:
        raise ValueError("Early-time window contains fewer than two samples.")
    gamma_eff = np.empty(eps_values.size)
    for k, eps in enumerate(eps_values):
        pvac = P_vac[k, mask]
        if np.any(pvac <= 0.0):
            raise ValueError(f"P_vac is non-positive for eps={eps} inside the fit window.")
        slope, _intercept = np.polyfit(t_fit, -np.log(pvac), 1)
        gamma_eff[k] = slope
        print(f"    eps={eps:g}: gamma_eff={slope:.6f}")
    _write_npz(
        out_path,
        eps=eps_values,
        gamma_eff=gamma_eff,
        N=np.int64(N),
        a=np.float64(a),
        t_over_a_min=np.float64(t_over_a_min),
        t_over_a_max=np.float64(t_over_a_max),
    )
    print(f"    saved {out_path}")
    return gamma_eff


if __name__ == "__main__":
    print(
        "NOTE: Table I and all QC/VQE/VQD/Trotter markers from the paper "
        "are not reproduced; this pipeline is exact diagonalization only."
    )
    N = 8
    a = 1.0
    m = 1.0
    g = 1.0
    data = ROOT / "data" / "ED"
    eps_quench = (0.5, 1.0, 1.5, 2.0)
    eps_fidelity = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
    times = np.linspace(0.0, 12.0, 121)
    states_path = data / "quench_states_n8.npz"

    print("[1/7] Chiral condensate & spectrum vs external field (N=8)...")
    condensate_and_spectrum(
        N=N,
        a=a,
        m=m,
        g=g,
        eps_values=np.linspace(0.0, 3.0, 81),
        out_path=data / "condensate_spectrum_n8.npz",
    )

    print("[2/7] Finite-size critical field scaling (N=8,10,12...)...")
    critical_field_scaling(
        N_values=(8, 10, 12, 14, 16),
        a=a,
        m=m,
        g=g,
        eps_values=np.linspace(0.0, 1.5, 30),
        out_path=data / "critical_field_scaling.npz",
    )

    print("[3/7] Quench sanity check (N=8)...")
    quench_conserved(
        N=N,
        a=a,
        m=m,
        g=g,
        eps_values=eps_quench,
        times=times,
        states_path=states_path,
        out_path=data / "quench_conserved_n8.npz",
    )

    print("[4/7] Spatial and site-resolved charge dynamics (N=8)...")
    charge_dynamics(
        N=N,
        a=a,
        m=m,
        g=g,
        eps_values=eps_quench,
        times=times,
        states_path=states_path,
        out_path=data / "charge_dynamics_n8.npz",
    )

    print("[5/7] Electric field energy dynamics (N=8)...")
    field_energy_dynamics(
        N=N,
        a=a,
        m=m,
        g=g,
        eps_values=eps_quench,
        times=times,
        states_path=states_path,
        out_path=data / "field_energy_n8.npz",
    )

    print("[6/7] Vacuum-state fidelity dynamics (N=8)...")
    vacuum_fidelity_dynamics(
        N=N,
        a=a,
        m=m,
        g=g,
        eps_values=eps_fidelity,
        times=times,
        out_path=data / "vacuum_fidelity_n8.npz",
    )

    print("[7/7] Effective decay rate vs field strength (N=8)...")
    effective_decay_rate(
        fid_path=data / "vacuum_fidelity_n8.npz",
        a=a,
        out_path=data / "decay_rate_n8.npz",
    )
