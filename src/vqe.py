"""VQE for the lattice Schwinger model, Section III.B of arXiv:2607.02894.

Eq. (18) is minimized with the RealAmplitudes trial state and SLSQP.
At ε = 0 several random initial parameters are tried and the lowest energy
is kept (Table I). The field scan then walks from ε = 0 to ε = 3. Each
later field starts from the previous parameters and from random restarts,
and the lowest energy is kept. Those states are the ψ0^VQE that Eq. (19)
penalizes against. Nothing here writes a file; scripts/run_vqe.py does that.
"""

import os
import warnings
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# Section III.B names RealAmplitudes. Qiskit 2.1 deprecated the class in favor
# of real_amplitudes(); the class is still that circuit, including its bounds.
warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"The class ``qiskit\.circuit\.library\.n_local\.real_amplitudes\.RealAmplitudes``.*",
)

import numpy as np
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import Statevector
from qiskit_algorithms.optimizers import SLSQP
from qiskit_algorithms.utils import algorithm_globals

from schwinger_model import build_hamiltonian, diagonalize


def ed_ground_pair(N, a, m, g, eps):
    """Ground energy and ground vector of build_hamiltonian() at one field.

    Returns evals[0] and evecs[:, 0] from diagonalize(H, k=1). At N = 8,
    a = m = g = 1, ε = 0, that eigenvalue is Table I's ED energy,
    −4.63805774. On the field scan it is the E0(ε) curve in Fig. 1(b).
    The vector is the state Table I's fidelity overlaps with.
    """
    H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
    evals, evecs = diagonalize(H, k=1)
    return H, float(evals[0]), np.asarray(evecs[:, 0]).reshape(-1)


def trial_ansatz(num_qubits, reps=None):
    """RealAmplitudes trial state ψ(θ) in Eq. (18).

    Section III.B names this ansatz and does not set repetitions or the
    entanglement pattern. With reps omitted this is the Qiskit class default:
    reps=3, entanglement='reverse_linear', parameter bounds (−π, π).
    qiskit_algorithms.VQE reads those bounds in validate_bounds and
    validate_initial_point. The statevector is in the same basis as evecs
    from diagonalize(), so Table I's overlap needs no reordering.
    """
    if reps is None:
        return RealAmplitudes(num_qubits)
    return RealAmplitudes(num_qubits, reps=reps)


def trial_statevector(ansatz, theta):
    """Amplitudes of ψ(θ) in Eq. (18).

    Overlap this vector with evecs[:, 0] from diagonalize() to get Table I's
    fidelity. RealAmplitudes amplitudes are real, and diagonalize() returns a
    real eigenvector of the same matrix, so the overlap is real up to a sign.
    """
    return np.asarray(Statevector(ansatz.assign_parameters(theta)).data).reshape(-1)


def energy_expectation(ansatz, theta, H):
    """Eq. (18), E(θ) = ⟨ψ(θ)|H|ψ(θ)⟩.

    H is the matrix from build_hamiltonian(), the operator diagonalize()
    factorizes. Variationally this lies at or above evals[0]. At the Table I
    point that floor is −4.63805774, and the published VQE energy
    −4.63766032 lies just above it.
    """
    psi = trial_statevector(ansatz, theta)
    return float(np.real(np.vdot(psi, H @ psi)))


def state_fidelity(psi, psi_ed):
    """Table I fidelity F = |⟨ψ0^ED|ψ0^VQE⟩|^2.

    psi_ed is evecs[:, 0] from diagonalize(build_hamiltonian(...), k=1).
    F = 1 means the VQE vector is that ED ground state up to a global phase.
    Table I quotes 99.9931% at ε = 0. Anything smaller is the wave-function
    error that goes with E(θ) sitting above evals[0].
    """
    return float(np.abs(np.vdot(psi_ed, psi)) ** 2)


def relative_energy_error(energy, energy_ed):
    """Table I relative error |E_VQE − E_ED| / |E_ED|.

    E_ED is evals[0] from diagonalize(). Table I quotes 8.569×10^−5 for this
    ratio against the ED energy −4.63805774. Eq. (18) is an upper bound, so
    the numerator is E_VQE − E_ED when the minimization has stayed above the
    exact ground energy.
    """
    return float(abs(energy - energy_ed) / abs(energy_ed))


def minimize_energy(ansatz, H, theta0, maxiter=1000, ftol=1e-6):
    """One SLSQP minimization of Eq. (18) from one initial parameter vector.

    Optimizer: qiskit_algorithms.optimizers.SLSQP with maxiter=1000, as in
    the Qiskit Algorithms tutorial "An Introduction to Algorithms using
    Qiskit". ftol stays at the class default, 1e-6, unless overridden.
    Bounds are RealAmplitudes.parameter_bounds. The returned energy is still
    only an upper bound on diagonalize(H) evals[0]; Section III.B repeats the
    minimization from new initial parameters because a single start can stop
    in a local minimum.
    """
    def objective(theta):
        return energy_expectation(ansatz, theta, H)

    return SLSQP(maxiter=maxiter, ftol=ftol).minimize(
        objective,
        np.asarray(theta0, dtype=float),
        bounds=ansatz.parameter_bounds,
    )


def random_parameters(ansatz, n_restarts, seed):
    """Initial parameters for the random restarts in Section III.B.

    The paper says multiple random restarts and does not say how they are
    drawn. qiskit_algorithms.VQE, through validate_initial_point, draws each
    parameter uniformly from the ansatz bounds when no initial point is
    given. RealAmplitudes sets those bounds to (−π, π). Each draw is one
    attempt to land in the basin closest to diagonalize() evals[0].
    """
    algorithm_globals.random_seed = seed
    low = np.array([bound[0] for bound in ansatz.parameter_bounds], dtype=float)
    high = np.array([bound[1] for bound in ansatz.parameter_bounds], dtype=float)
    return np.vstack(
        [algorithm_globals.random.uniform(low, high) for _ in range(n_restarts)]
    )


def _minimized_energy(job):
    """Eq. (18) energy and ψ(θ) from one SLSQP start.

    best_restart consumes the energy and keeps the lowest one. The state is
    the vector state_fidelity compares with evecs[:, 0].
    """
    ansatz, H, theta0, maxiter, ftol = job
    result = minimize_energy(ansatz, H, theta0, maxiter=maxiter, ftol=ftol)
    theta = np.asarray(result.x, dtype=float)
    psi = trial_statevector(ansatz, theta)
    energy = float(np.real(np.vdot(psi, H @ psi)))
    return energy, theta, psi, int(result.nfev)


def best_restart(ansatz, H, psi_ed, initial_points, maxiter=1000, ftol=1e-6):
    """Lowest Eq. (18) energy among the supplied SLSQP starts.

    Table I passes random initial parameters. field_scan passes those
    together with the warm start from the previous field. The kept vector
    is the upper bound on evals[0] that run_vqe.py stores. Its fidelity is
    state_fidelity against psi_ed = evecs[:, 0]. Per-start energies are
    returned so a discarded minimum can be compared with that eigenvalue.
    """
    jobs = [
        (ansatz, H, np.asarray(theta0, dtype=float), int(maxiter), float(ftol))
        for theta0 in initial_points
    ]
    if len(jobs) == 1:
        outcomes = [_minimized_energy(jobs[0])]
    else:
        workers = min(len(jobs), os.cpu_count() or 1)
        context = mp.get_context("fork")
        with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
            outcomes = list(pool.map(_minimized_energy, jobs))
    restart_energies = np.empty(len(outcomes))
    best = None
    for i, (energy, theta, psi, nfev) in enumerate(outcomes):
        fidelity = state_fidelity(psi, psi_ed)
        restart_energies[i] = energy
        print(
            f"    restart {i}: E={energy:.8f}  F={fidelity:.6f}  nfev={nfev}",
            flush=True,
        )
        if best is None or energy < best["energy"]:
            best = {"energy": energy, "theta": theta, "psi": psi, "fidelity": fidelity}
    return best, restart_energies


def field_scan(
    ansatz,
    theta_at_zero,
    eps_values,
    N,
    a,
    m,
    g,
    n_restarts,
    seed,
    maxiter=1000,
    ftol=1e-6,
):
    """Ground energy E0^VQE(ε) and state ψ0^VQE on the Fig. 1 grid.

    Eq. (18) at each stored ε. ε = 0 is the Table I vector and is not
    optimized again. At each later field the previous parameters and
    n_restarts random vectors are minimized, and the lowest energy is kept.
    The batch size is fixed. Fidelity to the exact ground state is recorded
    afterwards and does not open another batch. run_vqe.py writes these
    states. run_vqd.py loads them as the orthogonality reference in Eq. (19).
    """
    n_eps = len(eps_values)
    dim = 2 ** N
    energies = np.empty(n_eps)
    energies_ed = np.empty(n_eps)
    fidelities = np.empty(n_eps)
    thetas = np.empty((n_eps, ansatz.num_parameters))
    psis = np.empty((n_eps, dim), dtype=complex)
    theta = np.asarray(theta_at_zero, dtype=float)
    for i, eps in enumerate(eps_values):
        H, energy_ed, psi_ed = ed_ground_pair(N, a, m, g, eps)
        if i == 0:
            psi = trial_statevector(ansatz, theta)
            energy = float(np.real(np.vdot(psi, H @ psi)))
        else:
            # The warm start stays in the old vacuum across a first-order jump.
            # The lowest energy among that start and n_restarts random starts
            # is kept. Fidelity is not used to draw further starts.
            draws = random_parameters(ansatz, n_restarts, seed + i)
            initial_points = np.vstack([theta, draws])
            best, _restart_energies = best_restart(
                ansatz, H, psi_ed, initial_points, maxiter=maxiter, ftol=ftol
            )
            theta = best["theta"]
            psi = best["psi"]
            energy = best["energy"]
        energies[i] = energy
        energies_ed[i] = energy_ed
        fidelities[i] = state_fidelity(psi, psi_ed)
        thetas[i] = theta
        psis[i] = psi
        print(
            f"    eps={eps:.4f}: E_VQE={energy:.6f}  E0_ED={energy_ed:.6f}  "
            f"F={fidelities[i]:.4f}",
            flush=True,
        )
    return energies, energies_ed, fidelities, thetas, psis
