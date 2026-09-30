"""VQD first excited state of the lattice Schwinger model, Section III.B.

Eq. (19) is minimized with a RealAmplitudes trial state and SLSQP, the
same optimizer used for Eq. (18). The orthogonality reference ψ0^VQE is an
argument; this
module does not load it or run VQE. After the minimization, Eq. (20) is the
number compared with evals[1] from diagonalize(H, k=2), the E1(ε) curve in
Fig. 1(b). scripts/run_vqd.py loads the reference and writes the dataset.
"""

import os
import warnings
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"The class ``qiskit\.circuit\.library\.n_local\.real_amplitudes\.RealAmplitudes``.*",
)

import numpy as np
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.optimizers import SLSQP
from qiskit_algorithms.utils import algorithm_globals

from schwinger_model import build_hamiltonian, diagonalize


def trial_ansatz(num_qubits):
    """RealAmplitudes trial state ψ(θ) in Eq. (19).

    excited_scan minimizes Eq. (19) on this circuit. The vector it returns
    is the state whose overlap with evecs[:, 1] from diagonalize(H, k=2)
    is the Fig. 1 first-excited-state fidelity.
    """
    # Nine repetitions still leave the lowest Eq. (19) cost below fidelity
    # 0.99 near ε = 0.5, where the first excited state rearranges inside
    # the Q = −1 sector. Ten is where that minimum is evecs[:, 1].
    return RealAmplitudes(num_qubits, reps=10)


def _trial_amplitudes(theta, num_qubits):
    """Amplitudes of the RealAmplitudes state ψ(θ) in Eq. (19).

    trial_statevector returns this vector, and _penalized_cost_and_gradient
    differentiates the same amplitudes. RY layers with reverse-linear CX
    are the circuit from trial_ansatz.
    """
    theta = np.asarray(theta, dtype=float)
    n_layers = theta.size // num_qubits
    reps = n_layers - 1
    psi = np.zeros(1 << num_qubits, dtype=float)
    psi[0] = 1.0
    for layer in range(n_layers):
        for qubit in range(num_qubits):
            angle = 0.5 * theta[layer * num_qubits + qubit]
            cosine = np.cos(angle)
            sine = np.sin(angle)
            view = psi.reshape(-1, 2, 1 << qubit)
            amp0 = view[:, 0, :].copy()
            amp1 = view[:, 1, :].copy()
            view[:, 0, :] = cosine * amp0 - sine * amp1
            view[:, 1, :] = sine * amp0 + cosine * amp1
        if layer == reps:
            break
        for qubit in range(num_qubits - 2, -1, -1):
            view = psi.reshape(-1, 2, 2, 1 << qubit)
            held = view[:, 0, 1, :].copy()
            view[:, 0, 1, :] = view[:, 1, 1, :]
            view[:, 1, 1, :] = held
    return psi


def trial_statevector(ansatz, theta):
    """Amplitudes of the Eq. (19) trial state ψ(θ).

    This is the vector whose energy, after the penalty is dropped, is Eq. (20),
    and whose overlap with evecs[:, 1] from diagonalize(H, k=2) says whether
    that energy belongs to the ED first excited state or to a different vector
    at a similar energy.
    """
    return _trial_amplitudes(theta, ansatz.num_qubits)


def _penalized_cost_and_gradient(theta, num_qubits, hamiltonian, psi_reference, beta):
    """Eq. (19) and ∂C1/∂θ.

    minimize_excited gives both to SLSQP. The value is the penalized cost.
    Eq. (20) is the energy of the minimizer after this penalty is dropped.
    """
    theta = np.asarray(theta, dtype=float)
    psi = _trial_amplitudes(theta, num_qubits)
    final = psi.copy()
    h_psi = hamiltonian @ psi
    amplitude = np.vdot(psi_reference, psi)
    energy = float(psi @ h_psi)
    cost = energy + beta * float(np.abs(amplitude) ** 2)
    cotangent = 2.0 * h_psi + beta * 2.0 * np.real(psi_reference * amplitude)
    n_layers = theta.size // num_qubits
    reps = n_layers - 1
    gradient = np.empty(theta.size, dtype=float)
    index = theta.size - 1
    for layer in range(reps, -1, -1):
        for qubit in range(num_qubits - 1, -1, -1):
            angle = 0.5 * theta[index]
            cosine = np.cos(angle)
            sine = np.sin(angle)
            state = psi.reshape(-1, 2, 1 << qubit)
            derivative = cotangent.reshape(-1, 2, 1 << qubit)
            amp0 = state[:, 0, :].copy()
            amp1 = state[:, 1, :].copy()
            grad0 = derivative[:, 0, :].copy()
            grad1 = derivative[:, 1, :].copy()
            gradient[index] = 0.5 * np.sum(grad1 * amp0 - grad0 * amp1)
            derivative[:, 0, :] = cosine * grad0 + sine * grad1
            derivative[:, 1, :] = -sine * grad0 + cosine * grad1
            state[:, 0, :] = cosine * amp0 + sine * amp1
            state[:, 1, :] = -sine * amp0 + cosine * amp1
            index -= 1
        if layer == 0:
            break
        for qubit in range(num_qubits - 1):
            state = psi.reshape(-1, 2, 2, 1 << qubit)
            derivative = cotangent.reshape(-1, 2, 2, 1 << qubit)
            held = state[:, 0, 1, :].copy()
            state[:, 0, 1, :] = state[:, 1, 1, :]
            state[:, 1, 1, :] = held
            held = derivative[:, 0, 1, :].copy()
            derivative[:, 0, 1, :] = derivative[:, 1, 1, :]
            derivative[:, 1, 1, :] = held
    return final, cost, gradient


def penalty_coefficient(H):
    """β in Eq. (19).

    Section III.B leaves the value unspecified. qiskit_algorithms.VQD, when
    betas is not passed, sets β = 10 * Σ|c_i|, where c_i are the Pauli
    coefficients of the operator (see VQD.compute_eigenvalues). The operator
    here is SparsePauliOp.from_operator on the build_hamiltonian() matrix, so
    the coefficients belong to the same H that diagonalize() factorizes.
    """
    op = SparsePauliOp.from_operator(H.toarray())
    return float(10.0 * np.sum(np.abs(np.real(op.coeffs))))


def penalized_cost(ansatz, theta, H, psi_reference, beta):
    """Eq. (19), C1(θ) = ⟨ψ(θ)|H|ψ(θ)⟩ + β |⟨ψ0^VQE|ψ(θ)⟩|^2.

    H is build_hamiltonian() at this ε, and psi_reference is the VQE ground
    state at the same ε. The second term is what pushes the minimum off
    evals[0] from diagonalize(H, k=2) and toward evals[1]. The number
    returned here is the objective, not the excited-state energy; Eq. (20)
    drops the penalty after the minimization.
    """
    psi = trial_statevector(ansatz, theta)
    energy = float(np.real(np.vdot(psi, H @ psi)))
    overlap = float(np.abs(np.vdot(psi_reference, psi)) ** 2)
    return energy + beta * overlap


def excited_energy(psi, H):
    """Eq. (20), E1^VQD = ⟨ψ1^VQD|H|ψ1^VQD⟩.

    H is build_hamiltonian() at this ε. This is the expectation with the
    penalty removed, so it is the number placed against evals[1] from
    diagonalize(H, k=2). It is not bounded below by evals[1] as strictly as
    Eq. (18) is by evals[0], because the state was constrained to be
    orthogonal to the VQE ground state rather than to evecs[:, 0].
    """
    return float(np.real(np.vdot(psi, H @ psi)))


def minimize_excited(ansatz, H, psi_reference, beta, theta0):
    """One SLSQP minimization of Eq. (19) from one initial parameter vector.

    Same optimizer as the ground-state search: qiskit_algorithms.optimizers.SLSQP
    with maxiter=1000 (Qiskit Algorithms tutorial "An Introduction to
    Algorithms using Qiskit") and the class-default ftol. SLSQP receives
    ∂C1/∂θ from _penalized_cost_and_gradient. The state that is kept is the
    minimizer of Eq. (19). Its Eq. (20) energy is what gets compared with
    diagonalize(H, k=2) evals[1].
    """
    hamiltonian = np.real(np.asarray(H.toarray()))
    reference = np.asarray(psi_reference).reshape(-1)
    num_qubits = ansatz.num_qubits
    cache = {}

    def _evaluated(theta):
        key = np.asarray(theta, dtype=float).tobytes()
        if cache.get("key") == key:
            return cache["cost"], cache["gradient"]
        _psi, cost, gradient = _penalized_cost_and_gradient(
            theta, num_qubits, hamiltonian, reference, beta
        )
        cache["key"] = key
        cache["cost"] = cost
        cache["gradient"] = gradient
        return cost, gradient

    def objective(theta):
        return _evaluated(theta)[0]

    def jacobian(theta):
        return _evaluated(theta)[1]

    result = SLSQP(maxiter=1000).minimize(
        objective,
        np.asarray(theta0, dtype=float),
        jac=jacobian,
        bounds=ansatz.parameter_bounds,
    )
    theta = np.asarray(result.x, dtype=float)
    psi = trial_statevector(ansatz, theta)
    return {
        "theta": theta,
        "psi": psi,
        "cost": float(result.fun),
        "energy": excited_energy(psi, H),
        "overlap": float(np.abs(np.vdot(psi_reference, psi)) ** 2),
        "nfev": int(result.nfev),
    }


def random_parameters(ansatz, n_restarts, seed):
    """Initial parameters at the first field point.

    There is no previous field to warm-start from. The draw is the same one
    qiskit_algorithms.VQE uses for Eq. (18): uniform on
    RealAmplitudes.parameter_bounds, (−π, π). Each restart minimizes Eq. (19);
    the caller keeps the lowest penalized cost, then reads Eq. (20).
    """
    algorithm_globals.random_seed = seed
    low = np.array([bound[0] for bound in ansatz.parameter_bounds], dtype=float)
    high = np.array([bound[1] for bound in ansatz.parameter_bounds], dtype=float)
    return np.vstack(
        [algorithm_globals.random.uniform(low, high) for _ in range(n_restarts)]
    )


def ed_lowest_two(N, a, m, g, eps):
    """evals[0], evals[1], and evecs[:, 1] of build_hamiltonian() at this field.

    diagonalize(H, k=2) is the ED side of Fig. 1(b). evals[1] is the energy
    Eq. (20) is compared with. evecs[:, 1] is the vector the VQD fidelity is
    |⟨·|ψ1^VQD⟩|^2 against. evals[0] is returned so a VQD energy that fell
    back onto the ground state is visible next to the loaded VQE energy.
    """
    H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
    evals, evecs = diagonalize(H, k=2)
    return H, float(evals[0]), float(evals[1]), np.asarray(evecs[:, 1]).reshape(-1)


def _minimized_excited(job):
    """One SLSQP minimum of Eq. (19).

    excited_scan consumes the penalized cost and keeps the lowest one.
    Eq. (20) is the energy on that minimizer, with the penalty removed.
    """
    ansatz, H, psi_reference, beta, theta0 = job
    return minimize_excited(ansatz, H, psi_reference, beta, theta0)


def excited_scan(N, a, m, g, eps_values, psi_vqe, n_restarts, seed):
    """E1^VQD(ε) on a grid whose ground states are already known.

    psi_vqe[i] is ψ0^VQE at eps_values[i], the reference in Eq. (19), and
    it has to already be the ground state. At each field Eq. (19) is
    minimized from the previous parameters and from n_restarts random
    vectors, and the lowest penalized cost is kept. Further draws are
    tried while that state is not yet evecs[:, 1]. The first field has
    no previous parameters. The reported energy is Eq. (20), compared
    with evals[1] from diagonalize(H, k=2). run_vqd.py writes the scan.
    """
    eps_values = np.asarray(eps_values, dtype=float)
    psi_vqe = np.asarray(psi_vqe)
    n_eps = len(eps_values)
    dim = psi_vqe.shape[1]
    ansatz = trial_ansatz(N)
    energies = np.empty(n_eps)
    energies_ed = np.empty(n_eps)
    ground_ed = np.empty(n_eps)
    fidelities = np.empty(n_eps)
    overlaps = np.empty(n_eps)
    betas = np.empty(n_eps)
    thetas = np.empty((n_eps, ansatz.num_parameters))
    psis = np.empty((n_eps, dim), dtype=complex)
    theta = None
    for i, eps in enumerate(eps_values):
        H, e0, e1, psi_e1 = ed_lowest_two(N, a, m, g, eps)
        beta = penalty_coefficient(H)
        psi_ref = np.asarray(psi_vqe[i]).reshape(-1)
        # A warm start can stay orthogonal to the reference and still miss
        # the first excited state, so the lower penalized cost is kept.
        chosen = None
        fidelity = 0.0
        batch = 0
        n_drawn = 0
        while chosen is None or (fidelity < 0.99 and n_drawn < 48):
            if batch == 0 and i == 0:
                starts = random_parameters(ansatz, n_restarts, seed)
            elif batch == 0:
                draws = random_parameters(ansatz, n_restarts, seed + i)
                starts = np.vstack([theta, draws])
            else:
                starts = random_parameters(
                    ansatz, n_restarts, seed + 10007 * batch + i
                )
            batch += 1
            n_drawn += n_restarts
            jobs = [
                (ansatz, H, psi_ref, beta, np.asarray(th, dtype=float))
                for th in starts
            ]
            if len(jobs) == 1:
                candidates = [_minimized_excited(jobs[0])]
            else:
                workers = min(len(jobs), os.cpu_count() or 1)
                context = mp.get_context("fork")
                with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
                    candidates = list(pool.map(_minimized_excited, jobs))
            challenger = min(candidates, key=lambda candidate: candidate["cost"])
            if chosen is None or challenger["cost"] < chosen["cost"]:
                chosen = challenger
            fidelity = float(np.abs(np.vdot(psi_e1, chosen["psi"])) ** 2)
        theta = chosen["theta"]
        energies[i] = chosen["energy"]
        energies_ed[i] = e1
        ground_ed[i] = e0
        fidelities[i] = fidelity
        overlaps[i] = chosen["overlap"]
        betas[i] = beta
        thetas[i] = chosen["theta"]
        psis[i] = chosen["psi"]
        print(
            f"    eps={eps:.4f}: E1_VQD={chosen['energy']:.6f}  E1_ED={e1:.6f}  "
            f"F={fidelity:.4f}  overlap={chosen['overlap']:.3e}",
            flush=True,
        )
    return energies, energies_ed, ground_ed, fidelities, overlaps, betas, thetas, psis
