"""VQD first excited state of the lattice Schwinger model, Section III.B.

Eq. (19) is minimized with the same RealAmplitudes trial state and SLSQP
used for Eq. (18). The orthogonality reference ψ0^VQE is an argument; this
module does not load it or run VQE. After the minimization, Eq. (20) is the
number compared with evals[1] from diagonalize(H, k=2), the E1(ε) curve in
Fig. 1(b). scripts/run_vqd.py loads the reference and writes the dataset.
"""

import warnings

warnings.filterwarnings(
    "ignore",
    category=DeprecationWarning,
    message=r"The class ``qiskit\.circuit\.library\.n_local\.real_amplitudes\.RealAmplitudes``.*",
)

import numpy as np
from qiskit.circuit.library import RealAmplitudes
from qiskit.quantum_info import SparsePauliOp, Statevector
from qiskit_algorithms.optimizers import SLSQP
from qiskit_algorithms.utils import algorithm_globals

from schwinger_model import build_hamiltonian, diagonalize


def trial_ansatz(num_qubits, reps):
    """RealAmplitudes trial state ψ(θ) in Eq. (19).

    Same circuit Section III.B uses for Eq. (18). reps is the repetition
    count of the VQE ground states passed in (the Qiskit RealAmplitudes
    default unless that run used another). The statevector lines up with
    evecs from diagonalize(), so the fidelity against evecs[:, 1] is a
    direct overlap.
    """
    return RealAmplitudes(num_qubits, reps=reps)


def trial_statevector(ansatz, theta):
    """Amplitudes of the Eq. (19) trial state ψ(θ).

    This is the vector whose energy, after the penalty is dropped, is Eq. (20),
    and whose overlap with evecs[:, 1] from diagonalize(H, k=2) says whether
    that energy belongs to the ED first excited state or to a different vector
    at a similar energy.
    """
    return np.asarray(Statevector(ansatz.assign_parameters(theta)).data).reshape(-1)


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
    Algorithms using Qiskit") and the class-default ftol. The state that is
    kept is the minimizer of Eq. (19). Its Eq. (20) energy is what gets
    compared with diagonalize(H, k=2) evals[1].
    """
    def objective(theta):
        return penalized_cost(ansatz, theta, H, psi_reference, beta)

    result = SLSQP(maxiter=1000).minimize(
        objective,
        np.asarray(theta0, dtype=float),
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


def excited_scan(N, a, m, g, reps, eps_values, psi_vqe, n_restarts, seed):
    """E1^VQD(ε) on a grid whose ground states are already known.

    psi_vqe[i] is ψ0^VQE at eps_values[i], the reference in Eq. (19). At the
    first field, Eq. (19) is minimized from n_restarts random parameter
    vectors and the lowest penalized cost is kept. After that, the paragraph
    following Eq. (20) applies: the parameters at one field initialize the
    next. The reported energy is Eq. (20), compared with evals[1] from
    diagonalize(H, k=2) at that same field.
    """
    eps_values = np.asarray(eps_values, dtype=float)
    psi_vqe = np.asarray(psi_vqe)
    n_eps = len(eps_values)
    dim = psi_vqe.shape[1]
    ansatz = trial_ansatz(N, reps)
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
        if i == 0:
            draws = random_parameters(ansatz, n_restarts, seed)
            chosen = None
            for draw in draws:
                candidate = minimize_excited(ansatz, H, psi_ref, beta, draw)
                if chosen is None or candidate["cost"] < chosen["cost"]:
                    chosen = candidate
        else:
            chosen = minimize_excited(ansatz, H, psi_ref, beta, theta)
        theta = chosen["theta"]
        fidelity = float(np.abs(np.vdot(psi_e1, chosen["psi"])) ** 2)
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
