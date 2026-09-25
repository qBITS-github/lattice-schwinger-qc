"""Second-order Trotter quench of the lattice Schwinger model, Section III.C.

The initial state is the ε = 0 ground state from diagonalize(), evolved under
H(ε) by the Suzuki-Trotter step in Eq. (22) with a fixed Δt. The same initial
state and the same times are passed to evolve(), which is the exact propagator
Eq. (21). scripts/run_trotter.py chooses the fields and writes both trajectories.
"""

import numpy as np
from qiskit.circuit.library import PauliEvolutionGate
from qiskit.quantum_info import Operator, SparsePauliOp
from qiskit.synthesis import SuzukiTrotter

from schwinger_model import build_hamiltonian, diagonalize, evolve


def zero_field_ground_state(N, a, m, g):
    """ε = 0 ground state used as |ψ0⟩ in Eq. (21).

    evecs[:, 0] from diagonalize(build_hamiltonian(eps=0), k=1). This is the
    same vector scripts/run_exact_diagonalization.py feeds to evolve() for Figs. 3–6. It is the ED
    ground state, not the VQE state: the comparison below changes only the
    propagator, from evolve()'s exact exponential to Eq. (22).
    """
    H0 = build_hamiltonian(N, a=a, m=m, g=g, eps=0.0)
    _evals, evecs = diagonalize(H0, k=1)
    return np.asarray(evecs[:, 0]).reshape(-1)


def pauli_hamiltonian(N, a, m, g, eps):
    """Pauli operator passed to the Suzuki-Trotter synthesizer in Section III.C.

    SparsePauliOp.from_operator of the build_hamiltonian() matrix, so the
    operator is that matrix and not a second encoding of Eqs. (9)–(11).
    evolve() diagonalizes the same matrix. SuzukiTrotter walks the Pauli
    terms in the order from_operator returns (preserve_order=True, the class
    default).
    """
    H = build_hamiltonian(N, a=a, m=m, g=g, eps=float(eps))
    op = SparsePauliOp.from_operator(H.toarray())
    if not np.allclose(op.to_matrix(), H.toarray()):
        raise RuntimeError("Pauli operator does not match build_hamiltonian().")
    return H, op


def one_step_unitary(op, dt):
    """Matrix of one second-order step, Eq. (22) with r = 1.

    Section III.C: the full Pauli Hamiltonian is given to Suzuki-Trotter
    (order=2, reps=1) and the step is Δt. That is
    qiskit.synthesis.SuzukiTrotter(order=2, reps=1) synthesizing
    qiskit.circuit.library.PauliEvolutionGate, as in the Qiskit Algorithms
    tutorial "Quantum Real Time Evolution using Trotterization". One
    application approximates exp(−i H Δt). evolve() applies the exact
    exponential of the same H. They agree at t = 0 and separate by the
    O(Δt^2) global error of Eq. (22).
    """
    # SuzukiTrotter(order=2, reps=1) is the paper's synthesizer and the class
    # default; order is passed explicitly because Section III.C names it.
    synthesis = SuzukiTrotter(order=2, reps=1)
    gate = PauliEvolutionGate(op, time=dt, synthesis=synthesis)
    # Operator() of the gate itself is the exact exponential. synthesize()
    # is the product of one-Pauli rotations in Eq. (22).
    circuit = synthesis.synthesize(gate)
    return np.asarray(Operator(circuit).data)


def trotter_trajectory(psi0, unitary, n_steps):
    """States at t = 0, Δt, ..., n_steps·Δt by repeating the Eq. (22) step.

    Index k is the state evolve() returns at times[k] for the same psi0 and
    the same H. The k = 0 vector is psi0 itself, so the fidelity against
    evolve() starts at 1 and then measures the accumulated Trotter error.
    """
    psi = np.asarray(psi0, dtype=complex).reshape(-1)
    states = np.empty((n_steps + 1, psi.size), dtype=complex)
    states[0] = psi
    for k in range(n_steps):
        psi = unitary @ psi
        states[k + 1] = psi
    return states


def trajectory_fidelity(psi_trotter, psi_ed):
    """|⟨ψ_ED(t)|ψ_Trotter(t)⟩|^2 at each stored time.

    psi_ed is the array evolve() returns, shape (n_times, dim). This is 1 at
    t = 0. How far it falls by the last time is the ED-versus-Trotter
    comparison Figs. 3–6 are built from, before any observable is measured.
    """
    overlap = np.sum(np.conj(psi_ed) * psi_trotter, axis=-1)
    return np.abs(overlap) ** 2


def quench_compare(N, a, m, g, eps_values, times, dt):
    """Trotter and exact trajectories for each quench field.

    psi0 is the ε = 0 ED ground state. For each ε, one Eq. (22) step of
    length dt is repeated across `times`, and evolve(psi0, H(ε), times)
    produces the exact states at those same instants. The returned fidelity
    is |⟨ψ_exact(t)|ψ_Trotter(t)⟩|^2. `times` must be uniformly spaced by dt,
    which is the Δt = 0.1 grid of Section III.C when the caller passes
    linspace(0, 12, 121).
    """
    eps_values = np.asarray(eps_values, dtype=float)
    times = np.asarray(times, dtype=float)
    if times.size < 2 or not np.allclose(np.diff(times), dt):
        raise ValueError("times must be uniformly spaced by the Trotter step dt.")
    psi0 = zero_field_ground_state(N, a, m, g)
    n_steps = times.size - 1
    psi_trotter = np.empty((eps_values.size, times.size, psi0.size), dtype=complex)
    psi_ed = np.empty_like(psi_trotter)
    fidelity = np.empty((eps_values.size, times.size))
    for i, eps in enumerate(eps_values):
        H, op = pauli_hamiltonian(N, a, m, g, eps)
        print(f"    eps={eps:g}: building one Δt step", flush=True)
        unitary = one_step_unitary(op, dt)
        psi_trotter[i] = trotter_trajectory(psi0, unitary, n_steps)
        psi_ed[i] = evolve(psi0, H, times)
        fidelity[i] = trajectory_fidelity(psi_trotter[i], psi_ed[i])
        print(
            f"    eps={eps:g}: F(0)={fidelity[i, 0]:.6f}  "
            f"min F={fidelity[i].min():.6f}  F(T)={fidelity[i, -1]:.6f}",
            flush=True,
        )
    return psi0, psi_trotter, psi_ed, fidelity
