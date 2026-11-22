"""Schwinger operators built directly in a fixed-charge computational basis.

Site l=1 is the most significant bit, exactly as in schwinger_model.py.
Consequently it is Qiskit qubit N-l. No full 2**N matrix is built to project
onto a sector. H(eps) = H0 + eps*diag(V) + c*eps**2*I.
"""
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy import sparse
from scipy.linalg import eigh_tridiagonal
from scipy.sparse.linalg import eigsh


def charge_basis(N, Q=0):
    """Sorted bit strings and inverse map for integer charge Q (even N only)."""
    if not isinstance(N, (int, np.integer)) or N < 2 or N % 2:
        raise ValueError("Use an even number of sites, N >= 2")
    if Q is None or int(Q) != Q or not -N//2 <= Q <= N//2:
        raise ValueError("Q must be an integer between -N/2 and N/2")
    basis = np.array(sorted(sum(1 << b for b in bits)
                            for bits in combinations(range(N), N//2-int(Q))), dtype=np.int64)
    return basis, {int(state): i for i, state in enumerate(basis)}


def neutral_basis(N):
    return charge_basis(N, 0)


@dataclass
class Operators:
    N: int
    a: float
    m: float
    g: float
    Q: int
    basis: np.ndarray
    basis_index: dict
    H0: sparse.csr_matrix
    V: np.ndarray
    c: float
    q: np.ndarray
    occupation: np.ndarray

    def hamiltonian(self, eps=0.0, scalar=True):
        """Sparse physical Hamiltonian; scalar=False removes only the identity."""
        return self.H0 + sparse.diags(eps*self.V + (self.c*eps**2 if scalar else 0.0))


def build_sector_operators(N, a=1.0, m=1.0, g=1.0, Q=0, boundary="open"):
    """Build mass/electric diagonals and unequal-bit nearest-neighbor hopping.

    q and occupation have shape (sector dimension, number of sites).
    The cumulative charge L0 supplies the electric energy and field coupling.
    """
    if boundary != "open" or a <= 0 or g <= 0 or m < 0:
        raise ValueError("Require open boundaries, a,g > 0 and m >= 0")
    basis, inverse = charge_basis(N, Q)
    bits = (basis[:, None] >> np.arange(N-1, -1, -1)) & 1
    z = 1 - 2*bits
    stagger = (-1.)**np.arange(N)
    q = (z+stagger)/2
    L0 = np.cumsum(q, axis=1)[:, :-1]
    diagonal = (m/2)*(z@stagger) + (g*g*a/2)*np.sum(L0*L0, axis=1)
    rows, columns = [], []
    for l in range(N-1):
        left = np.flatnonzero((bits[:, l] == 0) & (bits[:, l+1] == 1))
        mask = (1 << (N-1-l)) | (1 << (N-2-l))
        right = np.array([inverse[int(basis[i]) ^ mask] for i in left], dtype=int)
        rows.extend(left); columns.extend(right)
        rows.extend(right); columns.extend(left)
    hopping = sparse.coo_matrix((np.full(len(rows), 1/(2*a)), (rows, columns)),
                                 shape=(len(basis), len(basis))).tocsr()
    return Operators(N, a, m, g, Q, basis, inverse,
                     hopping+sparse.diags(diagonal), g*g*a*np.sum(L0, axis=1),
                     g*g*a*(N-1)/2, q, q*q)


def sector_ground_state(ops, eps=0.0):
    """Lowest state in the requested sector, with a deterministic sparse start."""
    H = ops.hamiltonian(eps)
    if H.shape[0] <= 1024:
        energies, states = np.linalg.eigh(H.toarray())
        return float(energies[0]), states[:, 0].astype(complex)
    energies, states = eigsh(H, k=1, which="SA", tol=1e-12,
                            v0=np.linspace(1, 2, H.shape[0]))
    return float(energies[0]), states[:, 0].astype(complex)


def transition_spectrum(ops, eps=0.0):
    """All small-sector eigenstates and field-coupling weights.

    Degenerate eigenspaces are grouped by summed projector weight for response.
    Individual eigenvectors remain available for the state-population figure.
    Large sectors use lanczos_spectrum instead of a dense diagonalization.
    """
    H = ops.hamiltonian(eps)
    if H.shape[0] > 1024:
        raise ValueError("Use lanczos_spectrum above dimension 1024")
    energies, states = np.linalg.eigh(H.toarray())
    vacuum = states[:, 0]
    drive = ops.V*vacuum
    variance = float(np.vdot(drive, drive).real-abs(np.vdot(vacuum, drive))**2)
    gaps = energies-energies[0]
    weights = abs(states.conj().T@drive)**2
    weights[0] = 0
    keep = gaps > 1e-10
    grouped_gaps, grouped_weights = [], []
    for gap, weight in zip(gaps[keep], weights[keep]):
        if grouped_gaps and abs(gap-grouped_gaps[-1]) < 1e-9:
            grouped_weights[-1] += weight
        else:
            grouped_gaps.append(float(gap)); grouped_weights.append(float(weight))
    return {"gaps": np.array(grouped_gaps), "weights": np.array(grouped_weights),
            "energies": energies, "vectors": states, "variance": variance,
            "captured_weight": float(weights[keep].sum())}


def lanczos_spectrum(ops, eps=0.0, steps=128, ground=None):
    """Gauss/Lanczos measure of (V-<V>)|0>, including high-energy transitions.

    This is a quadrature measure, not individually identified eigenstates.
    Full reorthogonalization prevents duplicate Ritz weights. Repeat with
    more steps and compare the prediction for the actual waveform.
    """
    E0, vacuum = sector_ground_state(ops, eps) if ground is None else ground
    start = (ops.V-np.vdot(vacuum, ops.V*vacuum).real)*vacuum
    variance = float(np.vdot(start, start).real)
    if variance < 1e-28:
        return {"gaps": np.array([]), "weights": np.array([]), "variance": variance,
                "captured_weight": 0., "steps": 0}
    H = ops.hamiltonian(eps)
    vectors, alpha, beta = [], [], []
    v = start/np.sqrt(variance)
    for j in range(min(steps, len(v))):
        vectors.append(v)
        w = H@v-E0*v
        alpha.append(float(np.vdot(v, w).real))
        # Two passes keep the Krylov basis orthogonal near converged levels.
        basis = np.asarray(vectors)
        for _ in range(2):
            w -= (basis.conj()@w)@basis
        b = np.linalg.norm(w)
        if b < 1e-12 or j == min(steps, len(v))-1:
            break
        beta.append(b)
        v = w/b
    gaps, vectors = eigh_tridiagonal(np.array(alpha), np.array(beta))
    weights = variance*abs(vectors[0])**2
    return {"gaps": gaps, "weights": weights, "variance": variance,
            "captured_weight": float(weights.sum()), "steps": len(alpha)}
