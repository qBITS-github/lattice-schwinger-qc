"""
Lattice Schwinger model Hamiltonian (Chen, Cheng, Guo, arXiv:2607.02894),
built directly in the qubit/Pauli representation of Eqs. (8)-(11):

    H = H_kin + H_m + H_E

    H_kin = (1/4a) * sum_{n=1}^{N-1} (X_n X_{n+1} + Y_n Y_{n+1})
    H_m   = (m/2)  * sum_{n=1}^{N}   (-1)^{n+1} Z_n
    H_E   = (g^2 a/2) * sum_{n=1}^{N-1} [ eps + (1/2) sum_{l=1}^{n} (Z_l + (-1)^{l+1}) ]^2

Site indices n = 1..N in the paper <-> python indices idx = 0..N-1 (idx = n-1).
No Jordan-Wigner transformation needs to be done by hand -- the paper already
expresses H purely in terms of single- and two-qubit Pauli operators.

The paper's Eqs. (10)-(13) and (17) write the stagger as (-1)^n. With that
sign, a quench from the zero-field vacuum does not separate charge (Fig. 4
stays flat) and E_0(ε) misses Fig. 1. The sign (-1)^{n+1} is the stagger of a
zero-based site index, applied in the mass, the Gauss law, and every charge
observable. At ε = 0 it is equivalent to Z_n -> -Z_n, so the Table I energy
is unchanged. At ε > 0 it is the relative sign between the boundary field and
the staggered charge that reproduces the ED curves in Figs. 1, 4, 5, and 8.
"""

import numpy as np
from scipy import sparse

# --- single-qubit Pauli matrices (sparse) ---
I2 = sparse.identity(2, format="csr", dtype=complex)
SX = sparse.csr_matrix([[0, 1], [1, 0]], dtype=complex)
SY = sparse.csr_matrix([[0, -1j], [1j, 0]], dtype=complex)
SZ = sparse.csr_matrix([[1, 0], [0, -1]], dtype=complex)


def stagger(n):
    """(-1)^{n+1} for 1-based site n. See the module docstring."""
    return (-1) ** (n + 1)


def op_on_site(op, n, N):
    """Embed a single-qubit operator `op` on site n (0-indexed) of an N-qubit chain."""
    mats = [I2] * N
    mats[n] = op
    result = mats[0]
    for m in mats[1:]:
        result = sparse.kron(result, m, format="csr")
    return result


def build_hamiltonian(N, a=1.0, m=1.0, g=1.0, eps=0.0):
    """
    Build the N-qubit lattice Schwinger Hamiltonian as a sparse matrix
    of dimension 2^N x 2^N, exactly per Eqs. (9)-(11).
    """
    dim = 2 ** N
    identity_full = sparse.identity(dim, format="csr", dtype=complex)

    # cache single-site Z operators (idx = 0..N-1 <-> n = idx+1)
    Zs = [op_on_site(SZ, idx, N) for idx in range(N)]

    H = sparse.csr_matrix((dim, dim), dtype=complex)

    # --- H_kin: hopping term, n = 1..N-1 -> idx = 0..N-2 ---
    for idx in range(N - 1):
        Xn = op_on_site(SX, idx, N)
        Xn1 = op_on_site(SX, idx + 1, N)
        Yn = op_on_site(SY, idx, N)
        Yn1 = op_on_site(SY, idx + 1, N)
        H = H + (1.0 / (4 * a)) * (Xn @ Xn1 + Yn @ Yn1)

    # --- H_m: staggered mass term, n = 1..N -> idx = 0..N-1 ---
    for idx in range(N):
        n = idx + 1
        H = H + (m / 2.0) * stagger(n) * Zs[idx]

    # --- H_E: electric-field energy, n = 1..N-1 -> idx = 0..N-2 ---
    # L_n = eps + (1/2) * cumulative_{l=1}^{n} (Z_l + (-1)^{l+1})
    cum = sparse.csr_matrix((dim, dim), dtype=complex)
    for n in range(1, N):
        idx = n - 1  # this iteration adds the l = n term to the running sum
        cum = cum + 0.5 * (Zs[idx] + stagger(n) * identity_full)
        Ln = eps * identity_full + cum
        H = H + (g ** 2 * a / 2.0) * (Ln @ Ln)

    return H


def diagonalize(H, k=None):
    """
    Full or partial diagonalization.
    k=None  -> dense np.linalg.eigh, returns ALL eigenvalues/eigenvectors
               (fine for N <= ~12, dim <= 4096).
    k=int   -> sparse eigsh for the k lowest eigenpairs (use for larger N).
    """
    if k is None:
        Hd = H.toarray()
        assert np.allclose(Hd, Hd.conj().T, atol=1e-10), "H is not Hermitian!"
        evals, evecs = np.linalg.eigh(Hd)
    else:
        from scipy.sparse.linalg import eigsh
        evals, evecs = eigsh(H, k=k, which="SA")
        order = np.argsort(evals)
        evals, evecs = evals[order], evecs[:, order]
    return evals, evecs


def chiral_condensate(psi, N, a=1.0):
    """Gamma = (1/(2N a)) sum_n (-1)^{n+1} <Z_n>, Eq. (17) with the stagger of stagger()."""
    val = 0.0
    for idx in range(N):
        n = idx + 1
        Zn = op_on_site(SZ, idx, N)
        expval = np.real(psi.conj() @ (Zn @ psi))
        val += stagger(n) * expval
    return val / (2 * N * a)


def total_charge(psi, N):
    """Q_N = (1/2) sum_n (Z_n + (-1)^{n+1}), Eq. (12) with the stagger of stagger()."""
    val = 0.0
    for idx in range(N):
        n = idx + 1
        Zn = op_on_site(SZ, idx, N)
        expval = np.real(psi.conj() @ (Zn @ psi))
        val += expval + stagger(n)
    return val / 2


def local_charge(psi, N):
    """
    q_n(t) = (1/2)(<Z_n> + (-1)^{n+1}), Eq. (13) with the stagger of stagger() -- the staggered-site charge
    density, site by site. Returns array of length N (index idx = n-1).
    """
    q = np.zeros(N)
    for idx in range(N):
        n = idx + 1
        Zn = op_on_site(SZ, idx, N)
        expval = np.real(psi.conj() @ (Zn @ psi))
        q[idx] = 0.5 * (expval + stagger(n))
    return q


def spatial_point_charge(psi, N):
    """
    Q_i(t) = q_{2i-1}(t) + q_{2i}(t), Eq. (14): pairs adjacent staggered
    sites into physical spatial points. Returns array of length N//2.
    """
    q = local_charge(psi, N)
    return q[0::2] + q[1::2]


def vacuum_fidelity(psi_t, psi0):
    """
    P_vac(t) = |<psi0|psi(t)>|^2, Eq. (16).
    psi_t: either a single state vector, or an array of shape (n_times, dim)
    as returned by `evolve` -- vectorized over the time axis in that case.
    """
    overlap = np.asarray(psi_t) @ psi0.conj()
    return np.abs(overlap) ** 2


def field_energy(psi, N, a=1.0, g=1.0, eps=0.0):
    """H_E(t) expectation value, Eq. (15)/(11)."""
    dim = 2 ** N
    identity_full = sparse.identity(dim, format="csr", dtype=complex)
    Zs = [op_on_site(SZ, idx, N) for idx in range(N)]
    cum = sparse.csr_matrix((dim, dim), dtype=complex)
    total = 0.0
    for n in range(1, N):
        idx = n - 1
        cum = cum + 0.5 * (Zs[idx] + stagger(n) * identity_full)
        Ln = eps * identity_full + cum
        L2 = Ln @ Ln
        total += np.real(psi.conj() @ (L2 @ psi))
    return (g ** 2 * a / 2.0) * total


def evolve(psi0, H, times):
    """
    Exact real-time evolution psi(t) = exp(-i H t) psi0 at each t in `times`,
    via eigendecomposition of H (reused across all times -- efficient for
    scanning many t values, and this *is* the paper's ED benchmark curve).
    Returns array of shape (len(times), dim).
    """
    evals, evecs = diagonalize(H)  # H = evecs @ diag(evals) @ evecs^dagger
    c0 = evecs.conj().T @ psi0  # overlap coefficients in the energy eigenbasis
    psis = []
    for t in times:
        phase = np.exp(-1j * evals * t)
        psis.append(evecs @ (phase * c0))
    return np.array(psis)