# Building the Hamiltonian

The Pauli Hamiltonian is given in Eq. (8)-(11) of the paper. The Jordan-Wigner transformation has already been applied at this point, so we do not need to touch fermionic operators at all.

$$
H = H_{\mathrm{kin}} + H_m + H_E
$$

The kinetic term is:

$$
H_{\mathrm{kin}}
=
\frac{1}{4a}
\sum_{n=1}^{N-1}
\left(
\sigma_n^x \sigma_{n+1}^x
+
\sigma_n^y \sigma_{n+1}^y
\right)
$$

The mass term is:

$$
H_m
=
\frac{m}{2}
\sum_{n=1}^{N}
(-1)^n \sigma_n^z
$$

The electric-field term is:
$$
H_E
=
\frac{g^2 a}{2}
\sum_{n=1}^{N-1}
\left[
\epsilon
+
\frac{1}{2}
\sum_{l=1}^{n}
\left(
\sigma_l^z + (-1)^l
\right)
\right]^2
$$

## Physical Intuition

- In the basis used for exact diagonalization, every site is either occupied $(\sigma^{z} = +1)$ or empty $(\sigma^{z} = -1)$.
- $H_{\mathrm{kin}}$: The $XX + YY$ term is a hopping term — it moves fermions between neighboring sites. If the two sites differ, it swaps them; if they agree, it does nothing. It is non-diagonal.
- $H_m$: This is rest energy. It depends on which sites are occupied and on nothing else. It is diagonal.
- $H_E$: Gauss's law fixes the flux on every link. Link n carries the boundary value $\epsilon$ plus the net physical charge on sites 1 through n. It is diagonal.

So, $H_{\mathrm{kin}}$ is the only genuinely "quantum" component here since it forces superpositions.

## Setting up the simulation

We first create the Pauli spin matrices acting on one qubit.

- Identity matrix $I_2$ = 
$
\begin{bmatrix}
  1 & 0 \\
  0 & 1 
\end{bmatrix}
$

- Pauli-X (bit flip) $S_x$ = 
$
\begin{bmatrix}
  0 & 1 \\
  1 & 0 
\end{bmatrix}
$

- Pauli-Y (bit + phase flip) $S_y$ = 
$
\begin{bmatrix}
  0 & -i \\
  i & 0 
\end{bmatrix}
$

- Pauli-Z (phase flip) $S_z$ = 
$
\begin{bmatrix}
  1 & 0 \\
  0 & -1 
\end{bmatrix}
$

We then build the entire N-qubit operator i.e. the matrix representation of $\sigma_n^{\alpha}$ acting on the full $2^N$ Hilbert space.

- We start with a list of $N$ identities, then swap in the operator at position n, then take the tensor product of the whole list.
- This is equivalent to applying the operator only to the $n^{th}$ qubit while applying identity on the rest.

## Building the Hamiltonian

### 1. Kinetic term:
- Loop over indices $0$ to $N-2$ (one link at a time):
    - Apply $S_x$ to sites $n$ and $n+1$, apply $S_y$ to sites $n$ and $n+1$
    - Add $\frac{1}{4a}\left(\sigma_n^x\sigma_{n+1}^x + \sigma_n^y\sigma_{n+1}^y\right)$ to the running total
- Get $H_{\mathrm{kin}}$

### 2. Mass term:
- Loop over all N sites
- Multiply $\sigma_n^{z}$ by the scalar coefficient $\frac{m}{2}(-1)^n$
- Get $H_m$

### 3. Electric Field term:
- For every link $n$, construct an operator $L_n$, square it, then add it to the Hamiltonian
- We repeatedly accumulate terms on $L_n$ rather than recomputing from scratch at every $n$
- We then implement the outer sum to get $H_E$

## Diagonalization

We first check if the matrix is Hermitian ($H = H^\dagger$; for this real-valued $H$ that's equivalent to ordinary symmetry) and then use `np.linalg.eigh` to return all its eigenvalues/eigenvectors

---



