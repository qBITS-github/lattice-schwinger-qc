## Observables

The paper discusses the following observables:

1. Total charge

$$
Q_N
=
\frac{1}{2}
\sum_{n=1}^{N}
\left(
\sigma_n^z + (-1)^n
\right)
$$

2. Local charge density on the $n^{th}$ lattice site

$$
q_n(t)
=
\frac{1}{2}
\left(
\left\langle \sigma_n^z(t) \right\rangle
+
(-1)^n
\right)
$$

3. Local charge density on each spatial point

$$
Q_i(t)
=
q_{2i-1}(t) + q_{2i}(t)
$$

4. Total electric-field energy

$$
H_E(t)
=
\frac{g^2 a}{2}
\sum_{l=1}^{N-1}
\left\langle
E_l^2(t)
\right\rangle
$$

5. Vacuum-state fidelity

$$
P_{\mathrm{vac}}(t)
=
\left|
\left\langle \psi_0 \middle| \psi(t) \right\rangle
\right|^2
$$

6. Chiral condensate density

$$
\Gamma
=
\langle \bar{\psi}\psi\rangle
=
\frac{1}{2Na}
\sum_{n=1}^{N}
(-1)^n
\left\langle
\sigma_n^z
\right\rangle
$$


## Implementing the observables

**1. Chiral condensate**
Loop over each site $n$, compute the expectation value $\langle\sigma^z_n\rangle$ by sandwiching the operator between $\langle\psi|$ and $|\psi\rangle$, multiply by the staggering sign $(-1)^n$, and sum. Divide the total by $2Na$ at the end.

**2. Total charge**
Same site-by-site pattern as above: compute $\langle\sigma^z_n\rangle$ at each site, add $(-1)^n$, sum over all sites, divide by 2 at the end. Implemented independently rather than reusing another function.

**3. Local charge**
Also the same per-site pattern, but instead of summing everything into one number, each site's value $q_n = \tfrac12(\langle\sigma^z_n\rangle+(-1)^n)$ is kept separately and returned as an array. This is the building block the next observable uses.

**4. Spatial point charge**
Calls the local charge function above to get every $q_n$, then pairs up neighboring sites (one physical point = two adjacent staggered sites) by adding consecutive entries together.

**5. Electric field energy**
Reconstructs the same running cumulative sum $L_n$ that's used inside the Hamiltonian, but here it's evaluated as an expectation value $\langle\psi|L_n^2|\psi\rangle$ rather than left as an operator. Summed over all links with the prefactor applied.

**6. Vacuum-state fidelity**
No per-site loop at all — just the overlap $\langle\psi_0|\psi(t)\rangle$ between the fixed zero-field ground state and whatever state is passed in, squared in magnitude. Works on a single state or a whole batch of states at once, so it can be handed the entire output of the time evolution in one shot.

## Time evolution

Any quantum state evolves as $|\psi(t)\rangle = e^{-iHt}|\psi_0\rangle$. Writing $H$ in its own eigenbasis, $H|k\rangle = E_k|k\rangle$, this becomes

$$|\psi(t)\rangle = \sum_k c_k\, e^{-iE_k t}\,|k\rangle, \qquad c_k = \langle k|\psi_0\rangle$$

Each eigenstate just picks up a phase that spins at a rate set by its own energy, and the evolved state is the sum of all these spinning pieces.

- The code does exactly this in three steps.
    - First, diagonalize $H$ once to get all the energies $E_k$ and eigenstates $|k\rangle$. 
    - Second, project the initial state onto that eigenbasis to get the coefficients $c_k$. 
    - Third, for each requested time $t$, multiply each $c_k$ by its phase $e^{-iE_kt}$ and add the pieces back together in the original basis to get $|\psi(t)\rangle$.

Diagonalizing is the expensive part, so it's done only once at the start; each additional time point after that is just a cheap phase multiply and recombine, which is why sweeping through many times in a loop stays fast.

---
