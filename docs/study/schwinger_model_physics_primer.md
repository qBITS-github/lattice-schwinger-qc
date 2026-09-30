# Physics Primer: The Lattice Schwinger Model
### A from-first-principles derivation of Section II of Chen, Cheng & Guo (arXiv:2607.02894)

**Scope.** This document covers exactly what the roadmap specifies: enough physics to construct, understand, and debug the qubit Hamiltonian used in the paper (their Eqs. 1–17). It assumes only first-year undergraduate mathematics (linear algebra, ODEs, basic complex numbers) and physics (Newtonian mechanics, a first pass at Lagrangian/Hamiltonian mechanics, and Maxwell's equations at the level of "E and B fields, charge conservation"). Everything else — fields, Dirac fermions, gauge symmetry, second quantization, Jordan-Wigner — is built from scratch.

**How to read this.** Each of the five roadmap items is tagged **[Intuition]** or **[Mechanics]**. Intuition-tier sections are meant to be read once, understood conceptually, and moved past. The mechanics-tier section (Gauss's law elimination) is the one place where you should stop and reproduce the algebra by hand — it is flagged accordingly, with a full worked N=4 example, exactly as the roadmap prescribes.

---

## 0. Foundations: bridging first-year physics to field theory

You know Newtonian mechanics as: state = (position, momentum) for each particle, and dynamics from $F=ma$ or, better, from a Hamiltonian $H(q,p)$ via Hamilton's equations
$$
\dot q = \frac{\partial H}{\partial p}, \qquad \dot p = -\frac{\partial H}{\partial q}.
$$

**A field is just what you get when you promote "one particle, one position" to "infinitely many degrees of freedom, one at every point in space."** Instead of a single $q(t)$, you have $\phi(x,t)$ — a number (or operator) attached to *every point $x$*. Its dynamics still comes from a Hamiltonian, now written as an integral over space:
$$
H = \int dx\, \mathcal{H}(\phi, \pi, x),
$$
where $\mathcal{H}$ is the **Hamiltonian density** and $\pi$ is the momentum conjugate to $\phi$ at each point (exactly analogous to $p$ conjugate to $q$).

The paper's Eq. (2),
$$
H = \int dx \left[-i\bar\psi \gamma^1 \partial_1 \psi + m\bar\psi\psi + \tfrac{1}{2}E^2\right],
$$
is precisely this: a sum (integral) of a "kinetic energy density," a "mass/potential energy density," and a "field energy density" ($\tfrac12 E^2$, which you already know from $E\&M$ as the energy density stored in an electric field, in units where $\epsilon_0=1$).

**What's new relative to first-year physics:**

1. **The matter field $\psi$ is a *fermionic* quantum field**, not a classical number. Concretely: at each point (or, once discretized, each lattice site) there is a two-outcome quantum degree of freedom — "is there a fermion here or not" — obeying the Pauli exclusion principle. The cleanest way to encode "occupied / empty" is with **creation and annihilation operators** $\phi^\dagger, \phi$ satisfying the *anticommutation* relation
$$
\{\phi, \phi^\dagger\} \equiv \phi\phi^\dagger + \phi^\dagger\phi = 1, \qquad \phi^2 = (\phi^\dagger)^2 = 0.
$$
   The last relation *is* the Pauli exclusion principle in operator form: you cannot create the same fermion twice. The number operator $\hat n = \phi^\dagger\phi$ has eigenvalues $0$ (empty) or $1$ (occupied) only — check this yourself: $\hat n^2 = \phi^\dagger\phi\phi^\dagger\phi = \phi^\dagger(1-\phi^\dagger\phi)\phi = \phi^\dagger\phi - \phi^\dagger\phi^\dagger\phi\phi = \hat n$, so $\hat n(\hat n - 1)=0$.

2. **Pauli matrices as the language of two-level systems.** A single qubit is described by
$$
\sigma^x = \begin{pmatrix}0&1\\1&0\end{pmatrix},\quad
\sigma^y = \begin{pmatrix}0&-i\\i&0\end{pmatrix},\quad
\sigma^z = \begin{pmatrix}1&0\\0&-1\end{pmatrix},
$$
   with $\sigma^z|{\uparrow}\rangle = +|{\uparrow}\rangle$, $\sigma^z|{\downarrow}\rangle=-|{\downarrow}\rangle$, and raising/lowering operators $\sigma^\pm = (\sigma^x \pm i\sigma^y)/2$, which satisfy $\sigma^+|{\downarrow}\rangle = |{\uparrow}\rangle$, $\sigma^+|{\uparrow}\rangle=0$. These will turn out to be exactly what you need to represent "occupied/empty" once we get to Jordan-Wigner (Section 4).

3. **The electric field $E$ becomes, after discretization, a single number per link, $L_n$** (Section 2). This is the one genuinely new physical idea in the paper: because the theory lives in one spatial dimension, this field has no independent dynamics of its own (Section 3) — it's completely slaved to where the charges are.

With this vocabulary, the whole paper is: *(a)* put the continuum theory (Eq. 1–2) on a lattice, *(b)* use the fact that in 1D the gauge field carries no independent information to eliminate it algebraically, *(c)* map the remaining fermionic problem to qubits, *(d)* solve/evolve the resulting qubit Hamiltonian on a (simulated) quantum computer.

---

## 1. Gauge invariance — [Intuition]

### 1.1 Global symmetry, first (something you already half-know)

Consider a complex scalar-like field $\psi$ appearing in a Lagrangian only through $\bar\psi\psi$ and $\bar\psi\partial\psi$ type combinations (this is the structure of Eq. 1). If you multiply $\psi$ by a single, constant phase everywhere in space,
$$
\psi(x) \to e^{i\alpha}\psi(x), \qquad \alpha = \text{constant (the same at every point)},
$$
nothing observable changes — $\bar\psi\psi = \psi^\dagger\gamma^0\psi \to \psi^\dagger e^{-i\alpha}\gamma^0 e^{i\alpha}\psi = \bar\psi\psi$, unchanged. This is a **global U(1) symmetry** ("U(1)" just means "phases," i.e., the group of $1\times1$ unitary matrices $e^{i\alpha}$). By Noether's theorem — the general fact that every continuous symmetry of the dynamics implies a conserved quantity — this symmetry implies a conserved charge, which here is exactly the total electric charge (Eq. 12 in the paper, $Q_N$, is the lattice avatar of this).

### 1.2 Why *local* symmetry forces a new field into existence

Now ask a harder question: can you let the phase depend on position, $\alpha \to \alpha(x)$, and still have a symmetry? That is,
$$
\psi(x) \to e^{i\alpha(x)}\psi(x).
$$
Try it on the kinetic term $\bar\psi \gamma^1 \partial_1 \psi$. The derivative doesn't just see $\psi$; by the product rule it also picks up $\partial_1\alpha(x)$:
$$
\partial_1\left(e^{i\alpha(x)}\psi\right) = e^{i\alpha(x)}\left(\partial_1\psi + i(\partial_1\alpha)\psi\right).
$$
The extra term $i(\partial_1\alpha)\psi$ does **not** cancel — a *local* phase redefinition is no longer a symmetry of the free kinetic term, because "how fast $\psi$ is changing" now depends on how fast you *chose to twist your phase convention* at each point, which isn't physical.

The fix, first found by Weyl and is the origin of essentially every gauge theory in physics: introduce a new field $A_\mu$ (here, in 1+1D, just $A_1 \equiv A$) that transforms in exactly the compensating way,
$$
A_1(x) \to A_1(x) - \tfrac{1}{g}\partial_1\alpha(x),
$$
and replace the ordinary derivative by the **covariant derivative** $D_1 = \partial_1 + igA_1$ everywhere it appears (this is exactly the $D_\mu = \partial_\mu + igA_\mu$ in the paper's Eq. 1). Now check: under the combined transformation of $\psi$ *and* $A$,
$$
D_1\psi \to \left(\partial_1 + ig A_1 - i\partial_1\alpha\right)\left(e^{i\alpha}\psi\right) = e^{i\alpha}\big(\partial_1\psi + i(\partial_1\alpha)\psi + igA_1\psi - i(\partial_1\alpha)\psi\big) = e^{i\alpha}D_1\psi,
$$
i.e., $D_1\psi$ transforms exactly like $\psi$ itself (an overall phase), and the extra piece cancels *because $A$ absorbed it*. So: **a field that must exist purely to make a local symmetry consistent** is a gauge field. It is not optional decoration; it is forced into existence by the requirement of local phase freedom. This gauge field's field-strength ($F_{\mu\nu}$, giving the $-\tfrac14 F_{\mu\nu}F^{\mu\nu}$ term in Eq. 1) is exactly the electromagnetic field, and in 1+1D reduces to a single electric field $E$.

**Two-sentence summary (the roadmap's exercise):** A global symmetry costs nothing extra because you're applying the same transformation everywhere, so all the "bookkeeping" cancels identically. A local symmetry lets each point twist independently, so derivatives (which compare *neighboring* points) pick up an uncancelled leftover — you need a new field whose entire job is to "connect" neighboring points and absorb that leftover, and that field is the gauge field.

**Why this matters for the paper:** the gauge field $A$ (later $U_n$ on the lattice) is not an independent thing you get to choose freely — its whole existence is dictated by matter's need for local phase invariance. This is why, in Section 3, you'll be able to *eliminate* it entirely in favor of the matter configuration: it was never truly independent to begin with, at least not in one spatial dimension.

---

## 2. Lattice discretization — [Intuition]

### 2.1 From continuum to a grid

To simulate Eq. (1)–(2) on a computer (classical or quantum), space must become a finite set of points. Put $N$ sites at positions $x_n = na$, $n=1,\dots,N$, with lattice spacing $a$. The matter field becomes one operator $\phi_n$ per site; the gauge field, being naturally associated with the *link between* two sites (it's what you integrate along a path to compare phases at different points — a **Wilson line**), becomes one operator $U_n$ per link between sites $n$ and $n+1$, and correspondingly one electric field value $L_n$ per link.

```
matter:   φ_1    φ_2    φ_3    φ_4   ...   φ_N
sites:     ●------●------●------●-----...---●
links:        U_1,L_1  U_2,L_2  U_3,L_3      (N−1 links for open boundaries)
```

This directly explains the index ranges in Eq. (3): the hopping and electric-field sums run over $n=1,\dots,N-1$ (one term per **link**), while the mass sum runs over $n=1,\dots,N$ (one term per **site**).

### 2.2 Why staggered fermions, and the $(-1)^n$

A relativistic fermion (a Dirac fermion) in 1+1D has **two** components (particle-like and antiparticle-like), reflecting the two solutions of the Dirac equation (positive- and negative-energy solutions). Naively discretizing a two-component field on a single lattice with a naive derivative causes a notorious problem (**fermion doubling** — extra, unphysical low-energy fermion species appear; you don't need to derive this, just know it's why staggering exists). The **Kogut-Susskind staggered fermion** trick sidesteps it by using only **one** fermionic degree of freedom per site, but interpreting *even* and *odd* sites as the two different components of the original Dirac fermion — i.e., $N$ single-component sites together encode $N/2$ full two-component Dirac fermions.

This is the origin of the factor $(-1)^n$ that appears throughout:

- In the **mass term**, $m\sum_n (-1)^n \phi_n^\dagger\phi_n$ (Eq. 3): particle-type sites and antiparticle-type sites must contribute with *opposite sign* to the mass energy, because in the continuum theory a filled antiparticle state (a "hole" in the Dirac sea) has energy $+m$ relative to the vacuum only if the vacuum convention is set correctly — the alternating sign is exactly what makes the discretized theory reproduce a fermion **and** an antiparticle with the same mass $m$, rather than two fermions of different sign.
- In the **Gauss law**, the term $\tfrac{1-(-1)^n}{2}$ (Eq. 4) equals $1$ on odd sites and $0$ on even sites. This is the **background charge of the Dirac sea**: in the interacting vacuum, odd (antiparticle-type) sites are conventionally taken to be "filled" in the reference/vacuum configuration, so their contribution to the physical charge must have this offset subtracted off. You'll see exactly why this specific combination is needed once you do the Gauss-law derivation in Section 3.

**Roadmap's memorization check:** for $n=1,2,3,4,5,6,7,8$, $(-1)^n = -1,+1,-1,+1,-1,+1,-1,+1$. Odd $\Rightarrow -1$ (antiparticle-type, "background charge $+1$" needing subtraction), even $\Rightarrow +1$ (particle-type, "background charge $0$"). Internalizing this instantly is a genuine debugging tool: almost every sign bug in this kind of code is a staggering mismatch, and eyeballing "is $n$ odd or even, and does the sign in my code match" catches most of them without running anything.

`src/schwinger_model.py` writes the stagger as $(-1)^{n+1}$ on this same 1-based label. That is $(-1)$ raised to the Python index $n-1$. The extra $+1$ is only because Python indexes from zero. It does not change the physics.

---

## 3. Gauss's law elimination (Eqs. 4–6, 11) — [MECHANICS]

This is the one section where you should stop and do the algebra yourself before reading further. What follows is a complete worked version, matching the roadmap's four-step procedure, for $N=4$.

### 3.1 Why Gauss's law even applies here

In ordinary 3+1D electromagnetism, Gauss's law is $\nabla\cdot E = \rho$ — a *constraint*, not a dynamical equation, relating the electric field's spatial derivative to the local charge density. It has no time derivative in it, so it never tells you how $E$ evolves; it just restricts which field configurations are physically allowed at any instant.

In 1 spatial dimension, "$\nabla\cdot E$" degenerates to an ordinary derivative $\partial_x E$, and a first-order spatial derivative with no dynamics of its own means: **given the charge distribution and one boundary value, the entire field configuration is fixed algebraically.** There is no room left over for $E$ to have independent wave-like dynamics (that's a special feature of $D=1$ spatial dimension: a "field" with only one component and a first-order constraint has zero propagating degrees of freedom). This is the deep reason the paper can eliminate $L_n$ entirely — it was never carrying independent information.

### 3.2 Step 1 — read the recursion (Eq. 4)

$$
L_n - L_{n-1} = \phi_n^\dagger\phi_n - \frac{1-(-1)^n}{2}.
$$

Define the **background-subtracted charge on site $n$**:
$$
\rho_n \equiv \phi_n^\dagger\phi_n - \frac{1-(-1)^n}{2} =
\begin{cases}
\phi_n^\dagger\phi_n - 1, & n\ \text{odd} \\[2pt]
\phi_n^\dagger\phi_n, & n\ \text{even}.
\end{cases}
$$
Then Eq. (4) is simply $L_n = L_{n-1} + \rho_n$: **each link's field equals the previous link's field, plus the (background-subtracted) charge that just "passed" that site.** This is exactly discrete Gauss's law: the field only changes across a site if there's a charge sitting there.

### 3.3 Step 2 — solve the recursion by hand for $N=4$

With boundary condition $L_0 = \varepsilon$ (the imposed external field), unroll the recursion:
$$
\begin{aligned}
L_1 &= L_0 + \rho_1 = \varepsilon + \rho_1,\\
L_2 &= L_1 + \rho_2 = \varepsilon + \rho_1 + \rho_2,\\
L_3 &= L_2 + \rho_3 = \varepsilon + \rho_1+\rho_2+\rho_3.
\end{aligned}
$$
Written out fully in terms of occupation numbers ($\phi_n^\dagger\phi_n$):
$$
\begin{aligned}
L_1 &= \varepsilon + \phi_1^\dagger\phi_1 - 1,\\
L_2 &= \varepsilon + (\phi_1^\dagger\phi_1 - 1) + \phi_2^\dagger\phi_2,\\
L_3 &= \varepsilon + (\phi_1^\dagger\phi_1-1) + \phi_2^\dagger\phi_2 + (\phi_3^\dagger\phi_3-1)
     = \varepsilon - 2 + \phi_1^\dagger\phi_1+\phi_2^\dagger\phi_2+\phi_3^\dagger\phi_3.
\end{aligned}
$$
Compare term-by-term against the paper's general formula, Eq. (5):
$$
L_n = \varepsilon + \sum_{l=1}^n\left(\phi_l^\dagger\phi_l - \frac{1-(-1)^l}{2}\right).
$$
Setting $n=1,2,3$ in Eq. (5) reproduces exactly the three lines above — **this is the check**: if your hand-derived $L_1, L_2, L_3$ don't match Eq. (5) termwise, you've made an arithmetic slip in the unrolling and must find it before continuing (a very common one: dropping the $-1$ from an odd site, or applying it to an even site).

### 3.4 Step 3 — substitute into the electric-field energy, expand, and see the long-range structure

The third term of Eq. (3) is $\dfrac{g^2a}{2}\sum_{n=1}^{N-1} L_n^2$. For $N=4$ this is $L_1^2 + L_2^2 + L_3^2$ (times the prefactor). Take $L_2 = \varepsilon+\rho_1+\rho_2$ as a representative term and expand:
$$
L_2^2 = \varepsilon^2 + 2\varepsilon(\rho_1+\rho_2) + (\rho_1+\rho_2)^2 = \varepsilon^2 + 2\varepsilon\rho_1 + 2\varepsilon\rho_2 + \rho_1^2 + \rho_2^2 + 2\rho_1\rho_2.
$$
**The cross-term $2\rho_1\rho_2$ is the key structural feature.** It says site 1 and site 2 talk to each other *directly* in the Hamiltonian, even though they are not adjacent (well, they are here, but see $L_3^2$ below, which pairs sites 1 and 3). Doing the same for $L_3^2$:
$$
L_3^2 = \varepsilon^2 + 2\varepsilon(\rho_1+\rho_2+\rho_3) + \rho_1^2+\rho_2^2+\rho_3^2 + 2\rho_1\rho_2 + 2\rho_1\rho_3 + 2\rho_2\rho_3.
$$
Now $\rho_1\rho_3$ appears — sites 1 and 3, separated by an intervening site, interact directly. **This is exactly why the electric-field term becomes a long-range interaction (Eq. 6) after elimination**: squaring a *cumulative* sum always produces a cross-term between every pair of sites up to that point, regardless of physical separation. In the un-eliminated theory (Eq. 3), this term was a completely local, single-link energy $\tfrac{g^2a}{2}L_n^2$; the nonlocality is not extra physics — it's the bookkeeping cost of having substituted away the gauge field in favor of matter variables. Summing $L_1^2+L_2^2+L_3^2$ and collecting terms reproduces the structure of Eq. (6): a sum over $n$ of a squared cumulative charge sum $\left[\varepsilon+\sum_{l=1}^n \rho_l\right]^2$.

### 3.5 Step 4 — redo after Jordan-Wigner, check against Eq. (11)

Jordan-Wigner (worked in full in Section 4) gives the occupation number in terms of a spin variable:
$$
\phi_n^\dagger\phi_n = \frac{1+\sigma_n^z}{2}.
$$
(Sanity check on this formula alone: if $\sigma^z_n=+1$ ["spin up"], occupation is $1$ — site occupied; if $\sigma^z_n=-1$, occupation is $0$ — empty. Good, matches the convention where $\sigma^z=+1$ means "occupied.")

Substitute into $\rho_n$:
$$
\rho_n = \frac{1+\sigma_n^z}{2} - \frac{1-(-1)^n}{2} = \frac{\sigma_n^z + (-1)^n}{2}.
$$
This is a clean, useful intermediate result: **the background-subtracted charge is just $(\sigma^z_n + (-1)^n)/2$.** For $N=4$:
$$
\rho_1 = \frac{\sigma_1^z - 1}{2},\quad \rho_2=\frac{\sigma_2^z+1}{2},\quad \rho_3=\frac{\sigma_3^z-1}{2},\quad \rho_4=\frac{\sigma_4^z+1}{2}.
$$
So:
$$
\begin{aligned}
L_1 &= \varepsilon + \frac{\sigma_1^z-1}{2},\\
L_2 &= \varepsilon + \frac{\sigma_1^z+\sigma_2^z}{2}\quad\ (\text{the } -1,+1 \text{ cancel}),\\
L_3 &= \varepsilon + \frac{\sigma_1^z+\sigma_2^z+\sigma_3^z - 1}{2}.
\end{aligned}
$$
Compare against the paper's Eq. (11) general form,
$$
L_n \;\to\; \varepsilon + \frac12\sum_{l=1}^n\left(\sigma_l^z + (-1)^l\right),
$$
which for $n=1,2,3$ reproduces exactly the three lines above (do the sum of $(-1)^l$ inside each: for $n=1$, $\sum(-1)^l=-1$; for $n=2$, $-1+1=0$; for $n=3$, $-1+1-1=-1$ — matching the constants that survived above). $H_E$ in Eq. (11) is then just $\tfrac{g^2a}{2}\sum_n L_n^2$ with $L_n$ written this way — the same object as Section 3.4, now in qubit variables.

### 3.6 The pass/fail check, and the known traps

**The actual test (per the roadmap): implement `H_E(N)` in code, and check its $N=4$ output against the by-hand result above, *before* comparing to exact diagonalization.** If you skip straight to ED, a self-consistent-looking but wrong Hamiltonian (e.g., with a systematically flipped sign convention) can still produce plausible-looking spectra, making the bug much harder to isolate later.

Explicit traps, all of which you now have the derivation needed to check directly against:

| Trap | What it looks like | How to catch it |
|---|---|---|
| Off-by-one in the cumulative sum | Using $\sum_{l=1}^{n-1}$ instead of $\sum_{l=1}^{n}$ (or vice versa) | Your $L_1$ should contain exactly $\rho_1$ — a sum with the wrong bound will give $L_1 = \varepsilon$ (empty sum) or start from $\rho_2$ |
| Mixed stagger | Using $(-1)^l$ in one term and $(-1)^{l+1}$ in another | Keep a single exponent. In the code, $(-1)^{n+1}$ is the zero-based index and does not change the physics |
| $\varepsilon$ placed inside vs. outside the sum before squaring | $\left[\sum_l(\ldots+\varepsilon)\right]^2$ instead of $\left[\varepsilon + \sum_l(\ldots)\right]^2$ | These differ by a factor of $n$ multiplying $\varepsilon$ inside the bracket — compare the $\varepsilon^2$ coefficient of each $L_n^2$ term: it must always be exactly $1\cdot\varepsilon^2$, never $n^2\varepsilon^2$ |
| Forgetting the $\tfrac12$ from $\sigma^z\to$ occupation number | Using $\phi_n^\dagger\phi_n = \sigma_n^z$ instead of $(1+\sigma_n^z)/2$ | Check the $N=4$, $\varepsilon=0$, all-spins-down state: it should have $L_n = -n/2$-ish structure (specific numbers per the table above), not integer-valued $L_n$ |

**Exit criterion for this section (as stated in the roadmap):** you should now be able to write `H_E(N, ε)` for general $N$, defend every term in it by pointing at the by-hand $N=4$ derivation above, and extend it to a time-dependent $\varepsilon(t)$ (which only requires promoting $\varepsilon\to\varepsilon(t)$ in Eq. 11 and re-Trotterizing at each step — the *structure* of $H_E$ is unaffected, since Gauss's law is imposed at each instant of time, not just once).

---

## 4. Jordan-Wigner mapping — [Intuition, mostly review]

### 4.1 The core problem it solves

Fermionic operators anticommute *across all sites*, not just at the same site: $\{\phi_n,\phi_m^\dagger\}=0$ for $n\neq m$ (this is what makes them fermions rather than independent bosonic qubits — swapping two fermions must produce a minus sign, a direct consequence of the Pauli principle for identical particles). Qubit operators $\sigma_n^\pm$ at *different* sites, by contrast, simply **commute** ($[\sigma_n^-,\sigma_m^+]=0$ for $n\ne m$, since they act on different tensor-product factors). If you naively set $\phi_n \to \sigma_n^-$, you'd get a qubit Hamiltonian with the wrong exchange statistics — a bug that has nothing to do with lattice indices or signs and everything to do with this global anticommutation structure.

### 4.2 The fix: a string of $\sigma^z$'s as a "sign counter"

The paper's Eq. (7):
$$
\phi_n = \left(\prod_{l=1}^{n-1} i\sigma_l^z\right)\frac{\sigma_n^x - i\sigma_n^y}{2} = \left(\prod_{l=1}^{n-1} i\sigma_l^z\right)\sigma_n^-.
$$
The **string** $\prod_{l<n} i\sigma_l^z$ acts as a bookkeeping device: it multiplies by $-1$ for every occupied site to the left of $n$. Since $\sigma^z_l = +1$ (occupied) or $-1$ (empty), moving a fermion operator "past" an occupied site picks up exactly the $-1$ that anticommutation demands, while moving past an empty site picks up nothing — precisely reproducing $\{\phi_n,\phi_m^\dagger\}=0$ for $n\ne m$ using only mutually-commuting qubit operators. You don't need to re-derive this identity if you already trust Jordan-Wigner from a quantum-computing background; the only new content in Eq. (7) relative to the generic JW formula is that it's being applied directly to the **staggered** fermion (so "site $n$" already encodes the particle/antiparticle structure from Section 2, with no separate spin or flavor index to track).

### 4.3 Worked check: the 2-site hopping term

Apply Eq. (7) to sites 1 and 2. Since the string for $n=1$ is empty (product over $l=1$ to $0$),
$$
\phi_1 = \sigma_1^-, \qquad \phi_2 = i\sigma_1^z\,\sigma_2^-.
$$
Then
$$
\phi_1^\dagger\phi_2 = \sigma_1^+ \cdot i\sigma_1^z\sigma_2^- = i\,(\sigma_1^+\sigma_1^z)\,\sigma_2^-.
$$
Using $\sigma^+\sigma^z = -\sigma^+$ (check directly: $\sigma^+ = \begin{pmatrix}0&1\\0&0\end{pmatrix}$, $\sigma^z=\begin{pmatrix}1&0\\0&-1\end{pmatrix}$, so $\sigma^+\sigma^z = \begin{pmatrix}0&-1\\0&0\end{pmatrix} = -\sigma^+$):
$$
\phi_1^\dagger\phi_2 = -i\,\sigma_1^+\sigma_2^-.
$$
Adding the Hermitian conjugate,
$$
\phi_1^\dagger\phi_2 + \phi_2^\dagger\phi_1 = -i\sigma_1^+\sigma_2^- + i\sigma_1^-\sigma_2^+.
$$
Now use $\sigma^+ = (\sigma^x+i\sigma^y)/2$, $\sigma^- = (\sigma^x-i\sigma^y)/2$. Substituting and expanding (straightforward but slightly tedious algebra — the $\sigma^x\sigma^y$ and $\sigma^y\sigma^x$ cross-terms cancel between the two pieces) gives
$$
\phi_1^\dagger\phi_2+\phi_2^\dagger\phi_1 = \tfrac12\left(\sigma_1^x\sigma_2^x + \sigma_1^y\sigma_2^y\right),
$$
which is exactly the structure of $H_{\rm kin}$ in Eq. (9) (up to the overall $\tfrac{1}{2a}$ or $\tfrac{1}{4a}$ prefactor bookkeeping, and note that in the actual $N$-site Hamiltonian the $i\sigma^z$ string is present but cancels between $\phi_n^\dagger$ and $\phi_{n+1}$ for **nearest-neighbor** hopping specifically, which is *why the kinetic term stays local* even though the general JW map is nonlocal — a fact worth internalizing: JW strings only cause trouble for *non-nearest-neighbor* fermion bilinears, which is exactly what makes $H_E$, not $H_{\rm kin}$, the nonlocal term).

**Pass/fail:** you should be able to state, without looking it up, why the kinetic term becomes $\sigma^x\sigma^x+\sigma^y\sigma^y$ (an "XY" or "hopping" interaction) rather than a single Pauli string: it's because a fermion hop is a *superposition* of "was at $n$, now at $n+1$" and its conjugate, and $\sigma^+\sigma^-+\sigma^-\sigma^+$ (raise-then-lower and lower-then-raise) is exactly $\tfrac12(\sigma^x\sigma^x+\sigma^y\sigma^y)$ by the identity above.

---

## 5. Sign problem vs. tensor-network entanglement growth — [Intuition]

These are two **unrelated** obstructions that happen to both be invoked to justify "why quantum simulation" in the paper's introduction — conflating them is a common and important error to avoid.

### 5.1 The Monte Carlo sign problem (why the *Euclidean/path-integral* approach struggles)

Euclidean path-integral Monte Carlo evaluates expectation values by treating $e^{-S}$ (with $S$ the Euclidean action) as a probability weight and importance-sampling field configurations accordingly. This requires $e^{-S}$ to be real and non-negative. When the theory has a chemical potential, a topological ($\theta$) term, or — most relevantly for *this* paper — you try to compute **real-time** dynamics (where the weight is $e^{iS}$, not $e^{-S}$: purely oscillatory, never a probability), the weight is complex or oscillatory. Monte Carlo then computes a ratio of two quantities, each built from wildly cancelling positive and negative (or complex) contributions; the *signal* (the true expectation value) is exponentially small compared to the *statistical noise* from the cancellations, and the number of samples needed to resolve it grows exponentially with system volume or evolution time. **This is a statistical sampling failure**, not a Hilbert-space-size failure — it is about the impossibility of importance-sampling a non-positive measure.

### 5.2 Tensor networks / MPS / DMRG (why the *classical Hamiltonian* approach struggles specifically for this paper's setup)

A common oversimplification (which the roadmap explicitly flags) is "tensor networks fail because of exponential Hilbert space scaling." This is **not** generally true. A matrix product state (MPS) represents a 1D quantum state with a bond dimension $\chi$, and the classical resources needed scale with $\chi$, not with the raw Hilbert space dimension $2^N$. **Ground states of gapped, local 1D Hamiltonians obey an area law**: their entanglement entropy across any cut saturates to a constant as the system grows, so $\chi$ stays *bounded*, and MPS/DMRG solve such ground states efficiently even for very large $N$. This is precisely why DMRG is a gold-standard *ground-state* method for 1D lattice gauge theories, including the Schwinger model itself (see the Bañuls et al. and Byrnes et al. references in the paper).

What actually strains tensor networks is **real-time evolution after a quench** — exactly the scenario in this paper (Section IV.B). Following a sudden change in the Hamiltonian (switching on $\varepsilon$), entanglement entropy across a cut typically **grows linearly in time** (this is a fairly general feature of non-equilibrium dynamics away from the ground state, sometimes summarized as "a global quench generates a ballistically-spreading light cone of entanglement"). Representing this growing entanglement faithfully with an MPS requires the bond dimension $\chi$ to grow correspondingly — in the worst case, exponentially in the evolution time $t$ — even though $N$ itself might be modest. **This is a computational cost tied to *time*, not to system size.**

### 5.3 The distinguishing question (the roadmap's pass/fail check)

*"Would tensor networks also struggle with just the ground state of this same Hamiltonian?"* **No.** The ground state of the (static, $\varepsilon$-scan) Hamiltonian in Section IV.A of the paper is exactly the well-behaved, area-law case that DMRG handles efficiently — indeed, the paper itself uses VQE/ED for that part precisely because the system sizes are small enough that all methods agree, and the physics point of the paper is specifically the *dynamics* after the quench (Section IV.B), which is the genuinely hard regime for classical tensor-network methods and the actual motivation for the digital-quantum-simulation approach.

| Obstruction | Cause | Scales badly with |
|---|---|---|
| Monte Carlo sign problem | Complex/oscillatory path-integral weight breaks importance sampling | System volume *and/or* real evolution time |
| Tensor-network (MPS/DMRG) difficulty | Entanglement growth after a quench requires growing bond dimension | Evolution time specifically (not ground-state system size) |

---

## 6. Putting it all together: reading Eqs. (1)–(17) end to end

You now have every ingredient needed to read Section II of the paper as a single coherent derivation rather than a sequence of disconnected equations. Here is the full chain:

1. **Eq. (1)–(2):** Start with the continuum Schwinger model — 1+1D QED — written with a covariant derivative $D_\mu = \partial_\mu+igA_\mu$ (Section 1: this specific combination exists *because* of local U(1) gauge invariance) and choose temporal gauge $A_0=0$ so the Hamiltonian has the clean form "kinetic + mass + field energy."

2. **Eq. (3):** Discretize on $N$ staggered sites (Section 2). Fermion bilinears live on sites, the gauge link variable $U_n$ and electric field $L_n$ live on links; the $(-1)^n$ factors encode which sites represent particle-type vs. antiparticle-type degrees of freedom.

3. **Eq. (4)–(6):** Because 1+1D gauge fields carry no independent propagating degrees of freedom (Section 3.1), Gauss's law completely determines $L_n$ from the matter charge distribution plus one boundary value $\varepsilon$ (Section 3.2–3.4). Substituting this solution back into the electric-field-energy term converts the theory into a **purely fermionic**, but now **long-range interacting**, Hamiltonian (Eq. 6) — the price paid for eliminating the gauge field is that the remaining matter Hamiltonian is no longer nearest-neighbor.

4. **Eq. (7)–(11):** Jordan-Wigner (Section 4) converts this fermionic Hamiltonian into a qubit Hamiltonian implementable on a digital quantum computer. The **kinetic term stays local** ($\sigma^x\sigma^x+\sigma^y\sigma^y$, nearest-neighbor only — Section 4.3) because JW strings cancel for adjacent-site hopping; the **electric-field term stays long-range** (Eq. 11, Section 3.5) because it was already long-range *before* JW, and JW doesn't fix that — it just re-expresses the same cumulative-charge structure in Pauli variables.

5. **Eq. (12)–(17):** The physical observables of interest — total charge (a Noether-conserved quantity from Section 1.1, and a genuine sanity check on your Trotter circuit: it must stay exactly constant), local/spatial charge density, the residual electric-field energy, vacuum-state fidelity, and the chiral condensate — are all built from expectation values of $\sigma_n^z$ operators (or short strings thereof) in the time-evolved qubit state. None of these require any new physics beyond what's above; they're bookkeeping definitions applied to the state produced by evolving under the Hamiltonian just constructed.

**The physical story the whole apparatus is built to tell** (Sections IV–V of the paper, for context — not part of the roadmap's mechanics requirement, but the point of building all of the above): the theory has a genuine vacuum instability. In a strong enough background field $\varepsilon$, it becomes energetically favorable for the "vacuum" (the ground state at $\varepsilon=0$) to rearrange itself into a state with real particle-antiparticle pairs sitting near the open boundaries, screening the imposed field — the finite-lattice, finite-time echo of the continuum **Schwinger effect** (vacuum decay via pair production in a strong field). The static field-scan (chiral condensate, Fig. 1) shows this as sharp, finite-size level crossings; the real-time quench (Figs. 3–7) shows it dynamically as boundary charge separation, oscillatory energy exchange between the fermionic sector and the electric-field term, and an exponentially-decaying overlap with the original zero-field vacuum — exactly the phenomenon that Monte Carlo cannot easily access (Section 5.1) and that tensor networks find expensive specifically because it is a real-time quench (Section 5.2), which is the physics case for the paper's digital Trotterized approach in the first place.

---

## Summary table: exit criteria met

| # | Topic | Tier | Key result you should now be able to reproduce unaided |
|---|---|---|---|
| 1 | Gauge invariance | Intuition | Explain, without equations, why local (not global) symmetry forces a compensating field to exist |
| 2 | Lattice discretization | Intuition | State $(-1)^n$ for any $n$ instantly; know matter lives on sites, gauge field on links |
| 3 | Gauss's law elimination | **Mechanics** | Derive $L_n$ for $N=4$ from the recursion, expand $L_n^2$ to see the long-range cross-terms, and re-derive Eq. (11) in $\sigma^z$ variables from scratch |
| 4 | Jordan-Wigner | Intuition | Derive $\phi_1^\dagger\phi_2+\text{h.c.} = \tfrac12(\sigma_1^x\sigma_2^x+\sigma_1^y\sigma_2^y)$ by hand |
| 5 | Sign problem vs. entanglement growth | Intuition | Correctly answer: would tensor networks struggle with just the ground state? (No — it's the quench-driven entanglement growth in time that's expensive, a separate obstruction from the Monte Carlo sign problem) |

Once you can do all five without notes, the physics prerequisite for building and debugging `H_kin`, `H_m`, and `H_E` — and for understanding *why* the resulting simulation is physically interesting — is complete.
