# Definitions for the retained driven experiments

This is the mathematical specification for the eight driven figures. The
[intuitive walkthrough](results.md) explains their purpose and results.
The frozen baseline keeps its original implementation and output paths.

## Model, units and preparation

Use open boundaries, even $N$, $\hbar=1$, and the baseline's staggered-site
convention. The driven builder reproduces the baseline matrix restricted to
an integer-charge sector, as checked independently in the unit tests:

\[
H(\epsilon)=H_0+\epsilon V+c\epsilon^2I,\qquad
V=g^2a\sum_{n=1}^{N-1}L_n(0),\qquad c=\frac{g^2a(N-1)}2,
\]

\[
q_l=\frac{Z_l+(-1)^{l+1}}2,\qquad
Q=\sum_lq_l,\qquad L_n(0)=\sum_{l\le n}q_l,\qquad l=1,\ldots,N.
\]

The stored `V` is the diagonal of this operator. Dynamics use $Q=0$ and
start from the exact ground state of $H_0$ in that sector. The static spectrum
also diagonalizes other charge sectors for comparison. The $N=8$ neutral
sector has 70 states; $N=12,16$ have dimensions 924 and 12870.

All retained datasets have numerical $a=m=g=1$, equivalently $ag=m/g=1$.
Figures label energy in units of $g$ and time in units of $1/g$; all field
parameters $\epsilon$ are dimensionless. Plot recipes verify these parameters
rather than silently interpreting changed units as the same experiment.

## Three pulse protocols

1. **Two Sauter pulses:** nominal shapes
   $\epsilon_s\operatorname{sech}^2(t/T_s)+
   \epsilon_w\operatorname{sech}^2((t-t_d)/\tau)$.
   Each component is retained until its relative tail reaches
   $\eta=10^{-4}$, then multiplied by a quintic taper over one component
   width. Both the field and its derivative vanish at its support endpoints.
   Controls and convergence checks use these actual compact waveforms.
2. **Weak sine:** $\epsilon_w\sin(\omega t+\phi)
   \sin^2(\pi t/T_w)$ on $[0,T_w]$, zero outside. The coefficient
   $\epsilon_w$ is the carrier amplitude; the windowed waveform need not
   reach that value exactly. Both frequency sidebands are included in response.
3. **Plateau plus sine:** a quintic ramp up, a constant plateau containing the
   weak sine and buffers, then a quintic ramp down. Ramp durations are 10/g
   and 20/g. Fidelity with the lowest neutral state at the plateau background
   is measured immediately after ramp-up, before the probe begins.

The exact pulse classes, derivatives and support boundaries live in
`src/driven/pulses.py`. These protocols start and end at zero field. A ground state
at a static background in Figure 1 is not assumed to have been prepared by
any particular ramp.

## Endpoint and post-pulse observables

At the last support endpoint $t_f$, define

\[
P_{\rm exc}=1-|\langle0|\psi(t_f)\rangle|^2,\qquad
W=\langle H_0\rangle_{t_f}-E_0,
\]

\[
\Delta P_{\rm exc}=P_{s+w}-P_s-P_w,\qquad
\Delta W=W_{s+w}-W_s-W_w.
\]

The no-drive trajectory starts in the same exact vacuum and has zero
excitation up to numerical error. It is retained as a preparation/solver
control. We never renormalize states to hide drift. To reduce cancellation,
the code evaluates $P_{\rm exc}$ as the squared vacuum-orthogonal norm and
$W$ as $\langle\psi|(H_0-E_0)|\psi\rangle$.

The local occupation diagnostic is $o_l=q_l^2$. Each trajectory's recorded
`occupation` is $\langle o_l\rangle-\langle o_l\rangle_0$. After the field is
off it is averaged for $10/g$, using trapezoidal integration of the saved
output samples. For region $R$ the size figure displays

\[
\Delta\overline{o}_R=\frac1{|R|}\sum_{l\in R}
\big[\overline{\langle o_l\rangle}_{s+w}
-\overline{\langle o_l\rangle}_{s}
-\overline{\langle o_l\rangle}_{w}+\langle o_l\rangle_0\big].
\]

The vacuum term is added once because each of the three stored profiles was
already vacuum-subtracted. The central region is the middle four sites; the
edge region is the first two and last two. This quantity can be signed and
is not an asymptotic particle or pair count. No independent spatial error
bars were estimated for the displayed regional means.

For Figure 5, $p_n=|\langle n;0|\psi(t_f)\rangle|^2$ and
$\Delta p_n=p_n^{s+w}-p_n^s-p_n^w$. All 69 excited states are included;
$\sum_{n>0}\Delta p_n=\Delta P_{\rm exc}$ is checked when rendering.
Individual populations use the saved zero-field eigenbasis; they do not
specify transition histories. The plotted horizontal coordinate is
$(E_n-E_0)/g$.

## Response calculations and controls

The zero-background leading weak response uses the Fourier transform of the
implemented waveform:

\[
P_{\rm LR}=\sum_{n>0}|\langle n|V|0\rangle|^2
\left|\int\epsilon(t)e^{i(E_n-E_0)t}\,dt\right|^2.
\]

Small sectors use a full eigensystem. For the large-sector response checks,
Lanczos quadrature starts from $(V-\langle V\rangle)|0\rangle$, with full
reorthogonalization and a repeated calculation at twice the step count.
These quadrature nodes are not individually identified physical eigenstates.

The strong-pulse calculation expands the state as
$\psi(\lambda)=\psi_0+\lambda\psi_1+\lambda^2\psi_2+O(\lambda^3)$ about the
actual strong trajectory. The second-order observable includes interference
between $\psi_0$ and $\psi_2$. Its plotted comparison is the added response
$\delta P=P_{s+w}-P_s$, which differs from $\Delta P_{\rm exc}$ by $P_w$.
For weak polarity $\pm$,

\[
\delta P_{\rm even}=(P_++P_-)/2-P_s,\qquad
\delta P_{\rm odd}=(P_+-P_-)/2.
\]

Their sum equals $P_+-P_s$. In the fixed-budget control the actual compact
weak waveform is rescaled to hold $\int\epsilon_w(t)^2dt$ at its value for
$g\tau=1$, subject to the amplitude ceiling 0.04. This is an input-budget
convention, not a measurement of absorbed work. The equal-peak slow control
has strong width and amplitude $\epsilon_s+\epsilon_w$.

## Numerical checks and provenance

The primary method is chronological DOP853 with nominal `rtol=1e-10`,
`atol=1e-12`, `max_step=0.05/g`. Known pulse boundaries split the integration;
active short features impose smaller step caps. Output sampling and accepted
solver steps are stored separately. A scalar $c\epsilon^2I$ is omitted from
propagation because it changes only a global phase. Physical energy and
integrated source work include it:

\[
W_{\rm src}(t)=\int_{t_i}^{t}\dot\epsilon(u)
[\langle V\rangle_u+2c\epsilon(u)\langle\psi|\psi\rangle_u],du.
\]

We check $E(t)-E(t_i)-W_{\rm src}(t)$, norm drift, and post-pulse constancy
of $P_{\rm exc}$ and $W$. $Q=0$ is enforced by the basis, so it is not an
independent diagnostic of leakage. The builder is separately checked against
the unrestricted baseline Hamiltonian at small sizes.

Main runs repeat with half the step and tighter tolerances. The numerical
uncertainty estimate for each of `Y` and `W` is the absolute refinement
change with a nominal floor $10^{-10}$, augmented by the available tail and
independent-method differences. Constituent estimates are added for signed
nonadditivity. The operational excitation-resolution rule is
$|\Delta P_{\rm exc}|>\max(10^{-8},5u_{\Delta P})$. These are numerical
acceptance rules, not statistical intervals or proven bounds. Figure 8 uses
midpoint steps 0.02/g and 0.01/g and tail thresholds $10^{-5},10^{-6}$.

Original raw points, checksums, metadata and actual source archives remain
unchanged for the nine retained studies. Current recipes omit the unused
readout-duration analysis; original settings remain in historical archives.
Use a fresh suffix for reruns after source changes. See the
[running guide](guide.md) for the complete regeneration path.
