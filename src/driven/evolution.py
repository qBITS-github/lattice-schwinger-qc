"""Chronological propagation with independent adaptive and midpoint methods.

The scalar c*eps**2 is omitted from state propagation (global phase only).
Physical source work still includes it. Adaptive steps and output sampling
are separate. An observer lets large-sector runs stream observables without
retaining every state. No state is silently renormalized.
"""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.integrate import DOP853
from scipy.sparse.linalg import expm_multiply


@dataclass
class Trajectory:
    times: np.ndarray
    states: np.ndarray | None
    final_state: np.ndarray
    integration_times: np.ndarray
    work: np.ndarray
    observations: dict
    settings: dict


def sampling_grid(start, end, spacing, boundaries=()):
    """At most spacing apart, with exact pulse/readout boundaries inserted."""
    if end <= start or spacing <= 0:
        raise ValueError("Require increasing endpoints and positive spacing")
    grid = np.linspace(start, end, int(np.ceil((end-start)/spacing))+1)
    return np.unique(np.r_[grid, [x for x in boundaries if start <= x <= end]])


def step_cap(pulse, g=1.):
    """Initial pulse-resolution cap; convergence tests are still required."""
    return min([.05/g, *[x/10 for x in pulse.features]])


def propagate(ops, pulse, psi0, times, *, method="DOP853", rtol=1e-10,
              atol=1e-12, max_step=None, energy_shift=0., observer=None,
              store_states=True, resolution_scale=1.):
    """Evolve from times[0], accumulating the external source work.

    method='midpoint' uses two exponential half-steps for an independent
    chronological reference and Simpson work quadrature.
    An observer is called only at requested output times.
    """
    times = np.asarray(times, dtype=float)
    psi0 = np.asarray(psi0, dtype=complex)
    if (times.ndim != 1 or len(times) < 2 or not np.all(np.isfinite(times))
            or not np.all(np.diff(times) > 0) or psi0.shape != (len(ops.basis),)):
        raise ValueError("Invalid time grid or state dimension")
    if method not in ("DOP853", "midpoint"):
        raise ValueError("method must be DOP853 or midpoint")
    max_step = (.05 if method == "DOP853" else .02)/ops.g if max_step is None else max_step
    if max_step <= 0 or resolution_scale <= 0:
        raise ValueError("max_step must be positive")
    stored, work, observations = [], [], []

    def power(t, psi):
        return float(pulse.epsilon_dot(t)*(np.vdot(psi, ops.V*psi).real
                     +2*ops.c*pulse.epsilon(t)*np.vdot(psi, psi).real))

    def sample(t, psi, source_work):
        if store_states:
            stored.append(psi.copy())
        work.append(float(np.real(source_work)))
        if observer is not None:
            observations.append(observer(t, psi))

    sample(times[0], psi0, 0.)
    integration_times = [times[0]]
    if method == "DOP853":
        def rhs(t, state):
            psi = state[:-1]
            dpsi = -1j*(ops.H0@psi+(pulse.epsilon(t)*ops.V-energy_shift)*psi)
            return np.r_[dpsi, power(t, psi)]
        # Known support boundaries prevent stepping across an entire short
        # pulse. Only active components impose their fast-feature step cap.
        components = getattr(pulse,"components",(pulse,))
        boundaries = sorted(set([times[0],times[-1],
            *[t for t in pulse.boundaries if times[0]<t<times[-1]]]))
        state = np.r_[psi0,0j]
        next_output = 1
        segments=[]
        for left,right in zip(boundaries[:-1],boundaries[1:]):
            midpoint=(left+right)/2
            active=[p for p in components if p.interval[0]<midpoint<p.interval[1]]
            cap=min([max_step,*[resolution_scale*step_cap(p,ops.g) for p in active]])
            segments.append({"start":left,"end":right,"max_step":cap})
            solver = DOP853(rhs,left,state,right,rtol=rtol,atol=atol,max_step=cap)
            while solver.status == "running":
                solver.step()
                if solver.status == "failed": raise RuntimeError("DOP853 failed to advance")
                integration_times.append(solver.t)
                interpolant=solver.dense_output()
                while next_output<len(times) and times[next_output]<=solver.t:
                    y=interpolant(times[next_output])
                    sample(times[next_output],y[:-1],y[-1]); next_output+=1
            state=solver.y
        psi = state[:-1]
    else:
        psi = psi0.copy()
        t, source_work = times[0], 0.
        for output_time in times[1:]:
            for target in np.linspace(t, output_time, int(np.ceil((output_time-t)/max_step))+1)[1:]:
                h = target-t
                eps = float(pulse.epsilon(t+h/2))
                H = ops.H0+sparse.diags(eps*ops.V-energy_shift)
                A = -1j*h*H/2
                half = expm_multiply(A, psi, traceA=A.diagonal().sum())
                following = expm_multiply(A, half, traceA=A.diagonal().sum())
                source_work += h/6*(power(t, psi)+4*power(t+h/2, half)+power(target, following))
                psi = following; t = target
                integration_times.append(t)
            sample(output_time, psi, source_work)
    result = {key: np.asarray([row[key] for row in observations]) for key in observations[0]} if observations else {}
    return Trajectory(times, np.asarray(stored) if store_states else None, psi.copy(),
                      np.asarray(integration_times), np.asarray(work), result,
                      {"method": method, "rtol": rtol, "atol": atol,
                       "max_step": max_step, "energy_shift": energy_shift,
                       "resolution_scale": resolution_scale,
                       "scalar_phase_omitted": True,
                       "segments": segments if method == "DOP853" else None})
