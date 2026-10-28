"""Read the quantities used by the eight driven figures and their solver checks.

The observer is called by time_evolution; summarize turns its time series into
one endpoint record and the post-pulse spatial average used in the size plot.
"""
import numpy as np
from scipy.integrate import trapezoid


def make_observer(ops, vacuum, pulse):
    """Return (time, state) -> observables relative to the zero-field vacuum.

    Y is the squared norm outside the initial vacuum, not a particle count.
    W is the excess zero-field energy. Both use cancellation-resistant forms.
    Physical instantaneous energy includes the scalar field term omitted from
    propagation; this makes the independent source-work check meaningful.
    """
    occupation0 = abs(vacuum)**2 @ ops.occupation
    charge_squared = ops.q.sum(axis=1)**2
    E0 = float(np.vdot(vacuum, ops.H0@vacuum).real)

    def observe(t, psi):
        probability = abs(psi)**2
        norm = float(probability.sum())
        eps = float(pulse.epsilon(t))
        orthogonal = psi-vacuum*np.vdot(vacuum, psi)
        energy0 = np.vdot(psi, ops.H0@psi).real
        return {
            "Y": float(np.vdot(orthogonal, orthogonal).real),
            "W": float(np.vdot(psi, ops.H0@psi-E0*psi).real),
            "energy": float(energy0+eps*(probability@ops.V)+ops.c*eps**2*norm),
            "occupation": probability@ops.occupation-occupation0,
            "norm": norm, "Q2": float(probability@charge_squared),
        }
    return observe


def summarize(trajectory, pulse_end, N, post_duration=10.):
    """Endpoint results, post-pulse occupation average and conservation checks.

    Local occupation is vacuum-subtracted and may be negative. Its average is
    a spatial diagnostic; it does not define a count of produced pairs.
    """
    times, values = trajectory.times, trajectory.observations
    endpoint = int(np.argmin(abs(times-pulse_end)))
    if abs(times[endpoint]-pulse_end) > 1e-9:
        raise ValueError("Sampling must include pulse_end exactly")
    summary = {name: np.asarray(value[endpoint]).tolist() for name, value in values.items()}
    post = (times >= pulse_end) & (times <= pulse_end+post_duration+1e-10)
    if post_duration > 0 and post.sum() >= 2:
        summary["occupation_average"] = (trapezoid(values["occupation"][post], times[post], axis=0)
                                         /(times[post][-1]-times[post][0])).tolist()
    summary["W_per_site"] = summary["W"]/N
    summary["max_norm_error"] = float(np.max(abs(values["norm"]-1)))
    summary["max_charge_squared"] = float(np.max(abs(values["Q2"])))
    summary["max_work_residual"] = float(np.max(abs(values["energy"]-values["energy"][0]-trajectory.work)))
    summary["post_Y_drift"] = float(np.ptp(values["Y"][post]))
    summary["post_W_drift"] = float(np.ptp(values["W"][post]))
    return summary
