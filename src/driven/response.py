"""Finite-waveform weak response and the expansion about a strong pump."""
import numpy as np
from scipy.integrate import simpson, solve_ivp


def fourier(pulse, frequencies, spacing=None):
    """Integral epsilon(t) exp(i Omega t) dt of the implemented compact pulse.

    Both sine sidebands and every taper are included. Refine spacing to test
    quadrature independently of propagation; process frequency chunks to keep
    broadband Lanczos calculations small.
    """
    frequencies = np.atleast_1d(frequencies)
    start, end = pulse.interval
    spacing = min(min(pulse.features)/30, .02, 1/(30*max(1., np.max(abs(frequencies), initial=0)))) if spacing is None else spacing
    count = int(np.ceil((end-start)/spacing))
    count += count % 2
    times = np.linspace(start, end, count+1)
    signal = pulse.epsilon(times)
    result = np.empty(len(frequencies), dtype=complex)
    for i in range(0, len(frequencies), 32):
        result[i:i+32] = simpson(signal[:, None]*np.exp(1j*times[:, None]*frequencies[None, i:i+32]),
                                  x=times, axis=0)
    return result


def linear_response(spectrum, pulse):
    populations = spectrum["weights"]*abs(fourier(pulse, spectrum["gaps"]))**2
    return {"Y": float(populations.sum()), "W": float(spectrum["gaps"]@populations),
            "populations": populations}


def pump_expansion(ops, vacuum, strong, unit_weak, interval=None, rtol=1e-10,
                   atol=1e-12, max_step=.025):
    """Order-0,1,2 auxiliary states, with no factorial in their definition.

    psi(lambda)=psi0+lambda psi1+lambda**2 psi2. The psi2 interference term
    is necessary when the pump already excites a coherent superposition.
    The energy origin is shifted by the initial vacuum energy for efficiency.
    """
    interval = (min(strong.interval[0], unit_weak.interval[0]),
                max(strong.interval[1], unit_weak.interval[1])) if interval is None else interval
    E0 = float(np.vdot(vacuum, ops.H0@vacuum).real)
    def rhs(t, flat):
        state = flat.reshape(3, -1)
        derivative = np.asarray([ops.H0@v+(strong.epsilon(t)*ops.V-E0)*v for v in state])
        derivative[1:] += unit_weak.epsilon(t)*ops.V*state[:-1]
        return (-1j*derivative).ravel()
    result = solve_ivp(rhs, interval, np.r_[vacuum, np.zeros(2*len(vacuum), complex)],
                       method="DOP853", rtol=rtol, atol=atol, max_step=max_step)
    if not result.success:
        raise RuntimeError(result.message)
    states = result.y[:, -1].reshape(3, -1)
    overlaps = states@vacuum.conj()
    survival_coefficients = np.array([1-abs(overlaps[0])**2,
        -2*np.real(overlaps[0].conj()*overlaps[1]),
        -abs(overlaps[1])**2-2*np.real(overlaps[0].conj()*overlaps[2])])
    # Equivalent orthogonal formulation avoids subtractive cancellation at Y~0.
    perpendicular = states-overlaps[:, None]*vacuum
    coefficients = np.array([np.vdot(perpendicular[0], perpendicular[0]).real,
        2*np.vdot(perpendicular[0], perpendicular[1]).real,
        np.vdot(perpendicular[1], perpendicular[1]).real+2*np.vdot(perpendicular[0], perpendicular[2]).real])
    return {"coefficients": coefficients, "survival_coefficients": survival_coefficients,
            "states": states, "integration_times": result.t}


def squared_drive(pulse):
    """Input integral epsilon**2 dt; this is not absorbed work."""
    a, b = pulse.interval
    times = np.linspace(a, b, max(1001, int((b-a)/min(pulse.features)*50)|1))
    return float(simpson(pulse.epsilon(times)**2, x=times))
