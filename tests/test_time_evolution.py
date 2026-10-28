"""Time ordering, pulse calculus, work, weak response and pump expansion."""
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from driven.model import build_sector_operators, sector_ground_state, transition_spectrum
from driven.pulses import Constant, Sauter, Sine, Plateau, Sum
from driven.evolution import propagate, sampling_grid
from driven.observables import make_observer, summarize
from driven.response import linear_response, pump_expansion
from baseline.model import evolve


class EvolutionTests(unittest.TestCase):
    def setUp(self):
        self.ops = build_sector_operators(4)
        self.E0, self.vacuum = sector_ground_state(self.ops)

    def test_endpoints_and_derivatives(self):
        for pulse in (Sauter(.2, .4, .1), Sine(.1, 2.3, 3.), Plateau(.8, .4, .8)):
            for t in pulse.interval:
                self.assertEqual(float(pulse.epsilon(t)), 0.)
                self.assertEqual(float(pulse.epsilon_dot(t)), 0.)
            times = np.linspace(*pulse.interval, 61)[1:-1]
            np.testing.assert_allclose([pulse.epsilon(float(t)) for t in times],pulse.epsilon(times),atol=1e-14)
            np.testing.assert_allclose([pulse.epsilon_dot(float(t)) for t in times],pulse.epsilon_dot(times),atol=1e-14)
            finite = (pulse.epsilon(times+1e-6)-pulse.epsilon(times-1e-6))/2e-6
            np.testing.assert_allclose(pulse.epsilon_dot(times), finite, atol=1e-8)

    def test_short_pulse_is_resolved_between_output_times(self):
        pulse=Sine(2.,1000.,.002,start=.499)
        times=np.array([0.,1.])
        coarse=propagate(self.ops,pulse,self.vacuum,times,max_step=.05,energy_shift=self.E0)
        fine=propagate(self.ops,pulse,self.vacuum,times,max_step=.001,energy_shift=self.E0)
        np.testing.assert_allclose(coarse.final_state,fine.final_state,atol=1e-10)
        self.assertGreater(np.sum((coarse.integration_times>.499)&(coarse.integration_times<.501)),10)

    def test_zero_drive_and_old_constant_evolution(self):
        times = np.linspace(0, 2, 21)
        for eps in (0., .3):
            pulse = Constant(eps, 0, 2)
            result = propagate(self.ops, pulse, self.vacuum, times)
            reference = evolve(self.vacuum, self.ops.hamiltonian(eps), times)
            # The omitted eps**2 identity causes a known global phase.
            reference *= np.exp(1j*self.ops.c*eps**2*times[:, None])
            np.testing.assert_allclose(result.states, reference, atol=2e-10)

    def test_midpoint_convergence_work_and_stationarity(self):
        pulse = Sine(.1, 2.5, 3.)
        times = sampling_grid(0, 4, .1, [3.])
        observer = make_observer(self.ops, self.vacuum, pulse)
        adaptive = propagate(self.ops, pulse, self.vacuum, times, observer=observer)
        errors = []
        for h in (.04, .02, .01):
            result = propagate(self.ops, pulse, self.vacuum, times, method="midpoint", max_step=h)
            errors.append(np.linalg.norm(result.final_state-adaptive.final_state))
        self.assertLess(errors[-1], errors[0]/8)
        check = summarize(adaptive, 3., 4, 1.)
        self.assertLess(check["max_norm_error"], 1e-9)
        self.assertLess(check["max_work_residual"], 1e-9)
        self.assertLess(check["post_Y_drift"], 1e-9)
        self.assertLess(check["post_W_drift"], 1e-9)

    def test_weak_response_and_pump_second_order(self):
        spectrum = transition_spectrum(self.ops)
        pulse = Sine(.001, 2.5, 3.)
        times = np.array([0., 3.])
        observer = make_observer(self.ops, self.vacuum, pulse)
        run = propagate(self.ops, pulse, self.vacuum, times, observer=observer)
        prediction = linear_response(spectrum, pulse)
        self.assertLess(abs(run.observations["Y"][-1]/prediction["Y"]-1), .002)
        strong, weak = Sine(.15, 2., 3.), Sine(1., 2.5, 3.)
        expansion = pump_expansion(self.ops, self.vacuum, strong, weak)
        errors = []
        for amplitude in (.04, .02, .01):
            p = Sum((strong, Sine(amplitude, 2.5, 3.)))
            r = propagate(self.ops, p, self.vacuum, times,
                          observer=make_observer(self.ops, self.vacuum, p))
            predicted = np.polynomial.polynomial.polyval(amplitude, expansion["coefficients"])
            errors.append(abs(r.observations["Y"][-1]-predicted))
        self.assertLess(errors[-1], errors[0]/30)
        np.testing.assert_allclose(expansion["coefficients"], expansion["survival_coefficients"], atol=1e-8)
