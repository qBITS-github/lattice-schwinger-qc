"""Independent checks against the original Kronecker-product Hamiltonian."""
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"src"))
from driven.model import build_sector_operators, neutral_basis, sector_ground_state, transition_spectrum, lanczos_spectrum
from baseline.model import build_hamiltonian


class ModelTests(unittest.TestCase):
    def test_restriction_all_small_sectors_and_fields(self):
        for N in (2, 4, 6, 8):
            for Q in (-1, 0, 1):
                ops = build_sector_operators(N, a=.8, m=.7, g=1.2, Q=Q)
                np.testing.assert_allclose(ops.q.sum(axis=1), Q)
                for eps in (-.3, 0., .7):
                    H = build_hamiltonian(N, .8, .7, 1.2, eps)[ops.basis][:, ops.basis].toarray()
                    result = ops.hamiltonian(eps).toarray()
                    self.assertLess(np.linalg.norm(H-result)/max(1, np.linalg.norm(H)), 1e-12)
                    np.testing.assert_allclose(result, result.T)

    def test_dimensions_and_bit_order(self):
        for N, dim in ((8, 70), (12, 924), (16, 12870)):
            basis, inverse = neutral_basis(N)
            self.assertEqual(len(basis), dim)
            self.assertEqual(inverse[int(basis[-1])], dim-1)
        ops = build_sector_operators(4)
        i = ops.basis_index[0b1100]
        np.testing.assert_array_equal(ops.q[i], [0, -1, 1, 0])

    def test_all_n8_regression_anchors(self):
        ops = build_sector_operators(8)
        spectrum = transition_spectrum(ops)
        bright = spectrum["weights"] > 1e-8
        self.assertAlmostEqual(spectrum["energies"][0], -4.6380577409, places=9)
        self.assertAlmostEqual(spectrum["gaps"][bright][0], 2.6579940584, places=9)
        self.assertAlmostEqual(spectrum["weights"][bright][0], .2165681749, places=9)
        negative = build_sector_operators(8, Q=-1)
        crossing = brentq(lambda e: sector_ground_state(ops, e)[0]-sector_ground_state(negative, e)[0], .6, .8)
        self.assertAlmostEqual(crossing, .6923761274, places=9)
        at_crossing = transition_spectrum(ops, crossing)
        self.assertAlmostEqual(at_crossing["gaps"][0], 1.3512648343, places=9)
        self.assertAlmostEqual(spectrum["captured_weight"], spectrum["variance"], places=12)

    def test_lanczos_moments(self):
        ops = build_sector_operators(8)
        exact = transition_spectrum(ops)
        measure = lanczos_spectrum(ops, steps=40)
        for power in range(5):
            a = exact["weights"]@(exact["gaps"]**power)
            b = measure["weights"]@(measure["gaps"]**power)
            self.assertLess(abs(a-b)/max(1, abs(a)), 1e-10)
