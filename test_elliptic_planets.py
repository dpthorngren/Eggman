import numpy as np
from scipy.integrate import dblquad
import shapes
import unittest


def _integrateDebugTriangle(stack, params):
    '''Helper function that does a subset of triangle integrals for checking against IntegrateTriangle'''
    x0, y0, _, x1, y1, _, x2, y2, _ = stack[0]
    assert x0 == x1

    def top(x):
        return (x-x1)*(y2-y1)/(x2-x1) + y1
    
    def bot(x):
        return (x-x0)*(y2-y0)/(x2-x0) + y0

    def integrand(y, x):
        return shapes.integrand(x, y, params)
    
    # TODO: This function works but it leaks memory that sometimes causes a segfault on program close >.<
    return dblquad(integrand, x0, x2, bot, top)[0]


class TestTriangleIntegration(unittest.TestCase):
    def test_right(self):
        stack = np.zeros((12,9))

        # Flat function across the unit right triangle.
        params = np.array([1., 0., 0., 0., 0.])
        stack[0] = 0, 0, 0,   0, 1, 0,   1, 0, 0
        ref = _integrateDebugTriangle(stack, params)
        actual = shapes.integrateTriangle(stack, params, initEval=True)
        self.assertAlmostEqual(ref, actual, 6)

        # Quadratic function across an isosceles triangle
        params = np.array([-1., 5., 5., 0., 0.])
        stack[0] = 0, -1.5, 0, 0, 1.5, 0, 2, 1, 0
        ref = _integrateDebugTriangle(stack, params)
        actual = shapes.integrateTriangle(stack, params, initEval=True)
        self.assertAlmostEqual(ref, actual, 6)

        # Sinusoidal function across an irregular triangle
        params = np.array([0., 0., 0., 1., 1.5])
        stack[0] = .1, -1.5, 0, .1, 1., 0, 2.1, 1, 0
        ref = _integrateDebugTriangle(stack, params)
        actual = shapes.integrateTriangle(stack, params, initEval=True)
        self.assertAlmostEqual(ref, actual, 6)

        # Known function across an isosceles triangle
        params = np.array([1.5, 0., 0., 1., 0.])
        stack[0] = -.5, -2, 0, 0, 1., 0., 0.5, -2, 0.
        actual = shapes.integrateTriangle(stack, params, initEval=True)
        self.assertAlmostEqual(1.5*1.5, actual, 6)


if __name__ == "__main__":
    unittest.main()
