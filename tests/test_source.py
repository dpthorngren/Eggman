import numpy as np
from pytest import approx

from eggman.cy_eggman import get_brightness


def test_quadratic_darkening():
    # Exactly zero outside bounds
    assert get_brightness(0, [.1, .1], .9, .9) == 0.0
    assert get_brightness(0, [.1, .1], -.9, .85) == 0.0
    # Simple cases
    norm1 = (1. - .5/3 - .2/6) * np.pi
    for x, y in 0.5 * np.random.rand(100, 2):
        assert get_brightness(0., [0., 0.], x, y) == approx(1.0 / np.pi, 1e-12)
        nu = 1.0 - np.sqrt(1 - x*x - y*y)
        expect = (1.0 - .5*nu - .2*nu*nu) / norm1
        assert get_brightness(0., [0.5, 0.2], x, y) == approx(expect, 1e-12)


def test_nonlinear_darkening():
    # Exactly zero outside bounds
    assert get_brightness(1, [.5, .1, .05, .05], .9, .9) == 0.0
    assert get_brightness(1, [.5, .1, .05, .05], -.9, .85) == 0.0
    # Simple cases
    norm = (1 - .5/5. - .2/3. - 3.*.1/7. - .1/2.) * np.pi
    for x, y in 0.5 * np.random.rand(100, 2):
        assert get_brightness(1, [0., 0., 0., 0.], x, y) == approx(1.0 / np.pi, 1e-12)
        nu = np.sqrt(np.sqrt(1 - x*x - y*y))
        expect = (1. - .5 * (1.-nu) - .2 * (1. - nu**2) - .1 * (1. - nu**3) - .1 * (1. - nu**4)) / norm
        assert get_brightness(1., [0.5, 0.2, .1, .1], x, y) == approx(expect, 1e-12)
