import numpy as np
from pytest import approx

import eggman


def test_quadratic_darkening():
    # Exactly zero outside bounds
    source = eggman.LightSourceWrap('uniform', [1.], 'quadratic', [.1, .1])
    assert source.get_brightness_sphere(.9, .9) == 0.0
    assert source.get_brightness_sphere(-.9, .85) == 0.0
    # Simple cases
    norm1 = (1. - .5/3 - .2/6) * np.pi
    source1 = eggman.LightSourceWrap('uniform', [1.], 'quadratic', [0, 0])
    source2 = eggman.LightSourceWrap('uniform', [1.], 'quadratic', [.5, .2])
    for x, y in 0.5 * np.random.rand(100, 2):
        assert source1.get_brightness_sphere(x, y) == approx(1.0 / np.pi, 1e-12)
        nu = 1.0 - np.sqrt(1 - x*x - y*y)
        expect = (1.0 - .5*nu - .2*nu*nu) / norm1
        assert source2.get_brightness_sphere(x, y) == approx(expect, 1e-12)


def test_nonlinear_darkening():
    # Exactly zero outside bounds
    source = eggman.LightSourceWrap('uniform', [1.], 'nonlinear', [.5, .1, .05, .05])
    assert source.get_brightness_sphere(.9, .9) == 0.0
    assert source.get_brightness_sphere(-.9, .85) == 0.0
    # Simple cases
    norm = (1 - .5/5. - .2/3. - 3.*.1/7. - .1/2.) * np.pi
    source1 = eggman.LightSourceWrap('uniform', [1.], 'nonlinear', [0., 0., 0., 0.])
    source2 = eggman.LightSourceWrap('uniform', [1.], 'nonlinear', [.5, .2, .1, .1])
    for x, y in 0.5 * np.random.rand(100, 2):
        assert source1.get_brightness_sphere(x, y) == approx(1.0 / np.pi, 1e-12)
        nu = np.sqrt(np.sqrt(1 - x*x - y*y))
        expect = (1. - .5 * (1.-nu) - .2 * (1. - nu**2) - .1 * (1. - nu**3) - .1 * (1. - nu**4))
        assert source2.get_brightness_sphere(x, y) == approx(expect / norm, 1e-12)
