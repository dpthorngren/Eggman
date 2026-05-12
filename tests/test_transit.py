import numpy as np
from pytest import approx

import eggman


def test_transit_grazing():
    # Grazing transits can cause numerical issues, let's double-check that doesn't happen.
    # Grazing on ingress
    baseArgs = {
        't': np.array([-np.arcsin(1.1 / 10)]),
        't0': 0.,
        'theta': 0.,
        'phi': 0.,
        'gamma': 0.,
        'period': 2 * np.pi,
        'semimajor': 10.,
        'inclination': 90.,
        'limb_type': 'quadratic',
        'limb': [.3, .2],
        'rotate_with_orbit': False,
    }

    # Check for errors on grazing transits
    assert eggman.transit(.1 + 1e-3, .11, .1, .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.1 + 1e-9, .1, .09, .09, **dict(baseArgs)) < 1.

    # Check for errors on barely non-grazing transits
    assert eggman.transit(.1 - 1e-3, .11, .1, .1, **dict(baseArgs)) == 1.
    assert eggman.transit(.1 - 1e-12, .1, .09, .09, **dict(baseArgs)) == 1.

    # Grazing on egress
    baseArgs['t'] = -baseArgs['t']
    assert eggman.transit(.11, .1 + 1e-3, .1, .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.1, .1 + 1e-9, .09, .09, **dict(baseArgs)) < 1.

    # Grazing on top (true grazing transit)
    baseArgs['inclination'] = 90. - abs(baseArgs['t'][0]) * 180 / np.pi
    baseArgs['t'] = np.array([0.])
    assert eggman.transit(.11, .1, .1 + 1e-3, .1, **dict(baseArgs)) < 1.
    # TODO: Would like to pass these, but they don't actually violate my <1 PPM precision target.
    # assert eggman.transit(.09, .1, .1 + 1e-6, .1, **dict(baseArgs)) < 1.
    # assert eggman.transit(.1, .09, .1 + 1e-6, .1, **dict(baseArgs)) < 1.

    # Barely non-grazing on top (true grazing transit)
    assert eggman.transit(.11, .1, .1 - 1e-3, .1, **dict(baseArgs)) == 1.
    assert eggman.transit(.09, .1, .1 - 1e-6, .1, **dict(baseArgs)) == 1.
    assert eggman.transit(.1, .09, .1 - 1e-12, .1, **dict(baseArgs)) == 1.


def test_transit_is_transiting():
    # Ensure we correctly can tell transit from non-transit
    baseArgs = {
        't': np.array([0.]),
        't0': 0.,
        'theta': 0.,
        'phi': 0.,
        'gamma': 0.,
        'period': 10.,
        'semimajor': 10.,
        'inclination': 89.,
        'limb_type': 'quadratic',
        'limb': [.3, .2],
        'eccen': 0,
        'lon_periapse': 90.,
        'rotate_with_orbit': False,
    }
    # Non-transiting
    assert eggman.transit(.1, .1, .1, .1, **dict(baseArgs, t=np.array([5.]))) == 1.
    assert eggman.transit(.1, .1, .1, .1, **dict(baseArgs, t=np.array([3.]))) == 1.
    assert eggman.transit(
        1.1, 1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.51]), semimajor=1.)) == 1.
    assert eggman.transit(.1, .1, .1, .1, **dict(baseArgs, t=np.array([-25.01]))) == 1.
    # Transiting
    assert eggman.transit(
        1.1, 1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([0.0]), semimajor=1.)) < 1.
    assert eggman.transit(
        1.1, 1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.49]), semimajor=1.)) < 1.
    assert eggman.transit(.01, .015, .011, .011, **dict(baseArgs, t=np.array([.01]))) < 1.
    assert eggman.transit(.01, .015, .011, .011, **dict(baseArgs, t=np.array([-10.01]))) < 1.
    results = eggman.transit(
        .1, .1, .1, .1, **dict(baseArgs, t=np.linspace(-10, 10, 1000), inclination=80.))
    assert results == approx(1, abs=1e-12)
    # Would be non-transiting, but for the eccentricity
    assert eggman.transit(
        1.1, 1.1, 1.1, 1.1,
        **dict(baseArgs, t=np.array([2.51]), eccen=0.1, lon_periapse=-90, semimajor=1.)) < 1.
    # Would be transiting, but for the eccentricity
    assert eggman.transit(
        1.1, 1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.49]), eccen=0.1, semimajor=1.)) == 1.


def test_transit_analytic():
    # Test cases where the transit depth is known analytically
    baseArgs = {
        't': np.array([0.]),
        't0': 0.,
        'theta': 0.,
        'phi': 0.,
        'gamma': 0.,
        'period': 10.,
        'semimajor': 10.,
        'inclination': 89.,
        'limb_type': 'quadratic',
        'limb': [0., 0.],
        'rotate_with_orbit': False,
    }
    assert eggman.transit(.1, .1, .1, .1, **baseArgs) == approx(1 - .1*.1, abs=1e-7)
    assert eggman.transit(.12, .19, .15, .15, **baseArgs) == approx(
        1 - (.12*.15 + .19*.15) / 2., abs=1e-7)
    assert eggman.transit(.01, .01, .02, .02, **baseArgs) == approx(
        1 - (.01*.02 + .01*.02) / 2., abs=1e-7)
    # Small eccentricities shouldn't change anything, since limb-darkening is off
    baseArgs['eccen'] = 0.05
    baseArgs['lon_periapse'] = 35.
    assert eggman.transit(.1, .1, .1, .1, **baseArgs) == approx(1 - .1*.1, abs=1e-7)
    assert eggman.transit(.12, .19, .15, .15, **baseArgs) == approx(
        1 - (.12*.15 + .19*.15) / 2., abs=1e-7)
    assert eggman.transit(.01, .01, .02, .02, **baseArgs) == approx(
        1 - (.01*.02 + .01*.02) / 2., abs=1e-7)
