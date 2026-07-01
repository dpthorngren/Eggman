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
    assert eggman.transit(.1 + 1e-3, .09, -1, .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.1 + 1e-9, .1, .09, .09, **dict(baseArgs)) < 1.

    # Check for errors on barely non-grazing transits
    assert eggman.transit(.1 - 1e-3, .11, .1, .1, **dict(baseArgs)) == 1.
    assert eggman.transit(.1 - 1e-3, .11, -1, .1, **dict(baseArgs)) == 1.
    assert eggman.transit(.1 - 1e-12, .1, .09, .09, **dict(baseArgs)) == 1.

    # Grazing on egress
    baseArgs['t'] = -baseArgs['t']
    assert eggman.transit(.11, .1 + 1e-3, .1, .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.1, .1 + 1e-9, .09, .09, **dict(baseArgs)) < 1.

    # Grazing on top (true grazing transit)
    baseArgs['inclination'] = 90. - abs(baseArgs['t'][0]) * 180 / np.pi
    baseArgs['t'] = np.array([0.])
    assert eggman.transit(.11, .1, .1 + 1e-3, .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.1 + 1e-3, .1 - 1e-3, -1., .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.09, .1, .1 + 1e-6, .1, **dict(baseArgs)) < 1.
    assert eggman.transit(.1, .09, .1 + 1e-6, .1, **dict(baseArgs)) < 1.

    # Barely non-grazing on top (true grazing transit)
    assert eggman.transit(.11, .1, .1 - 1e-3, .1, **dict(baseArgs)) == 1.
    assert eggman.transit(.1 - 1e-3, .1 - 1e-3, -1, .1, **dict(baseArgs)) == 1.
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
    assert eggman.transit(.1, .2, -1, .1, **baseArgs) == approx(1 - (.1*.1 + .2*.2) / 2., abs=1e-7)
    assert eggman.transit(.12, .19, .15, .15, **baseArgs) == approx(
        1 - (.12*.15 + .19*.15) / 2., abs=1e-7)
    assert eggman.transit(.01, .01, .02, .02, **baseArgs) == approx(
        1 - (.01*.02 + .01*.02) / 2., abs=1e-7)
    # Small eccentricities shouldn't change anything, since limb-darkening is off
    baseArgs['eccen'] = 0.05
    baseArgs['lon_periapse'] = 35.
    assert eggman.transit(.1, .1, .1, .1, **baseArgs) == approx(1 - .1*.1, abs=1e-7)
    assert eggman.transit(.1, .2, -1, .1, **baseArgs) == approx(1 - (.1*.1 + .2*.2) / 2., abs=1e-7)
    assert eggman.transit(.12, .19, .15, .15, **baseArgs) == approx(
        1 - (.12*.15 + .19*.15) / 2., abs=1e-7)
    assert eggman.transit(.01, .01, .02, .02, **baseArgs) == approx(
        1 - (.01*.02 + .01*.02) / 2., abs=1e-7)


def test_transit_numerical():
    # Test cases where the transit depth was calculated from other codes
    baseArgs = {
        't': np.array([.001, .01]),
        't0': 0.,
        'theta': 0.,
        'phi': 0.,
        'gamma': 0.,
        'period': 1.,
        'semimajor': 15.,
        'inclination': 90.,
        'limb_type': 'quadratic',
        'limb': [0.1, 0.3],
        'rotate_with_orbit': False,
    }

    # Symmetric case, comparing with batman and catwoman
    result = eggman.transit(.12, .12, .12, .12, **baseArgs)
    assert result[0] == approx(0.984304107347, abs=1e-7)
    assert result[1] == approx(0.989985655450, abs=1e-7)

    # Symmetric case, comparing with batman and catwoman
    result = eggman.transit(.1, .1, .1, .10, **baseArgs)
    assert result[0] == approx(0.989098764152, abs=1e-7)
    assert result[1] == approx(0.992627976697, abs=1e-7)

    # Asymmetric discontinuous poles case 1, comparing with catwoman
    result = eggman.transit(.11, .1, -1, .1, **baseArgs)
    assert result[0] == approx(0.987955283022, abs=1e-7)
    assert result[1] == approx(0.992339221135, abs=1e-7)

    # Asymmetric discontinuous poles case 2, comparing with catwoman
    # Note the lower precision of tests.  I think this is catwoman's fault, but I'm not certain.
    baseArgs['inclination'] = 89
    result = eggman.transit(.11, .09, -1, .1, **baseArgs)
    assert result[0] == approx(0.989036885966, abs=1e-6)
    assert result[1] == approx(0.995337628970, abs=1e-6)
