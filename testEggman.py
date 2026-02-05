import numpy as np
from pytest import approx

import eggman


def testGrazing():
    # Grazing transits can cause numerical issues, let's double-check that doesn't happen.
    # Grazing on ingress
    baseArgs = {
        't': np.array([-np.arcsin(1.1 / 10)]),
        't0': 0.,
        'period': 2 * np.pi,
        'semimajor': 10.,
        'inclination': 90.,
        'limbType': 'quadratic',
        'limb': [.3, .2]
    }

    # Check for errors on grazing transits
    assert eggman.asymmetricTransit(.1 + 1e-3, .11, .1, **dict(baseArgs)) < 1.
    assert eggman.asymmetricTransit(.1 + 1e-6, .09, -1, **dict(baseArgs)) < 1.
    assert eggman.asymmetricTransit(.1 + 1e-9, .1, .09, **dict(baseArgs)) < 1.

    # Check for errors on barely non-grazing trants
    assert eggman.asymmetricTransit(.1 - 1e-3, .11, .1, **dict(baseArgs)) == 1.
    assert eggman.asymmetricTransit(.1 - 1e-6, .09, -1, **dict(baseArgs)) == 1.
    assert eggman.asymmetricTransit(.1 - 1e-12, .1, .09, **dict(baseArgs)) == 1.

    # Grazing on egress
    baseArgs['t'] = -baseArgs['t']
    assert eggman.asymmetricTransit(.11, .1 + 1e-3, .1, **dict(baseArgs)) < 1.
    assert eggman.asymmetricTransit(.09, .1 + 1e-6, -1, **dict(baseArgs)) < 1.
    assert eggman.asymmetricTransit(.1, .1 + 1e-9, .09, **dict(baseArgs)) < 1.

    # Grazing on top (true grazing transit)
    baseArgs['inclination'] = 90. - abs(baseArgs['t'][0]) * 180 / np.pi
    baseArgs['t'] = np.array([0.])
    assert eggman.asymmetricTransit(.11, .1, .1 + 1e-3, **dict(baseArgs)) < 1.
    assert eggman.asymmetricTransit(.09, .1, .1 + 1e-6, **dict(baseArgs)) < 1.
    assert eggman.asymmetricTransit(.1, .09, .1 + 1e-6, **dict(baseArgs)) < 1.

    # Barely non-grazing on top (true grazing transit)
    assert eggman.asymmetricTransit(.11, .1, .1 - 1e-3, **dict(baseArgs)) == 1.
    assert eggman.asymmetricTransit(.09, .1, .1 - 1e-6, **dict(baseArgs)) == 1.
    assert eggman.asymmetricTransit(.1, .09, .1 - 1e-12, **dict(baseArgs)) == 1.


def testIsTransiting():
    # Ensure we correctly can tell transit from non-transit
    baseArgs = {
        't': np.array([0.]),
        't0': 0.,
        'period': 10.,
        'semimajor': 10.,
        'inclination': 89.,
        'limbType': 'quadratic',
        'limb': [.3, .2]
    }
    # Non-transiting
    assert eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([5.]))) == 1.
    assert eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([3.]))) == 1.
    assert eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.51]), semimajor=1.)) == 1.
    assert eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([-25.01]))) == 1.
    # Transiting
    assert eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.49]), semimajor=1.)) < 1.
    assert eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=np.array([.01]))) < 1.
    assert eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=np.array([-10.01]))) < 1.
    results = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.linspace(-10, 10, 1000), inclination=80.))
    assert results == approx(1, abs=1e-12)


def testAsymmetricAnalytic():
    # Test cases where the transit depth is known analytically
    baseArgs = {
        't': np.array([0.]),
        't0': 0.,
        'period': 10.,
        'semimajor': 10.,
        'inclination': 89.,
        'limbType': 'quadratic',
        'limb': [0., 0.]
    }
    assert eggman.asymmetricTransit(.1, .1, .1, **baseArgs) == approx(1 - .1*.1, abs=1e-6)
    assert eggman.asymmetricTransit(.1, .2, -1, **baseArgs) == approx(1 - (.1*.1 + .2*.2) / 2., abs=1e-6)
    assert eggman.asymmetricTransit(.12, .19, .15, **baseArgs) == approx(1 - (.12*.15 + .19*.15) / 2., abs=1e-6)
    assert eggman.asymmetricTransit(.01, .01, .02, **baseArgs) == approx(1 - (.01*.02 + .01*.02) / 2., abs=1e-6)


def testAsymmetricNumerical():
    baseArgs = {
        't': np.array([0.001, .01]),
        't0': 0.,
        'period': 1.,
        'semimajor': 15.,
        'inclination': 90.,
        'limbType': 'quadratic',
        'limb': [.1, .3]
    }
    # Test cases where the transit depth was calculated from other codes.
    result = eggman.asymmetricTransit(.11, .1, -1, **baseArgs)
    assert result[0] == approx(0.9879552837718931, abs=1e-7)
    assert result[1] == approx(0.9923392211352866, abs=1e-7)
