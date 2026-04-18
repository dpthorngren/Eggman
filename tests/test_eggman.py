import numpy as np
from pytest import approx

import eggman


def testKepler():
    '''Tests that eggman._solve_kepler successfully solves Kepler's equation.'''
    # Known values
    for eccen in np.random.rand(100):
        assert eggman._solve_kepler(0, eccen) == approx(0., 1e-12)
        assert eggman._solve_kepler(np.pi, eccen) == approx(np.pi, 1e-12)

    # Check that other solutions actually solve the equation
    for ma in np.random.rand(100) * 50 - 10:
        assert eggman._solve_kepler(ma, 0) == approx(ma, 1e-12)
        for eccen in np.random.rand(1000):
            if eccen == 1:
                continue
            ea = eggman._solve_kepler(ma, eccen)
            assert ea - eccen * np.sin(ea) == approx(ma, 1e-12)


def testOrbitPositions():
    # Simple cases with known results
    for eccen, lon_periapse, inc in np.random.rand(100, 3):
        # t=0 is mid-transit by construction, so x=0
        assert eggman.orbit_to_position(0., 10, 15, eccen, inc, lon_periapse)[0] == approx(0., abs=1e-9)
        # t=period/2 is mid-secondary-transit if lon_periapse=+/090, so x=0
        assert eggman.orbit_to_position(5.0, 10, 15, eccen, inc, 90.)[0] == approx(0., abs=1e-9)
        assert eggman.orbit_to_position(5.0, 10, 15, eccen, inc, -90.)[0] == approx(0., abs=1e-9)
    # Circular orbit checks
    assert eggman.orbit_to_position(0., 4., 1., 0.0, 90., 0.) == approx([0., 0., 1.], abs=1e-10)
    assert eggman.orbit_to_position(0.5, 4., 1., 0.0, 90., 0.) == approx([np.sqrt(.5), 0., np.sqrt(.5)], abs=1e-10)
    for inc, t in np.random.rand(100, 2) * np.pi:
        theta = t * 2 * np.pi / 4.
        pos = [np.sin(theta), np.cos(inc) * np.cos(theta), np.sin(inc) * np.cos(theta)]
        assert eggman.orbit_to_position(t, 4., 1., 0.0, inc * 180 / np.pi, 0.) == approx(pos, abs=1e-9)
    # Checking the signs of the output for correct orientation
    pos = eggman.orbit_to_position(0., 2.3, 13.5, 0.53, 82., 13.)
    assert pos[0] == approx(0, abs=1e-10)
    assert pos[1] > 0
    assert pos[2] > 0
    pos = eggman.orbit_to_position(.13, 2.3, 13.5, 0.53, 105., 12.)
    assert pos[0] > 0
    assert pos[1] < 0
    assert pos[2] > 0
    pos = eggman.orbit_to_position(1.17, 2.3, 13.5, 0.53, 85., 95.)
    assert pos[0] < 0
    assert pos[1] < 0
    assert pos[2] < 0
    # Test cases from Batman (https://github.com/lkreidberg/batman, batman-package 2.5.3 on pip)
    # Periastron and apoastron positions
    t_peri = -0.3569835522205464
    pos = eggman.orbit_to_position(t_peri, 5.3, 10., 0.3, 87., 45.)
    assert np.linalg.norm(pos) == approx(10. * (1.-.3), 1e-9)
    t_apo = 5.3/2. - 0.3569835522205464
    pos = eggman.orbit_to_position(t_apo, 5.3, 10., 0.3, 87., 45.)
    assert np.linalg.norm(pos) == approx(10. * (1.+.3), 1e-9)
    # Distance on-sky from star
    # eggman.orbit_to_position args: t, period, semimajor, eccen, inclination, lon_periapse
    # batman.rsky args: t, t_conj, period, semimajor, inclination, eccen, lon_periapse, transittype, nthreads
    pos = eggman.orbit_to_position(2.1, 5., 15., .3, 89., 35.)
    # batman._rsky._rsky(np.array([2.1]), 0., 5., 15., 89.*np.pi/180., .3, 35*np.pi/180., 0, 1)
    assert np.linalg.norm(pos[:2]) == approx(15.80092098, abs=1e-6)
    pos = eggman.orbit_to_position(-1.2, 4., 13., .3, 89., 63.)
    # batman._rsky._rsky(np.array([-1.2]), 0., 4., 13., 89.*np.pi/180., .3, 63*np.pi/180., 0, 1)
    assert np.linalg.norm(pos[:2]) == approx(6.59137875, abs=1e-6)
    pos = eggman.orbit_to_position(-1.2, 4., 13., .3, 85., 63.)
    # batman._rsky._rsky(np.array([-1.2]), 0., 4., 13., 85.*np.pi/180., .3, 63*np.pi/180., 0, 1)
    assert np.linalg.norm(pos[:2]) == approx(6.68045597, abs=1e-6)


# def testGrazing():
#     # Grazing transits can cause numerical issues, let's double-check that doesn't happen.
#     # Grazing on ingress
#     baseArgs = {
#         't': np.array([-np.arcsin(1.1 / 10)]),
#         't0': 0.,
#         'period': 2 * np.pi,
#         'semimajor': 10.,
#         'inclination': 90.,
#         'limbType': 'quadratic',
#         'limb': [.3, .2]
#     }
#
#     # Check for errors on grazing transits
#     assert eggman.asymmetricTransit(.1 + 1e-3, .11, .1, **dict(baseArgs)) < 1.
#     assert eggman.asymmetricTransit(.1 + 1e-6, .09, -1, **dict(baseArgs)) < 1.
#     assert eggman.asymmetricTransit(.1 + 1e-9, .1, .09, **dict(baseArgs)) < 1.
#
#     # Check for errors on barely non-grazing trants
#     assert eggman.asymmetricTransit(.1 - 1e-3, .11, .1, **dict(baseArgs)) == 1.
#     assert eggman.asymmetricTransit(.1 - 1e-6, .09, -1, **dict(baseArgs)) == 1.
#     assert eggman.asymmetricTransit(.1 - 1e-12, .1, .09, **dict(baseArgs)) == 1.
#
#     # Grazing on egress
#     baseArgs['t'] = -baseArgs['t']
#     assert eggman.asymmetricTransit(.11, .1 + 1e-3, .1, **dict(baseArgs)) < 1.
#     assert eggman.asymmetricTransit(.09, .1 + 1e-6, -1, **dict(baseArgs)) < 1.
#     assert eggman.asymmetricTransit(.1, .1 + 1e-9, .09, **dict(baseArgs)) < 1.
#
#     # Grazing on top (true grazing transit)
#     baseArgs['inclination'] = 90. - abs(baseArgs['t'][0]) * 180 / np.pi
#     baseArgs['t'] = np.array([0.])
#     assert eggman.asymmetricTransit(.11, .1, .1 + 1e-3, **dict(baseArgs)) < 1.
#     assert eggman.asymmetricTransit(.09, .1, .1 + 1e-6, **dict(baseArgs)) < 1.
#     assert eggman.asymmetricTransit(.1, .09, .1 + 1e-6, **dict(baseArgs)) < 1.
#
#     # Barely non-grazing on top (true grazing transit)
#     assert eggman.asymmetricTransit(.11, .1, .1 - 1e-3, **dict(baseArgs)) == 1.
#     assert eggman.asymmetricTransit(.09, .1, .1 - 1e-6, **dict(baseArgs)) == 1.
#     assert eggman.asymmetricTransit(.1, .09, .1 - 1e-12, **dict(baseArgs)) == 1.
#
#
# def testIsTransiting():
#     # Ensure we correctly can tell transit from non-transit
#     baseArgs = {
#         't': np.array([0.]),
#         't0': 0.,
#         'period': 10.,
#         'semimajor': 10.,
#         'inclination': 89.,
#         'limbType': 'quadratic',
#         'limb': [.3, .2],
#         'eccen': 0,
#         'lonPeriapse': 90.
#     }
#     # Non-transiting
#     assert eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([5.]))) == 1.
#     assert eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([3.]))) == 1.
#     assert eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.51]), semimajor=1.)) == 1.
#     assert eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([-25.01]))) == 1.
#     # Transiting
#     assert eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.49]), semimajor=1.)) < 1.
#     assert eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=np.array([.01]))) < 1.
#     assert eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=np.array([-10.01]))) < 1.
#     results = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.linspace(-10, 10, 1000), inclination=80.))
#     assert results == approx(1, abs=1e-12)
#     # Would be non-transiting, but for the eccentricity
#     assert eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.51]), eccen=0.1, lonPeriapse=-90, semimajor=1.)) < 1.
#     # Would be transiting, but for the eccentricity
#     assert eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.49]), eccen=0.1, semimajor=1.)) == 1.
#
#
# def testAsymmetricAnalytic():
#     # Test cases where the transit depth is known analytically
#     baseArgs = {
#         't': np.array([0.]),
#         't0': 0.,
#         'period': 10.,
#         'semimajor': 10.,
#         'inclination': 89.,
#         'limbType': 'quadratic',
#         'limb': [0., 0.]
#     }
#     assert eggman.asymmetricTransit(.1, .1, .1, **baseArgs) == approx(1 - .1*.1, abs=1e-7)
#     assert eggman.asymmetricTransit(.1, .2, -1, **baseArgs) == approx(1 - (.1*.1 + .2*.2) / 2., abs=1e-7)
#     assert eggman.asymmetricTransit(.12, .19, .15, **baseArgs) == approx(1 - (.12*.15 + .19*.15) / 2., abs=1e-7)
#     assert eggman.asymmetricTransit(.01, .01, .02, **baseArgs) == approx(1 - (.01*.02 + .01*.02) / 2., abs=1e-7)
#     # Small eccentricities shouldn't change anything, since limb-darkening is off
#     baseArgs['eccen'] = 0.05
#     baseArgs['lonPeriapse'] = 35.
#     assert eggman.asymmetricTransit(.1, .1, .1, **baseArgs) == approx(1 - .1*.1, abs=1e-7)
#     assert eggman.asymmetricTransit(.1, .2, -1, **baseArgs) == approx(1 - (.1*.1 + .2*.2) / 2., abs=1e-7)
#     assert eggman.asymmetricTransit(.12, .19, .15, **baseArgs) == approx(1 - (.12*.15 + .19*.15) / 2., abs=1e-7)
#     assert eggman.asymmetricTransit(.01, .01, .02, **baseArgs) == approx(1 - (.01*.02 + .01*.02) / 2., abs=1e-7)
#
#
# def testAsymmetricNumerical():
#     baseArgs = {
#         't': np.array([0.001, .01]),
#         't0': 0.,
#         'period': 1.,
#         'semimajor': 15.,
#         'inclination': 90.,
#         'limbType': 'quadratic',
#         'limb': [.1, .3]
#     }
#     # Test cases where the transit depth was calculated from other codes.
#     result = eggman.asymmetricTransit(.11, .1, -1, **baseArgs)
#     assert result[0] == approx(0.9879552837718931, abs=1e-7)
#     assert result[1] == approx(0.9923392211352866, abs=1e-7)
