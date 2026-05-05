import numpy as np
from eggman.cy_eggman import _solve_kepler
from pytest import approx

import eggman


def test_kepler():
    '''Tests that eggman._solve_kepler successfully solves Kepler's equation.'''
    # Known values
    for eccen in np.random.rand(100):
        assert _solve_kepler(0, eccen) == approx(0., 1e-12)
        assert _solve_kepler(np.pi, eccen) == approx(np.pi, 1e-12)

    # Check that other solutions actually solve the equation
    for ma in np.random.rand(100) * 50 - 10:
        assert _solve_kepler(ma, 0) == approx(ma, 1e-12)
        for eccen in np.random.rand(1000):
            if eccen == 1:
                continue
            ea = _solve_kepler(ma, eccen)
            assert ea - eccen * np.sin(ea) == approx(ma, 1e-12)


def test_signs():
    # Checking that the planet is in the correct part of the orbit to catch sign issues.
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


def test_orbit_positions():
    # Simple cases with known results
    for eccen, lon_periapse, inc in np.random.rand(100, 3):
        # t=0 is mid-transit by construction, so x=0
        assert eggman.orbit_to_position(0., 10, 15, eccen, inc, lon_periapse)[0] == approx(
            0., abs=1e-9)
        # t=period/2 is mid-secondary-transit if lon_periapse=+/090, so x=0
        assert eggman.orbit_to_position(5.0, 10, 15, eccen, inc, 90.)[0] == approx(0., abs=1e-9)
        assert eggman.orbit_to_position(5.0, 10, 15, eccen, inc, -90.)[0] == approx(0., abs=1e-9)
    # Circular orbit checks
    assert eggman.orbit_to_position(0., 4., 1., 0.0, 90., 0.) == approx([0., 0., 1.], abs=1e-10)
    assert eggman.orbit_to_position(0.5, 4., 1., 0.0, 90.,
                                    0.) == approx([np.sqrt(.5), 0., np.sqrt(.5)], abs=1e-10)
    for inc, t in np.random.rand(100, 2) * np.pi:
        theta = t * 2 * np.pi / 4.
        pos = [np.sin(theta), np.cos(inc) * np.cos(theta), np.sin(inc) * np.cos(theta)]
        assert eggman.orbit_to_position(t, 4., 1., 0.0, inc * 180 / np.pi, 0.) == approx(
            pos, abs=1e-9)


def test_numerical():
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
