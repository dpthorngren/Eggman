import numpy as np
from pytest import approx

from eggman import Orbit


def test_kepler():
    '''Tests that eggman._solve_kepler successfully solves Kepler's equation.'''
    # Known values
    for eccen in np.random.rand(100):
        assert Orbit.solve_kepler(0, eccen) == approx(0., 1e-12)
        assert Orbit.solve_kepler(np.pi, eccen) == approx(np.pi, 1e-12)

    # Check that other solutions actually solve the equation
    for ma in np.random.rand(100) * 50 - 10:
        assert Orbit.solve_kepler(ma, 0) == approx(ma, 1e-12)
        for eccen in np.random.rand(1000):
            if eccen == 1:
                continue
            ea = Orbit.solve_kepler(ma, eccen)
            assert ea - eccen * np.sin(ea) == approx(ma, 1e-12)


def test_signs():
    # Checking that the planet is in the correct part of the orbit to catch sign issues.
    orb = Orbit(2.3, 0.0, 13.5, 0.53, 82., 13.)
    pos = orb.get_position(0.)
    assert pos[0] == approx(0, abs=1e-10)
    assert pos[1] > 0
    assert pos[2] > 0
    orb = Orbit(2.3, 0.0, 13.5, 0.53, 105., 12.)
    pos = orb.get_position(.13)
    assert pos[0] > 0
    assert pos[1] < 0
    assert pos[2] > 0
    orb = Orbit(2.3, 0.0, 13.5, 0.53, 85., 95.)
    pos = orb.get_position(1.17)
    assert pos[0] < 0
    assert pos[1] < 0
    assert pos[2] < 0


def test_orbit_positions():
    # Simple cases with known results
    for eccen, lon_periapse, inc in np.random.rand(100, 3):
        # t=0 is mid-transit by construction, so x=0
        orb = Orbit(10, 0., 15, eccen, inc, lon_periapse)
        assert orb.get_position(0.)[0] == approx(0., abs=1e-9)
        # t=period/2 is mid-secondary-transit if lon_periapse=+/90, so x=0
        orb = Orbit(10, 0., 15, eccen, inc, 90.)
        assert orb.get_position(5.0)[0] == approx(0., abs=1e-9)
        orb = Orbit(10, 0., 15, eccen, inc, -90.)
        assert orb.get_position(5.0)[0] == approx(0., abs=1e-9)
        # Shifting t0 and t by the same amount has no effect
        shift = 10 * np.random.rand()
        orb = Orbit(10, 0., 15, eccen, inc)
        orb2 = Orbit(10, shift, 15, eccen, inc)
        assert orb.get_position(0.) == approx(orb2.get_position(shift), abs=1e-12)
    # Circular orbit checks
    orb = Orbit(4., 0.5, 1., 0.0, 90., 0.)
    assert orb.get_position(0.5) == approx([0., 0., 1.], abs=1e-10)
    orb = Orbit(4., 0., 1., 0.0, 90., 0.)
    assert orb.get_position(0.5) == approx([np.sqrt(.5), 0., np.sqrt(.5)], abs=1e-10)
    for inc, t in np.random.rand(100, 2) * np.pi:
        theta = t * 2 * np.pi / 4.
        pos = [np.sin(theta), np.cos(inc) * np.cos(theta), np.sin(inc) * np.cos(theta)]
        orb = Orbit(4., 0., 1., 0.0, inc * 180 / np.pi, 0.)
        assert orb.get_position(t) == approx(pos, abs=1e-9)


def test_numerical():
    # Test cases from Batman (https://github.com/lkreidberg/batman, batman-package 2.5.3 on pip)
    # Periastron and apoastron positions
    t_peri = -0.3569835522205464
    orb = Orbit(5.3, 0., 10., 0.3, 87., 45.)
    pos = orb.get_position(t_peri)
    assert np.linalg.norm(pos) == approx(10. * (1.-.3), 1e-9)
    t_apo = 5.3/2. - 0.3569835522205464
    orb = Orbit(5.3, 0., 10., 0.3, 87., 45.)
    assert np.linalg.norm(orb.get_position(t_apo)) == approx(10. * (1.+.3), 1e-9)
    # Distance on-sky from star
    # eggman.OrbitWrap args: period, t0, semimajor, eccen, inclination, lon_periapse
    # batman.rsky args: t, t_conj, period, semimajor, inclination, eccen, lon_periapse, transittype, nthreads
    # batman._rsky._rsky(np.array([2.1]), 0., 5., 15., 89.*np.pi/180., .3, 35*np.pi/180., 0, 1)
    orb = Orbit(5., 0., 15., .3, 89., 35.)
    assert np.linalg.norm(orb.get_position(2.1)[:2]) == approx(15.80092098, abs=1e-6)
    orb = Orbit(4., 0., 13., .3, 89., 63.)
    pos = orb.get_position(-1.2)
    # batman._rsky._rsky(np.array([-1.2]), 0., 4., 13., 89.*np.pi/180., .3, 63*np.pi/180., 0, 1)
    assert np.linalg.norm(pos[:2]) == approx(6.59137875, abs=1e-6)
    orb = Orbit(4., 0., 13., .3, 85., 63.)
    # batman._rsky._rsky(np.array([-1.2]), 0., 4., 13., 85.*np.pi/180., .3, 63*np.pi/180., 0, 1)
    assert np.linalg.norm(orb.get_position(-1.2)[:2]) == approx(6.68045597, abs=1e-6)
