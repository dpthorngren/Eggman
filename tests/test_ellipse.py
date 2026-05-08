import numpy as np
from pytest import approx

from eggman.cy_eggman import EllipseWrap
from eggman.utils import rotation_matrix


def test_properties():
    ci = np.cos(10.)
    si = np.sin(10.)
    e1 = np.array([2. * ci, 2. * si, 0.])
    e2 = np.array([-.75 * si, .75 * ci, 0.])
    ell = EllipseWrap(e1, e2)
    assert np.linalg.norm(ell.e1) == approx(2.0)
    assert np.linalg.norm(ell.e2) == approx(.75)
    assert ell.x_size == approx(np.sqrt(e1[0]**2 + e2[0]**2))
    assert ell.y_size == approx(np.sqrt(e1[1]**2 + e2[1]**2))


def test_rotation_constructor():
    rot = rotation_matrix(37.3, 0)
    ell = EllipseWrap.create_from_rot_radii(2.1, 1.5, rot)
    assert ell.x_size == approx(2.1)
    assert ell.y_size == approx(1.5 * np.cos(37.3 * np.pi / 180))
    rot = rotation_matrix(53.1, 1)
    ell = EllipseWrap.create_from_rot_radii(2.1, 1.5, rot)
    assert ell.x_size == approx(2.1 * np.cos(53.1 * np.pi / 180))
    assert ell.y_size == approx(1.5)
    rot = rotation_matrix(53.1, 1) @ rotation_matrix(37.3, 0)
    ell = EllipseWrap.create_from_rot_radii(3.0, 1.1, rot)
    assert 3.0 * np.cos(53.1 * np.pi / 180) < ell.x_size < 3
    assert ell.y_size == approx(1.1 * np.cos(37.3 * np.pi / 180))
    rot = rotation_matrix(53.1, 2)
    ell = EllipseWrap.create_from_rot_radii(0.5, 0.9, rot)
    assert ell.e1[2] == 0
    assert ell.e2[2] == 0
    assert 0.5 < ell.x_size < 0.9
    assert 0.5 < ell.y_size < 0.9


def test_circle():
    for angle in np.random.rand(20):
        rot = rotation_matrix(angle, 2)
        ell = EllipseWrap.create_from_rot_radii(1.5, 1.5, rot)
        assert np.linalg.norm(ell.e1) == approx(1.5)
        assert np.linalg.norm(ell.e2) == approx(1.5)
        assert ell.x_size == approx(1.5)
        assert ell.y_size == approx(1.5)
        ell = EllipseWrap.create_from_rot_radii(0.9, 0.9, rot)
        for x, y in np.random.rand(200, 2):
            expect = x*x + y*y < .9**2
            assert ell.line_intersects(x, y) == expect
            expect = (x,y) if expect else .9*np.array([x, y]) / np.sqrt(x*x+y*y)
            assert ell.nearest_to_line(x, y)[:2] == approx(expect, abs=1e-6)
