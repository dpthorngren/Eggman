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


def test_circles():
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
            expect = (x, y) if expect else .9 * np.array([x, y]) / np.sqrt(x*x + y*y)
            assert ell.nearest_to_line(x, y)[:2] == approx(expect, abs=1e-6)


def test_ellipses():
    for angle, r1, r2 in np.random.rand(20, 3):
        angle *= 360.
        r1 = 2*r1 + .2
        r2 = 2*r2 + .2
        rot = rotation_matrix(angle, 2)
        assert rot.T @ rot == approx(np.eye(3), rel=1e-12)
        ell = EllipseWrap.create_from_rot_radii(r1, r2, rot)
        assert np.linalg.norm(ell.e1) == approx(r1, rel=1e-9)
        assert np.linalg.norm(ell.e2) == approx(r2, rel=1e-9)
        # Check reverse transform
        assert rot.T @ np.row_stack(ell.e1) == approx([r1, 0., 0.])
        assert rot.T @ np.row_stack(ell.e2) == approx([0., r2, 0.])
        # Area is conserved under rotation
        area = 2 * np.pi * np.linalg.norm(ell.e1) * np.linalg.norm(ell.e2)
        assert area == approx(2 * np.pi * r1 * r2, rel=1e-9)
        # Line intersection and nearest-point check
        for x, y in 6 * np.random.rand(200, 2) - 3:
            xp, yp, zp = np.ravel(rot.T @ np.row_stack([x, y, 0.]))
            assert zp == 0.
            assert xp**2 + yp**2 == approx(x**2 + y**2, 1e-12)
            dist = (xp / r1)**2 + (yp / r2)**2
            assert ell.line_intersects(x, y) == (dist < 1)
            nearest = ell.nearest_to_line(x, y)
            if dist < 1:
                assert nearest == approx([x, y, 0.], 1e-9)
            else:
                # Point should be on the ellipse boundary
                assert nearest[2] == approx(0., rel=1e-12)
                nearestp = np.ravel(rot.T @ np.row_stack(nearest))
                assert nearestp[2] == approx(0., rel=1e-12)
                dist = (nearestp[0] / r1)**2 + (nearestp[1] / r2)**2
                assert dist == approx(1, rel=1e-12)
                # Point should be closer to the target than the origin
                displacement = np.array([x - nearest[0], y - nearest[1]])
                assert np.linalg.norm(displacement) < np.linalg.norm([x, y])
                # Normal vector at nearest point should point from it to the target point
                displacement = np.array([xp - nearestp[0], yp - nearestp[1]])
                normal_vec = np.array([nearestp[0] / (r1**2), nearestp[1] / (r2**2)])
                normal_vec /= np.linalg.norm(normal_vec)
                displacement /= np.linalg.norm(displacement)
                assert normal_vec == approx(displacement, rel=1e-9)
