import numpy as np
from pytest import approx

from eggman.cy_eggman import BiellipsoidWrap
from eggman.utils import rotation_matrix


def _decompose_position_(x, y, z):
    '''Determines radius, sin and cosine of the azimuthal angle, and sin and cosine of
    the inclination angle consistent with the given position. These are only orbital
    elements if the orbit is circular. There are multiple solutions if z=0, assumes
    i=90 in that case. Used for creating valid test cases.'''
    rad = np.sqrt(x*x + y*y + z*z)
    st = x / rad if rad > 0 else 0.
    ct = np.sign(z) * np.sqrt(1 - st*st)
    rad = np.sqrt(y*y + z*z)
    ci = np.sign(ct) * y / rad if rad > 0 else 0.
    si = np.sqrt(1 - ci*ci)
    return rad, ct, st, ci, si


def test_biellipsoid_rotation():
    # Check that the rotation matrix produced by initialization is good.
    bell = BiellipsoidWrap(0., 0., 0., 0., 0., 0., 1., 1., 1., 1.)
    assert bell.rot.flat == approx([1., 0., 0., 0., 1., 0., 0., 0., 1.])
    # Make sure theta rotates counterclockwise along view vector
    bell = BiellipsoidWrap(0., 0., 0., 15., 0., 0., 1., 1., 1., 1.)
    assert bell.forward == approx([np.cos(15 * np.pi / 180), np.sin(15 * np.pi / 180), 0.])
    # Make sure phi rotates counterclockwise around up vector (looking down)
    bell = BiellipsoidWrap(0., 0., 0., 0., 25., 0., 1., 1., 1., 1.)
    assert bell.rot[:, 0] == approx([np.cos(25 * np.pi / 180), 0., -np.sin(25 * np.pi / 180)])
    # Make sure gamma rotates counterclockwise around the biellipse forward vector
    bell = BiellipsoidWrap(0., 0., 0., 0., 0., 35., 1., 1., 1., 1.)
    assert bell.forward == approx([1., 0., 0.])
    assert bell.up == approx([0., np.cos(35 * np.pi / 180), np.sin(35 * np.pi / 180)])
    # Valid rotation matrix must have R @ R.T == I
    for angles in 179. * np.random.rand(100, 3) - 90.:
        bell = BiellipsoidWrap(0., 0., 0., angles[0], angles[1], angles[2], 1., 1., 1., 1.)
        assert bell.rot.T @ bell.rot == approx(np.eye(3))
    # Compare with rotation matrix helper
    bell = BiellipsoidWrap(3., 5., -2., 30., 10., 0., 3., 2., 1., 1.5)
    expect = rotation_matrix(10, 1) @ rotation_matrix(30, 2)
    assert bell.rot == approx(expect, abs=1e-12)
    bell = BiellipsoidWrap(3., 5., -2., 105., 67.3, 13.2, 3., 2., 1., 1.5)
    expect = rotation_matrix(67.3, 1) @ rotation_matrix(105, 2) @ rotation_matrix(13.2, 0)
    assert bell.rot == approx(expect, abs=1e-12)


def test_biellipsoid_orbit_rotation():
    # Valid rotation matrix must have R @ R.T == I
    for angles in 179. * np.random.rand(100, 3) - 90.:
        pos = np.random.rand(3)
        ci = _decompose_position_(*pos)[3]
        bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 1., ci)
        assert bell.rot.T @ bell.rot == approx(np.eye(3))
    # Check that rotations go in the right direction
    pos = [1.5, 0., 0.]
    angles = [0., 0., 0.]
    bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 2., 0.)
    assert bell.x_bounds() == approx([-0.5, 3.5])
    assert bell.y_bounds() == approx([-1, 1])
    pos = [-1.5, 0., 0.]
    bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 2., 1.)
    assert bell.x_bounds() == approx([-3.5, .5])
    assert bell.y_bounds() == approx([-1, 1])
    # Check that rotations combine properly
    pos = [1.5, 0., 0.]
    angles = [0., 0., 0.]
    bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 2., 0.)
    assert bell.x_bounds() == approx([-0.5, 3.5])
    assert bell.y_bounds() == approx([-1, 1])
    pos = [0., 0., 3.]
    angles = [0., 180., 0.]
    bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 2., 0.)
    assert bell.rot == approx(np.diag([-1, 1, -1]))
    # Verify no rotation at orbital origin (0, 0, 1)
    pos = [0., 0., 2.]
    angles = [0., 0., 0.]
    bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 2., -2)
    assert bell.rot == approx(np.eye(3))
    bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 2., 0.)
    assert bell.rot == approx(np.eye(3))
    # Compare with rotation matrix helper
    bell = BiellipsoidWrap(3., 0., 3., 30., 10., 0., 3., 2., 1., 1.5, 0.)
    expect = rotation_matrix(10, 1) @ rotation_matrix(30, 2) @ rotation_matrix(45., 1)
    assert bell.rot == approx(expect, abs=1e-12)
    # Check that the correct axis points towards the origin
    for pos in np.random.rand(100, 3) * 2 - 1:
        pos /= np.linalg.norm(pos)
        angles = [0., 0., 0.]
        ci = _decompose_position_(*pos)[3]
        bell = BiellipsoidWrap(*pos, *angles, 1., 1., 1., 1., ci)
        assert bell.side == approx(pos, 1e-12)
        assert bell.up.dot(pos) == approx(0, abs=1e-12)


def test_biellipsoid_bounds():
    # Check bounds (axis-aligned maximum x and y extents)
    bell = BiellipsoidWrap(0., 0., 0., 0., 0, 0, 1., 2., 3., 4.)
    assert bell.x_bounds() == approx([-2, 1])
    assert bell.y_bounds() == approx([-3, 3])
    bell = BiellipsoidWrap(0., 0., 0., 15., 0, 0, 1., 2., 1., 1.)
    xmin, xmax = bell.x_bounds()
    assert -2 < xmin < -2 * np.cos(15 * np.pi / 180)
    assert xmax == approx(1.0)
    ymin, ymax = bell.y_bounds()
    assert -2.0 < ymin < -1.0
    assert ymax == approx(1.0)


def test_transforms():
    for i in range(20):
        angles = 179. * np.random.rand(3) - 90.
        radii = 1.5 * np.random.rand(4) + .5
        position = 10 * np.random.rand(3) - 5.
        bell = BiellipsoidWrap(*position, *angles, *radii)
        # Check origin transforms first
        assert bell.world_to_sphere(position) == approx([0., 0., 0.])
        assert bell.world_to_aligned(position) == approx([0., 0., 0.])
        assert bell.sphere_to_world([0., 0., 0.]) == approx(position)
        assert bell.aligned_to_world([0., 0., 0.]) == approx(position)
        assert bell.sphere_to_aligned([0., 0., 0.]) == approx([0., 0., 0.])
        assert bell.aligned_to_sphere([0., 0., 0.]) == approx([0., 0., 0.])
        # Check unit vectors
        for v, r in zip(np.eye(3), bell.rot):
            assert bell.world_to_aligned(position + v) == approx(r)
        # Check random points
        for i in range(20):
            loc = 10 * np.random.rand(3) - 5.
            assert bell.aligned_to_world(loc) == approx(bell.rot @ loc + position)
            loc_sphere = bell.world_to_sphere(loc)
            assert bell.sphere_to_world(loc_sphere) == approx(loc)
            assert bell.aligned_to_world(bell.sphere_to_aligned(loc_sphere)) == approx(loc)
            assert bell.aligned_to_sphere(bell.world_to_aligned(loc)) == approx(loc_sphere)


def test_biellipsoid_point_visible():
    bell = BiellipsoidWrap(0., 0., 0., 30., 1.5, 0, 1.33, 1.33, 1.33, 1.33)
    assert not bell.is_visible(0, 0, 0)
    assert not bell.is_visible(0., 0., -10.0)
    assert bell.is_visible(1., 1., 0.)
    assert not bell.is_visible(0.7, 0.7, 0)
    assert bell.is_visible(0., 0., 10.0)
    bell = BiellipsoidWrap(5., 3., 0., 30., 1.5, 0, .15, .1, .1, .1)
    assert not bell.is_visible(5, 3, 0)
    assert not bell.is_visible(5., 3., -10.0)
    assert not bell.is_visible(5.02, 2.98, 0.)
    assert not bell.is_visible(5.07, 3.07, 0)
    assert bell.is_visible(5., 3., 10.0)


def test_biellipsoid_intersect():
    bell = BiellipsoidWrap(0., 0., 0., 0, 0, 0, 1.33, 1.33, 1.33, 1.33)
    hit = bell.line_project(.5, .5)
    assert hit[0] == 0.5
    assert hit[1] == 0.5
    assert np.linalg.norm(hit / 1.33) == approx(1.0)
    # Changing the far ellipsoid shouldn't matter
    bell = BiellipsoidWrap(0., 0., 0., 0, 0, 0, 1.2, 0.5, 1.2, 1.2)
    hit = bell.line_project(.5, .5)
    assert hit[0] == 0.5
    assert hit[1] == 0.5
    assert np.linalg.norm(hit / 1.2) == approx(1.0)
    for i in range(100):
        loc = 10 * np.random.rand(3) - 5
        angles = 179. * np.random.rand(3) - 90. * np.array([0., 0., 0.])
        radii = 1.5 * np.random.rand(4) + .5
        bell = BiellipsoidWrap(*loc, *angles, *radii)
        for pos in np.random.rand(100, 2) * 4 - 2 + loc[:2]:
            hit = bell.line_project(pos[0], pos[1])
            if bell.line_intersects(pos[0], pos[1]):
                assert hit[:2] == approx(pos[:2])
                pos_rot = np.ravel(bell.rot.T @ np.row_stack(hit - loc))
                pos_rot[0] /= bell.r_forward if pos_rot[0] > 0 else bell.r_back
                pos_rot[1] /= bell.r_up
                pos_rot[2] /= bell.r_side
                assert np.linalg.norm(pos_rot) == approx(1)
            else:
                assert np.all(np.isnan(hit))


def test_biellipsoid_slice_ylimits():
    # Check that the ylimits for a given x are correctly calculated
    # Circle
    bell = BiellipsoidWrap(0., 0., 0., 0, 0, 0, 1.33, 1.33, 1.33, 1.33)
    bounds = bell.slice_ylimits(0.)
    assert bounds[0] == approx(-1.33, 1e-12)
    assert bounds[1] == approx(1.33, 1e-12)
    bell = BiellipsoidWrap(0., 0., 0., 0, 0, 0, 1., 1., 1., 1.)
    bounds = bell.slice_ylimits(-.5)
    assert bounds[0] == approx(-np.sqrt(1 - .25), 1e-12)
    assert bounds[1] == approx(np.sqrt(1 - .25), 1e-12)
    bounds = bell.slice_ylimits(.5)
    assert bounds[0] == approx(-np.sqrt(1 - .25), 1e-12)
    assert bounds[1] == approx(np.sqrt(1 - .25), 1e-12)

    # Simple Ellipses
    bell = BiellipsoidWrap(0., 0., 0., 90., 0, 0, 1., 2., 2., 1.)
    assert bell.slice_ylimits(0.) == approx([-2., 1.])
    bell = BiellipsoidWrap(0., 0., 0., 0., 0, 0, 1., 1., 2., 1.)
    bounds = bell.slice_ylimits(-0.5)
    assert bounds == approx([-2 * np.sqrt(1 - .25), 2 * np.sqrt(1 - .25)], 1e-12)
    bell = BiellipsoidWrap(0., 0., 0., 90., 0, 0, 1., 2., 1., 1.)
    assert bell.slice_ylimits(0.) == approx([-2., 1.])
    bell = BiellipsoidWrap(0., 0., 0., 90., 0, 0, 1., 2., 1., 1.)
    bounds = bell.slice_ylimits(-0.2)
    assert bounds[0] == approx(-2 * np.sqrt(1 - .2*.2), 1e-12)
    assert bounds[1] == approx(np.sqrt(1 - .2*.2), 1e-12)
