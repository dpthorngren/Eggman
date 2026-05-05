import numpy as np
from pytest import approx

from eggman.cy_eggman import _biellipsoid_slice_ylimits, biellipsoid_dump


def test_biellipsoid_rotation():
    # Check that the rotation matrix produced by initialization is good.
    bell = biellipsoid_dump(0., 0., 0., 0., 0., 0., 1., 1., 1., 1.)
    assert bell['rot'] == approx([1., 0., 0., 0., 1., 0., 0., 0., 1.])
    # Make sure theta rotates counterclockwise along view vector
    bell = biellipsoid_dump(0., 0., 0., 15., 0., 0., 1., 1., 1., 1.)
    assert bell['rot'][::3] == approx([np.cos(15 * np.pi / 180), np.sin(15 * np.pi / 180), 0.])
    # Make sure phi rotates counterclockwise around up vector (looking down)
    bell = biellipsoid_dump(0., 0., 0., 0., 25., 0., 1., 1., 1., 1.)
    assert bell['rot'][::3] == approx([np.cos(25 * np.pi / 180), 0., -np.sin(25 * np.pi / 180)])
    # Make sure gamma rotates counterclockwise around the biellipse forward vector
    bell = biellipsoid_dump(0., 0., 0., 0., 0., 35., 1., 1., 1., 1.)
    assert bell['rot'][::3] == approx([1., 0., 0.])
    assert bell['rot'][1::3] == approx([0., np.cos(35 * np.pi / 180), np.sin(35 * np.pi / 180)])
    # Valid rotation matrix must have R @ R.T == I
    for angles in 179. * np.random.rand(100, 3) - 90.:
        bell = biellipsoid_dump(0., 0., 0., angles[0], angles[1], angles[2], 1., 1., 1., 1.)
        rot = np.array(bell['rot']).reshape(3, 3)
        assert rot.T @ rot == approx(np.eye(3))


def test_biellipsoid_bounds():
    # Check bounds (axis-aligned maximum x and y extents)
    bell = biellipsoid_dump(0., 0., 0., 0., 0, 0, 1., 2., 3., 4.)
    assert bell['x_bounds'] == approx(dict(min=-2, max=1))
    assert bell['y_bounds'] == approx(dict(min=-3, max=3))
    bell = biellipsoid_dump(0., 0., 0., 15., 0, 0, 1., 2., 1., 1.)
    assert -2 < bell['x_bounds']['min'] < -2 * np.cos(15 * np.pi / 180)
    assert bell['x_bounds']['max'] == approx(1.0)
    assert -2.0 < bell['y_bounds']['min'] < -1.0
    assert bell['y_bounds']['max'] == approx(1.0)


def test_biellipsoid_slice_ylimits():
    # Check that the ylimits for a given x are correctly calculated
    # Circle
    bounds = _biellipsoid_slice_ylimits(0., 0, 0, 0, 1.33, 1.33, 1.33, 1.33)
    assert bounds[0] == approx(-1.33, 1e-12)
    assert bounds[1] == approx(1.33, 1e-12)
    bounds = _biellipsoid_slice_ylimits(-.5, 0, 0, 0, 1., 1., 1., 1.)
    assert bounds[0] == approx(-np.sqrt(1 - .25), 1e-12)
    assert bounds[1] == approx(np.sqrt(1 - .25), 1e-12)
    bounds = _biellipsoid_slice_ylimits(.5, 0, 0, 0, 1., 1., 1., 1.)
    assert bounds[0] == approx(-np.sqrt(1 - .25), 1e-12)
    assert bounds[1] == approx(np.sqrt(1 - .25), 1e-12)

    # Simple Ellipses
    bounds = _biellipsoid_slice_ylimits(0., 90., 0, 0, 1., 2., 2., 1.)
    assert bounds == approx([-2., 1.])
    bounds = _biellipsoid_slice_ylimits(-.5, 0., 0, 0, 1., 1., 2., 1.)
    assert bounds == approx([-2 * np.sqrt(1 - .25), 2 * np.sqrt(1 - .25)], 1e-12)
    bounds = _biellipsoid_slice_ylimits(0., 90., 0, 0, 1., 2., 1., 1.)
    assert bounds == approx([-2., 1.])
    bounds = _biellipsoid_slice_ylimits(-.2, 90., 0, 0, 1., 2., 1., 1.)
    assert bounds[0] == approx(-2 * np.sqrt(1 - .2*.2), 1e-12)
    assert bounds[1] == approx(np.sqrt(1 - .2*.2), 1e-12)
