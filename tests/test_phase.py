import numpy as np
from pytest import approx, raises

import eggman


def test_phase_trivial():
    p = eggman.PhaseIntegratorWrap()
    p.add_star("quadratic", [0.2, 0.1])
    p.add_planet(.1, .11, .08, 1.3, 0., 0., 0., 0., 8., 5., 89., 0.)
    assert p.get_n_objects() == 2
    p.set_time(0.)

    # Check parent star struct
    orbit, shape, light, rotate = p[0]
    assert not rotate
    assert orbit.semimajor == approx(0., abs=1e-12)
    assert shape.position == approx([0., 0., 0.], abs=1e-6)
    assert shape.radii == approx([1., 1., 1., 1.], abs=1e-6)
    assert light.source_type == "uniform"
    assert light.limb_type == "quadratic"
    assert light.limb_params == approx([.2, .1], abs=1e-12)
    assert light.get_brightness_sphere(0., 0.) > (1 / np.pi)
    assert light.get_brightness_sphere(0.99, 0.) < (1 / np.pi)

    # Check planet structs
    orbit, shape, light, rotate = p[1]
    assert rotate
    assert orbit.semimajor == approx(5., abs=1e-12)
    assert shape.position == approx(
        [0., 5 * np.cos(89 * np.pi / 180), 5 * np.sin(89 * np.pi / 180)], abs=1e-6)
    assert shape.radii == approx([.1, .11, .08, 1.3], abs=1e-6)
    assert light.source_type == "none"
    assert light.get_brightness(1., 0., 0.) == 0.
    with raises(IndexError):
        p[2]

    # Test single integrations
    p.set_time(0.)
    assert p.integrate_single(1) == 0.
    assert p.integrate_single(0) < 1
    p.set_time(2.)
    assert p.integrate_single(1) == 0.
    assert p.integrate_single(0) == approx(1., rel=1e-7)

    # Test phase integral
    times = np.linspace(0, 8, 5)
    result = p.phase_curve_integral(times)
    assert len(result) == len(times)
    assert result[0] < 1.
    assert result[1] == approx(1., rel=1e-7)
    assert result[2] == approx(1., rel=1e-7)
    assert result[3] == approx(1., rel=1e-7)
    assert result[0] == approx(result[-1], rel=1e-7)
