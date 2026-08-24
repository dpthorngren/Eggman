import numpy as np
from pytest import approx, raises

import eggman


def test_phase_trivial():
    p = eggman.PlanetSystem()
    p.add_star("quadratic_limb", [0.2, 0.1])
    p.add_planet(.1, .11, .08, 1.3, 8., 5., inclination=89.)
    assert p.get_n_objects() == 2
    p.set_time(0.)

    # Check parent star struct
    orbit, shape, light, rotate = p[0]
    assert not rotate
    assert orbit.semimajor == approx(0., abs=1e-12)
    assert shape.position == approx([0., 0., 0.], abs=1e-6)
    assert shape.radii == approx([1., 1., 1., 1.], abs=1e-6)
    assert light.source_type == "quadratic_limb"
    assert light.source_params == approx([1 / np.pi, .2, .1], abs=1e-12)
    assert light.get_brightness_sphere(0., 0.) > (1 / np.pi)
    assert light.get_brightness_sphere(0.99, 0.) < (1 / np.pi)

    # Check planet structs
    orbit, shape, light, rotate = p[1]
    assert rotate
    assert orbit.semimajor == approx(5., abs=1e-12)
    assert shape.position == approx(
        [0., 5 * np.cos(89 * np.pi / 180), 5 * np.sin(89 * np.pi / 180)], abs=1e-6)
    assert shape.radii == approx([.1, .11, .08, 1.3], abs=1e-6)
    assert light.source_type == "no_emission"
    bell = eggman.Shape(0., 0., 0., 0., 0., 0., 1., 1., 1., 1.)
    assert light.get_brightness(1., 0., bell) == 0.
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


def test_rings():
    # Lambertian star, face-on-ring exact check
    p = eggman.PlanetSystem()
    p.add_star("Lambertian", [1.0 / np.pi])
    p.add_planet(.1, .11, .08, 1.3, 8., 5.)
    p.add_ring(.16, .12, 1)
    assert p[0][1].area == approx(np.pi, abs=1e-12)
    assert p[1][1].area == approx(np.pi * .08 * (.1+.11) / 2.0, abs=1e-12)
    assert p[2][1].area == approx(np.pi * (.16**2 - .12**2), abs=1e-12)
    assert p.get_n_objects() == 3
    t = np.array([1e-5, 4.0])
    depths = p.phase_curve_integral(t)
    assert depths[0] == approx(1 - (p[1][1].area + p[2][1].area) / p[0][1].area, abs=1e-7)
    assert depths[1] == approx(1, abs=1e-6)

    # Quadratic limb plausibility check
    p = eggman.PlanetSystem()
    p.add_star("quadratic_limb", [0.2, 0.1])
    p.add_planet(.1, .11, .08, 1.3, 8., 5., inclination=89.)
    assert p.get_n_objects() == 2
    t = np.array([0.1, 4.0])
    depths = p.phase_curve_integral(t)
    assert depths[0] < 1
    assert depths[1] == approx(1, abs=1e-6)
    p.add_ring(.16, .12, 1, gamma=60)
    depths2 = p.phase_curve_integral(t)
    assert depths2[0] < depths[0]
    assert depths2[1] == approx(depths[1], abs=1e-6)
