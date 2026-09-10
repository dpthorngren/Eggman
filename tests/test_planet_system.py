import numpy as np
from pytest import approx, raises

import eggman


def test_system_integration_trivial():
    p = eggman.PlanetSystem(max_steps=400)
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

    # Test integral
    times = np.linspace(0, 8, 5)
    result = p.integrate(times)
    assert len(result) == len(times)
    assert result[0] < 1.
    assert result[1] == approx(1., rel=1e-7)
    assert result[2] == approx(1., rel=1e-7)
    assert result[3] == approx(1., rel=1e-7)
    assert result[0] == approx(result[-1], rel=1e-7)


def test_phase_curve():
    p = eggman.PlanetSystem()
    p.add_star("quadratic_limb", [0.2, 0.1])
    planet_light = eggman.LightSource('emission_map', [100, 100])
    result = planet_light.set_emission_map(lambda x, y, z: 1e-2 + 1e-3*z*z)
    for i, res in enumerate(result):
        loc = planet_light.get_emission_location(i)
        assert res == approx(1e-2 + 1e-3 * loc[2]**2)
    p.add_planet(.12, .11, .09, .115, 2., 5., source=planet_light)
    times = np.linspace(0, 2, 50)
    result = p.integrate(times)
    assert np.all(np.isfinite(result))


def test_secondary_transits():
    p = eggman.PlanetSystem(max_steps=200)
    p.add_star("quadratic_limb", [0.2, 0.1])
    planet_light = eggman.LightSource('lambertian', [.01])
    p.add_planet(.1, .1, .1, .1, 8, 4., inclination=87.1, theta=23, source=planet_light)
    assert p.get_n_objects() == 2

    times = np.linspace(2, 6, 1000)
    result = p.integrate(times)
    l_planet = np.pi * .1**2 * .01
    for t, r, in zip(times, result):
        if abs(t - 4.0) > .35:
            # Out of transit
            assert r == approx(1 + l_planet)
        elif abs(t - 4.0) < .28:
            # Behind star
            assert r == approx(1.0)
        else:
            # Secondary ingress / egress
            assert r < 1. + l_planet, f"{t-4=}, {r=}"
            assert r > 1., f"{t-4=}, {r=}"


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
    depths = p.integrate(t)
    assert depths[0] == approx(1 - (p[1][1].area + p[2][1].area) / p[0][1].area, abs=1e-7)
    assert depths[1] == approx(1, abs=1e-6)

    # Quadratic limb plausibility check
    p = eggman.PlanetSystem()
    p.add_star("quadratic_limb", [0.2, 0.1])
    p.add_planet(.1, .11, .08, 1.3, 8., 5., inclination=89.)
    assert p.get_n_objects() == 2
    t = np.array([0.1, 4.0])
    depths = p.integrate(t)
    assert depths[0] < 1
    assert depths[1] == approx(1, abs=1e-6)
    p.add_ring(.16, .12, 1, gamma=60)
    depths2 = p.integrate(t)
    assert depths2[0] < depths[0]
    assert depths2[1] == approx(depths[1], abs=1e-6)
