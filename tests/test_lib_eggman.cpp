#include "light_source.hpp"
#include "math_utils.hpp"
#include "orbit.hpp"
#include "planet_system.hpp"
#include "shape.hpp"
#include "transit_integral.hpp"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>

using namespace std;

#define DRAND(a, b) (drand48() * (b - a) + a)
#define ANNOUNCE_TEST() cout << "Running " << __func__ << "." << endl;
#define TEST_ASSERT(A, OP, B, out)                                                                 \
    if (!((A)OP(B))) {                                                                             \
        cout << "ERROR in " << __func__ << " at line " << __LINE__ << endl;                        \
        cout << "  Assertion \"" << #A << #OP << #B << "\" (" << (A) << #OP << (B)                 \
             << ") is false." << endl;                                                             \
        out += 1;                                                                                  \
    }
#define TEST_APPROX(A, B, atol, out)                                                               \
    if (fabs((A) - (B)) > atol) {                                                                  \
        cout << "  ERROR in " << __func__ << " at line " << __LINE__ << ":" << endl;               \
        cout << "    Assertion \"" << #A << "~=" << #B << "\" (" << (A) << "~=" << (B)             \
             << ") is false." << endl;                                                             \
        out += 1;                                                                                  \
    }

int test_biellipsoid() {
    ANNOUNCE_TEST();
    int errors = 0;
    Shape shp = Shape(.15, .1, .1, .1);
    TEST_APPROX(shp.x_bounds().min, -.1, 1e-12, errors);
    TEST_APPROX(shp.x_bounds().max, .15, 1e-12, errors);
    TEST_APPROX(shp.y_bounds().min, -.1, 1e-12, errors);
    TEST_APPROX(shp.y_bounds().max, .1, 1e-12, errors);
    shp.set_position((Vec3){1.0, 0.1, 0.3});
    TEST_APPROX(shp.x_bounds().min, 1 - .1, 1e-12, errors);
    TEST_APPROX(shp.x_bounds().max, 1 + .15, 1e-12, errors);
    TEST_APPROX(shp.y_bounds().min, .1 + -.1, 1e-12, errors);
    TEST_APPROX(shp.y_bounds().max, .1 + .1, 1e-12, errors);
    TEST_ASSERT(shp.shape_type, ==, Biellipsoid, errors);

    // Test visibility
    shp = Shape(1.33, 1.33, 1.33, 1.33);
    TEST_ASSERT(shp.is_visible({0, 0, 0}), ==, false, errors);
    TEST_ASSERT(shp.is_visible({0, 0, -10}), ==, false, errors);
    TEST_ASSERT(shp.is_visible({1., 1., 0.}), ==, true, errors);
    TEST_ASSERT(shp.is_visible({.7, .7, 0.}), ==, false, errors);
    TEST_ASSERT(shp.is_visible({0, 0, 10}), ==, true, errors);

    // Test line projection
    shp = Shape(1., 1., 1., 1.);
    TEST_ASSERT(shp.shape_type, ==, Sphere, errors);
    Vec3 loc;
    double mu;
    bool hit = shp.raycast(0, 0, &mu, &loc);
    TEST_ASSERT(hit, ==, true, errors);
    TEST_APPROX(loc.x, 0., 1e-9, errors);
    TEST_APPROX(loc.y, 0., 1e-9, errors);
    TEST_APPROX(loc.z, 1., 1e-9, errors);
    TEST_APPROX(mu, 1.0, 1e-9, errors);
    hit = shp.raycast(0.5, 0, &mu, &loc);
    TEST_ASSERT(hit, ==, true, errors);
    TEST_APPROX(loc.x, 0.5, 1e-9, errors);
    TEST_APPROX(loc.y, 0., 1e-9, errors);
    TEST_APPROX(loc.z, sqrt(1 - .5 * .5), 1e-9, errors);
    shp.raycast(0.5, 0., &mu, &loc);
    TEST_APPROX(0.5, loc.x, 1e-9, errors);
    TEST_APPROX(0, loc.y, 1e-9, errors);
    TEST_APPROX(mu, loc.z, 1e-9, errors);
    TEST_APPROX(mu, cos(M_PI / 6), 1e-9, errors);
    shp.raycast(1 - 1e-12, 0., &mu);
    TEST_APPROX(mu, 0., 1e-4, errors);
    // Aligned but no longer spherical
    shp.set_radii(1.5, 1.4, 1.0, 0.9);
    hit = shp.raycast(0.5, 0.0, &mu, &loc);
    TEST_ASSERT(hit, ==, true, errors);
    TEST_APPROX(loc.x, .5 / 1.5, 1e-9, errors);
    TEST_APPROX(loc.y, 0., 1e-9, errors);
    TEST_ASSERT(loc.z, >, 0, errors);
    hit = shp.raycast(-0.5, 0.0, &mu, &loc);
    TEST_ASSERT(hit, ==, true, errors);
    TEST_APPROX(loc.x, -.5 / 1.4, 1e-9, errors);
    TEST_APPROX(loc.y, 0.0, 1e-9, errors);
    TEST_ASSERT(loc.z, >, 0, errors);
    TEST_ASSERT(mu, >, 0, errors);
    TEST_ASSERT(mu, <, 1, errors);
    hit = shp.raycast(0.0, 0.0, &mu, &loc);
    TEST_APPROX(mu, 1.0, 1e-9, errors);

    // Test area calculation
    shp = Shape(1., 1.5, 1.33, 1.);
    shp.set_position({1., 1., 1.});
    TEST_APPROX(shp.get_area(), M_PI * 1.33 * (1.5 + 1.) / 2.0, 1e-9, errors);
    shp.set_rotation(0., 90., 0.);
    TEST_APPROX(shp.get_area(), M_PI * 1.33, 1e-9, errors);
    for (int i = 0; i < 10; i++) {
        shp.set_rotation(i * 180. / 10., 0., 0.);
        TEST_APPROX(shp.get_area(), M_PI * 1.33 * (1.5 + 1.) / 2.0, 1e-9, errors);
    }

    // Test forward/backward determination
    shp = Shape(2., 1.5, 1., 0.9);
    shp.set_position({4., 3., 2.});
    shp.set_rotation(0., 10., 0.);
    TEST_ASSERT(shp.is_forward_2d(3, 3.), ==, false, errors);
    TEST_ASSERT(shp.is_forward_2d(4, 3.), ==, false, errors);
    TEST_ASSERT(shp.is_forward_2d(5.3, 3.), ==, true, errors);
    shp.set_rotation(0., -10., 0.);
    TEST_ASSERT(shp.is_forward_2d(3, 3.2), ==, false, errors);
    TEST_ASSERT(shp.is_forward_2d(4, 3.), ==, true, errors);
    TEST_ASSERT(shp.is_forward_2d(5., 2.95), ==, true, errors);
    shp.set_rotation(89., -10., 0.);
    TEST_ASSERT(shp.is_forward_2d(4., 2.), ==, false, errors);
    TEST_ASSERT(shp.is_forward_2d(4., 4.), ==, true, errors);
    shp.set_rotation(30., 0., 0.);
    TEST_ASSERT(shp.is_forward_2d(3, 3.), ==, false, errors);
    TEST_ASSERT(shp.is_forward_2d(4.1, 3.1), ==, true, errors);

    // Test raycasts against random shapes
    double x, y;
    bool forward;
    for (int i = 0; i < 100; i++) {
        shp = Shape(DRAND(.5, 2), DRAND(.5, 2), DRAND(.5, 2), DRAND(.5, 2));
        shp.set_position({DRAND(-5, 5), DRAND(-5, 5), DRAND(-5, 5)});
        shp.set_rotation(DRAND(0, 179), DRAND(0, 179), DRAND(0, 179));
        for (int j = 0; j < 100; j++) {
            x = DRAND(-2, 2) + shp.position.x;
            y = DRAND(-2, 2) + shp.position.y;
            hit = shp.raycast(x, y, &mu, &loc);
            TEST_ASSERT(hit, ==, shp.line_intersects(x, y), errors);
        }
        if (hit) {
            TEST_APPROX(loc.x * loc.x + loc.y * loc.y + loc.z * loc.z, 1.0, 1e-12, errors);
            forward = loc.x > 0.;
            loc = shp.sphere_to_world(loc);
            TEST_APPROX(loc.x, x, 1e-12, errors);
            TEST_APPROX(loc.y, y, 1e-12, errors);
            TEST_ASSERT(forward, ==, shp.is_forward_2d(loc.x, loc.y), errors);
            TEST_ASSERT(forward, ==, shp.is_forward(loc), errors);
        }
    }

    return errors;
}

int test_rings() {
    ANNOUNCE_TEST();
    int errors = 0;
    Shape shp = Shape(1.0, 0.7, 0., 0.);
    Bounds bottom, top;
    top = shp.slice_ylimits(-0.3, &bottom, 0);
    TEST_APPROX(bottom.min, -sqrt(1. - 0.09), 1e-9, errors);
    TEST_APPROX(bottom.max, -sqrt(.7 * .7 - 0.09), 1e-9, errors);
    TEST_APPROX(top.min, sqrt(.7 * .7 - 0.09), 1e-9, errors);
    TEST_APPROX(top.max, sqrt(1. - 0.09), 1e-9, errors);
    TEST_ASSERT(shp.raycast(0., 0., nullptr, nullptr), ==, false, errors);
    TEST_ASSERT(shp.raycast(0.3, 0.23, nullptr, nullptr), ==, false, errors);
    TEST_ASSERT(shp.raycast(10., 30., nullptr, nullptr), ==, false, errors);
    TEST_ASSERT(shp.raycast(0.9, 0.1, nullptr, nullptr), ==, true, errors);
    TEST_ASSERT(shp.raycast(-0.1, 0.9, nullptr, nullptr), ==, true, errors);
    TEST_ASSERT(shp.raycast(-.8 / sqrt(2.), .8 / sqrt(2.), nullptr, nullptr), ==, true, errors);
    TEST_APPROX(shp.get_area(), M_PI * (1 - .7 * .7), 1e-9, errors);

    shp.set_rotation(0., 0., 60.);
    top = shp.slice_ylimits(-0.3, &bottom, 0);
    TEST_APPROX(bottom.min, -0.5 * sqrt(1. - 0.09), 1e-9, errors);
    TEST_APPROX(bottom.max, -0.5 * sqrt(.7 * .7 - 0.09), 1e-9, errors);
    TEST_APPROX(top.min, 0.5 * sqrt(.7 * .7 - 0.09), 1e-9, errors);
    TEST_APPROX(top.max, 0.5 * sqrt(1. - 0.09), 1e-9, errors);
    top = shp.slice_ylimits(-0.3, &bottom, -1);
    TEST_APPROX(bottom.min, -0.5 * sqrt(1. - 0.09), 1e-9, errors);
    TEST_APPROX(bottom.max, -0.5 * sqrt(.7 * .7 - 0.09), 1e-9, errors);
    TEST_APPROX(top.min, 0., 1e-9, errors);
    TEST_APPROX(top.max, 0., 1e-9, errors);
    top = shp.slice_ylimits(-0.3, &bottom, 1);
    TEST_APPROX(bottom.min, 0., 1e-9, errors);
    TEST_APPROX(bottom.max, 0., 1e-9, errors);
    TEST_APPROX(top.min, 0.5 * sqrt(.7 * .7 - 0.09), 1e-9, errors);
    TEST_APPROX(top.max, 0.5 * sqrt(1. - 0.09), 1e-9, errors);
    TEST_ASSERT(shp.raycast(0., 0.8, nullptr, nullptr), ==, false, errors);
    TEST_APPROX(shp.get_area(), 0.5 * M_PI * (1 - .7 * .7), 1e-9, errors);
    return errors;
}

int test_orbital_position() {
    ANNOUNCE_TEST();
    int errors = 0;
    Shape shp = Shape(.2, .18, .15, .23);
    Orbit orb = Orbit(2., 0., 5., 0., 85., 90.);
    shp.position_from_orbit(0.0, orb);
    TEST_APPROX(shp.position.x, 0., 1e-12, errors);
    TEST_APPROX(shp.position.y, 5. * cos(85 * M_PI / 180), 1e-12, errors);
    TEST_APPROX(shp.position.z, 5. * sin(85 * M_PI / 180), 1e-12, errors);
    shp.position_from_orbit(0.5, orb);
    TEST_APPROX(shp.position.x, 5., 1e-12, errors);
    TEST_APPROX(shp.position.y, 0., 1e-12, errors);
    TEST_APPROX(shp.position.z, 0., 1e-12, errors);
    return errors;
}

int test_transit() {
    ANNOUNCE_TEST();
    int errors = 0;
    Orbit orb = Orbit(2., 0., 5., 0.01, 85., 90.);
    double source_params[MAX_SOURCE_PARAMS] = {1.0 / M_PI, 0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0,        0.0, 0.0, 0.0, 0.0, 0.0};
    LightSource star = LightSource(QuadraticLimb, source_params);

    const int n_times = 100;
    double outputs[n_times];
    double times[n_times];
    for (int i = 0; i < n_times; i++) {
        times[i] = i * orb.get_period() / n_times;
        outputs[i] = -1.0;
    }

    // No limb darkening, sphere
    transit_integral(times, outputs, n_times, orb, star, 0., 0., 0., .1, .1, .1, .1, true);
    for (int i = 0; i < n_times; i++) { // General bounds
        TEST_ASSERT(outputs[i], <=, 1., errors);
        TEST_ASSERT(outputs[i], >=, 0., errors);
    }
    // Out of transit
    TEST_APPROX(outputs[n_times / 4], 1.0, 1e-12, errors);
    TEST_APPROX(outputs[n_times / 2], 1.0, 1e-12, errors);
    // In-transit (pi omitted because star brightness is 1/pi)
    TEST_APPROX(outputs[0], 1.0 - .1 * .1, 1e-7, errors);
    TEST_ASSERT(outputs[1], <, 1.0, errors);
    TEST_ASSERT(outputs[n_times - 2], <, 1.0, errors);

    // No limb darkening, oblate spheroid
    orb = Orbit(2., 0., 5., 0.01, 90., 90.);
    transit_integral(times, outputs, n_times, orb, star, 0., 0., 0., .11, .11, .08, .11, true);
    // Out of transit
    TEST_APPROX(outputs[n_times / 4], 1.0, 1e-12, errors);
    TEST_APPROX(outputs[n_times / 2], 1.0, 1e-12, errors);
    // In-transit (pi omitted because star brightness is 1/pi)
    TEST_APPROX(outputs[0], 1.0 - .11 * .08, 1e-7, errors);
    TEST_ASSERT(outputs[1], <, 1.0, errors);
    TEST_ASSERT(outputs[n_times - 2], <, 1.0, errors);

    // No limb darkening, asymmetric transit
    orb = Orbit(2., 0., 5., 0.01, 90., 90.);
    transit_integral(times, outputs, n_times, orb, star, 0., 0., 0., .11, .09, -1., .1, false);
    // Out of transit
    TEST_APPROX(outputs[n_times / 4], 1.0, 1e-12, errors);
    TEST_APPROX(outputs[n_times / 2], 1.0, 1e-12, errors);
    // In-transit (pi omitted because star brightness is 1/pi)
    TEST_APPROX(outputs[0], 1.0 - (.11 * .11 + .09 * .09) / 2.0, 1e-7, errors);
    TEST_ASSERT(outputs[1], <, 1.0, errors);
    TEST_ASSERT(outputs[n_times - 2], <, 1.0, errors);

    // Limb Darkening, symmetric transit (ref from catwoman)
    orb = Orbit(1., 0., 15., 0, 90., 90.);
    source_params[1] = .1;
    source_params[2] = .3;
    double times2[2] = {.001, .01};
    star = LightSource(QuadraticLimb, source_params);
    transit_integral(times2, outputs, 2, orb, star, 0., 0., 0., .1, .1, .1, .1, true);
    TEST_APPROX(outputs[0], 0.989098764152, 1e-7, errors);
    TEST_APPROX(outputs[1], 0.992627976697, 1e-7, errors);

    // Limb Darkening, asymmetric transit (ref from catwoman)
    star = LightSource(QuadraticLimb, source_params);
    orb = Orbit(1., 0., 15., 0, 90., 90.);
    transit_integral(times2, outputs, 2, orb, star, 0., 0., 0., .11, .1, -1, .1, true);
    TEST_APPROX(outputs[0], 0.987955283022, 1e-7, errors);
    TEST_APPROX(outputs[1], 0.992339221135, 1e-7, errors);

    // Limb Darkening, asymmetric transit, slightly inclined (ref from catwoman)
    star = LightSource(QuadraticLimb, source_params);
    orb = Orbit(1., 0., 15., 0, 89., 90.);
    transit_integral(times2, outputs, 2, orb, star, 0., 0., 0., .11, .09, -1, .1, true);
    // Note the reduced precision.  I *think* this is Catwoman's fault, but it's hard to tell.
    TEST_APPROX(outputs[0], 0.989036885966, 1e-6, errors);
    TEST_APPROX(outputs[1], 0.995337628970, 1e-6, errors);
    return errors;
}

int test_planetary_system() {
    ANNOUNCE_TEST();
    int errors = 0;
    PlanetSystem p;

    // Add the star
    Orbit orb = Orbit();
    double source_params[MAX_SOURCE_PARAMS] = {1.0 / M_PI, 0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0,        0.0, 0.0, 0.0, 0.0, 0.0};
    Shape shp = Shape();
    LightSource star = LightSource(QuadraticLimb, source_params);
    TEST_APPROX(star.limb_norm, 1.0, 1e-9, errors)
    TEST_APPROX(star.get_brightness(0., 0., shp), 1 / M_PI, 1e-9, errors)
    p.add_object(orb, shp, star);
    TEST_ASSERT(p.get_n_objects(), ==, 1, errors);
    TEST_APPROX(p.integrate_single(0), 1.0, 1e-7, errors);

    // Add the planet
    orb = Orbit(10., 0., 5., 0.00, 89., 90.);
    shp = Shape(.12, .09, .08, .08);
    source_params[0] = 1e-6;
    LightSource planet = LightSource(Lambertian, source_params);
    TEST_APPROX(star.limb_norm, 1.0, 1e-9, errors);
    TEST_APPROX(planet.get_brightness(0., 0., shp), 1e-6, 1e-9, errors)
    p.add_object(orb, shp, planet);
    TEST_ASSERT(p.get_n_objects(), ==, 2, errors);

    // Test sources individually
    p.set_time(2.5);
    TEST_APPROX(p.shapes[0].position.x, 0.0, 1e-9, errors);
    TEST_APPROX(p.shapes[0].position.y, 0.0, 1e-9, errors);
    TEST_APPROX(p.shapes[0].position.z, 0.0, 1e-9, errors);
    TEST_APPROX(p.shapes[1].position.x, 5.0, 1e-9, errors);
    TEST_APPROX(p.shapes[1].position.y, 0.0, 1e-9, errors);
    TEST_APPROX(p.shapes[1].position.z, 0.0, 1e-9, errors);
    double r = p.integrate_single(0);
    TEST_APPROX(r, 1.0, 1e-7, errors);
    r = p.integrate_single(1);
    TEST_APPROX(r, .12 * .09 * 1e-6, 1e-7, errors);

    // Test a full phase curve
    const int n_times = 4;
    double time[n_times];
    double result[n_times];
    for (int i = 0; i < n_times; i++) {
        time[i] = 10. * i / n_times;
        result[i] = -1.;
    }
    p.phase_curve_integral(time, result, n_times);
    for (int i = 0; i < n_times; i++) {
        TEST_ASSERT(result[i], >, 0., errors);
        TEST_ASSERT(result[i], <, 2., errors);
    }
    TEST_ASSERT(result[0], <, 1, errors);
    TEST_ASSERT(result[1], >, 1, errors);
    TEST_APPROX(result[2], 1, 1e-6, errors);
    TEST_ASSERT(result[3], >, 1, errors);
    TEST_APPROX(result[3], 1 + .12 * .09 * 1e-6, 1e-6, errors);
    TEST_APPROX(result[2], result[1], 1e-6, errors);

    // Test star with quadratic limb-darkening
    source_params[1] = 0.2;
    source_params[2] = 0.1;
    source_params[0] = 1 / M_PI;
    shp = Shape();
    star = LightSource(QuadraticLimb, source_params);
    double expect = 1. - .2 / 3 - .1 / 6.;
    TEST_APPROX(star.limb_norm, expect, 1e-9, errors);
    double nu = 1.0 - sqrt(1 - 0.5 * 0.5);
    expect = (1.0 - .2 * nu - .1 * nu * nu) / (expect * M_PI);
    TEST_APPROX(star.get_brightness_sphere(0.5, 0.), expect, 1e-9, errors);
    p.clear_objects();
    p.add_object(Orbit(), shp, star);
    TEST_ASSERT(p.get_n_objects(), ==, 1, errors);
    TEST_APPROX(p.integrate_single(0), 1.0, 1e-7, errors);

    // Test full day-night phase curve
    p.clear_objects();
    star = LightSource(QuadraticLimb, source_params);
    shp = Shape();
    p.add_object(Orbit(), shp, star);
    source_params[0] = 1e-1;
    source_params[1] = 1e-2;
    orb = Orbit(10., 0., 5., 0.00, 89., 90.);
    planet = LightSource(DayNight, source_params);
    shp = Shape(0.1, 0.1, 0.1, 0.1);
    p.add_object(orb, shp, planet, true);
    TEST_ASSERT(p.get_n_objects(), ==, 2, errors);
    p.phase_curve_integral(time, result, n_times);
    TEST_ASSERT(result[0], <, 1, errors);
    TEST_APPROX(result[1], 1 + M_PI * .1 * .1 * 1.1e-1 / 2., 1e-9, errors);
    TEST_APPROX(result[2], 1, 1e-9, errors);
    TEST_APPROX(result[3], result[1], 1e-9, errors);

    // Test the bounds sorting mechanism.
    Bounds b[5] = {
        {-1, 1.}, {-3, -0.5}, {-.1, 0.1}, {-4.0, -0.7}, {0.9, 1.5},
    };
    int n_bounds = process_bounds(b, 5);
    TEST_ASSERT(n_bounds, ==, 2, errors);
    TEST_APPROX(b[0].min, -0.5, 1e-12, errors);
    TEST_APPROX(b[0].max, -0.1, 1e-12, errors);
    TEST_APPROX(b[1].min, 0.1, 1e-12, errors);
    TEST_APPROX(b[1].max, 0.9, 1e-12, errors);

    // Test the inverse bounds sorting
    b[0] = {-1, 1.};
    b[1] = {-3, -0.5};
    b[2] = {-.1, 0.1};
    b[3] = {-4.0, -0.7};
    b[4] = {0.9, 1.5};
    n_bounds = process_bounds(b, 5, true);
    TEST_ASSERT(n_bounds, ==, 3, errors);
    TEST_APPROX(b[0].min, -1, 1e-12, errors);
    TEST_APPROX(b[0].max, -0.5, 1e-12, errors);
    TEST_APPROX(b[1].min, -0.1, 1e-12, errors);
    TEST_APPROX(b[1].max, 0.1, 1e-12, errors);
    TEST_APPROX(b[2].min, 0.9, 1e-12, errors);
    TEST_APPROX(b[2].max, 1.0, 1e-12, errors);
    // Again but without occluders on the edges
    b[0] = {-1, 1.};
    b[1] = {-0.5, 0.1};
    b[2] = {0.3, 0.4};
    n_bounds = process_bounds(b, 3, true);
    TEST_ASSERT(n_bounds, ==, 2, errors);
    TEST_APPROX(b[0].min, -.5, 1e-12, errors);
    TEST_APPROX(b[0].max, 0.1, 1e-12, errors);
    TEST_APPROX(b[1].min, 0.3, 1e-12, errors);
    TEST_APPROX(b[1].max, 0.4, 1e-12, errors);
    return errors;

    // Test a ring system by comparison with planets
    double d_inner, d_outer, d_ring;
    p.clear_objects();
    source_params[0] = 1 / M_PI;
    source_params[1] = 0.22;
    source_params[2] = 0.17;
    star = LightSource(QuadraticLimb, source_params);
    p.add_object(Orbit(), Shape(), star);
    orb = Orbit(10., 0., 5., 0.00, 89., 90.);
    LightSource nosource = LightSource(NoEmission, source_params);
    p.add_object(orb, Shape(.2 * cos(M_PI * 70), .2, .2, .2), nosource);
    p.set_time(.1);
    d_inner = p.integrate_single(0);
    p.shapes[1].set_radii(.05 * cos(M_PI * 70), .05, .05, .05);
    d_outer = p.integrate_single(0);
    TEST_ASSERT(d_outer, >, d_inner, errors);
    Shape ring = Shape(.2, .05, 0., 0.);
    ring.set_rotation(0., 0., 60.);
    p.clear_objects();
    p.add_object(Orbit(), Shape(), star);
    p.add_object(orb, ring, nosource);
    d_ring = p.integrate_single(0);
    TEST_APPROX(d_ring, 1 - ((1 - d_inner) - (1 - d_outer)), 1e-6, errors);
}

int test_light_source() {
    ANNOUNCE_TEST();
    int errors = 0;
    double source_params[MAX_SOURCE_PARAMS] = {1.0 / M_PI, 0.01, 0.1, 0., 0., 0.,
                                               0.,         0.,   0.,  0., 0., 0.};
    // No source should always return 0
    Shape shp = Shape();
    LightSource s = LightSource(NoEmission, source_params);
    TEST_APPROX(s.get_brightness_sphere(0., 0.), 0., 1e-9, errors);
    TEST_APPROX(s.get_brightness_sphere(0.5, 0.5), 0., 1e-9, errors);
    TEST_APPROX(s.get_brightness(0.9, 0.8, shp), 0., 1e-9, errors);
    // Lambertian emitter of total brightness 1
    double expect = 1 / M_PI;
    s = LightSource(Lambertian, source_params);
    TEST_APPROX(s.get_brightness_sphere(0., 0.), expect, 1e-9, errors);
    TEST_APPROX(s.get_brightness_sphere(0.5, 0.5), expect, 1e-9, errors);
    TEST_APPROX(s.get_brightness(0.9, 0.8, shp), 0., 1e-9, errors);
    // Standard star (brightness 1, quadratic limb darkening)
    s = LightSource(QuadraticLimb, source_params);
    // Small params ~= lambertian
    TEST_APPROX(s.get_brightness_sphere(0.5, 0.), expect, 1e-1, errors);
    // Center is brighter, limb is darker
    double center = s.get_brightness_sphere(0., 0.);
    double edge = s.get_brightness_sphere(0.9, 0.);
    TEST_ASSERT(center, >, expect, errors);
    TEST_ASSERT(edge, <, expect, errors);
    // Now with more reasonable limb params
    source_params[1] = 0.2;
    source_params[2] = 0.1;
    s = LightSource(QuadraticLimb, source_params);
    TEST_ASSERT(s.get_brightness_sphere(0., 0.), >, center, errors);
    TEST_ASSERT(s.get_brightness_sphere(0.9, 0.), <, edge, errors);
    double norm = (1. - .2 / 3 - .1 / 6);
    TEST_APPROX(s.limb_norm, norm, 1e-9, errors);
    expect = 1.0 - sqrt(1 - 0.7 * 0.7);
    expect = (1.0 - .2 * expect - .1 * expect * expect) / (norm * M_PI);
    TEST_APPROX(s.get_brightness_sphere(0.7, 0.0), expect, 1e-9, errors)
    // Testing day-night emission type
    source_params[0] = 1.5;
    source_params[1] = 0.3;
    s = LightSource(DayNight, source_params);
    TEST_APPROX(s.get_brightness_sphere(0., 0.), 0.3, 1e-12, errors);
    TEST_APPROX(s.get_integrated_brightness(shp), M_PI * 0.3, 1e-9, errors);
    shp.set_rotation(0., 45., 0.);
    TEST_APPROX(s.get_integrated_brightness(shp), M_PI_2 * (1.8 - 1.2 * sqrt(0.5)), 1e-9, errors);
    TEST_APPROX(s.get_brightness(-.9, 0, shp), 1.5, 1e-12, errors);
    TEST_APPROX(s.get_brightness(.9, 0, shp), 0.3, 1e-12, errors);
    TEST_APPROX(s.get_brightness(5, 5, shp), 0., 1e-12, errors);
    TEST_APPROX(s.get_brightness(-.707, -.707, shp), 1.5, 1e-12, errors);
    TEST_APPROX(s.get_brightness(.707, -.707, shp), 0.3, 1e-12, errors);
    shp.set_rotation(0., 90., 0.);
    TEST_APPROX(s.get_integrated_brightness(shp), M_PI * 0.9, 1e-9, errors);
    shp.set_rotation(0., 180., 0.);
    TEST_APPROX(s.get_integrated_brightness(shp), M_PI * 1.5, 1e-9, errors);
    return errors;
}

int main() {
    // Note: This tester is not as rigorous as the Pytest tests.  It is intended mainly to catch
    // memory management errors and distinguish between C++ errors and Cython wrapping errors.
    cout << "===== Running C++ Tester =====" << endl;
    int seed = time(NULL);
    cout << "Random seed: " << seed << endl;
    srand48(seed);
    int errors = 0;
    errors += test_biellipsoid();
    errors += test_rings();
    errors += test_orbital_position();
    errors += test_transit();
    errors += test_planetary_system();
    errors += test_light_source();
    if (errors == 0) {
        cout << "All tests passed." << endl << endl;
    } else {
        cout << "ERROR: " << errors << " test(s) failed." << endl << endl;
    }
    return 0;
}
