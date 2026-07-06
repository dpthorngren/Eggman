#include "biellipsoid.hpp"
#include "light_source.hpp"
#include "orbit.hpp"
#include "phase_curve_integral.hpp"
#include "transit_integral.hpp"
#include <iostream>

using namespace std;

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
    Biellipsoid bell = Biellipsoid(.15, .1, .1, .1);
    bell.update_derived();
    TEST_APPROX(bell.x_bounds().min, -.1, 1e-12, errors);
    TEST_APPROX(bell.x_bounds().max, .15, 1e-12, errors);
    TEST_APPROX(bell.y_bounds().min, -.1, 1e-12, errors);
    TEST_APPROX(bell.y_bounds().max, .1, 1e-12, errors);
    bell.set_position((Vec3){1.0, 0.1, 0.3});
    TEST_APPROX(bell.x_bounds().min, 1 - .1, 1e-12, errors);
    TEST_APPROX(bell.x_bounds().max, 1 + .15, 1e-12, errors);
    TEST_APPROX(bell.y_bounds().min, .1 + -.1, 1e-12, errors);
    TEST_APPROX(bell.y_bounds().max, .1 + .1, 1e-12, errors);

    // Test line projection
    bell = Biellipsoid(1., 1., 1., 1.);
    bell.update_derived();
    Vec3 loc;
    double mu;
    bool hit = bell.raycast(0, 0, &mu, &loc);
    TEST_ASSERT(hit, ==, true, errors);
    TEST_APPROX(loc.x, 0., 1e-9, errors);
    TEST_APPROX(loc.y, 0., 1e-9, errors);
    TEST_APPROX(loc.z, 1., 1e-9, errors);
    TEST_APPROX(mu, 1.0, 1e-9, errors);
    hit = bell.raycast(0.5, 0, &mu, &loc);
    TEST_ASSERT(hit, ==, true, errors);
    TEST_APPROX(loc.x, 0.5, 1e-9, errors);
    TEST_APPROX(loc.y, 0., 1e-9, errors);
    TEST_APPROX(loc.z, sqrt(1 - .5 * .5), 1e-9, errors);
    bell.raycast(0.5, 0., &mu, &loc);
    TEST_APPROX(mu, loc.z, 1e-9, errors);
    TEST_APPROX(mu, cos(M_PI / 6), 1e-9, errors);
    bell.raycast(1 - 1e-12, 0., &mu);
    TEST_APPROX(mu, 0., 1e-4, errors);

    // Test area calculation
    bell = Biellipsoid(1., 1.5, 1.33, 1.);
    bell.set_position({1., 1., 1.});
    TEST_APPROX(bell.get_area(), M_PI * 1.33 * (1.5 + 1.) / 2.0, 1e-9, errors);
    bell.set_rotation(0., 90., 0.);
    TEST_APPROX(bell.get_area(), M_PI * 1.33, 1e-9, errors);
    for (int i = 0; i < 10; i++) {
        bell.set_rotation(i * 180. / 10., 0., 0.);
        TEST_APPROX(bell.get_area(), M_PI * 1.33 * (1.5 + 1.) / 2.0, 1e-9, errors);
    }
    return errors;
}

int test_orbital_position() {
    ANNOUNCE_TEST();
    int errors = 0;
    Biellipsoid bell = Biellipsoid(.2, .18, .15, .23);
    Orbit orb = Orbit(2., 0., 5., 0., 85., 90.);
    bell.position_from_orbit(0.0, orb);
    TEST_APPROX(bell.position.x, 0., 1e-12, errors);
    TEST_APPROX(bell.position.y, 5. * cos(85 * M_PI / 180), 1e-12, errors);
    TEST_APPROX(bell.position.z, 5. * sin(85 * M_PI / 180), 1e-12, errors);
    bell.position_from_orbit(0.5, orb);
    TEST_APPROX(bell.position.x, 5., 1e-12, errors);
    TEST_APPROX(bell.position.y, 0., 1e-12, errors);
    TEST_APPROX(bell.position.z, 0., 1e-12, errors);
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

int test_phase_curve() {
    ANNOUNCE_TEST();
    int errors = 0;
    PhaseIntegrator p;

    // Add the star
    Orbit orb = Orbit();
    double source_params[MAX_SOURCE_PARAMS] = {1.0 / M_PI, 0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0,        0.0, 0.0, 0.0, 0.0, 0.0};
    Biellipsoid bell = Biellipsoid();
    LightSource star = LightSource(QuadraticLimb, source_params);
    TEST_APPROX(star.limb_norm, 1.0, 1e-9, errors)
    TEST_APPROX(star.get_brightness(0., 0., bell), 1 / M_PI, 1e-9, errors)
    p.add_object(orb, bell, star);
    TEST_ASSERT(p.get_n_objects(), ==, 1, errors);
    TEST_APPROX(p.integrate_single(0), 1.0, 1e-7, errors);

    // Add the planet
    orb = Orbit(10., 0., 5., 0.00, 89., 90.);
    bell = Biellipsoid(.12, .09, .08, .08);
    source_params[0] = 1e-6;
    LightSource planet = LightSource(Lambertian, source_params);
    TEST_APPROX(star.limb_norm, 1.0, 1e-9, errors);
    TEST_APPROX(planet.get_brightness(0., 0., bell), 1e-6, 1e-9, errors)
    p.add_object(orb, bell, planet);
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
    bell = Biellipsoid();
    star = LightSource(QuadraticLimb, source_params);
    double expect = 1. - .2 / 3 - .1 / 6.;
    TEST_APPROX(star.limb_norm, expect, 1e-9, errors);
    double nu = 1.0 - sqrt(1 - 0.5 * 0.5);
    expect = (1.0 - .2 * nu - .1 * nu * nu) / (expect * M_PI);
    TEST_APPROX(star.get_brightness_sphere(0.5, 0.), expect, 1e-9, errors);
    p.clear_objects();
    p.add_object(Orbit(), bell, star);
    TEST_ASSERT(p.get_n_objects(), ==, 1, errors);
    TEST_APPROX(p.integrate_single(0), 1.0, 1e-7, errors);

    return errors;
}

int test_light_source() {
    ANNOUNCE_TEST();
    int errors = 0;
    double source_params[MAX_SOURCE_PARAMS] = {1.0 / M_PI, 0.01, 0.1, 0., 0., 0.,
                                               0.,         0.,   0.,  0., 0., 0.};
    // No source should always return 0
    Biellipsoid bell = Biellipsoid();
    LightSource s = LightSource(None, source_params);
    TEST_APPROX(s.get_brightness_sphere(0., 0.), 0., 1e-9, errors);
    TEST_APPROX(s.get_brightness_sphere(0.5, 0.5), 0., 1e-9, errors);
    TEST_APPROX(s.get_brightness(0.9, 0.8, bell), 0., 1e-9, errors);
    // Lambertian emitter of total brightness 1
    double expect = 1 / M_PI;
    s = LightSource(Lambertian, source_params);
    TEST_APPROX(s.get_brightness_sphere(0., 0.), expect, 1e-9, errors);
    TEST_APPROX(s.get_brightness_sphere(0.5, 0.5), expect, 1e-9, errors);
    TEST_APPROX(s.get_brightness(0.9, 0.8, bell), expect, 1e-9, errors);
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
    return errors;
}

int main() {
    // Note: This tester is not as rigorous as the Pytest tests.  It is intended mainly to catch
    // memory management errors and distinguish between C++ errors and Cython wrapping errors.
    cout << "===== Running C++ Tester =====" << endl;
    int errors = 0;
    errors += test_biellipsoid();
    errors += test_orbital_position();
    errors += test_transit();
    errors += test_phase_curve();
    errors += test_light_source();
    if (errors == 0) {
        cout << "All tests passed." << endl << endl;
    } else {
        cout << "ERROR: " << errors << " test(s) failed." << endl << endl;
    }
    return errors;
}
