#include "biellipsoid.hpp"
#include "light_source.hpp"
#include "orbit.hpp"
#include "phase_curve_integral.hpp"
#include "transit_integral.hpp"
#import <cmath>
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
        cout << "ERROR in " << __func__ << " at line " << __LINE__ << endl;                        \
        cout << "  Assertion \"" << #A << "~=" << #B << "\" (" << (A) << "~=" << (B)               \
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
    double source_params[MAX_SOURCE_PARAMS] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double limb_params[MAX_LIMB_PARAMS] = {0.0, 0.0, 0.0, 0.0};
    LightSource star = LightSource(1, source_params, 1, limb_params);

    const int n_times = 100;
    double outputs[n_times];
    double times[n_times];
    for (int i = 0; i < n_times; i++) {
        times[i] = i * orb.get_period() / n_times;
        outputs[i] = -1.0;
    }

    transit_integral(times, outputs, n_times, orb, star, 0., 0., 0., .1, .09, .08, .07, true);

    // General bounds
    for (int i = 0; i < n_times; i++) {
        TEST_ASSERT(outputs[i], <=, 1., errors);
        TEST_ASSERT(outputs[i], >=, 0., errors);
    }
    // Out of transit
    TEST_APPROX(outputs[n_times / 4], 1.0, 1e-12, errors);
    TEST_APPROX(outputs[n_times / 2], 1.0, 1e-12, errors);
    // In-transit
    TEST_ASSERT(outputs[0], <, 1.0, errors);
    TEST_ASSERT(outputs[1], <, 1.0, errors);
    TEST_ASSERT(outputs[n_times - 2], <, 1.0, errors);
    return errors;
}

int test_phase_curve() {
    ANNOUNCE_TEST();
    int errors = 0;
    PhaseIntegrator p;

    // Add the star
    Orbit orb = Orbit();
    double source_params[MAX_SOURCE_PARAMS] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double limb_params[MAX_LIMB_PARAMS] = {0., 0., 0.0, 0.0};
    Biellipsoid bell = Biellipsoid();
    bell.update_derived();
    LightSource star = LightSource(1, source_params, 1, limb_params);
    TEST_APPROX(star.limb_norm, M_PI, 1e-9, errors)
    TEST_APPROX(star.get_brightness(0., 0., 0.), 1 / M_PI, 1e-9, errors)
    p.add_object(orb, bell, star);
    TEST_ASSERT(p.get_n_objects(), ==, 1, errors);

    // Add the planet
    orb = Orbit(10., 0., 5., 0.00, 89., 90.);
    bell = Biellipsoid(.12, .09, .08, .08);
    bell.update_derived();
    source_params[0] = 1e-6;
    LightSource planet = LightSource(1, source_params, 0, limb_params);
    TEST_APPROX(star.limb_norm, M_PI, 1e-9, errors);
    TEST_APPROX(planet.get_brightness(0., 0., 0.), 1e-6 / M_PI, 1e-9, errors)
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
    if (errors == 0) {
        cout << "All tests passed." << endl << endl;
    } else {
        cout << "ERROR: " << errors << " test(s) failed." << endl << endl;
    }
    return errors;
}
