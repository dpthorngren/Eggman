#include "biellipsoid.hpp"
#include "light_source.hpp"
#include "orbit.hpp"
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
    bell.position_from_orbit(0.0, &orb);
    TEST_APPROX(bell.position.x, 0., 1e-12, errors);
    TEST_APPROX(bell.position.y, 5. * cos(85 * M_PI / 180), 1e-12, errors);
    TEST_APPROX(bell.position.z, 5. * sin(85 * M_PI / 180), 1e-12, errors);
    bell.position_from_orbit(0.5, &orb);
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

    transit_integral(times, outputs, n_times, &orb, &star, 0., 0., 0., .1, .09, .08, .07, true);

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

int main() {
    // Note: This tester is not as rigorous as the Pytest tests.  It is intended mainly to catch
    // memory management errors and distinguish between C++ errors and Cython wrapping errors.
    cout << "===== Running C++ Tester =====" << endl;
    int errors = 0;
    errors += test_biellipsoid();
    errors += test_orbital_position();
    errors += test_transit();
    if (errors == 0) {
        cout << "All tests passed." << endl << endl;
    } else {
        cout << "ERROR: " << errors << " test(s) failed." << endl << endl;
    }
    return errors;
}
