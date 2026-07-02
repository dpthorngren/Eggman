#include "biellipsoid.hpp"
#include "light_source.hpp"
#include "orbit.hpp"
#include "phase_curve_integral.hpp"


int main() {
    PhaseIntegrator p;

    // Add the star
    Orbit orb = Orbit();
    double source_params[MAX_SOURCE_PARAMS] = {1.0 / M_PI, 0.0, 0.0, 0.0, 0.0, 0.0,
                                               0.0,        0.0, 0.0, 0.0, 0.0, 0.0};
    double limb_params[MAX_LIMB_PARAMS] = {0.2, 0.1, 0., 0.};
    Biellipsoid bell = Biellipsoid();
    LightSource star = LightSource(1, source_params, 1, limb_params);
    p.add_object(orb, bell, star);

    // Add the planet
    orb = Orbit(10., 0., 5., 0.00, 89., 90.);
    bell = Biellipsoid(.12, .09, .08, .08);
    source_params[0] = 1e-3;
    limb_params[0] = 0.;
    limb_params[1] = 0.;
    LightSource planet = LightSource(1, source_params, 0, limb_params);
    p.add_object(orb, bell, planet);

    // Test a full phase curve
    const int n_times = 1000;
    double time[n_times];
    double result[n_times];
    for (int i = 0; i < n_times; i++) {
        time[i] = orb.get_period() * i / n_times;
        result[i] = -1.;
    }
    p.phase_curve_integral(time, result, n_times);
    printf("%f\n", result[0]);

    return 0;
}
