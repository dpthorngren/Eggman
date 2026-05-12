#include "transit_integral.hpp"


double transit_integrand(double y, void *params) {
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    return g->emitter.get_brightness(g->x, y);
}


double transit_inner_integral(double x, void *params) {
    int code;
    double result, err;
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    g->x = x;

    // Get y bounds and check for no overlap between planet and star (at this x)
    Bounds ylim = g->bell.slice_ylimits(x);
    if (ylim.min >= ylim.max) {
        return 0.;
    }

    code = gsl_integration_qag(
        g->integrand, ylim.min, ylim.max, 1e-7, 1e-7, 100, 1, g->work, &result, &err
    );
    return (code == 0) ? result : NAN;
}


void transit_integral(
    double *times, double *outputs, int n, Orbit *orb, LightSource *emitter, double theta,
    double phi, double gamma, double r_forward, double r_back, double r_up, double r_side,
    bool rotate_with_orbit
) {
    int i;
    int code = 0;
    Vec3 position = {0., 0., 0.};
    Vec3 nearest;
    Biellipsoid bell = Biellipsoid(r_forward, r_back, r_up, r_side);
    bell.set_rotation(theta, phi, gamma);
    double result, err;
    Bounds x_lim, y_lim;
    double x_min, x_max;

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();

    // Prepare the inner (y) integral variables
    gsl_integration_workspace *workspaceInner = gsl_integration_workspace_alloc(100);
    gsl_function integInner;
    TransitIntegralParams g = {*emitter, bell, 0., workspaceInner, &integInner};
    integInner.function = &transit_integrand;
    integInner.params = &g;
    g.integrand = &integInner;

    // Now prepare the outer (x) integral variables
    gsl_integration_workspace *workspaceOuter = gsl_integration_workspace_alloc(100);
    gsl_function integOuter;
    integOuter.function = &transit_inner_integral;
    integOuter.params = &g;

    for (i = 0; i < n; i++) {
        position = orb->get_position(times[i]);
        if (position.z < 0) {
            outputs[i] = 1.0;
            continue;
        }
        g.bell.set_position(position);
        if (rotate_with_orbit) {
            g.bell.set_rotation(theta, phi, gamma, orb->get_cos_inc());
        }

        // Get the integration bounds
        x_lim = g.bell.x_bounds();
        y_lim = g.bell.x_bounds();
        x_min = fmax(x_lim.min, -1.);
        x_max = fmin(x_lim.max, 1.);

        // Quick bounding-box check to skip most non-transits
        if ((x_min >= x_max) || (y_lim.min > 1.) || (y_lim.max < -1)) {
            outputs[i] = 1.0;
            continue;
        }

        // Split the integral around the nearest point to help integrator find non-zero areas
        nearest = g.bell.nearest_to_line(0., 0.);
        if (nearest.z > 1.) {
            outputs[i] = 1.0;
            continue;
        }
        code = gsl_integration_qag(
            &integOuter, x_min, nearest.x, 1e-7, 1e-5, 100, 1, workspaceOuter, &result, &err
        );
        outputs[i] = (code != 0) ? NAN : 1 - result;
        code = gsl_integration_qag(
            &integOuter, nearest.x, x_max, 1e-7, 1e-5, 100, 1, workspaceOuter, &result, &err
        );
        outputs[i] = (code != 0) ? NAN : outputs[i] - result;
    }

    // Cleanup
    gsl_integration_workspace_free(workspaceInner);
    gsl_integration_workspace_free(workspaceOuter);
    return;
}
