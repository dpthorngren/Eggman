#include "eggman.h"


double transit3d_integrand(double y, void *params) {
    Transit3dIntegralParams *g = (Transit3dIntegralParams *)params;
    return get_source_brightness(&g->emitter, g->x, y);
}


double transit3d_inner_integral(double x, void *params) {
    int code;
    double result, err;
    Transit3dIntegralParams *g = (Transit3dIntegralParams *)params;
    g->x = x;

    // Get y bounds and check for no overlap between planet and star (at this x)
    Bounds ylim = get_ylimits(&g->bell, x);
    if (ylim.min >= ylim.max) {
        return 0.;
    }

    code = gsl_integration_qag(
        g->integrand, ylim.min, ylim.max, 1e-7, 1e-7, 100, 1, g->work, &result, &err
    );
    return (code == 0) ? result : NAN;
}


void transit3d_integral(
    double *times, double *outputs, int n, Orbit *orb, LightSource *emitter, double theta,
    double phi, double gamma, double r_forward, double r_back, double r_up, double r_side
) {
    int i;
    int code = 0;
    Vec3 position = {0., 0., 0.};
    Vec3 nearest;
    Biellipsoid bell =
        create_biellipsoid(position, theta, phi, gamma, r_forward, r_back, r_up, r_side);
    double result, err;
    double x_min, x_mid, x_max;

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();

    // Prepare the inner (y) integral variables
    gsl_integration_workspace *workspaceInner = gsl_integration_workspace_alloc(100);
    gsl_function integInner;
    Transit3dIntegralParams g = {*emitter, bell, 0., workspaceInner, &integInner};
    integInner.function = &transit3d_integrand;
    integInner.params = &g;
    g.integrand = &integInner;

    // Now prepare the outer (x) integral variables
    gsl_integration_workspace *workspaceOuter = gsl_integration_workspace_alloc(100);
    gsl_function integOuter;
    integOuter.function = &transit3d_inner_integral;
    integOuter.params = &g;

    for (i = 0; i < n; i++) {
        position = get_position(orb, times[i]);
        if (position.z < 0) {
            outputs[i] = 1.0;
            continue;
        }
        // TODO: Reorient the biellipsoid rather than just repositioning it.
        set_position(&g.bell, position);

        // Get the integration bounds
        x_min = fmax(g.bell.xbounds.min, -1.);
        x_max = fmin(g.bell.xbounds.max, 1.);

        // Quick bounding-box check to skip most non-transits
        if ((x_min >= x_max) || (g.bell.ybounds.min > 1.) || (g.bell.ybounds.max < -1)) {
            outputs[i] = 1.0;
            continue;
        }

        // Split the integral around the nearest point to help integrator find non-zero areas
        nearest = nearest_point_3d(&g.bell);
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
