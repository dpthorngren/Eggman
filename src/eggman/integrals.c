#include "eggman.h"


Vec3 nearest_point(double xe, double ye, double a, double b) {
    // Returns the x, y, and squared dist of the nearest point to the origin on the
    // ellipse, unless it contains the origin, in which case returns [0, 0, -1]
    int i;
    double axe = fabs(xe);
    double aye = fabs(ye);
    if ((axe / a) * (axe / a) + (aye / b) * (aye / b) < 1.) {
        return (Vec3){0., 0., -1};
    }

    double x, y, d;
    double x0 = axe - a;
    double x1 = axe;
    double y0 = aye;
    double y1 = aye - b;
    double d0 = x0 * x0 + y0 * y0;
    double d1 = x1 * x1 + y1 * y1;
    // Midpoint method
    for (i = 0; i < 15; i++) {
        x = (x0 + x1) / 2.0;
        y = aye - b * sqrt(1 - ((axe - x) / a) * ((axe - x) / a));
        d = x * x + y * y;
        if (d1 < d0) {
            x0 = x;
            y0 = y;
            d0 = d;
        } else {
            x1 = x;
            y1 = y;
            d1 = d;
        }
    }
    return (Vec3){x0 * SIGN(xe), y0 * SIGN(ye), d0};
}


double transit2d_integrand(double y, void *params) {
    Transit2dIntegralParams *g = (Transit2dIntegralParams *)params;
    return get_source_brightness(&g->emitter, g->x, y);
}


double transit2d_inner_integral(double x, void *params) {
    int code;
    double result, err;
    Transit2dIntegralParams *g = (Transit2dIntegralParams *)params;
    g->x = x;

    // Which side of the biellipse are we on?
    double r_side = (x < g->xe) ? g->r_back : g->r_forward;
    // If r_up is negative, use the discontinuous two-circle setup from Catwoman
    double r_up = (g->r_up < 0.) ? r_side : g->r_up;
    // Calculate the bounds of integration in y
    double ys = (x - g->xe) / r_side;
    ys = sqrt(1 - ys * ys);
    double y0 = fmax(g->ye - r_up * ys, -sqrt(1 - x * x));
    double y1 = fmin(g->ye + r_up * ys, sqrt(1 - x * x));

    // Check for no overlap between planet and star (at this x)
    if (y0 >= y1) {
        return 0.;
    }

    code = gsl_integration_qag(g->integrand, y0, y1, 1e-7, 1e-7, 100, 1, g->work, &result, &err);
    if (code != 0) {
        return NAN;
    }
    return result;
}


void transit2d_integral(
    double *times, double *outputs, int n, Orbit *orb, LightSource *emitter, double r_forward,
    double r_back, double r_up
) {
    int i;
    int code = 0;
    Vec3 position;
    Vec3 nearest;
    double result, err;
    double x_min, x_max;
    double r_side, r_top;
    double r_up_bound = (r_up < 0.) ? fmax(r_forward, r_back) : r_up;

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();

    // Prepare the inner (y) integral variables
    gsl_integration_workspace *workspaceInner = gsl_integration_workspace_alloc(100);
    gsl_function integInner;
    Transit2dIntegralParams g = {*emitter,       0.,         0., r_forward, r_back, r_up, 0.,
                                 workspaceInner, &integInner};
    integInner.function = &transit2d_integrand;
    integInner.params = &g;
    g.integrand = &integInner;

    // Now prepare the outer (x) integral variables
    gsl_integration_workspace *workspaceOuter = gsl_integration_workspace_alloc(100);
    gsl_function integOuter;
    integOuter.function = &transit2d_inner_integral;
    integOuter.params = &g;

    for (i = 0; i < n; i++) {
        position = get_position(orb, times[i]);
        if (position.z < 0) {
            outputs[i] = 1.0;
            continue;
        }
        g.xe = position.x;
        g.ye = position.y;

        // Compute limits of integration on the x axis
        x_min = fmax(g.xe - r_back, -1.);
        x_max = fmin(g.xe + r_forward, 1.);

        // Quick bounding-box check to skip most non-transits
        if ((x_min >= x_max) || (g.ye - r_up_bound > 1.) || (g.ye + r_up_bound < -1.)) {
            outputs[i] = 1.0;
            continue;
        }

        r_side = (g.xe > 0) ? r_back : r_forward;
        nearest = nearest_point(g.xe, g.ye, r_side, (r_up > 0) ? r_up : r_side);
        if ((r_up > 0) && (nearest.z >= 1)) {
            outputs[i] = 1.0;
            continue;
        }

        code = gsl_integration_qag(
            &integOuter, x_min, nearest.x, 1e-7, 1e-7, 100, 1, workspaceOuter, &result, &err
        );
        outputs[i] = (code != 0) ? NAN : 1 - result;
        code = gsl_integration_qag(
            &integOuter, nearest.x, x_max, 1e-7, 1e-7, 100, 1, workspaceOuter, &result, &err
        );
        outputs[i] = (code != 0) ? NAN : outputs[i] - result;
    }

    // Cleanup
    gsl_integration_workspace_free(workspaceInner);
    gsl_integration_workspace_free(workspaceOuter);
    return;
}
