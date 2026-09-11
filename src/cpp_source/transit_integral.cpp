#include "transit_integral.hpp"
#include "math_utils.hpp"


double transit_integrand(double y, void *params) {
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    double mu = sqrt(fmax(1 - g->x * g->x - y * y, 0));
    double nu;
    if (g->limb[3] < 0) {
        nu = 1 - mu;
        return 1. - g->limb[0] * nu - g->limb[1] * nu * nu;
    } else {
        nu = sqrt(mu);
        return 1. - g->limb[0] * (1. - nu) - g->limb[1] * (1. - mu) - g->limb[2] * (1. - mu * nu) -
               g->limb[3] * (1. - mu * mu);
    }
}


double transit_inner_integral(double x, void *params) {
    int code;
    double result, err;
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    g->x = x;

    // Get y bounds of integral: overlap between planet and star
    double y_star = sqrt(1 - x * x);
    double y_planet = (x - g->xe) / g->a;
    y_planet = g->b * sqrt(1 - y_planet * y_planet);
    Bounds ylim = {fmax(-y_star, g->ye - y_planet), fmin(y_star, g->ye + y_planet)};

    // Check for no overlap between star and planet at x
    if (ylim.min >= ylim.max) {
        return 0.;
    }

    code = gsl_integration_qag(
        g->integrand, ylim.min, ylim.max, .1 * g->atol, .1 * g->rtol, g->max_steps, 1, g->work,
        &result, &err
    );
    if (integration_failed(code, result, err, g->atol, g->rtol)) {
        return NAN;
    }
    return result;
}


void transit_integral(
    double *times, double *outputs, int n, const Orbit &orb, double r_forward, double r_back,
    double r_up, double limb0, double limb1, double limb2, double limb3, double theta, double atol,
    double rtol, int max_steps
) {
    int code = 0;
    Vec3 loc;
    double result, err;
    Bounds xb, yb;

    double st = sin(theta);
    double ct = cos(theta);
    double limb_norm;
    if (limb3 < 0) {
        limb_norm = M_PI * (1. - limb0 / 3. - limb1 / 6.);
    } else {
        limb_norm = M_PI * (1 - limb0 / 5. - limb1 / 3. - 3. * limb2 / 7. - limb3 / 2.);
    }

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();

    // Prepare the inner (y) integral variables
    gsl_integration_workspace *workspaceInner = gsl_integration_workspace_alloc(max_steps);
    gsl_function integInner;
    TransitIntegralParams g = {
        r_forward,   r_back, 0.,        0., {limb0, limb1, limb2, limb3},
        atol,        rtol,   max_steps, 0., workspaceInner,
        &integInner,
    };
    integInner.function = &transit_integrand;
    integInner.params = &g;

    // Now prepare the outer (x) integral variables
    gsl_integration_workspace *workspaceOuter = gsl_integration_workspace_alloc(max_steps);
    gsl_function integOuter;
    integOuter.function = &transit_inner_integral;
    integOuter.params = &g;

    for (int i = 0; i < n; i++) {
        loc = orb.get_position(times[i]);
        if (loc.z < 0) {
            outputs[i] = 1.0;
            continue;
        }

        // Rotate planet so it's axes are aligned with the coordinate system
        g.xe = ct * loc.x - st * loc.y;
        g.ye = st * loc.x + ct * loc.y;

        // Get planet bounding box
        xb = {g.xe - r_back, g.xe + r_forward};
        if (r_up < 0) {
            yb = {g.ye - fmax(r_back, r_forward), g.ye + fmax(r_back, r_forward)};
        } else {
            yb = {g.ye - r_up, g.ye + r_up};
        }

        // Clip to stellar bounding box and check for trivial non-transits
        xb.min = fmax(xb.min, -1.);
        xb.max = fmin(xb.max, 1.);
        if ((xb.min >= xb.max) || (yb.min > 1.) || (yb.max < -1)) {
            outputs[i] = 1.0;
            continue;
        }

        // Integrate the back side of the planet
        outputs[i] = 1.0;
        if (g.xe > xb.min) {
            g.a = r_back;
            g.b = r_up < 0 ? r_back : r_up;
            code = gsl_integration_qag(
                &integOuter, xb.min, g.xe, .1 * atol, .1 * rtol, max_steps, 1, workspaceOuter,
                &result, &err
            );
            if (integration_failed(code, result, err, atol, rtol)) {
                outputs[i] = NAN;
                continue;
            }
            outputs[i] -= result / limb_norm;
        }

        // Integrate the forward side of the planet
        if (g.xe < xb.max) {
            g.a = r_forward;
            g.b = r_up < 0 ? r_forward : r_up;
            code = gsl_integration_qag(
                &integOuter, g.xe, xb.max, .1 * atol, .1 * rtol, max_steps, 1, workspaceOuter,
                &result, &err
            );
            if (integration_failed(code, result, err, atol, rtol)) {
                outputs[i] = NAN;
                continue;
            }
            outputs[i] -= result / limb_norm;
        }
    }

    // Cleanup
    gsl_integration_workspace_free(workspaceInner);
    gsl_integration_workspace_free(workspaceOuter);
    return;
}
