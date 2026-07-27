#include "transit_integral.hpp"


double transit_integrand(double y, void *params) {
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    return g->emitter.get_brightness_sphere(g->x, y);
}


double transit_inner_integral(double x, void *params) {
    int code;
    double result, err;
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    g->x = x;

    // Get y bounds and check for no overlap between planet and star (at this x)
    double y_bound = sqrt(1 - x * x);
    Bounds ylim = g->shp.slice_ylimits(x);
    ylim.min = fmax(ylim.min, -y_bound);
    ylim.max = fmin(ylim.max, y_bound);
    if (ylim.min >= ylim.max) {
        return 0.;
    }

    code = gsl_integration_qag(
        g->integrand, ylim.min, ylim.max, 1e-8, 1e-8, 100, 1, g->work, &result, &err
    );
    if (integration_failed(code, result, err, 1e-9, 1e-7)) {
        return NAN;
    }
    return result;
}


void transit_integral(
    double *times, double *outputs, int n, const Orbit &orb, const LightSource &emitter,
    double theta, double phi, double gamma, double r_forward, double r_back, double r_up,
    double r_side, bool rotate_with_orbit, double atol, double rtol
) {
    int code = 0;
    Vec3 split_point;
    double result, err;
    Bounds x_lim, y_lim;
    double x_min, x_max, x;
    double st = sin(theta);
    double ct = cos(theta);
    bool discontinuous_pole = false;

    if (r_up < 0) {
        r_up = r_back;
        discontinuous_pole = true;
    }
    Shape shp = Shape(r_forward, r_back, r_up, r_side);
    if (!(discontinuous_pole || rotate_with_orbit)) {
        shp.set_rotation(theta, phi, gamma);
    }

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();

    // Prepare the inner (y) integral variables
    gsl_integration_workspace *workspaceInner = gsl_integration_workspace_alloc(100);
    gsl_function integInner;
    TransitIntegralParams g = {emitter, shp, 0., atol, rtol, workspaceInner, &integInner};
    integInner.function = &transit_integrand;
    integInner.params = &g;
    g.integrand = &integInner;

    // Now prepare the outer (x) integral variables
    gsl_integration_workspace *workspaceOuter = gsl_integration_workspace_alloc(100);
    gsl_function integOuter;
    integOuter.function = &transit_inner_integral;
    integOuter.params = &g;

    for (int i = 0; i < n; i++) {
        g.shp.position_from_orbit(times[i], orb, rotate_with_orbit && (!discontinuous_pole));
        if (g.shp.position.z < 0) {
            outputs[i] = 1.0;
            continue;
        }

        // Rotate and get planet bounding box
        if (discontinuous_pole) {
            x = ct * g.shp.position.x + st * g.shp.position.y;
            g.shp.position.y = -st * g.shp.position.x + ct * g.shp.position.y;
            g.shp.position.x = x;
            x_lim = {x - r_back, x + r_forward};
            y_lim = {
                g.shp.position.y - fmax(r_back, r_forward),
                g.shp.position.y + fmax(r_back, r_forward)
            };
        } else {
            if (rotate_with_orbit) {
                g.shp.set_rotation(theta, phi, gamma, orb.get_cos_inc());
            }
            x_lim = g.shp.x_bounds();
            y_lim = g.shp.y_bounds();
        }

        // Quick bounding-box check to skip most non-transits
        x_min = fmax(x_lim.min, -1.);
        x_max = fmin(x_lim.max, 1.);
        if ((x_min >= x_max) || (y_lim.min > 1.) || (y_lim.max < -1)) {
            outputs[i] = 1.0;
            continue;
        }

        if (discontinuous_pole) {
            // Split the integral around the middle of the planet to allow for a discontinuous pole
            split_point = g.shp.position;
            g.shp.set_radii(r_forward, r_back, r_back, r_side);
        } else {
            // Split the integral around the nearest point to help integrator find non-zero areas
            split_point = g.shp.nearest_to_line(0., 0.);
            if (split_point.z > 1.) {
                outputs[i] = 1.0;
                continue;
            }
        }
        code = gsl_integration_qag(
            &integOuter, x_min, split_point.x, .1 * atol, .1 * rtol, 100, 1, workspaceOuter,
            &result, &err
        );
        if (integration_failed(code, result, err, atol, rtol)) {
            outputs[i] = NAN;
            continue;
        }
        outputs[i] = 1 - result;
        if (discontinuous_pole) {
            g.shp.set_radii(r_forward, r_back, r_forward, r_side);
        }
        code = gsl_integration_qag(
            &integOuter, split_point.x, x_max, .1 * atol, .1 * rtol, 100, 1, workspaceOuter,
            &result, &err
        );
        if (integration_failed(code, result, err, atol, rtol)) {
            outputs[i] = NAN;
            continue;
        }
        outputs[i] -= result;
    }

    // Cleanup
    gsl_integration_workspace_free(workspaceInner);
    gsl_integration_workspace_free(workspaceOuter);
    return;
}
