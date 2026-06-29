#include "transit_integral.hpp"
#include <gsl/gsl_errno.h>
#include <iostream>


double transit_integrand(double y, void *params) {
    TransitIntegralParams *g = (TransitIntegralParams *)params;
    double mu = sqrt(g->x * g->x + y * y);
    return g->emitter.get_brightness(mu, 0., 0.);
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
        g->integrand, ylim.min, ylim.max, 1e-5, 1e-7, 100, 1, g->work, &result, &err
    );
    if ((code != 0) &&
        (((code != GSL_EMAXITER) && (code != GSL_EROUND)) || (err > 1e-9 + 1e-6 * fabs(result)))) {
        std::cout << "INTEGRATION ERROR (Inner integral) " << code << ": " << gsl_strerror(code)
                  << "Output = " << result << ", err = " << err << std::endl;
        return NAN;
    }
    return result;
}


void transit_integral(
    double *times, double *outputs, int n, const Orbit &orb, const LightSource &emitter,
    double theta, double phi, double gamma, double r_forward, double r_back, double r_up,
    double r_side, bool rotate_with_orbit
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
    Biellipsoid bell = Biellipsoid(r_forward, r_back, r_up, r_side);
    if (!(discontinuous_pole || rotate_with_orbit)) {
        bell.set_rotation(theta, phi, gamma);
    }

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();

    // Prepare the inner (y) integral variables
    gsl_integration_workspace *workspaceInner = gsl_integration_workspace_alloc(100);
    gsl_function integInner;
    TransitIntegralParams g = {emitter, bell, 0., workspaceInner, &integInner};
    integInner.function = &transit_integrand;
    integInner.params = &g;
    g.integrand = &integInner;

    // Now prepare the outer (x) integral variables
    gsl_integration_workspace *workspaceOuter = gsl_integration_workspace_alloc(100);
    gsl_function integOuter;
    integOuter.function = &transit_inner_integral;
    integOuter.params = &g;

    for (int i = 0; i < n; i++) {
        g.bell.position_from_orbit(times[i], orb, rotate_with_orbit && (!discontinuous_pole));
        if (g.bell.position.z < 0) {
            outputs[i] = 1.0;
            continue;
        }

        // Rotate and get planet bounding box
        if (discontinuous_pole) {
            x = ct * g.bell.position.x + st * g.bell.position.y;
            g.bell.position.y = -st * g.bell.position.x + ct * g.bell.position.y;
            g.bell.position.x = x;
            x_lim = (Bounds){x - r_back, x + r_forward};
            y_lim = (Bounds){
                g.bell.position.y - fmax(r_back, r_forward),
                g.bell.position.y + fmax(r_back, r_forward)
            };
        } else {
            if (rotate_with_orbit) {
                g.bell.set_rotation(theta, phi, gamma, orb.get_cos_inc());
            }
            x_lim = g.bell.x_bounds();
            y_lim = g.bell.y_bounds();
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
            split_point = g.bell.position;
            g.bell.set_radii(r_forward, r_back, r_back, r_side);
        } else {
            // Split the integral around the nearest point to help integrator find non-zero areas
            split_point = g.bell.nearest_to_line(0., 0.);
            if (split_point.z > 1.) {
                outputs[i] = 1.0;
                continue;
            }
        }
        code = gsl_integration_qag(
            &integOuter, x_min, split_point.x, 1e-6, 1e-9, 100, 1, workspaceOuter, &result, &err
        );
        if ((code != 0) && (((code != GSL_EMAXITER) && (code != GSL_EROUND)) ||
                            (err > 1e-9 + 1e-6 * fabs(result)))) {
            std::cout << "INTEGRATION ERROR (Outer integral 1) at i, t = " << i << ", " << times[i]
                      << code << ": " << gsl_strerror(code) << "Output = " << result
                      << ", err = " << err << std::endl;
            outputs[i] = NAN;
        }
        outputs[i] = 1 - result / M_PI;
        if (discontinuous_pole) {
            g.bell.set_radii(r_forward, r_back, r_forward, r_side);
        }
        code = gsl_integration_qag(
            &integOuter, split_point.x, x_max, 1e-6, 1e-9, 100, 1, workspaceOuter, &result, &err
        );
        if ((code != 0) && (((code != GSL_EMAXITER) && (code != GSL_EROUND)) ||
                            (err > 1e-9 + 1e-6 * fabs(result)))) {
            std::cout << "INTEGRATION ERROR (Outer integral 2) at i, t = " << i << ", " << times[i]
                      << code << ": " << gsl_strerror(code) << "Output = " << result
                      << ", err = " << err << std::endl;
            outputs[i] = NAN;
        }
        outputs[i] -= result / M_PI;
    }

    // Cleanup
    gsl_integration_workspace_free(workspaceInner);
    gsl_integration_workspace_free(workspaceOuter);
    return;
}
