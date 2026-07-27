#ifndef TRANSIT_INTEGRAL_HPP
#define TRANSIT_INTEGRAL_HPP

#include "shape.hpp"
#include "light_source.hpp"
#include "orbit.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

typedef struct {
    LightSource emitter; // The star being transited
    Shape shp;    // The planet doing the transiting
    double x;            // The x to evaluate the inner integral at (for this step)
    double atol;
    double rtol;
    // GSL integration working variables
    gsl_integration_workspace *work;
    gsl_function *integrand;
} TransitIntegralParams;

double transit_integrand(double y, void *params);
double transit_inner_integral(double x, void *params);
void transit_integral(
    double *times, double *outputs, int n, const Orbit &orb, const LightSource &emitter, double theta,
    double phi, double gamma, double r_forward, double r_back, double r_up, double r_side,
    bool rotate_with_orbit, double atol = 1e-6, double rtol = 1e-3
);

#endif
