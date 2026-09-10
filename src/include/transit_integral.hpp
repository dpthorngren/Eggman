#ifndef TRANSIT_INTEGRAL_HPP
#define TRANSIT_INTEGRAL_HPP

#include "orbit.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

typedef struct {
    // Integrand settings
    double a;  // Radius along the x axis
    double b;  // Radius along the y axis
    double xe; // x position of the ellipse
    double ye; // y position of the ellipse
    double limb[4];

    // Integration settings
    double atol;
    double rtol;
    int max_steps;

    // Working variables
    double x; // The x to evaluate the inner integral at (for this step)

    // GSL integration working variables
    gsl_integration_workspace *work;
    gsl_function *integrand;
} TransitIntegralParams;

double transit_integrand(double y, void *params);
double transit_inner_integral(double x, void *params);

void transit_integral(
    double *times, double *outputs, int n, const Orbit &orb, double r_forward, double r_back,
    double r_up, double limb0, double limb1, double limb2, double limb3, double theta = 0.,
    double atol = 1e-6, double rtol = 1e-3, int max_steps = 100
);

#endif
