#ifndef PHASE_CURVE_INTEGRAL_HPP
#define PHASE_CURVE_INTEGRAL_HPP

#include "biellipsoid.hpp"
#include "light_source.hpp"
#include "orbit.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define MAX_PHASE_OBJECTS 8

double phase_curve_integrand(double y, void *params);
double phase_curve_inner_integral(double y, void *params);

class PhaseIntegrator {
  private:
    int n_objects; // Total number of objects considered (< MAX_PHASE_OBJECTS)
    double x;      // The x to evaluate the inner integral at (for this step)
    int i_target;  // Which object is currently being integrated
    bool relevant[MAX_PHASE_OBJECTS];

    // GSL integration working variables
    gsl_integration_workspace *workspaceOuter;
    gsl_integration_workspace *workspaceInner;
    gsl_function integInner;
    gsl_function integOuter;

  public:
    Orbit orbits[MAX_PHASE_OBJECTS];
    Biellipsoid shapes[MAX_PHASE_OBJECTS];
    LightSource lights[MAX_PHASE_OBJECTS];
    bool rotate_with_orbit[MAX_PHASE_OBJECTS];
    Bounds xlim[MAX_PHASE_OBJECTS];
    Bounds ylim[MAX_PHASE_OBJECTS];

    PhaseIntegrator();
    PhaseIntegrator(PhaseIntegrator &p);
    ~PhaseIntegrator();

    int add_object(const Orbit &orb, const Biellipsoid &bell, const LightSource &source, bool rot_with_orbit = false);
    int get_n_objects() const;
    void clear_objects();
    void set_time(double t);
    double integrate_single(int i);
    void phase_curve_integral(double *times, double *outputs, int n);

    friend double phase_curve_integrand(double y, void *params);
    friend double phase_curve_inner_integral(double x, void *params);
};

#endif
