#ifndef PHASE_CURVE_INTEGRAL_HPP
#define PHASE_CURVE_INTEGRAL_HPP

#include "light_source.hpp"
#include "orbit.hpp"
#include "shape.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define MAX_PHASE_OBJECTS 8

double phase_curve_integrand(double y, void *params);
double phase_curve_inner_integral(double y, void *params);
int process_bounds(Bounds *b, int n_relevant);

class PlanetSystem {
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
    int parent_indices[MAX_PHASE_OBJECTS];
    Shape shapes[MAX_PHASE_OBJECTS];
    LightSource lights[MAX_PHASE_OBJECTS];
    bool rotate_with_orbit[MAX_PHASE_OBJECTS];
    Bounds xlim[MAX_PHASE_OBJECTS];
    Bounds ylim[MAX_PHASE_OBJECTS];
    double atol;
    double rtol;

    PlanetSystem(double atol = 1e-6, double rtol = 1e-3);
    PlanetSystem(PlanetSystem &p);
    ~PlanetSystem();

    int add_object(
        const Orbit &orb, const Shape &bell, const LightSource &source,
        bool rot_with_orbit = false, int parent_index = -1
    );
    int get_n_objects() const;
    void clear_objects();
    void set_time(double t);
    double integrate_single(int i);
    void phase_curve_integral(double *times, double *outputs, int n);

    friend double phase_curve_integrand(double y, void *params);
    friend double phase_curve_inner_integral(double x, void *params);
};

#endif
