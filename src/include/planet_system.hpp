#ifndef PLANET_SYSTEM_HPP
#define PLANET_SYSTEM_HPP

#include "light_source.hpp"
#include "orbit.hpp"
#include "shape.hpp"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#define MAX_SYSTEM_OBJECTS 8

double emission_integrand(double y, void *params);
double emission_outer_integral(double y, void *params);
int process_bounds(Bounds *b, int n_relevant, bool invert = false);

class PlanetSystem {
  private:
    int max_steps; // The max steps used in the GSL integrator
    int n_objects; // Total number of objects considered (< MAX_SYSTEM_OBJECTS)

    // Working variables during the integration process
    //// The x to evaluate the inner integral at (for this step)
    double x;
    // Which object is currently being integrated
    int i_target;
    // Whether each object is relevant (occluding) during integration of i_target
    bool relevant[MAX_SYSTEM_OBJECTS];
    double invert_integral;

    // GSL integration working variables
    gsl_integration_workspace *workspaceOuter;
    gsl_integration_workspace *workspaceInner;
    gsl_function integInner;
    gsl_function integOuter;

  public:
    Orbit orbits[MAX_SYSTEM_OBJECTS];
    int parent_indices[MAX_SYSTEM_OBJECTS];
    Shape shapes[MAX_SYSTEM_OBJECTS];
    LightSource lights[MAX_SYSTEM_OBJECTS];
    bool rotate_with_orbit[MAX_SYSTEM_OBJECTS];
    Bounds xlim[MAX_SYSTEM_OBJECTS];
    Bounds ylim[MAX_SYSTEM_OBJECTS];
    double atol;
    double rtol;

    PlanetSystem(double atol = 1e-6, double rtol = 1e-3, int max_steps = 100);
    PlanetSystem(PlanetSystem &p);
    PlanetSystem &operator=(const PlanetSystem &other);
    ~PlanetSystem();

    int add_object(
        const Orbit &orb, const Shape &bell, const LightSource &source, bool rot_with_orbit = false,
        int parent_index = -1
    );
    int get_n_objects() const;
    void clear_objects();
    void set_time(double t);
    double integrate_single(int i);
    void integrate(double *times, double *outputs, int n);

    friend double emission_integrand(double y, void *params);
    friend double emission_outer_integral(double x, void *params);
};

#endif
