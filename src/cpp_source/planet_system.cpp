#include "planet_system.hpp"
#include "math_utils.hpp"
#include <cmath>

double phase_curve_integrand(double y, void *params) {
    PlanetSystem *p = (PlanetSystem *)params;
    return p->lights[p->i_target].get_brightness(p->x, y, p->shapes[p->i_target]);
}

int process_bounds(Bounds *b, int n_relevant) {
    // Iterate through potential occluders, modifying the integration bounds stack as-needed
    // Bound stack 'b' should be at least length n_relevant and contain the overall range
    // in position [0], followed by the bounds of the potential occluders.
    int i, j, k;
    Bounds occ;
    int n_bounds = 1;
    for (i = 1; i < n_relevant; i++) {
        occ = b[i];
        j = 0;
        while (j < n_bounds) {
            if (occ.min < b[j].min) {
                if (occ.max < b[j].min) {
                    // No overlap -> skip
                } else if (occ.max < b[j].max) {
                    // Upper overlap -> truncate
                    b[j].min = occ.max;
                } else {
                    // Fully occluded -> remove entry
                    for (int k = j; k < n_bounds - 1; k++) {
                        b[k] = b[k + 1];
                    }
                    n_bounds -= 1;
                    continue;
                }
            } else if (occ.max < b[j].max) {
                // Contained within bounds -> split integration area
                for (k = n_bounds - 1; k >= j; k--) {
                    b[k + 1] = b[k];
                }
                b[j].max = occ.min;
                b[j + 1].min = occ.max;
                j += 1;
                n_bounds += 1;
            } else if (occ.min < b[j].max) {
                // Lower overlap -> truncate
                b[j].max = occ.min;
            }
            j += 1;
        }
    }
    return n_bounds;
}


double phase_curve_inner_integral(double x, void *params) {
    int i;
    PlanetSystem *p = (PlanetSystem *)params;
    p->x = x;

    // Determine the range(s) to integrate over
    // Max size of the bounds "stack" is low, so using an array.
    int n_bounds = 1;
    Bounds b[MAX_PHASE_OBJECTS + 1];
    // Initially, plan to integrate over the entire target
    b[0] = p->shapes[p->i_target].slice_ylimits(x);
    if (b[0].min >= b[0].max) {
        return 0.;
    }
    // Get the bounds of occluding objects for processing
    for (i = 0; i < p->n_objects; i++) {
        if (p->relevant[i]) {
            b[n_bounds] = p->shapes[i].slice_ylimits(x);
            if (b[n_bounds].min < b[n_bounds].max) {
                n_bounds += 1;
            }
        }
    }
    n_bounds = process_bounds(b, n_bounds);

    // Conduct the integration for each range identified.
    int code;
    double total = 0;
    double result, err;
    for (i = 0; i < n_bounds; i++) {
        code = gsl_integration_qag(
            &p->integInner, b[i].min, b[i].max, .1 * p->atol, .1 * p->rtol, 100, 1,
            p->workspaceInner, &result, &err
        );
        if (integration_failed(code, result, err, p->atol, p->rtol)) {
            return NAN;
        }
        total += result;
    }
    return total;
}


PlanetSystem::PlanetSystem(double atol, double rtol) {
    x = 0.;
    i_target = 0;
    n_objects = 0;
    this->atol = atol;
    this->rtol = rtol;

    // Prepare the inner (y) integral variables
    workspaceInner = gsl_integration_workspace_alloc(100);
    workspaceOuter = gsl_integration_workspace_alloc(100);
    integInner.function = &phase_curve_integrand;
    integInner.params = this;

    // Now prepare the outer (x) integral variables
    integOuter.function = &phase_curve_inner_integral;
    integOuter.params = this;

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();
};

PlanetSystem::PlanetSystem(PlanetSystem &p) {
    x = p.x;
    i_target = p.i_target;
    n_objects = p.n_objects;
    atol = p.atol;
    rtol = p.rtol;
    workspaceInner = gsl_integration_workspace_alloc(100);
    workspaceOuter = gsl_integration_workspace_alloc(100);
    integInner.function = &phase_curve_integrand;
    integInner.params = this;
    integOuter.function = &phase_curve_inner_integral;
    integOuter.params = this;
    for (int i = 0; i < n_objects; i++) {
        orbits[i] = p.orbits[i];
        shapes[i] = p.shapes[i];
        lights[i] = p.lights[i];
        rotate_with_orbit[i] = p.rotate_with_orbit[i];
        xlim[i] = Bounds();
        ylim[i] = Bounds();
    }
}

PlanetSystem::~PlanetSystem() {
    if (workspaceInner != nullptr) {
        gsl_integration_workspace_free(workspaceInner);
        workspaceInner = nullptr;
    }
    if (workspaceOuter != nullptr) {
        gsl_integration_workspace_free(workspaceOuter);
        workspaceOuter = nullptr;
    }
}

int PlanetSystem::add_object(
    const Orbit &orb, const Shape &shp, const LightSource &source, bool rot_with_orbit,
    int parent_index
) {
    // Note: These are all data-only structs, so these are copy operations
    orbits[n_objects] = orb;
    parent_indices[n_objects] = parent_index;
    shapes[n_objects] = shp;
    lights[n_objects] = source;
    rotate_with_orbit[n_objects] = rot_with_orbit;
    xlim[n_objects] = shapes[n_objects].x_bounds();
    ylim[n_objects] = shapes[n_objects].y_bounds();
    n_objects += 1;
    return 0;
}

void PlanetSystem::clear_objects() { n_objects = 0; }

void PlanetSystem::set_time(double t) {
    Vec3 origin;
    for (int i = 0; i < n_objects; i++) {
        origin = parent_indices[i] > 0 ? shapes[parent_indices[i]].position : (Vec3){0., 0., 0.};
        shapes[i].position_from_orbit(t, orbits[i], rotate_with_orbit[i], origin);
        xlim[i] = shapes[i].x_bounds();
        ylim[i] = shapes[i].y_bounds();
    }
}

double PlanetSystem::integrate_single(int it) {
    // Skip non-emitting objects
    if (lights[it].stype == None) {
        return 0.;
    }
    // Determine which objects might occlude the target object
    double result, err;
    i_target = it;
    double xmin = xlim[i_target].min;
    double xmax = xlim[i_target].max;

    // Detect trivially non-overlapping objects and mark them to be skipped
    int n_occluders = 0;
    for (int i = 0; i < n_objects; i++) {
        relevant[i] =
            // Objects cannot occlude themselves
            ((i != i_target) &&
             // Objects behind the target cannot occlude it
             (shapes[i].position.z > shapes[i_target].position.z) &&
             // Objects that don't overlap in x with the target cannot occlude it
             (xlim[i].max > xmin) && (xlim[i].min < xmax) &&
             // Objects that don't overlap in y with the target cannot occlude it
             (ylim[i].max > ylim[i_target].min) && (ylim[i].min < ylim[i_target].max));
        n_occluders += relevant[i];
    }
    if (n_occluders == 0) {
        result = lights[it].get_integrated_brightness(shapes[it]);
        if (!isnan(result)) {
            return result;
        }
    }
    int code = gsl_integration_qag(
        &integOuter, xmin, xmax, .1 * atol, .1 * rtol, 100, 1, workspaceOuter, &result, &err
    );
    if (integration_failed(code, result, err, atol, rtol)) {
        return NAN;
    }
    return result;
}

int PlanetSystem::get_n_objects() const { return n_objects; }

void PlanetSystem::phase_curve_integral(double *times, double *outputs, int n) {
    // TODO: Adjust to separately get the brightness of each component
    int i, j;
    double result = 0.;

    for (i = 0; i < n; i++) {
        set_time(times[i]);

        // Integrate each object in sequence
        result = 0.;
        for (j = 0; j < n_objects; j++) {
            result += integrate_single(j);
        }
        outputs[i] = result;
    }
    return;
}
