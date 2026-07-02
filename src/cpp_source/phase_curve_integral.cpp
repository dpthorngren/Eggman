#include "phase_curve_integral.hpp"
#include "math_utils.hpp"
#include <cmath>

double phase_curve_integrand(double y, void *params) {
    PhaseIntegrator *p = (PhaseIntegrator *)params;
    Vec3 pos = p->shapes[p->i_target].line_project(p->x, y, true);
    return p->lights[p->i_target].get_brightness(pos.x, pos.y, pos.z);
}

double phase_curve_inner_integral(double x, void *params) {
    PhaseIntegrator *p = (PhaseIntegrator *)params;
    p->x = x;
    // Determine the range(s) to integrate over
    // Max size of the bounds "stack" is low, so using an array.
    Bounds b[MAX_PHASE_OBJECTS + 1];
    // Initially, plan to integrate over the entire target
    b[0] = p->shapes[p->i_target].slice_ylimits(x);
    int n = 1;
    if (b[0].min >= b[0].max) {
        return 0.;
    }

    // Iterate through potential occluders, modifying the integration bounds stack as-needed
    int i, j, k;
    Bounds occ;
    for (i = 0; i < p->n_objects; i++) {
        // Skip occluders flagged as irrelevant in integrate_single(...)
        if (!p->relevant[i]) {
            continue;
        }

        occ = p->shapes[i].slice_ylimits(x);
        j = 0;
        while (j < n) {
            if (occ.min < b[j].min) {
                if (occ.max < b[j].min) {
                    // No overlap -> skip
                } else if (occ.max < b[j].max) {
                    // Upper overlap -> truncate
                    b[i].min = occ.max;
                } else {
                    // Fully occluded -> remove entry
                    for (int k = j; k < n - 1; k++) {
                        b[k] = b[k + 1];
                    }
                    n -= 1;
                    continue;
                }
            } else if (occ.max < b[j].max) {
                // Contained within bounds -> split integration area
                for (k = j; k < n; k++) {
                    b[k + 1] = b[k];
                }
                b[j].max = occ.min;
                b[j + 1].min = occ.max;
                j += 1;
                n += 1;
            } else if (occ.min < b[j].max) {
                // Lower overlap -> truncate
                b[i].max = occ.min;
            }
            j += 1;
        }
    }

    // Conduct the integration for each range identified.
    int code;
    double total = 0;
    double result, err;
    for (i = 0; i < n; i++) {
        code = gsl_integration_qag(
            &p->integInner, b[i].min, b[i].max, 1e-8, 1e-8, 100, 1, p->workspaceInner, &result, &err
        );
        if (integration_failed(code, result, err, 1e-9, 1e-7)) {
            return NAN;
        }
        total += result;
    }
    return total;
}


PhaseIntegrator::PhaseIntegrator() {
    x = 0.;
    i_target = 0;
    n_objects = 0;

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

PhaseIntegrator::PhaseIntegrator(PhaseIntegrator &p) {
    x = p.x;
    i_target = p.i_target;
    n_objects = p.n_objects;
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

PhaseIntegrator::~PhaseIntegrator() {
    if (workspaceInner != nullptr) {
        gsl_integration_workspace_free(workspaceInner);
        workspaceInner = nullptr;
    }
    if (workspaceOuter != nullptr) {
        gsl_integration_workspace_free(workspaceOuter);
        workspaceOuter = nullptr;
    }
}

int PhaseIntegrator::add_object(
    const Orbit &orb, const Biellipsoid &bell, const LightSource &source, bool rot_with_orbit
) {
    // Note: These are all data-only structs, so these are copy operations
    orbits[n_objects] = orb;
    shapes[n_objects] = bell;
    lights[n_objects] = source;
    rotate_with_orbit[n_objects] = rot_with_orbit;
    xlim[n_objects] = shapes[n_objects].x_bounds();
    ylim[n_objects] = shapes[n_objects].y_bounds();
    n_objects += 1;
    return 0;
}

void PhaseIntegrator::clear_objects() { n_objects = 0; }

void PhaseIntegrator::set_time(double t) {
    for (int i = 0; i < n_objects; i++) {
        shapes[i].position_from_orbit(t, orbits[i]);
        xlim[i] = shapes[i].x_bounds();
        ylim[i] = shapes[i].y_bounds();
    }
}

double PhaseIntegrator::integrate_single(int it) {
    // Skip non-emitting objects
    if (lights[it].source_type == 0) {
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
        if (!isnan(result)){
            return result;
        }
    }
    int code = gsl_integration_qag(
        &integOuter, xmin, xmax, 1e-8, 1e-8, 100, 1, workspaceOuter, &result, &err
    );
    if (integration_failed(code, result, err, 1e-9, 1e-7)) {
        return NAN;
    }
    return result;
}

int PhaseIntegrator::get_n_objects() const { return n_objects; }

void PhaseIntegrator::phase_curve_integral(double *times, double *outputs, int n) {
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
