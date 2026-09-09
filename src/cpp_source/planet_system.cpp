#include "planet_system.hpp"
#include "math_utils.hpp"
#include <cmath>

double emission_integrand(double y, void *params) {
    PlanetSystem *p = (PlanetSystem *)params;
    return p->lights[p->i_target].get_brightness(p->x, y, p->shapes[p->i_target]);
}

int process_bounds(Bounds *b, int n_relevant, bool invert) {
    // Iterate through potential occluders, modifying the integration bounds stack as-needed
    // Bound stack 'b' should be at least length n_relevant and contain the overall range
    // in position [0], followed by the bounds of the potential occluders.
    int i, j, k;
    Bounds occ;
    Bounds overall = b[0];
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
    double last_max = b[0].max;
    double this_max;
    if (invert) {
        i = 1; // Read position
        j = 0; // Write position <= i
        // Unoccluded range reaches lower bound?
        if (b[0].min > overall.min) {
            b[j++] = {overall.min, b[0].min};
        }
        while (i < n_bounds) {
            this_max = b[i].max;
            b[j++] = {last_max, b[i++].min};
            last_max = this_max;
        }
        // Unoccluded range reaches upper bound?
        if (last_max < overall.max) {
            b[j++] = {last_max, overall.max};
        }
        n_bounds = j;
    }
    return n_bounds;
}


double emission_outer_integral(double x, void *params) {
    int i;
    PlanetSystem *p = (PlanetSystem *)params;
    p->x = x;
    int zcut;

    // Determine the range(s) to integrate over
    // Max size of the bounds "stack" is low, so using an array.
    int n_bounds = 1;
    Bounds b[2 * MAX_SYSTEM_OBJECTS + 1];
    // Initially, plan to integrate over the entire target
    b[0] = p->shapes[p->i_target].slice_ylimits(x);
    if (b[0].min >= b[0].max) {
        return 0.;
    }
    // Get the bounds of occluding objects for processing
    for (i = 0; i < p->n_objects; i++) {
        if (p->relevant[i]) {
            if (p->shapes[i].shape_type == Ring) {
                b[n_bounds + 1] = {0., 0.};
                zcut = p->parent_indices[i] == p->i_target ? 1 : 0;
                b[n_bounds] = p->shapes[i].slice_ylimits(x, &(b[n_bounds + 1]), zcut);
                if (b[n_bounds].min < b[n_bounds].max) {
                    n_bounds += 1;
                }
            } else {
                b[n_bounds] = p->shapes[i].slice_ylimits(x);
            }
            if (b[n_bounds].min < b[n_bounds].max) {
                n_bounds += 1;
            }
        }
    }
    n_bounds = process_bounds(b, n_bounds, p->invert_integral);

    // Conduct the integration for each range identified.
    int code;
    double total = 0;
    double result, err;
    for (i = 0; i < n_bounds; i++) {
        code = gsl_integration_qag(
            &p->integInner, b[i].min, b[i].max, .1 * p->atol, .1 * p->rtol, p->max_steps, 1,
            p->workspaceInner, &result, &err
        );
        if (integration_failed(code, result, err, p->atol, p->rtol)) {
            return NAN;
        }
        total += result;
    }
    return total;
}


PlanetSystem::PlanetSystem(double atol, double rtol, int max_steps) {
    x = 0.;
    invert_integral = false;
    i_target = 0;
    n_objects = 0;
    this->atol = atol;
    this->rtol = rtol;

    // Prepare the inner (y) integral variables
    max_steps = CLAMP(max_steps, 10, 5000);
    this->max_steps = max_steps;
    workspaceInner = gsl_integration_workspace_alloc(max_steps);
    workspaceOuter = gsl_integration_workspace_alloc(max_steps);
    integInner.function = &emission_integrand;
    integInner.params = this;

    // Now prepare the outer (x) integral variables
    integOuter.function = &emission_outer_integral;
    integOuter.params = this;

    // Do not crash the program due to lack of precision
    gsl_set_error_handler_off();
};

PlanetSystem::PlanetSystem(PlanetSystem &p) {
    x = p.x;
    invert_integral = p.invert_integral;
    i_target = p.i_target;
    n_objects = p.n_objects;
    atol = p.atol;
    rtol = p.rtol;
    max_steps = p.max_steps;
    workspaceInner = gsl_integration_workspace_alloc(max_steps);
    workspaceOuter = gsl_integration_workspace_alloc(max_steps);
    integInner.function = &emission_integrand;
    integInner.params = this;
    integOuter.function = &emission_outer_integral;
    integOuter.params = this;
    for (int i = 0; i < n_objects; i++) {
        orbits[i] = p.orbits[i];
        shapes[i] = p.shapes[i];
        lights[i] = p.lights[i];
        rotate_with_orbit[i] = p.rotate_with_orbit[i];
        xlim[i] = p.xlim[i];
        ylim[i] = p.ylim[i];
        lum_cache[i] = p.lum_cache[i];
    }
}

PlanetSystem &PlanetSystem::operator=(const PlanetSystem &p) {
    if (&p == this) {
        return *this;
    }
    x = p.x;
    invert_integral = p.invert_integral;
    i_target = p.i_target;
    n_objects = p.n_objects;
    atol = p.atol;
    rtol = p.rtol;
    if (max_steps != p.max_steps || workspaceInner == nullptr || workspaceOuter == nullptr) {
        max_steps = p.max_steps;
        if (workspaceInner != nullptr) {
            gsl_integration_workspace_free(workspaceInner);
            workspaceInner = nullptr;
        }
        if (workspaceOuter != nullptr) {
            gsl_integration_workspace_free(workspaceOuter);
            workspaceOuter = nullptr;
        }
        workspaceInner = gsl_integration_workspace_alloc(max_steps);
        workspaceOuter = gsl_integration_workspace_alloc(max_steps);
    }
    integInner.function = &emission_integrand;
    integInner.params = this;
    integOuter.function = &emission_outer_integral;
    integOuter.params = this;
    for (int i = 0; i < n_objects; i++) {
        orbits[i] = p.orbits[i];
        shapes[i] = p.shapes[i];
        lights[i] = p.lights[i];
        rotate_with_orbit[i] = p.rotate_with_orbit[i];
        xlim[i] = p.xlim[i];
        ylim[i] = p.ylim[i];
        lum_cache[i] = p.lum_cache[i];
    }
    return *this;
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
    orbits[n_objects] = orb;
    parent_indices[n_objects] = parent_index;
    shapes[n_objects] = shp;
    lights[n_objects] = source;
    rotate_with_orbit[n_objects] = rot_with_orbit;
    xlim[n_objects] = shapes[n_objects].x_bounds();
    ylim[n_objects] = shapes[n_objects].y_bounds();
    lum_cache[n_objects] = NAN;
    n_objects += 1;
    return 0;
}

void PlanetSystem::clear_objects() { n_objects = 0; }

void PlanetSystem::set_time(double t) {
    Vec3 origin;
    for (int i = 0; i < n_objects; i++) {
        origin = parent_indices[i] >= 0 ? shapes[parent_indices[i]].position : (Vec3){0., 0., 0.};
        shapes[i].position_from_orbit(t, orbits[i], rotate_with_orbit[i], origin);
        xlim[i] = shapes[i].x_bounds();
        ylim[i] = shapes[i].y_bounds();
    }
}

double PlanetSystem::integrate_unoccluded_single(int it, bool may_integrate) {
    // Skip non-emitting objects
    if (lights[it].stype == NoEmission) {
        return 0.;
    }
    // Used previously cached result, if available
    if (isfinite(lum_cache[it])) {
        return lum_cache[it];
    }
    // Use closed-form brightness if available
    double result = lights[it].get_integrated_brightness(shapes[it]);
    if (!isnan(result)) {
        return result;
    }
    if (!may_integrate) {
        return NAN;
    }

    // Temporarily mark all objects as irrelevant (not occluding)
    bool prev_relevant[MAX_SYSTEM_OBJECTS];
    for (int i = 0; i < n_objects; i++) {
        prev_relevant[i] = relevant[i];
        relevant[i] = false;
    }

    // Conduct the integral
    i_target = it;
    double err;
    int code = gsl_integration_qag(
        &integOuter, xlim[it].min, xlim[it].max, .1 * atol, .1 * rtol, max_steps, 1, workspaceOuter,
        &result, &err
    );
    if (integration_failed(code, result, err, atol, rtol)) {
        return NAN;
    }

    // May cache result if object doesn't rotate
    if (!rotate_with_orbit[it]) {
        lum_cache[it] = result;
    }

    // Restore relevancy array to prrevious state
    for (int i = 0; i < n_objects; i++) {
        relevant[i] = prev_relevant[i];
    }

    return result;
}

double PlanetSystem::integrate_single(int it) {
    // Skip non-emitting objects
    if (lights[it].stype == NoEmission) {
        return 0.;
    }
    // Determine which objects might occlude the target object
    double result, err, baseline_flux;
    double area = 0.;
    i_target = it;
    invert_integral = false;
    double xmin = xlim[it].min;
    double xmax = xlim[it].max;

    // Detect trivially non-overlapping objects and mark them to be skipped
    int n_occluders = 0;
    for (int i = 0; i < n_objects; i++) {
        relevant[i] =
            // Objects cannot occlude themselves
            ((i != it) &&
             // Objects behind the target cannot occlude it
             (shapes[i].position.z > shapes[it].position.z) &&
             // Objects that don't overlap in x with the target cannot occlude it
             (xlim[i].max > xmin) && (xlim[i].min < xmax) &&
             // Objects that don't overlap in y with the target cannot occlude it
             (ylim[i].max > ylim[it].min) && (ylim[i].min < ylim[it].max));
        if (relevant[i]) {
            n_occluders += 1;
            area += shapes[i].get_area();
            // If object is fully occluded, no need to integrate
            if (shapes[i].fully_contains(shapes[it])) {
                return 0.;
            }
        }
    }

    if (n_occluders == 0) {
        return integrate_unoccluded_single(it, true);
    }

    // If the occluder's area is small and the total brightness is available,
    // integrate only the occluded area, and subtract from the total.
    if (area < 0.5 * shapes[it].get_area()) {
        invert_integral = true;
        baseline_flux = integrate_unoccluded_single(it, true);
        // Since we're integrating occluded area, recalculate x bounds.
        xmin = INFINITY;
        xmax = -INFINITY;
        for (int i = 0; i < n_objects; i++) {
            if (relevant[i]) {
                xmin = fmin(xlim[i].min, xmin);
                xmax = fmax(xlim[i].max, xmax);
            }
        }
        xmin = fmax(xmin, xlim[it].min);
        xmax = fmin(xmax, xlim[it].max);
    }

    // Conduct the integral
    int code = gsl_integration_qag(
        &integOuter, xmin, xmax, .1 * atol, .1 * rtol, max_steps, 1, workspaceOuter, &result, &err
    );
    if (integration_failed(code, result, err, atol, rtol)) {
        return NAN;
    }
    if (invert_integral) {
        return baseline_flux - result;
    }
    return result;
}

int PlanetSystem::get_n_objects() const { return n_objects; }

void PlanetSystem::integrate(double *times, double *outputs, int n) {
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

void PlanetSystem::reset_cache() {
    for (int i = 0; i < n_objects; i++) {
        lum_cache[i] = NAN;
    }
    return;
}
