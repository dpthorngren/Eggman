#include "light_source.hpp"
#include "math_utils.hpp"
#include <cmath>

LightSource::LightSource() {
    stype = NoEmission;
    for (int i = 0; i < MAX_SOURCE_PARAMS; i++) {
        params[i] = 0.;
    }
    limb_norm = 1.0;
    m = 0;
    n = 0;
}

LightSource::LightSource(SourceType type, double *params) {
    this->stype = type;
    for (int i = 0; i < MAX_SOURCE_PARAMS; i++) {
        this->params[i] = params[i];
    }
    // Normalization to apply to limb darkening calculations (see Mandel & Agol 2002)
    switch (type) {
    case QuadraticLimb:
        limb_norm = 1 - params[1] / 3.;
        limb_norm -= params[2] / 6.;
        break;
    case NonLinearLimb:
        limb_norm = 1 - params[1] / 5.;
        limb_norm -= params[2] / 3.;
        limb_norm -= params[3] * 3 / 7.;
        limb_norm -= params[4] / 2.;
        break;
    default:
        limb_norm = 1.0;
        break;
    }

    // Set up emission map if needed
    if (stype == EmissionMap) {
        n = (int)fmax(params[0], 1);
        m = (int)fmax(params[1], 1);
        emission_map.reserve(n * m + 2);
        for (int i = 0; i < n * m + 2; i++) {
            emission_map.push_back(0.);
        }
    } else {
        n = 0;
        m = 0;
    }
}

double LightSource::get_brightness(double x, double y, const Shape &shp) const {
    double mu, nu, limb_coeff;
    Vec3 hit;
    switch (stype) {
    case NoEmission:
        return 0.;
    case Lambertian:
        return shp.line_intersects(x, y) ? params[0] : 0.;
    case QuadraticLimb:
        if (!shp.raycast(x, y, &mu, nullptr)) {
            return 0; // Point doesn't intersect biellipsoid
        }
        nu = 1 - mu;
        limb_coeff = 1 - params[1] * nu - params[2] * nu * nu;
        return params[0] * limb_coeff / limb_norm;
    case NonLinearLimb:
        if (!shp.raycast(x, y, &mu, nullptr)) {
            return 0; // Point doesn't intersect biellipsoid
        }
        nu = sqrt(mu);
        limb_coeff = 1.0 - params[1] * (1 - nu);
        limb_coeff -= params[2] * (1 - mu);
        limb_coeff -= params[3] * (1 - mu * nu);
        limb_coeff -= params[4] * (1 - mu * mu);
        return params[0] * limb_coeff / limb_norm;
    case DayNight:
        if (!shp.raycast(x, y, nullptr, &hit)) {
            return 0; // Point doesn't intersect biellipsoid
        }
        return hit.z < 0 ? params[0] : params[1];
    case EmissionMap:
        if (!shp.raycast(x, y, nullptr, &hit)) {
            return 0; // Point doesn't intersect biellipsoid
        }
        return interp_emission(hit);
    default:
        return NAN;
    }
}


double LightSource::get_brightness_sphere(double x, double y) const {
    double mu, nu, limb_coeff;
    double rsq = x * x + y * y;
    if (rsq > 1.) {
        return 0.; // Point is not on the sphere
    }
    switch (stype) {
    case NoEmission:
        return 0.;
    case Lambertian:
        return params[0];
    case QuadraticLimb:
        nu = 1 - CLAMP(sqrt(1 - rsq), 0., 1.0);
        limb_coeff = 1 - params[1] * nu - params[2] * nu * nu;
        return params[0] * limb_coeff / limb_norm;
    case NonLinearLimb:
        mu = CLAMP(sqrt(1 - rsq), 0., 1.0);
        nu = sqrt(mu);
        limb_coeff = 1.0 - params[1] * (1 - nu);
        limb_coeff -= params[2] * (1 - mu);
        limb_coeff -= params[3] * (1 - mu * nu);
        limb_coeff -= params[4] * (1 - mu * mu);
        return params[0] * limb_coeff / limb_norm;
    case DayNight:
        // On an unrotated unit sphere, only the night side is visible
        return params[1];
    default:
        return NAN;
    }
}

double LightSource::get_integrated_brightness(const Shape &shp) const {
    double result;
    switch (stype) {
    case NoEmission:
        return 0.;
    case Lambertian:
        return params[0] * shp.get_area();
    case QuadraticLimb:
        return shp.shape_type == Sphere ? params[0] * shp.get_area() : NAN;
    case NonLinearLimb:
        return shp.shape_type == Sphere ? params[0] * shp.get_area() : NAN;
    case DayNight:
        if (shp.shape_type == Sphere) {
            result = params[0] + params[1] + (params[1] - params[0]) * shp.rot.zz;
            return result * shp.get_area() / 2.0;
        }
        return NAN;
    default:
        return NAN;
    }
}

double LightSource::get_emission_point(int i) {
    if ((i < emission_map.size()) && (i >= 0)) {
        return emission_map[i];
    }
    return NAN;
}

int LightSource::set_emission_point(int i, double value) {
    if ((i < emission_map.size()) && (i >= 0)) {
        emission_map[i] = value;
        return 0;
    }
    return 1;
}

Vec3 LightSource::get_emission_location(int i) {
    if ((i > n * m + 1) || (i < 0)) {
        return {0., 0., 0.};
    } else if (i == n * m) {
        return {0., 0., -1.};
    } else if (i == n * m + 1) {
        return {0., 0., 1.};
    }
    // Yes the integer roundoff is intended.
    int j = i / n;
    i = i % n;
    Vec3 result = {(double)i, (double)j, 0.};
    gridmap_to_sphere(result, n, m);
    return result;
}
