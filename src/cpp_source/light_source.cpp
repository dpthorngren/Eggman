#include "light_source.hpp"
#include "math_utils.hpp"
#include <cmath>

LightSource::LightSource() {
    stype = None;
    for (int i = 0; i < MAX_SOURCE_PARAMS; i++) {
        params[i] = 0.;
    }
    limb_norm = 1.0;
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
}

double LightSource::get_brightness(double x, double y, const Shape &shp) const {
    double mu, nu, limb_coeff;
    switch (stype) {
    case None:
        return 0.;
    case Lambertian:
        return params[0];
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
    case None:
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
    default:
        return NAN;
    }
}

double LightSource::get_integrated_brightness(const Shape &shp) const {
    switch (stype) {
    case None:
        return 0.;
    case Lambertian:
        return params[0] * shp.get_area();
    case QuadraticLimb:
        return shp.is_sphere ? params[0] * shp.get_area() : NAN;
    case NonLinearLimb:
        return shp.is_sphere ? params[0] * shp.get_area() : NAN;
    default:
        return NAN;
    }
}
