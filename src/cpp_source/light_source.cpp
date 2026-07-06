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

double LightSource::get_brightness(double mu, double sin_lat, double lon) const {
    // See light_source.hpp for the type codes used here
    double base, limb_coeff, nu, sqrtmu;

    mu = CLAMP(mu, 0., 1.);
    switch (stype) {
    case None:
        return 0.;
    case Lambertian:
        return params[0];
    case QuadraticLimb:
        nu = 1 - mu;
        limb_coeff = 1 - params[1] * nu - params[2] * nu * nu;
        return params[0] * limb_coeff / limb_norm;
    case NonLinearLimb:
        sqrtmu = sqrt(mu);
        limb_coeff = 1.0 - params[1] * (1 - sqrtmu);
        limb_coeff -= params[2] * (1 - mu);
        limb_coeff -= params[3] * (1 - mu * sqrtmu);
        limb_coeff -= params[4] * (1 - mu * mu);
        return params[0] * limb_coeff / limb_norm;
    default:
        return NAN;
    }
}

double LightSource::get_brightness_sphere(double x, double y) const {
    double rsq = x * x + y * y;
    if (rsq > 1.) {
        // Point is not on the sphere
        return 0.;
    }
    double mu = sqrt(1 - rsq);
    // In view coordinates, the unit sphere has:
    // x, y, z = [cos(lon) cos(lat), sin(lat), sin(lon) cos(lat)]
    // x^2 + y^2 + z^2 = 1
    // tan(lon) = z / x = +/- sqrt(1 - x^2 - y^2) / x
    // Negative case faces viewer, so:
    double lon = atan2(-sqrt(1 - rsq), x);
    // For the unit sphere, y is the sin(latitude)
    return get_brightness(mu, y, lon);
}

double LightSource::get_integrated_brightness(Biellipsoid &bell) {
    switch (stype) {
    case None:
        return 0.;
    case Lambertian:
        return params[0] * bell.get_area();
    case QuadraticLimb:
        return params[0] * bell.get_area();
    case NonLinearLimb:
        return params[0] * bell.get_area();
    default:
        return NAN;
    }
}
