#include "light_source.hpp"
#include "math_utils.hpp"
#include <cmath>

LightSource::LightSource() {
    limb_type = 0;
    source_type = 0;
    for (int i = 0; i < MAX_LIMB_PARAMS; i++) {
        limb_params[i] = 0.;
    }
    for (int i = 0; i < MAX_SOURCE_PARAMS; i++) {
        source_params[i] = 0.;
    }
    limb_norm = 1.0;
}

LightSource::LightSource(
    int source_type, double *source_params, int limb_type, double *limb_params
) {
    this->source_type = source_type;
    for (int i = 0; i < MAX_SOURCE_PARAMS; i++) {
        this->source_params[i] = source_params[i];
    }

    this->limb_type = limb_type;
    for (int i = 0; i < MAX_LIMB_PARAMS; i++) {
        this->limb_params[i] = limb_params[i];
    }

    // Normalization to apply to limb darkening calculations (see Mandel & Agol 2002)
    switch (limb_type) {
    case 0:
        limb_norm = 1.0;
        break;
    case 1:
        limb_norm = 1 - limb_params[0] / 3.;
        limb_norm -= limb_params[1] / 6.;
        break;
    case 2:
        limb_norm = 1 - limb_params[0] / 5.;
        limb_norm -= limb_params[1] / 3.;
        limb_norm -= limb_params[2] * 3 / 7.;
        limb_norm -= limb_params[3] / 2.;
        break;
    default:
        limb_norm = NAN;
        break;
    }
}

double LightSource::get_brightness(double mu, double sin_lat, double lon) const {
    // See light_source.hpp for the type codes used here
    double base, limb_coeff, nu, sqrtmu;

    switch (source_type) {
    case 0: // No emission
        return 0.;
    case 1: // Uniform emission
        base = source_params[0];
        break;
    case 2: // Dayside, nightside and pole brightness
        // TODO (placeholder)
        base = 1.0 / M_PI;
        break;
    default:
        return NAN;
    }

    mu = CLAMP(mu, 0., 1.);
    switch (limb_type) {
    case 0: // No limb darkening
        limb_coeff = 1.0;
        break;
    case 1: // Quadratic darkening
        nu = 1 - mu;
        limb_coeff = 1 - limb_params[0] * nu - limb_params[1] * nu * nu;
        break;
    case 2: // Non-linear darkening
        sqrtmu = sqrt(mu);
        limb_coeff = 1.0 - limb_params[0] * (1 - sqrtmu);
        limb_coeff -= limb_params[1] * (1 - mu);
        limb_coeff -= limb_params[2] * (1 - mu * sqrtmu);
        limb_coeff -= limb_params[3] * (1 - mu * mu);
        break;
    default:
        return NAN;
    }
    return base * limb_coeff / limb_norm;
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
    switch (source_type) {
    case 0: // No emission
        return 0.;
    case 1: // Flat emission
        return source_params[0] * bell.get_area();
    }
    return NAN;
}
