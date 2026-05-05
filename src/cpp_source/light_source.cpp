#include "light_source.hpp"
#include <cmath>

double LightSource::get_brightness(double x, double y) {
    // Temp helper variables
    double result, nu, sqrtnu;
    // Squared normalized radius from the source center
    double rsq = x * x + y * y;

    // Not within source disk.
    if (rsq > 1) {
        return 0.;
    }

    // Quadratic limb-darkening on a sphere (see Mandel & Agol 2002)
    if (source_type == 0) {
        nu = 1. - sqrt(fmax(1. - rsq, 0.));
        result = 1 - limb[0] * nu - limb[1] * nu * nu;
        return result / normalization;
    }
    // Non-linear limb-darkening on a sphere (see Mandel & Agol 2002)
    if (source_type == 1) {
        nu = sqrt(fmax(1 - rsq, 0.));
        sqrtnu = sqrt(nu);
        result = 1.0 - limb[0] * (1 - sqrtnu);
        result -= limb[1] * (1 - nu);
        result -= limb[2] * (1 - nu * sqrtnu);
        result -= limb[3] * (1 - nu * nu);
        return result / normalization;
    }
    // Invalid source type.
    return NAN;
}

LightSource::LightSource() {
    source_type = 0;
    limb[0] = 0.;
    limb[1] = 0.;
    limb[2] = 0.;
    limb[3] = 0.;
    normalization = M_PI;
}

LightSource::LightSource(int type, double limb0, double limb1, double limb2, double limb3) {
    // Source type should be 0 or 1 and limb params must be a length-4 array
    // Total brightness of the source (see Mandel & Agol 2002)
    this->source_type = type;
    limb[0] = limb0;
    limb[1] = limb1;
    limb[2] = limb2;
    limb[3] = limb3;
    double result = 1.;
    if (source_type == 0) {
        result -= limb[0] / 3.;
        result -= limb[1] / 6.;
        normalization = result * M_PI;
    } else if (source_type == 1) {
        result -= limb[0] / 5.;
        result -= limb[1] / 3.;
        result -= limb[2] * 3 / 7.;
        result -= limb[3] / 2.;
        normalization = result * M_PI;
    }
}
