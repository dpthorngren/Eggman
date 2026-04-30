#include "eggman.h"

double get_source_brightness(LightSource *source, double x, double y) {
    // Temp helper variables
    double result, nu, sqrtnu;
    // Squared normalized radius from the source center
    double rsq = x * x + y * y;

    // Not within source disk.
    if (rsq > 1) {
        return 0.;
    }

    // Quadratic limb-darkening on a sphere (see Mandel & Agol 2002)
    if (source->source_type == 0) {
        nu = 1. - sqrt(fmax(1. - rsq, 0.));
        result = 1 - source->limb_params[0] * nu - source->limb_params[1] * nu * nu;
        return result / source->normalization;
    }
    // Non-linear limb-darkening on a sphere (see Mandel & Agol 2002)
    if (source->source_type == 1) {
        nu = sqrt(fmax(1 - rsq, 0.));
        sqrtnu = sqrt(nu);
        result = 1.0 - source->limb_params[0] * (1 - sqrtnu);
        result -= source->limb_params[1] * (1 - nu);
        result -= source->limb_params[2] * (1 - nu * sqrtnu);
        result -= source->limb_params[3] * (1 - nu * nu);
        return result / source->normalization;
    }
    // Invalid source type.
    return NAN;
}

LightSource create_light_source(int source_type, double limb_params[4]) {
    // Source type should be 0 or 1 and limb params must be a length-4 array
    // Total brightness of the source (see Mandel & Agol 2002)
    LightSource s = {source_type,    limb_params[0], limb_params[1],
                     limb_params[2], limb_params[3], 0.};
    double result = 1.;
    if (source_type == 0) {
        result -= limb_params[0] / 3.;
        result -= limb_params[1] / 6.;
        s.normalization = result * M_PI;
    } else if (source_type == 1) {
        result -= limb_params[0] / 5.;
        result -= limb_params[1] / 3.;
        result -= limb_params[2] * 3 / 7.;
        result -= limb_params[3] / 2.;
        s.normalization = result * M_PI;
    }
    return s;
}
