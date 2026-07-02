#ifndef LIGHT_SOURCE_HPP
#define LIGHT_SOURCE_HPP

#include "biellipsoid.hpp"
#define MAX_SOURCE_PARAMS 12
#define MAX_LIMB_PARAMS 4

/* Describes a light source in terms of its emission map and limb darkening.
 *
 * Base source type code, number of parameters, and notes:
 * - 0, 0, No emission (also bypasses limb darkening)
 * - 1, 1, Uniform emission at a given constant brightness
 * - 2, 4, Dayside, nightside, and pole brightness, with transition size
 *
 * Limb larkening type code, number of parameters, and notes:
 * - 0, 0, Lambertian (appears uniformly bright as mu terms cancel)
 * - 1, 2, Quadratic formula of Mandel & Agol 2002
 * - 2, 4, Non-linear formula of Mandel & Agol 2002
 */
class LightSource {
  public:
    int source_type;
    double source_params[MAX_SOURCE_PARAMS];
    int limb_type;
    double limb_params[MAX_LIMB_PARAMS];
    double limb_norm;

    LightSource();
    LightSource(int source_type, double *source, int limb_type, double *limb_params);

    double get_brightness(double mu, double sin_lat, double lon) const;
    // Wrapper for get_brightness on a unit sphere given x and y coordinates
    double get_brightness_sphere(double x, double y) const;
    // Total luminosity if the object is completely unobscured
    // Returns NAN if this cannot be computed without an actual integration
    double get_integrated_brightness(Biellipsoid &bell);
};

#endif
