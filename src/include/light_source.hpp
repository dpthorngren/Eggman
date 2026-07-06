#ifndef LIGHT_SOURCE_HPP
#define LIGHT_SOURCE_HPP

#include "biellipsoid.hpp"
#define MAX_SOURCE_PARAMS 12

typedef enum {
    None,          // 0 params, no emission
    Lambertian,    // 1 params, uniform brightness regardless of angle
    QuadraticLimb, // 3 params, quadratic formula of Mandel & Agol 2002
    NonLinearLimb, // 5 params, non-linear formula of Mandel & Agol 2002
    DayNight,      // 3 params
} SourceType;


class LightSource {
  public:
    SourceType stype;
    double params[MAX_SOURCE_PARAMS];
    double limb_norm;

    LightSource();
    LightSource(SourceType type, double *params);

    double get_brightness(double mu, double sin_lat, double lon) const;
    // Wrapper for get_brightness on a unit sphere given x and y coordinates
    double get_brightness_sphere(double x, double y) const;
    // Total luminosity if the object is completely unobscured
    // Returns NAN if this cannot be computed without an actual integration
    double get_integrated_brightness(Biellipsoid &bell);
};

#endif
