#ifndef LIGHT_SOURCE_HPP
#define LIGHT_SOURCE_HPP

#include "shape.hpp"
#define MAX_SOURCE_PARAMS 12

typedef enum {
    NoEmission,    // 0 params, no emission
    Lambertian,    // 1 params, uniform brightness regardless of angle
    QuadraticLimb, // 3 params, quadratic formula of Mandel & Agol 2002
    NonLinearLimb, // 5 params, non-linear formula of Mandel & Agol 2002
    DayNight,      // 2 params
} SourceType;


class LightSource {
  public:
    SourceType stype;
    double params[MAX_SOURCE_PARAMS];
    double limb_norm;

    LightSource();
    LightSource(SourceType type, double *params);

    // Brightness of this source for the point on Biellipsoid shp through x, y (view space)
    double get_brightness(double x, double y, const Shape &shp) const;
    // Wrapper for get_brightness on a unit sphere given x and y coordinates
    double get_brightness_sphere(double x, double y) const;
    // Total luminosity if the object is completely unobscured
    // Returns NAN if this cannot be computed without an actual integration
    double get_integrated_brightness(const Shape &shp) const;
};

#endif
