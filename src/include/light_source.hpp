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
    EmissionMap,   // 2 params
} SourceType;


class LightSource {
  private:
    double *emission_map;
    int n;
    int m;

  public:
    SourceType stype;
    double params[MAX_SOURCE_PARAMS];
    double limb_norm;

    LightSource();
    LightSource(SourceType type, double *params);
    LightSource(const LightSource &other);
    LightSource &operator=(const LightSource &other);
    ~LightSource();

    // Brightness of this source for the point on Biellipsoid shp through x, y (view space)
    double get_brightness(double x, double y, const Shape &shp) const;
    // Wrapper for get_brightness on a unit sphere given x and y coordinates
    double get_brightness_sphere(double x, double y) const;
    // Total luminosity if the object is completely unobscured
    // Returns NAN if this cannot be computed without an actual integration
    double get_integrated_brightness(const Shape &shp) const;

    // Emission Map Accessors
    inline int get_n() const { return n; }
    inline int get_m() const { return m; }
    inline int get_map_size() const { return n * m + 2; }
    double get_emission_point(int i);
    int set_emission_point(int i, double value);
    Vec3 get_emission_location(int i);
    double interp_emission(Vec3 loc) const { return interp_gridmap(loc, n, m, emission_map); }
};

#endif
