#ifndef ORBIT_HPP
#define ORBIT_HPP

#include "math_utils.hpp"

// Solves the kepler equation for a given mean anomaly and eccentricity
double solve_kepler(double mean_anomaly, double eccen);

class Orbit {
  private:
    double period;    // The orbital period of the planet (any consistent units)
    double semimajor; // The semimajor axis, in stellar radius units.
    double eccen;     // The orbital eccentricity
    double s_inc;     // The sin of the inclination (90 degrees is an edge-on orbit)
    double c_inc;     // The cos of the inclination (90 degrees is an edge-on orbit)
    double s_tap;     // Sin of the true anomaly at periapse
    double c_tap;     // Cos of the true anomaly at periapse
    double t_p;       // Time of periapse passage (relative to t0, same units as the period)

  public:
    Orbit();
    Orbit(
        double period, double t0, double semimajor, double eccen, double inclination,
        double lon_periapse
    );
    // Get the position of the orbiting object in view space at time t
    Vec3 get_position(double t);

    // Trivial Accessors
    double get_period();
    double get_semimajor();
    double get_eccen();
    double get_cos_inc();
};

#endif
