#include "orbit.hpp"
#include <cmath>

double solve_kepler(double mean_anomaly, double eccen) {
    double e_anom = mean_anomaly;
    double e_new = mean_anomaly;
    int i = 0;
    if (eccen > 0.95) {
        // Direct iteration
        for (i = 0; i < 50; i++) {
            e_new = mean_anomaly + eccen * sin(e_anom);
            if (fabs(e_new - e_anom) < 1e-12)
                break;
            e_anom = e_new;
        }
    }
    if (eccen > 0.) {
        // Newton-Raphson
        for (i = 0; i < 50; i++) {
            e_new -= (e_anom - eccen * sin(e_anom) - mean_anomaly) / (1 - eccen * cos(e_anom));
            if (fabs(e_new - e_anom) < 1e-12)
                break;
            e_anom = e_new;
        }
    }
    return e_new;
}

Orbit::Orbit() {
    period = 1.;
    semimajor = 0.;
    eccen = 0.;
    s_inc = 1.;
    c_inc = 0.;
    s_tap = 0.;
    c_tap = 1.;
    t_p = 0.;
}

Orbit::Orbit(
    double period, double t0, double semimajor, double eccen, double inclination,
    double lon_periapse
) {
    this->period = period;
    this->semimajor = semimajor;
    this->eccen = eccen;
    // Precompute inclination sin and cosine
    inclination *= M_PI / 180.;
    s_inc = sin(inclination);
    c_inc = cos(inclination);
    // Precompute periapse data
    double ta_p = (lon_periapse - 90.) * M_PI / 180.;
    s_tap = sin(ta_p);
    c_tap = cos(ta_p);
    double ea_p = atan2(sqrt(1 - eccen * eccen) * s_tap, eccen + c_tap);
    double ma_p = ea_p - eccen * sin(ea_p);
    t_p = ma_p * period / (2 * M_PI) + t0;
}

Vec3 Orbit::get_position(double t) {
    // Skip objects locked at the origin (convenient for other objects)
    if (semimajor == 0.) {
        return {0., 0., 0.};
    }
    // Solve Kepler's equation in the rotated and inclined frame, relative to
    // inferior conjunction
    double ma = 2 * M_PI * (t - t_p) / period;
    double ea = solve_kepler(ma, eccen);
    // Position in the rotated frame where lon_periapse = 0, inclination = 0
    double x_rot = semimajor * sqrt(1 - eccen * eccen) * sin(ea);
    double y_rot = semimajor * (cos(ea) - eccen);
    // Transform back to frame where lon_periapse != 0 (still inclination = 0)
    double x_inc = x_rot * c_tap + y_rot * s_tap;
    double y_inc = x_rot * -s_tap + y_rot * c_tap;
    // Transform back to the view frame (x = right, y = up, z = towards observer)
    // centered on star
    return {x_inc, y_inc * c_inc, y_inc * s_inc};
}

// Trivial Accessors
double Orbit::get_period() { return period; }
double Orbit::get_semimajor() { return semimajor; }
double Orbit::get_eccen() { return eccen; }
double Orbit::get_cos_inc() { return c_inc; }
