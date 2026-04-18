cimport cython
from libc cimport math
from libc.math cimport M_PI as pi, NAN as nan


cdef extern from "orbit.c":
    ctypedef struct Orbit:
        pass
    void init_orbit(Orbit* orb, double period, double semimajor, double eccen, double inclination, double lon_periapse)
    double solve_kepler(double mean_anomaly, double eccen)
    void get_position(Orbit* orb, double t, double* x, double* y, double* z)


cdef extern from "biellipsoid.c":
    ctypedef struct Biellipsoid:
        double[3] position
        double[9] rot
        double[3] f1
        double[3] f2
        double[3] b1
        double[3] b2
        double[4] bounds
    void init_biellipsoid(Biellipsoid* f, double x, double y, double z, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)


cpdef (double, double, double) orbit_to_position(double t, double period, double semimajor, double eccen, double inclination, double lon_periapse) noexcept
