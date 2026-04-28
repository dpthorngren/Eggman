cimport cython
from libc cimport math
from libc.math cimport M_PI as pi, NAN as nan


cdef extern from "orbit.c":
    ctypedef struct Orbit:
        pass
    Orbit create_orbit(double period, double semimajor, double eccen, double inclination, double lon_periapse)
    double solve_kepler(double mean_anomaly, double eccen)
    void get_position(Orbit* orb, double t, double* x, double* y, double* z)


cdef extern from "biellipsoid.c":
    ctypedef struct Mat3:
        double xx
        double xy
        double xz
        double yx
        double yy
        double yz
        double zx
        double zy
        double zz
    ctypedef struct Vec3:
        double x
        double y
        double z
    ctypedef struct Bounds:
        double min
        double max
    ctypedef struct Biellipsoid:
        Vec3 position
        double[4] radii
        Mat3 rot
        Vec3 f_limb
        Vec3 b_limb
        Vec3 break_offset
        Vec3 f1
        Vec3 f2
        Vec3 b1
        Vec3 b2
        Bounds xbounds
        Bounds ybounds
    Biellipsoid create_biellipsoid(double x, double y, double z, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)
    Bounds get_bounds(Biellipsoid* bell, int axis)
    Bounds get_ylimits(Biellipsoid* bell, double x)

cpdef dict biellipsoid_dump(double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)

cpdef (double, double, double) orbit_to_position(double t, double period, double semimajor, double eccen, double inclination, double lon_periapse) noexcept

cpdef (double, double) _biellipsoid_ylimits(double x, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)
