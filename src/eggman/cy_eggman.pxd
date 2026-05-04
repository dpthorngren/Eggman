cimport cython
from libc cimport math
from libc.math cimport M_PI as pi, NAN as nan


cdef extern from "eggman.h":
    ctypedef struct Orbit:
        pass
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
    ctypedef struct LightSource:
        int source_type
        double limb_params[4]


# Orbit functions and wrappers
cdef extern from "orbit.c":
    Orbit create_orbit(double period, double semimajor, double eccen, double inclination, double lon_periapse)
    double solve_kepler(double mean_anomaly, double eccen)
    Vec3 get_position(Orbit* orb, double t)

cpdef (double, double, double) orbit_to_position(double t, double period, double semimajor, double eccen, double inclination, double lon_periapse) noexcept


# Biellipsoid functions and wrappers
cdef extern from "biellipsoid.c":
    Biellipsoid create_biellipsoid(Vec3 position, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)
    Bounds get_bounds(Biellipsoid* bell, int axis)
    Bounds get_ylimits(Biellipsoid* bell, double x)

cpdef dict biellipsoid_dump(double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)

cpdef (double, double) _biellipsoid_ylimits(double x, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)


# Light source functions and wrappers
cdef extern from "light_source.c":
    LightSource create_light_source(int source_type, double limb_params[4])
    double get_source_brightness(LightSource* source, double x, double y)

cpdef double get_brightness(int source_type, list[double] limb_params, double x, double y)

# 2-d Integration function and wrapper
cdef extern from "transit2d.c":
    void transit2d_integral(double *times, double *out_depths, int n, Orbit *orb, LightSource* emitter, double r_forward, double r_back, double r_up)

cpdef object transit2d(double rMorning, double rEvening, double rPole, double[::1] t, double t0, double period, double semimajor, double inclination, str limbType, object limb, double eccen=?, double lonPeriapse=?)

# 3-d Integration function and wrapper
cdef extern from "transit3d.c":
    void transit3d_integral(double *times, double *outputs, int n, Orbit *orb, LightSource *emitter, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side)

cpdef object transit3d(double r_forward, double r_back, double rPole, double rSide, double theta, double phi, double gamma, double[::1] t, double t0, double period, double semimajor, double inclination, str limbType, object limb, double eccen=?, double lon_periapse=?)
