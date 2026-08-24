cimport cython
from libc cimport math
from libc.math cimport M_PI as pi, NAN as nan
import numpy.typing

ctypedef double[::1] Array1d_f64
ctypedef double[:,:] Array2d_f64

cdef extern from "math_utils.hpp":
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


# ===== Orbit class and wrapper =====
cdef extern from "orbit.cpp":
    double solve_kepler(double mean_anomaly, double eccen)
    cdef cppclass COrbit "Orbit":
        COrbit() except +
        COrbit(double, double, double, double, double, double) except +
        Vec3 get_position(double t) except +
        double get_period()
        double get_semimajor()
        double get_eccen()
        double get_cos_inc()


cdef class Orbit:
    cdef COrbit corbit
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly


# ===== Shape class and wrapper =====
cdef extern from "ellipse.cpp":
    cdef cppclass CEllipse "Ellipse":
        Vec3 e1
        Vec3 e2
        double det
        double x_size
        double y_size
        CEllipse()
        CEllipse(Vec3 e1, Vec3 e2)
        void get_ybounds(double x, Vec3 &out_min, Vec3 &out_max)
        bint line_intersects(double x, double y, Vec3* out)
        Vec3 nearest_to_line(double xt, double yt)


cdef class Ellipse:
    cdef CEllipse cell


cdef extern from "shape.cpp":
    cdef cppclass CShape "Shape":
        Vec3 position
        double r_forward
        double r_back
        double r_up
        double r_side
        Mat3 rot
        CEllipse f_limb
        CEllipse b_limb
        CEllipse joint

        CShape()
        CShape(double r_forward, double r_backward, double r_up, double r_side)
        void set_rotation(double theta, double phi, double gamma, double ci)
        void set_position(Vec3 new_position)
        void set_radii(double r_forward, double r_back, double r_up, double r_side)
        void position_from_orbit(double t, const COrbit &orb, bint rotate_with_orbit, Vec3 origin)
        void update_derived()
        Vec3 forward_vector()
        bint is_forward(Vec3 loc)
        bint is_forward_2d(double x, double y, bint local)
        bint is_visible(Vec3 loc)
        Bounds x_bounds()
        Bounds y_bounds()
        Bounds slice_ylimits(double x, Bounds* out2=nullptr, int zcut=0)
        bint line_intersects (double x, double y)
        bint raycast(double x, double y, double *mu_out, Vec3 *hit_out)
        double get_area()
        Vec3 world_to_aligned(Vec3 loc)
        Vec3 world_to_sphere(Vec3 loc)
        Vec3 aligned_to_world(Vec3 loc)
        Vec3 aligned_to_sphere(Vec3 loc)
        Vec3 sphere_to_world(Vec3 loc)
        Vec3 sphere_to_aligned(Vec3 loc)


cdef class Shape:
    cdef CShape cshape
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly


# ===== Light Source class and wrapper =====
cdef extern from "light_source.cpp":
    const int MAX_SOURCE_PARAMS
    ctypedef enum SourceType:
        NoEmission
        Lambertian
        QuadraticLimb
        NonLinearLimb
        DayNight
    cdef cppclass CLightSource "LightSource":
        SourceType stype
        double params[MAX_SOURCE_PARAMS]
        double limb_norm

        CLightSource()
        CLightSource(int source_type, double *params)

        double get_brightness(double x, double y, const CShape &bell)
        double get_brightness_sphere(double x, double y)


ctypedef double[MAX_SOURCE_PARAMS] Array_SourceParams

cdef class LightSource:
    cdef CLightSource csource
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly


# ===== Transit integration function and wrapper =====
cdef extern from "transit_integral.cpp":
    void transit_integral(double *times, double *outputs, int n, const COrbit &orb, const CLightSource &emitter, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side, bint rotate_with_orbit, double atol, double rtol)

cpdef object transit(double r_forward, double r_back, double r_up, double r_side, double theta, double phi, double gamma, double[::1] t, double t0, double period, double semimajor, double inclination, str limbType, object limb, double eccen=?, double lon_periapse=?, bint rotate_with_orbit=?, double atol=?, double rtol=?)


# ===== Phase Curve Class and Wrapper =====
cdef extern from "planet_system.cpp":
    const int MAX_PHASE_OBJECTS
    cdef cppclass CPlanetSystem "PlanetSystem":
        COrbit orbits[MAX_PHASE_OBJECTS]
        CShape shapes[MAX_PHASE_OBJECTS]
        CLightSource lights[MAX_PHASE_OBJECTS]
        bint rotate_with_orbit[MAX_PHASE_OBJECTS]
        Bounds xlim[MAX_PHASE_OBJECTS]
        Bounds ylim[MAX_PHASE_OBJECTS]

        CPlanetSystem()
        CPlanetSystem(CPlanetSystem &p)
        CPlanetSystem(double atol, double rtol)

        int add_object(const COrbit &orb, const CShape &bell, const CLightSource &source, bint rot_with_orbit, int parent_index)
        int get_n_objects() const
        void clear_objects()
        void set_time(double t)
        double integrate_single(int i)
        double rtol
        double atol
        void phase_curve_integral(double *times, double *outputs, int n)

cdef class PlanetSystem:
    cdef CPlanetSystem cps
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly
