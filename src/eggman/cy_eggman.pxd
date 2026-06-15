cimport cython
from libc cimport math
from libc.math cimport M_PI as pi, NAN as nan


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
    cdef cppclass Orbit:
        Orbit() except +
        Orbit(double, double, double, double, double, double) except +
        Vec3 get_position(double t) except +
        double get_period()
        double get_semimajor()
        double get_eccen()
        double get_cos_inc()


cdef class OrbitWrap:
    cdef Orbit orb
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly


# ===== Biellipsoid class and wrapper =====
cdef extern from "ellipse.cpp":
    cdef cppclass Ellipse:
        Vec3 e1
        Vec3 e2
        double x_size
        double y_size
        Ellipse()
        Ellipse(Vec3 e1, Vec3 e2)
        void get_ybounds(double x, Vec3 *out_min, Vec3 *out_max)
        bint line_intersects(double x, double y, Vec3* out)
        Vec3 nearest_to_line(double xt, double yt)


cdef class EllipseWrap:
    cdef Ellipse ell


cdef extern from "biellipsoid.cpp":
    cdef cppclass Biellipsoid:
        Vec3 position
        double r_forward
        double r_back
        double r_up
        double r_side
        Mat3 rot
        Ellipse f_limb
        Ellipse b_limb

        Biellipsoid()
        Biellipsoid(double r_forward, double r_backward, double r_up, double r_side)
        void set_rotation(double theta, double phi, double gamma, double ci)
        void set_position(Vec3 new_position)
        void set_radii(double r_forward, double r_back, double r_up, double r_side)
        void update_derived()
        Vec3 forward_vector()
        bint is_forward(Vec3 loc)
        bint is_visible(Vec3 loc)
        Bounds x_bounds()
        Bounds y_bounds()
        Bounds slice_ylimits(double x)
        bint line_intersects (double x, double y)
        Vec3 line_project (double x, double y)
        Vec3 world_to_aligned(Vec3 loc)
        Vec3 world_to_sphere(Vec3 loc)
        Vec3 aligned_to_world(Vec3 loc)
        Vec3 aligned_to_sphere(Vec3 loc)
        Vec3 sphere_to_world(Vec3 loc)
        Vec3 sphere_to_aligned(Vec3 loc)


cdef class BiellipsoidWrap:
    cdef Biellipsoid bell
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly


# ===== Light Source class and wrapper =====
cdef extern from "light_source.cpp":
    const int MAX_SOURCE_PARAMS
    const int MAX_LIMB_PARAMS
    cdef cppclass LightSource:
        int source_type
        double source_params[MAX_SOURCE_PARAMS]
        int limb_type
        double limb_params[MAX_LIMB_PARAMS]
        double limb_norm

        LightSource()
        LightSource(int source_type, double *source, int limb_type, double *limb_params)

        double get_brightness(double mu, double cos_lat, double lon)
        double get_brightness_sphere(double x, double y)
        bint uses_mu()
        bint uses_latlon()


cdef class LightSourceWrap:
    cdef LightSource source
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly


# ===== Transit integration function and wrapper =====
cdef extern from "transit_integral.cpp":
    void transit_integral(double *times, double *outputs, int n, const Orbit &orb, const LightSource &emitter, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side, bint rotate_with_orbit)

cpdef object transit(double r_forward, double r_back, double r_up, double r_side, double theta, double phi, double gamma, double[::1] t, double t0, double period, double semimajor, double inclination, str limbType, object limb, double eccen=?, double lon_periapse=?, bint rotate_with_orbit=?)


# ===== Phase Curve Class and Wrapper =====
cdef extern from "phase_curve_integral.cpp":
    const int MAX_PHASE_OBJECTS
    cdef cppclass PhaseIntegrator:
        Orbit orbits[MAX_PHASE_OBJECTS]
        Biellipsoid shapes[MAX_PHASE_OBJECTS]
        LightSource lights[MAX_PHASE_OBJECTS]
        bint rotate_with_orbit[MAX_PHASE_OBJECTS]
        Bounds xlim[MAX_PHASE_OBJECTS]
        Bounds ylim[MAX_PHASE_OBJECTS]

        PhaseIntegrator()
        PhaseIntegrator(PhaseIntegrator &p)

        int add_object(const Orbit &orb, const Biellipsoid &bell, const LightSource &source, bint rot_with_orbit)
        int get_n_objects() const
        void clear_objects()
        void set_time(double t)
        double integrate_single(int i)
        void phase_curve_integral(double *times, double *outputs, int n)

cdef class PhaseIntegratorWrap:
    cdef PhaseIntegrator pci
    # All wrapper functions are implemented in Python, so aren't declared here.
    # For use from Cython, use the C++ class directly
