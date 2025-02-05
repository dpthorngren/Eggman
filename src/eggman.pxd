cimport cython
from libc.math cimport sin, cos, sqrt, atan2, acos, log, exp
from libc.math cimport M_PI as pi, NAN as nan


cdef extern from "SpiceUsr.h":
    ctypedef double SpiceDouble
    ctypedef bint SpiceBoolean
    ctypedef struct SpicePlane:
        pass
    ctypedef struct SpiceEllipse:
        pass
    

cdef extern from "SpiceZfc.h":
    # Initializes a SpicePlane given a normal vector and an offset from the origin
    cdef void nvc2pl_c(SpiceDouble* normal, SpiceDouble offset, SpicePlane* outPlane)
    # Gets the ellipse defined by the intersection of an ellipsoid and a plane
    cdef void inedpl_c(SpiceDouble a, SpiceDouble b, SpiceDouble c, SpicePlane* plane, SpiceEllipse* limbEllipse, SpiceBoolean* found)
    # Projects an ellipse onto the given plane
    cdef void pjelpl_c(SpiceEllipse* inEllipse, SpicePlane* projPlane, SpiceEllipse* outEllipse)
    # Gets the center and constructing vectors of an ellipse
    cdef void el2cgv_c(SpiceEllipse* inEllipse, SpiceDouble* center, SpiceDouble* vec1, SpiceDouble* vec2)
    # Gets the major and minor axis vectors from the constructing vectors.
    cdef void saelgv_c(SpiceDouble* vec1, SpiceDouble* vec2, SpiceDouble* majorAxis, SpiceDouble* minorAxis)


cdef extern from "gsl/gsl_math.h":
    ctypedef struct gsl_function:
        double (*function)(double x, void *params)
        void * params


cdef extern from "gsl/gsl_roots.h":
    ctypedef struct gsl_root_fsolver:
        pass
    ctypedef struct gsl_root_fsolver_type:
        pass

    const gsl_root_fsolver_type* gsl_root_fsolver_brent

    gsl_root_fsolver* gsl_root_fsolver_alloc(const gsl_root_fsolver_type* T)
    void gsl_root_fsolver_free(gsl_root_fsolver *s)
    int gsl_root_fsolver_set(gsl_root_fsolver *s, gsl_function *f, double x_lower, double x_upper)
    int gsl_root_fsolver_iterate(gsl_root_fsolver *s)
    double gsl_root_fsolver_root(const gsl_root_fsolver *s)
    double gsl_root_fsolver_x_lower(const gsl_root_fsolver *s)
    double gsl_root_fsolver_x_upper(const gsl_root_fsolver *s)

cdef extern from "gsl/gsl_integration.h":
    ctypedef struct gsl_integration_workspace:
        size_t size

    gsl_integration_workspace *gsl_integration_workspace_alloc(size_t n)
    void gsl_integration_workspace_free(gsl_integration_workspace *w)
    int gsl_integration_qag(const gsl_function *f, double a, double b, double epsabs, double epsrel, size_t limit, int key, gsl_integration_workspace *workspace, double *result, double *abserr)

cdef extern from "gsl/gsl_errno.h":
    ctypedef struct gsl_error_handler_t:
        pass
    gsl_error_handler_t* gsl_set_error_handler(gsl_error_handler_t *new_handler)
    gsl_error_handler_t* gsl_set_error_handler_off()


cdef struct IntegralParams:
    double a
    double b
    double xe
    double ye
    double tNear
    double tFar
    double rsq
    int limbType
    double limb[5]
    gsl_root_fsolver* solver
    gsl_function* func

cdef struct BruteIntegralParams:
    double a
    double b
    double xe
    double ye
    double x
    int limbType
    double limb[5]
    gsl_integration_workspace* work
    gsl_function* integrand


cpdef (double, double, double, double) orbitGeometry(double a, double b, double c, double semimajor, double theta, double phi)

cpdef double transitDepth(double a, double b, double c, double semimajor, double theta, double phi, double[:] limb)

cpdef double transitIntegral(double a, double b, double xe, double ye, double[:] limb, int preferBrute=?)
