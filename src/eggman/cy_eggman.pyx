# distutils: language = c++

import numpy as np

include "orbit_wrap.pyx"
include "biellipsoid_wrap.pyx"
include "light_source_wrap.pyx"
include "ellipse_wrap.pyx"
include "phase_wrap.pyx"


cpdef object transit(double r_forward, double r_back, double r_up, double r_side, double theta, double phi, double gamma, double[::1] t, double t0, double period, double semimajor, double inclination, str limb_type, object limb, double eccen=0, double lon_periapse=90, bint rotate_with_orbit=True):
    '''Wrapper for the c function transit3d_integral'''
    cdef Orbit orb = Orbit(period, t0, semimajor, eccen, inclination, lon_periapse)
    cdef LightSource emitter = LightSourceWrap("uniform", [1./np.pi], limb_type, limb).source
    assert (r_forward > 0) and (r_back > 0) and (r_side > 0), "Radii must be positive."
    if r_up < 0:
        assert phi == 0 and gamma == 0, "Nonzero phi and gamma are not supported for catwoman emulation mode (r_up < 0)."
        assert not rotate_with_orbit, "Rotation with orbit is not supported for catwoman emulation mode (r_up < 0)."
    results = np.full((len(t),), np.nan)
    cdef double[:] results_view = results
    transit_integral(&(t[0]), &(results_view[0]), len(t), orb, emitter, theta, phi, gamma, r_forward, r_back, r_up, r_side, rotate_with_orbit)
    return results
