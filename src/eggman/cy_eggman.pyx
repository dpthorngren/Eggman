# distutils: language = c++

import numpy as np

include "orbit_wrap.pyx"
include "biellipsoid_wrap.pyx"


cpdef object transit(double r_forward, double r_back, double r_pole, double r_side, double theta, double phi, double gamma, double[::1] t, double t0, double period, double semimajor, double inclination, str limb_type, object limb, double eccen=0, double lon_periapse=90):
    '''Wrapper for the c function transit3d_integral'''
    cdef Orbit orb = Orbit(period, t0, semimajor, eccen, inclination, lon_periapse)
    cdef LightSource emitter = light_source(limb_type, limb)
    assert (r_forward > 0) and (r_back > 0) and (r_side > 0), "Radii must be positive."
    assert r_pole > 0, "This function doesn't support catwoman-style transit emulation, use transit2d instead."

    results = np.full((len(t),), np.nan)
    cdef double[:] results_view = results
    transit_integral(&(t[0]), &(results_view[0]), len(t), &orb, &emitter, theta, phi, gamma, r_forward, r_back, r_pole, r_side)
    return results
