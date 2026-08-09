# distutils: language = c++

import numpy as np

include "orbit.py"
include "shape.py"
include "light_source.py"
include "ellipse.py"
include "planet_system.py"


cpdef object transit(double r_forward, double r_back, double r_up, double r_side, double theta, double phi, double gamma, double[::1] t, double t0, double period, double semimajor, double inclination, str source_type, object source_params, double eccen=0, double lon_periapse=90, bint rotate_with_orbit=True, double atol=1e-6, double rtol=1e-3):
    '''Wrapper for the c function transit3d_integral'''
    cdef COrbit orb = COrbit(period, t0, semimajor, eccen, inclination, lon_periapse)
    cdef CLightSource emitter = LightSource(source_type, source_params).csource
    assert (r_forward > 0) and (r_back > 0) and (r_side > 0), "Radii must be positive."
    if r_up < 0:
        assert phi == 0 and gamma == 0, "Nonzero phi and gamma are not supported for catwoman emulation mode (r_up < 0)."
        assert not rotate_with_orbit, "Rotation with orbit is not supported for catwoman emulation mode (r_up < 0)."
    results = np.full((len(t),), np.nan)
    cdef double[:] results_view = results
    transit_integral(&(t[0]), &(results_view[0]), len(t), orb, emitter, theta, phi, gamma, r_forward, r_back, r_up, r_side, rotate_with_orbit, atol, rtol)
    return results
