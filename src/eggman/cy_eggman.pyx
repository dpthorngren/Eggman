# import numpy as np


cpdef double _solve_kepler(double mean_anomaly, double eccen):
    '''Wrapper for c function solve_kepler, for testing purposes.'''
    return solve_kepler(mean_anomaly, eccen)


cpdef (double, double, double) orbit_to_position(double t, double period, double semimajor, double eccen, double inclination, double lon_periapse) noexcept:
    '''Calculates position in the view frame (x right, y up, z toward observer) given t
        and orbital elements.  t=0 is fixed as mid-transit.'''
    cdef Orbit orb = Orbit()
    init_orbit(&orb, period, semimajor, eccen, inclination, lon_periapse)
    cdef double x, y, z
    get_position(&orb, t, &x, &y, &z)
    return x, y, z
