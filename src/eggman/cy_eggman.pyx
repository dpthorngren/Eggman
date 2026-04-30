import numpy as np


cpdef dict biellipsoid_dump(double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side):
    '''Wrapper for the c function init_biellipsoid for visualization and testing purposes.'''
    cdef Vec3 position = Vec3(0., 0., 0.)
    cdef Biellipsoid bell = create_biellipsoid(position, theta, phi, gamma, r_forward, r_back, r_up, r_side)
    result = dict(bell)
    for key in ['position', 'rot', 'f_limb', 'b_limb', 'break_offset', 'f1', 'f2', 'b1', 'b2']:
        result[key] = np.array([result[key][k] for k in sorted(result[key].keys())])
    result['radii'] = np.array(result['radii'])
    return result


cpdef (double, double) _biellipsoid_ylimits(double x, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side):
    '''Wrapper for c functions init_biellipsoid and get_ybounds for testing purposes.'''
    cdef Vec3 position = Vec3(0., 0., 0.)
    cdef Biellipsoid bell = create_biellipsoid(position, theta, phi, gamma, r_forward, r_back, r_up, r_side)
    cdef Bounds b = get_ylimits(&bell, x)
    return b.min, b.max


cpdef double _solve_kepler(double mean_anomaly, double eccen):
    '''Wrapper for c function solve_kepler, for testing purposes.'''
    return solve_kepler(mean_anomaly, eccen)


cpdef (double, double, double) orbit_to_position(double t, double period, double semimajor, double eccen, double inclination, double lon_periapse) noexcept:
    '''Calculates position in the view frame (x right, y up, z toward observer) given t
        and orbital elements.  t=0 is fixed as mid-transit.'''
    cdef Orbit orb = create_orbit(period, semimajor, eccen, inclination, lon_periapse)
    cdef Vec3 position = get_position(&orb, t)
    return position.x, position.y, position.z
