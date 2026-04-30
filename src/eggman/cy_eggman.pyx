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


cpdef double get_brightness(int source_type, list[double] limb_params, double x, double y):
    assert source_type in [0, 1]
    assert len(limb_params) == [2, 4][source_type]
    if source_type == 0:
        limb_params += [0., 0.]
    cdef double[4] lp = [limb_params[0], limb_params[1], limb_params[2], limb_params[3]]
    cdef LightSource source = create_light_source(source_type, lp)
    return get_source_brightness(&source, x, y)
