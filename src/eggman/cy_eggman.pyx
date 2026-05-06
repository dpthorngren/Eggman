# distutils: language = c++

import numpy as np


def _struct_to_arr_(s: dict):
    return np.array([s[k] for k in sorted(s.keys())])


cpdef _solve_kepler(double mean_anomaly, double eccen):
    '''Wrapper for the c++ function solve_kepler, for testing purposes.'''
    return solve_kepler(mean_anomaly, eccen)


cpdef (double, double, double) orbit_to_position(double t, double period, double semimajor, double eccen, double inclination, double lon_periapse) noexcept:
    '''Calculates position in the view frame (x right, y up, z toward observer) given t
        and orbital elements.  t=0 is fixed as mid-transit.'''
    cdef Orbit orb = Orbit(period, 0., semimajor, eccen, inclination, lon_periapse)
    cdef Vec3 position = orb.get_position(t)
    return position.x, position.y, position.z


cpdef dict biellipsoid_dump(double x, double y, double z, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side, double ci=-2):
    '''Creates a biellipsoid and returns its contents for visualization and testing purposes.'''
    cdef Biellipsoid bell = Biellipsoid(r_forward, r_back, r_up, r_side)
    cdef Vec3 position = Vec3(x, y, z)
    bell.set_position(position)
    bell.set_rotation(theta, phi, gamma, ci)
    return dict(
        position=_struct_to_arr_(bell.position),
        radii=np.array([bell.r_forward, bell.r_back, bell.r_up, bell.r_side]),
        rot=_struct_to_arr_(bell.rot),
        x_bounds=bell.x_bounds(),
        y_bounds=bell.y_bounds(),
        f1=_struct_to_arr_(bell.f_limb.e1),
        f2=_struct_to_arr_(bell.f_limb.e2),
        b1=_struct_to_arr_(bell.b_limb.e1),
        b2=_struct_to_arr_(bell.b_limb.e2),
    )


cpdef (double, double) _biellipsoid_slice_ylimits(double x, double theta, double phi, double gamma, double r_forward, double r_back, double r_up, double r_side, double ci=-2):
    '''Creates a biellipsoid and returns the results of slice_ylimits at the given x.'''
    cdef Biellipsoid bell = Biellipsoid(r_forward, r_back, r_up, r_side)
    bell.set_rotation(theta, phi, gamma, ci)
    cdef Bounds b = bell.slice_ylimits(x)
    return b.min, b.max


cdef LightSource light_source(str limbType, object limb):
    '''Wrapper for the LightSource constructor, adjusted to take Python arguments.'''
    cdef int limb_code = -1
    limbType = limbType.lower().strip()
    if limbType == "quadratic":
        assert len(limb) == 2, "Error: Quadratic limb darkening takes exactly two parameters."
        limb_code = 0
        limb = [limb[0], limb[1], 0., 0.]
    elif limbType == "nonlinear":
        assert len(limb) == 4, "Error: Nonlinear limb darkening takes exactly two parameters."
        limb_code = 1
    else:
        raise ValueError("limbType not recognized, must be one of quadratic, nonlinear.")
    cdef LightSource emitter = LightSource(limb_code, limb[0], limb[1], limb[2], limb[3])
    return emitter


cpdef double _get_brightness(str limb_type, list[double] limb_params, double x, double y):
    cdef LightSource source = light_source(limb_type, limb_params)
    return source.get_brightness(x, y)


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
