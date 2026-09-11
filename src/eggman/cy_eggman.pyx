# distutils: language = c++

import numpy as np

include "orbit.py"
include "shape.py"
include "light_source.py"
include "ellipse.py"
include "planet_system.py"


cpdef object asymmetricTransit(double rMorning, double rEvening, double rPole, double[:] t, double t0, double period, double semimajor, double inclination, str limbType, object limb, double eccen=0, double lonPeriapse=90., double theta=0., double atol=1e-6, double rtol=1e-3, int max_steps = 200):
    '''Calculates the transit of a piecewise-elliptical planet.  Assumes the same projected shape regardless of its position
    in the orbit.  Uses the same model as catwoman (two spheres split down the middle) if rPole is negative.

    Parameters:
        rMorning        The radius of the planet at the morning (right) side of the planet at the equator relative to
                            the stellar radius.
        rEvening        The radius of the planet at the evening (left) side of the planet at the equator relative top
                            the stellar radius.
        rPole           The radius of the planet at the poles (top and bottom) relative to the stellar radius.  If -1,
                            rPole is set to rMorning one morning side and rEvening on the evening side.
        t               A 1-d Numpy array of observation times to be simulated (must be the same unit as t0 and period).
        t0              The time of a mid-transit (needn't be the observed one).
        period          The orbital period of the planet.
        semimajor       The semimajor axis of the planet relative to the stellar radius.
        inclination     The inclination of the orbit in degrees (near 90 for transiting planets).
        limbType        The type of limb darkening to use: 'quadratic' or 'nonlinear'.
        limb            An iterable of limb darkening parameters, length 2 for quadratic and 4 for nonlinear.
        eccen           The orbital eccentricity of the planet.  Must be 0 <= e < 1, defaults to zero.
        lonPeriapse     The longitude of the periapse of the planet's orbit, in degrees.  The default of 90 means the pariapse occurs at mid-transit.
        theta           The angle to rotate the planet by on-sky in the counter-clockwise direction.
        atol            The absolute tolerance for the integration results.
        rtol            The relative tolerance for the integration results.
        max_steps       The maximum number of steps to use in integration; if this is exceeded, an error is printed and nan is returned.

    Returns:
        The relative flux from the star at the time given, so 1 if the planet is out of transit.
    '''
    cdef COrbit orb = COrbit(period, t0, semimajor, eccen, inclination, lonPeriapse)
    assert (rMorning > 0) and (rEvening > 0), "Radii must be positive."
    cdef double limbParams[4]
    if limbType == "quadratic":
        limbParams[0] = limb[0]
        limbParams[1] = limb[1]
        limbParams[2] = -1
        limbParams[3] = -1
    elif limbType == "nonlinear":
        limbParams[0] = limb[0]
        limbParams[1] = limb[1]
        limbParams[2] = limb[2]
        limbParams[3] = limb[3]
    else:
        raise ValueError("Limb type must be 'quadratic' or 'nonlinear'")
    results = np.full((len(t),), np.nan)
    cdef double[:] results_view = results
    transit_integral(&(t[0]), &(results_view[0]), len(t), orb, rMorning, rEvening, rPole, limbParams[0], limbParams[1], limbParams[2], limbParams[3], theta, atol, rtol, max_steps)
    return results
