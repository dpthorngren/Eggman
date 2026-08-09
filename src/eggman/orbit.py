# flake8: noqa: E302
import cython
import numpy as np
from cython.cimports.eggman import cy_eggman as cye  # type: ignore


@cython.cclass
class Orbit:
    '''Describes an orbit so that the position at a given time may be calculated.

    This class wraps the C++ class Orbit, plus solve_kepler as a static function.'''
    corbit: cye.COrbit

    @property
    def period(self) -> float:
        '''The orbital period.'''
        return self.corbit.get_period()

    @property
    def semimajor(self) -> float:
        '''The orbital semimajor axis.'''
        return self.corbit.get_semimajor()

    @property
    def eccen(self) -> float:
        '''The orbital eccentricity.'''
        return self.corbit.get_eccen()

    @property
    def cos_inc(self) -> float:
        '''The cosine of the inclination, which is useful for orienting objects along their orbits.'''
        return self.corbit.get_cos_inc()

    def __init__(
            self,
            period: float,
            t0: float,
            semimajor: float,
            eccen: float = 0,
            inclination: float = 90,
            lon_periapse: float = 90):
        '''Initialize an orbit object with the given orbital elements.

        Args:
            period: The orbital period.  Any unit may be used so long as it is used consistently.
            t0: The epoch time, when the orbit passes x=0 in the positive direction; same units as period.
            semimajor: The orbital semimajor axis.  Any unit may be used, but this sets the output units of
                the get_position function.  Usually users will wish to use units of stellar radius.
            eccen: The orbital eccentricity, it must be positive and less than 1, as this object cannot
                work with unbound orbits.
            inclination: The orbital inclination in degrees.  90 is an edge-on orbit, so transiting planets
                will typically be near this. i=0 is a clockwise orbit, starting from [0, a, 0] at t=0.
            lon_periapse: The longitude of the orbit's periapse.  A value of 90 places the periapse at
                [-a*(1-e), 0., 0.].
        '''
        self.corbit = cye.COrbit(period, t0, semimajor, eccen, inclination, lon_periapse)

    def get_position(self, t: float):
        '''Calculates the position of the orbiting object at a given time.

        Args:
            t: The time to get the orbital position for.  The units must be the same as those used for
                `period` and `t0` on initialization.

        Returns:
            The position of the orbiting object at time t, as a one-dimensional numpy array of size 3.
        '''
        pos = self.corbit.get_position(t)
        return np.array([pos.x, pos.y, pos.z])

    @staticmethod
    def solve_kepler(mean_anomaly: float, eccen: float):
        '''Solves the kepler equation for the true anomaly given the mean anomaly and eccentricity.'''
        return cye.solve_kepler(mean_anomaly, eccen)
