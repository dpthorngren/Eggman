
cdef class OrbitWrap:
    '''A wrapper class for the C++ class Orbit, plus solve_kepler as a static function.'''

    @property
    def period(self) -> double:
        return self.orb.get_period()

    @property
    def semimajor(self) -> double:
        return self.orb.get_semimajor()

    @property
    def eccen(self) -> double:
        return self.orb.get_eccen()

    def __init__(self, double period, double t0, double semimajor, double eccen=0,
                 double inclination=90, double lon_periapse=90):
        self.orb = Orbit(period, t0, semimajor, eccen, inclination, lon_periapse)

    def get_position(self, double t):
        cdef Vec3 pos = self.orb.get_position(t)
        return np.array([pos.x, pos.y, pos.z])

    @staticmethod
    def solve_kepler(double mean_anomaly, double eccen):
        return solve_kepler(mean_anomaly, eccen)
