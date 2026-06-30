
cdef class PhaseIntegratorWrap:

    def add_object(self, OrbitWrap orbit, BiellipsoidWrap bell, LightSourceWrap source, bint rotate_with_orbit=False):
        cdef Orbit c_orb = orbit.orb
        cdef Biellipsoid c_bell = bell.bell
        cdef LightSource c_source = source.source
        self.pci.add_object(c_orb, c_bell, c_source, rotate_with_orbit)

    def add_planet(self, double r_forward, double r_back, double r_up, double r_side, double theta, double phi, double gamma, double t0, double period, double semimajor, double inclination, double eccen=0, double lon_periapse=90., bint rotate_with_orbit=True, LightSourceWrap source=None):
        cdef Orbit orb = Orbit(period, t0, semimajor, eccen, inclination, lon_periapse)
        cdef Biellipsoid bell = Biellipsoid(r_forward, r_back, r_up, r_side)
        bell.set_rotation(theta, phi, gamma, orb.get_cos_inc())
        cdef LightSource csource = LightSource()
        if source is not None:
            csource = source.source
        self.pci.add_object(orb, bell, csource, rotate_with_orbit)

    def add_star(self, str limb_type, list[double] limb_params, double brightness=1./np.pi):
        cdef Orbit orb = Orbit()
        cdef Biellipsoid bell = Biellipsoid()
        cdef LightSource source = LightSourceWrap("uniform", [brightness], limb_type, limb_params).source
        self.pci.add_object(orb, bell, source, False)

    def set_time(self, double t):
        self.pci.set_time(t)

    def integrate_single(self, int i):
        return self.pci.integrate_single(i)

    def phase_curve_integral(self, double[::1] times):
        results = np.full((len(times),), np.nan)
        cdef double[:] results_view = results
        self.pci.phase_curve_integral(&(times[0]), &(results_view[0]), len(times))
        return results

    def get_n_objects(self):
        return self.pci.get_n_objects()

    def clear_objects(self):
        return self.pci.clear_objects()

    def __getitem__(self, int i):
        if i < 0:
            i = self.pci.get_n_objects() - i
        if i >= self.pci.get_n_objects() or i < 0:
            raise IndexError
        cdef OrbitWrap orbit = OrbitWrap(5., 0., 0.)
        orbit.orb = self.pci.orbits[i]
        cdef BiellipsoidWrap bell = BiellipsoidWrap(0., 0., 0., 0., 0., 0., 0., 0., 0., 0.)
        bell.bell = self.pci.shapes[i]
        cdef LightSourceWrap source = LightSourceWrap("none", [], "lambertian", [])
        source.source = self.pci.lights[i]
        return (orbit, bell, source, self.pci.rotate_with_orbit[i])
