
cdef class PlanetSystem:

    def __init__(self, atol=1e-6, rtol=1e-3):
        '''Initialize a blank PhaseIntegratorWrap object, setting the integration tolerances.

        Args:
            atol: the absolute tolerance for the integrators.
            rtol: the tolerance for the integrators relative to the result.'''
        self.cps.atol = atol
        self.cps.rtol = rtol

    def add_object(self, Orbit orbit, Shape shape, LightSource source, bint rotate_with_orbit=False, int parent_index=-1):
        '''Adds an object of the given specifications to the system.

        Args:
            orbit: An Orbit object describing the orbit of the object being added.
            shape: A Shape object describing the shape of the object being added.
            source: A LightSource object describing the emission of the object (may be non-emitting).
            rotate_with_orbit: whether or not the object rotates as it moves through its orbit (tidal locking).
            parent_index: The index of a previously added object that the new object should orbit (like a moon); if negative,
                the new object orbits the origin.
        '''
        self.cps.add_object(orbit.corbit, shape.cshape, source.csource, rotate_with_orbit, parent_index)

    def add_planet(self, double r_forward, double r_back, double r_up, double r_side, double period, double semimajor,
                   double t0=0, double eccen=0, double inclination=90, double lon_periapse=90., double theta=0,
                   double phi=0, double gamma=0, bint rotate_with_orbit=True, LightSource source=None,
                   int parent_index=-1):
        '''Adds a planet of the given specifications to the system.  This is a wrapper for add_object, whose arguments are passed through to Shape and Orbit.

        Args:
            r_forward: The forward radius of the object -- along the +x direction for an identity rotation matrix.
            r_back: The forward radius of the object -- along the -x direction for an identity rotation matrix.
            r_up: The upward radius of the object -- along the y axis for an identity rotation matrix.
            r_side: The upward radius of the object -- along the z axis for an identity rotation matrix.
            period: The orbital period.  Any unit may be used so long as it is used consistently.
            semimajor: The orbital semimajor axis.  Any unit may be used, but this sets the output units of
                the get_position function.  Usually users will wish to use units of stellar radius.
            t0: The epoch time, when the orbit passes x=0 in the positive direction; same units as period.
            eccen: The orbital eccentricity, it must be positive and less than 1, as this object cannot
                work with unbound orbits.
            inclination: The orbital inclination in degrees.  90 is an edge-on orbit, so transiting planets
                will typically be near this. i=0 is a clockwise orbit, starting from [0, a, 0] at t=0.
            lon_periapse: The longitude of the orbit's periapse.  A value of 90 places the periapse at
                [-a*(1-e), 0., 0.].
            theta: The counter-clockwise rotation of the biellipsoid around the z axis.
            phi: The counter-clockwise rotation of the biellipsoid around the y axis.
            gamma: The counter-clockwise rotation of the biellipsoid around the x axis.
            rotate_with_orbit: whether or not the object rotates as it moves through its orbit (tidal locking).
            source: A LightSource object describing the emission of the object (if None, the planet is not-emitting).
            parent_index: The index of a previously added object that the new object should orbit (like a moon); if negative,
                the new object orbits the origin.
        '''
        cdef COrbit orb = COrbit(period, t0, semimajor, eccen, inclination, lon_periapse)
        cdef CShape shape = CShape(r_forward, r_back, r_up, r_side)
        shape.position_from_orbit(0., orb, rotate_with_orbit)
        shape.set_rotation(theta, phi, gamma, orb.get_cos_inc())
        cdef CLightSource csource = CLightSource()
        if source is not None:
            csource = source.csource
        self.cps.add_object(orb, shape, csource, rotate_with_orbit, parent_index)

    def add_star(self, str limb_type, list[double] limb_params):
        '''Adds a star to the system.  This is a wrapper for add_object, and places a unit-sphere at the origin.

        Args:
            limb_type: A string indicating the type of limb darkening of the star.  Should be one of "lambertian" (flat),
                "quadratic_limb", or "nonlinear_limb".
            limb_params: A list of floats representing the parameters of the limb darkening, whose length depends on
                the limb type (0, 2, and 4 respectively from the above list).  If provided, an additional parameter at
                the start sets the overall brightness (1/pi results in a total luminosity of 1, and so is the default).
            '''
        cdef CLightSource source = LightSource(limb_type, limb_params).csource
        self.cps.add_object(COrbit(), CShape(), source, False, -1)

    def set_time(self, double t):
        '''Sets the time of the system, moving the objects into position based on their orbits.  This is automatically
        called in phase_curve_integral but can be useful for debugging and with plot_objects.

        Args:
            t: the time, in the same units as the period and epoch time t0 specified in the object orbits.'''
        self.cps.set_time(t)

    def integrate_single(self, int i):
        '''Get the brightness of a single target at the current time (set by set_time, defaulted to 0.)

        Args:
            i: Which object to get the brightness of -- these are indexed increasing from 0 in order of insertion.'''
        return self.cps.integrate_single(i)

    def phase_curve_integral(self, double[::1] times):
        '''Get the brightness of all objects in the system, added together, at the given times.  Once objects are added,
        this is the standard function to call for getting a light curve.

        Args:
            times: A numpy array of times to calculate for, in the same units as periods and epoch time t0 specified in the orbits.

        Returns:
            A Numpy array of brightnesses at the specified times.'''
        results = np.full((len(times),), np.nan)
        cdef double[:] results_view = results
        self.cps.phase_curve_integral(&(times[0]), &(results_view[0]), len(times))
        return results

    def get_n_objects(self):
        '''Returns the number of objects that have been successfully inserted so far.'''
        return self.cps.get_n_objects()

    def clear_objects(self):
        '''Deletes all the objects so far added, starting over from scratch.  This is slightly faster than creating a new object.'''
        return self.cps.clear_objects()

    def __getitem__(self, int i):
        '''Get an object previously added to the system.

        Args:
            i: The index of the object, which were set in order of insertion starting from 0.

        Returns:
            A tuple of an OrbitWrap, Shape, LightSource, and bool object that describing the i'th object -- the final
                bool is whether the object rotates with its orbit.
            '''
        if i < 0:
            i = self.cps.get_n_objects() - i
        if i >= self.cps.get_n_objects() or i < 0:
            raise IndexError
        cdef Orbit orbit = Orbit(5., 0., 0.)
        orbit.corbit = self.cps.orbits[i]
        cdef Shape shape = Shape(0., 0., 0., 0.)
        shape.cshape = self.cps.shapes[i]
        cdef LightSource source = LightSource("none", [])
        source.csource = self.cps.lights[i]
        return (orbit, shape, source, self.cps.rotate_with_orbit[i])

    def plot_objects(self, res=200, axis=None, arg_list=list(), **args):
        '''Create a simple plot of the system using matplotlib, at the current time previously set by set_time (defaulting to 0).

        If the object is emissive, it is plotted with pyplot.pcolormesh, else it is plotted with pyplot.fill_between.

        Args:
            res: The resolution to plot the objects at, as an integer.
            axis: The matplotlib axis object to plot on.  If left as the default None, new axes will be created.
            arg_list: A list of dictionaries containing keyword arguments to pass to pyplot.pcolormesh
                for each object in the system in order.  These settings override global settings from **args.
            **args: Additional keyword arguments to be passed to pyplot.pcolormesh for all objects in the system,
                unless overridden by arg_list.'''
        for i in range(self.get_n_objects()):
            _, shape, source, _ = self[i]
            arguments = args.copy()
            if arg_list:
                arguments.update(arg_list[i])
            arguments['zorder'] = shape.position[2]
            source.plot_brightness(shape, res, axis, **arguments)
