import typing

import cython
import numpy as np
# Cython cimports (LSP cannot process)
from cython.cimports.eggman import cy_eggman as cye  # type: ignore

if not cython.compiled and typing.TYPE_CHECKING:
    # Allow editor LSP to find existing compiled code...
    from eggman import LightSource, Orbit, Shape
else:  # ...but provide the real cimports during compilation
    from cython.cimports.eggman.cy_eggman import LightSource, Orbit, Shape  # type: ignore


@cython.cclass
class PlanetSystem:
    cps: cye.CPlanetSystem

    def __init__(self, atol=1e-6, rtol=1e-3):
        '''Initialize a PlanetSystem with the integration tolerances.

        Args:
            atol: the absolute tolerance for the integrators.
            rtol: the tolerance for the integrators relative to the result.'''
        self.cps.atol = atol
        self.cps.rtol = rtol

    def add_object(
            self,
            orbit: Orbit,
            shape: Shape,
            source: LightSource,
            rotate_with_orbit: cython.bint = False,
            parent_index: int = -1):
        '''Adds an object of the given specifications to the system.

        Args:
            orbit: An Orbit object describing the orbit of the object being added.
            shape: A Shape object describing the shape of the object being added.
            source: A LightSource object describing the emission of the object.
            rotate_with_orbit: whether or not the object rotates as it moves through its orbit
                (tidally locked to the parent or origin).
            parent_index: The index of a previously added object that the new object should orbit
                (like a moon); if negative, the new object orbits the origin.
        '''
        self.cps.add_object(
            orbit.corbit, shape.cshape, source.csource, rotate_with_orbit, parent_index)

    def add_planet(
            self,
            r_forward: float,
            r_back: float,
            r_up: float,
            r_side: float,
            period: float,
            semimajor: float,
            t0: float = 0,
            eccen: float = 0,
            inclination: float = 90,
            lon_periapse: float = 90.,
            theta: float = 0,
            phi: float = 0,
            gamma: float = 0,
            source: LightSource = None,
            rotate_with_orbit: cython.bint = True,
            parent_index: int = -1):
        '''Adds a planet of the given specifications to the system.  This is a wrapper for
            add_object, whose arguments are passed through to Shape and Orbit.

        Args:
            r_forward: The forward radius of the object -- along the +x direction for an identity
                rotation matrix.
            r_back: The forward radius of the object -- along the -x direction for an identity
                rotation matrix.
            r_up: The upward radius of the object -- along the y axis for an identity rotation
                matrix.
            r_side: The upward radius of the object -- along the z axis for an identity rotation
                matrix.
            period: The orbital period.  Any unit may be used so long as it is used consistently.
            semimajor: The orbital semimajor axis.  Any unit may be used, but this sets the output
                units of the get_position function.  Usually users will wish to use units of
                stellar radius.
            t0: The epoch time, when the orbit passes x=0 in the positive direction; must have the
                same units as period.
            eccen: The orbital eccentricity, it must be positive and less than 1, as this object
                cannot work with unbound orbits.
            inclination: The orbital inclination in degrees.  90 is an edge-on orbit, so transiting
                planets will typically be near this. i=0 is a clockwise orbit, starting from
                [0, a, 0] at t=0.
            lon_periapse: The longitude of the orbit's periapse.  A value of 90 places the periapse
                at [-a*(1-e), 0., 0.].
            theta: The counter-clockwise rotation of the biellipsoid around the z axis.
            phi: The counter-clockwise rotation of the biellipsoid around the y axis.
            gamma: The counter-clockwise rotation of the biellipsoid around the x axis.
            source: A LightSource object describing the emission of the object.
            rotate_with_orbit: whether or not the object rotates as it moves through its orbit
                (tidally locked to the parent or origin).
            parent_index: The index of a previously added object that the new object should orbit
                (like a moon); if negative, the new object orbits the origin.'''
        cython.declare(orb=cye.COrbit, shape=cye.CShape, csource=cye.CLightSource, origin=cye.Vec3)
        orb = cye.COrbit(period, t0, semimajor, eccen, inclination, lon_periapse)
        shape = cye.CShape(r_forward, r_back, r_up, r_side)
        origin = cye.Vec3(0., 0., 0.)
        if parent_index >= 0:
            origin = self.cps.shapes[parent_index].position
        shape.position_from_orbit(0., orb, rotate_with_orbit, origin)
        shape.set_rotation(theta, phi, gamma, orb.get_cos_inc())
        csource = cye.CLightSource()
        if source is not None:
            csource = source.csource
        self.cps.add_object(orb, shape, csource, rotate_with_orbit, parent_index)

    def add_star(self, limb_type: str, limb_params: list[float]):
        '''Adds a star to the system.  This is a wrapper for add_object, and places a unit-sphere
            at the origin.

        Args:
            limb_type: A string indicating the type of limb darkening of the star.  Should be one
                of "lambertian" (flat), "quadratic_limb", or "nonlinear_limb".
            limb_params: A list of floats representing the parameters of the limb darkening, whose
                length depends on the limb type (0, 2, and 4 respectively from the above list).  If
                provided, an additional parameter at the start sets the overall brightness
                (1/pi results in a total luminosity of 1, and so is the default).'''
        source: cye.CLightSource = LightSource(limb_type, limb_params).csource
        self.cps.add_object(cye.COrbit(), cye.CShape(), source, False, -1)

    def add_ring(
            self,
            r_outer: float,
            r_inner: float,
            parent_index: int,
            theta: float = 0,
            phi: float = 0,
            gamma: float = 0):
        '''Adds a ring around the given object. With rotations set to zero, the ring is face-on.
            This is a wrapper for add_object, whose arguments are passed through to Shape and Orbit.

        Args:
            r_outer: The outer radius of the ring.
            r_inner: The inner radius of the ring, that is, within this radius the ring will not
                block any light.  This should be large enough that it is not touching the
                parent object if the latter is emissive.
            theta: The counter-clockwise rotation of the ring around the z axis.
            phi: The counter-clockwise rotation of the ring around the y axis.
            gamma: The counter-clockwise rotation of the ring around the x axis.
            parent_index: The index of the object to place a ring around.'''
        cython.declare(orb=cye.COrbit, shape=cye.CShape, csource=cye.CLightSource, origin=cye.Vec3)
        orb = cye.COrbit()
        shape = cye.CShape(r_outer, r_inner, 0, 0.)
        origin = cye.Vec3(0., 0., 0.)
        if parent_index >= 0:
            origin = self.cps.shapes[parent_index].position
        shape.position_from_orbit(0., orb, False, origin)
        shape.set_rotation(theta, phi, gamma, orb.get_cos_inc())
        csource = cye.CLightSource()
        self.cps.add_object(orb, shape, csource, False, parent_index)

    def set_time(self, t: float) -> None:
        '''Sets the time of the system, moving the objects into position based on their orbits.
            This is automatically called in phase_curve_integral but can be useful for debugging
            and with plot_objects.

        Args:
            t: the time, in the same units as the period and epoch time t0 specified in the Orbit
                objects.'''
        self.cps.set_time(t)

    def integrate_single(self, i: int) -> float:
        '''Get the brightness of a single target at the current time (set by set_time,
            defaulted to 0.)

        Args:
            i: Which object to get the brightness of -- these are indexed increasing from 0 in
                order of insertion.'''
        return self.cps.integrate_single(i)

    def phase_curve_integral(self, times: cye.Array1d_f64):
        '''Get the brightness of all objects in the system, added together, at the given times.
            This is the standard function to call for getting a light curve.

        Args:
            times: A numpy array of times to calculate for, in the same units as periods and epoch
                time t0 specified in the orbits.

        Returns:
            A Numpy array of brightnesses at the specified times.'''
        results = np.full((len(times),), np.nan)
        results_view: cye.Array1d_f64 = results
        self.cps.phase_curve_integral(
            cython.address(times[0]), cython.address(results_view[0]), len(times))
        return results

    def get_n_objects(self) -> int:
        '''Returns the number of objects that have been successfully inserted so far.'''
        return self.cps.get_n_objects()

    def clear_objects(self) -> None:
        '''Deletes all the objects so far added, starting over from scratch.  This is slightly
            faster than creating a new object.'''
        return self.cps.clear_objects()

    def __getitem__(self, i: int) -> tuple[Orbit, Shape, LightSource, bool]:
        '''Get an object previously added to the system.

        Args:
            i: The index of the object, which were set in order of insertion starting from 0.

        Returns:
            A tuple of an OrbitWrap, Shape, LightSource, and bool object that describing the i'th
                object -- the final bool is whether the object rotates with its orbit.
            '''
        if i < 0:
            i = self.cps.get_n_objects() - i
        if i >= self.cps.get_n_objects() or i < 0:
            raise IndexError
        orbit: Orbit = Orbit(5., 0., 0.)
        orbit.corbit = self.cps.orbits[i]
        shape: Shape = Shape(0., 0., 0., 0.)
        shape.cshape = self.cps.shapes[i]
        source: LightSource = LightSource("no_emission", [])
        source.csource = self.cps.lights[i]
        return (orbit, shape, source, self.cps.rotate_with_orbit[i])

    def plot_objects(self, res=200, axis=None, arg_list=list(), **args):
        '''Create a simple plot of the system using matplotlib, at the current time (set by
            set_time or 0 by default), using matplotlib.pcolormesh.

        Args:
            res: The resolution to plot the objects at, as an integer.
            axis: The matplotlib axis object to plot on.  If left as the default None, new axes
                will be created.
            arg_list: A list of dictionaries containing keyword arguments to pass to
                pyplot.pcolormesh. for each object in the system in order.  These settings override
                global settings from **args.
            **args: Additional keyword arguments to be passed to pyplot.pcolormesh for all objects
                in the system, unless overridden by arg_list.'''
        for i in range(self.get_n_objects()):
            _, shape, source, _ = self[i]
            arguments = args.copy()
            if arg_list:
                arguments.update(arg_list[i])
            if shape.r_up == 0:
                arguments['zorder'] = shape.position[2] + .5
                source.plot_brightness(shape, res, axis, ring_front=True, **arguments)
                arguments['zorder'] = shape.position[2] - .5
                source.plot_brightness(shape, res, axis, ring_front=False, **arguments)
            else:
                arguments['zorder'] = shape.position[2]
                source.plot_brightness(shape, res, axis, **arguments)
