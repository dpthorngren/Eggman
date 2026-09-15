General System Modelling
=========================

For more general calculations than standard transits, the :class:`eggman.PlanetSystem` class is used.  The basic framework is to initialize the class, add a star (or stars), then add any planets, moons, or rings. Objects are added with one of the following functions:

- :func:`eggman.PlanetSystem.add_object`
- :func:`eggman.PlanetSystem.add_star`
- :func:`eggman.PlanetSystem.add_planet`
- :func:`eggman.PlanetSystem.add_ring`

Although ``add_object`` is completely general, the other functions are easier to use and recommended unless you need to do something very unusual.  The ``add_planet`` function in particular has a large number of optional arguments you can use for setting the orbit, rotation, emission, etc. -- you can click through to the API docs for each function.

Transits
--------------------
We'll start with a simple example of a non-spherical object transiting a star with quadratic limb darkening.  The ``add_star`` function always defines the star as a unit sphere and varies only the limb darkening format and parameters used.  The ``add_planet`` function is much more elaborate -- the first four arguments are the radius of the planet in the prograde (morning), retrograde evening), polar, and radial (away from the star) directions. Then the ``period`` and ``semimajor`` axis arguments are required.  Additional arguments allow you to refine the orbit (eccentricity, inclination, etc), but we'll keep things simple here.

.. code:: python

   import eggman

    f = eggman.PlanetSystem()
    f.add_star("quadratic_limb", [0.2, 0.1])
    f.add_planet(.115, .11, .11, .12, period=10., semimajor=5., inclination=86.)

    times = np.linspace(0, 10, 1000)
    flux = f.integrate(times)

In this case, we've specified a planet with a slightly larger radius on the morning side than evening, and an even larger radius towards and away from the star -- this could occur for very close-orbiting hot Jupiters like WASP-12 b.  The total brightness of the system is obtained via the ``integrate`` function, which takes an array of times.

Moons and Rings
--------------------
Adding moons and rings around planets is a relatively simple modification from a simple transit.  Moons are achieved using the ``parent_index`` parameter, which causes them to orbit the n'th object added to the system rather than the origin.  The remaining arguments are the same, and note that length scales are in the same units (where the parent star has a default radius of 1), not relative to the parent object.

.. code:: python

    f = eggman.PlanetSystem()
    f.add_star("quadratic_limb", [0.2, 0.1])
    f.add_planet(.115, .11, .11, .12, 50., 300., inclination=89.)
    # The moon is given a parent object to orbit with parent_index
    f.add_planet(.05, .05, .05, .05, period=3.2, semimajor=.5, inclination=86., parent_index=1)
    # The ring also needs a parent_index, but has only outer and inner radii
    f.add_ring(0.15, 0.13, gamma=70, parent_index=1)

For good measure we have added a ring onto the planet as well.  This has fewer parameters, with the outer and inner radius being the most important.  The rotation system is the same as for planets, but because rings are azimuthally symmetrical, there is more than one way to achieve the same geometry.  The easiest setup is to use theta to set the angle where the ring is closest to the observer (counter-clockwise from bottom), gamma sets the inclination (0 is face-on, 90 is edge-on), and phi is left at 0.  Finally, as with the moon, you must specify the parent object index.

Phase Curves and Eclipses
-------------------------
Evaluating phase curves requires specifying the brightness of the planet across its surface.  This is done using the :class:`eggman.LightSource` class, and has the following options for the type:

- ``no_emission`` completely dark, the same as not specifying a planet's light source.
- ``lambertian`` uniform emission of light across the surface.
- ``quadratic_limb`` limb darkening using the quadratic functional form, usually only used for stars.
- ``nonlinear_limb`` limb darkening using the non-linear functional form, usually only used for stars.
- ``day_night`` a planet with one brightness across its day side and another on the night side.
- ``emission map`` interpolates the brightness on a grid of values provided by the user.  Extremely general, but much slower.

It is important to note for the latter two examples that rotating the planet rotates the source as well, so e.g. setting ``gamma=90`` in ``add_planet`` results in a planet whose "day side" points backwards along its orbit!  This isn't usually what you want, although it can make for an efficient approximation to a hot spot offset.

Here is an example using the day-night brightness map on a spherical planet.  The light source took two parameters -- the day and night side brightnesses -- and was specified as an argument to the ``add_planet`` function.

.. code:: python

    f = eggman.PlanetSystem()
    f.add_star("quadratic_limb", [0.2, 0.1])
    source = eggman.LightSource("day_night", [0.01, 0.001])
    f.add_planet(.1, .1, .1, .1, 10., 5., source=source)

For a more general emission function, we can instead use the ``emission_map`` source type.  This is extremely flexible but significantly slower out-of-transit compared to other functions, as the brightness at a given rotation isn't available in closed form.  This can be mitigated by integrating at fewer time steps when out of transit/eclipse.  The source takes two arguments: the latitudinal and longitudinal grid resolutions, which must be at least 2 each.  The definition requires an extra step as well, calling :func:`set_emission_map` with a function that takes a position on the unit sphere and returns a brightness at that point.

.. code:: python

    f = eggman.PlanetSystem()
    f.add_star("quadratic_limb", [0.2, 0.1])
    source = eggman.LightSource("emission_map", [500, 100])
    source.set_emission_map(lambda x, y, z: 1. - .3*z*z)
    f.add_planet(.1, .1, .1, .1, 10., 5., source=source)

For simplicity, we have defined a source function using a lambda function hat only depends on z, but you can define a full python function with a significantly more elaborate formula if you wish.
