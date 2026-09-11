General System Modelling
=========================

Docs are still in progress, sorry.  In the meantime see the PlanetSystem API docs. Example code:


.. code:: python

   import eggman

    f = eggman.PlanetSystem()
    f.add_star("quadratic_limb", [0.2, 0.1])
    f.add_planet(.15, .15, .15, .15, period=10., semimajor=5., inclination=86.)

    times = np.linspace(0, 10, 1000)
    flux = f.integrate(times)
