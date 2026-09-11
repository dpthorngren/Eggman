Basic Transit Modelling
=========================

.. image:: img/eggmanComparison.png
   :width: 30%
   :align: right

Transits are the most common use-case for Eggman, so it provides a specialized function that is more efficient than the general integrator class. As shown in the figure, the planet is represented as two ellipses joined along a common axis, meaning that the pole, morning, and evening radii are set independently. The ``theta`` parameter rotates the object on-sky counter-clockwise.

A second mode can be used by setting the polar radius to a negative number.  In this case, the plant's radius is discontinuous at the pole, as is done by `catwoman <https://github.com/KathrynJones1/catwoman>`__.

The orbit is specified by the ``t0`` (time of mid-transit), ``period``, ``semimajor``, ``inclination``, ``eccen``, and ``lonPeriapse`` arguments.  Units of distance are relative to stellar radius and angles are in degrees. Time units may be anything so long as they are used consistently.

The star's limb darkening is specified by the ``limbType`` and ``limb`` arguments; the former is either ``quadratic`` or ``nonlinear`` and the latter is a list of parameters (2 and 4 respectively).  These refer to the following formulas from `Mandel & Agol (2002) <https://ui.adsabs.harvard.edu/abs/2002ApJ...580L.171M/abstract>`__ respectively:

:math:`I_\mathrm{q}(\mu)/N(\gamma_0, \gamma_1) = 1 - \gamma_0 (1-\mu) - \gamma_1 (1-\mu)^2`
:math:`I_\mathrm{nl}(\mu)/N(\gamma_0, \gamma_1, \gamma_2, \gamma_3) = 1 - \gamma_0 (1 - \sqrt{\mu}) - \gamma_1 (1-\mu) - \gamma_2 (1-\mu^{3/2}) - \gamma_3 (1-\mu^2)`

The functions `N(...)` normalize the star to a total brightness of 1.

.. autofunction:: eggman.asymmetricTransit
