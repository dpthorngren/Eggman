# Eggman
Code for calculating the geometry and transit depths of piecewise ellipsoidal objects.

## Installation
Eggman relies on three packages that the user must install themselves before compiling: Cython, NASA's CSPICE library, and the GNU Scientific Library (GSL).  See below for instructions on how to do this.  The latter two are C libraries so the Python setup script cannot retrieve them itself.

Once the prerequisites are installed, eggman can be compiled and installed from the terminal with `pip install .` from the eggman directory.

### Cython
Cython can be installed in pip with `pip install cython`.

### CSPICE
The [CSPICE library](https://naif.jpl.nasa.gov/naif/toolkit.html) can be installed by navigating the terminal to the eggman directory and running the included install script with `csh getCSPICE.csh`.  This requires csh be installed, which it is by default on Mac and can be acquired from the system package manager on Linux (sorry, I'd have used Bash if CSPICE didn't use csh already).  As of right now, other (Anaconda or system) installations of CSPICE do not link correctly with eggman.

### GSL
The GNU Scientific Library (GSL) can be installed through essentially any system package manager: `anaconda::gsl` for [Anaconda](https://anaconda.org/anaconda/gsl), `gsl` [for MacPorts](https://ports.macports.org/port/gsl/), `gsl` for [Homebrew](https://formulae.brew.sh/formula/gsl)and `libgsl-dev` for Linux using apt-get.  It can also be installed directly from the [GSL website](https://www.gnu.org/software/gsl/); just make sure you install it such that the compiler can locate it.

## Usage
For now the only production ready function is `asymmetricTransit`, which calculates the transit of a piecewise-ellipsoidal planet consisting of two half-ellipses attached at the location of the planet (diagram TODO).  They have the same polar radius *unless* the polar radius is set to a negative number, which eggman interprets to mean that the two half-circle model of [catwoman](https://github.com/KathrynJones1/catwoman) is to be used, such that the polar radius is `rMorning` for the morning side and `rEvening` for the evening side.
```python
import numpy as np
import eggman

eggman.asymmetricTransit(
    rMorning=.12,     # Relative to stellar radius
    rEvening=.1,
    rPole=.11,
    t=0.001,          # In any unit, so long as it's the same as t0 and period.
    t0=0,
    period=1.,
    semimajor=10.,    # Relative to the stellar radius
    inclination=89.,  # In degrees
    limbType=0,
    limb=np.array([.3, .2, 0., 0., 0.])
)
```

## Limb Darkening Settings
The limb darkening of the star is set by the `limbType` and `limb` parameters.  `limbType` is an integer specifying the darkening formula to use and `limb` is the parameters of that formula.  The latter must *always* be a length-5 numpy array, regardless of how many parameters are actually used.  The limb darkening options are (note the zero-indexing):

0. Quadratic: $I(\mu)/N = 1 - \gamma_0 (1-\mu) - \gamma_1 (1-\mu)^2$
1. Non-linear: $I(\mu)/N = 1 - \gamma_0 (1 - \sqrt{\mu}) - \gamma_1 (1-\mu) - \gamma_2 (1-\mu^{3/2}) - \gamma_3 (1-\mu^2)$

where $N$ is a normalization factor, $\gamma$ are model parameters (specified in `limb`), $\mu$ is the cosine of the angle between the viewer-star vector and the star-surface-point vector.  See [Mandel & Agol (2002)](https://ui.adsabs.harvard.edu/abs/2002ApJ...580L.171M/abstract) for more information.
