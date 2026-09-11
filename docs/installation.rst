Installation
======================
Eggman relies on the GNU Scientific Library (GSL), which the user must install before Eggman. This is a C library so the Python setup script cannot retrieve them itself. It can be installed through essentially any system package manager: ``anaconda::gsl`` for `Anaconda <https://anaconda.org/anaconda/gsl>`__, ``gsl`` for `MacPorts <https://ports.macports.org/port/gsl/>`__, ``gsl`` for `Homebrew <https://formulae.brew.sh/formula/gsl>`__, and ``libgsl-dev`` for Linux using apt-get.  It can also be installed directly from the `GSL website <https://www.gnu.org/software/gsl/>`__; just make sure you install it such that the compiler can locate it.

You may install Eggman directly from Github using pip:

.. code-block:: console

    pip install git+https://github.com/dpthorngren/Eggman#egg=eggman


Alternatively, if you'd like to download the git repository somewhere specific, you can go to that directory and use:

.. code-block:: console

    git clone git@github.com:dpthorngren/Eggman.git
    cd Eggman
    pip install .

If you'd like to run the tests, install the test dependencies with ``pip install .[develop]``, then run them with ``pytest`` from the Eggman directory. 
