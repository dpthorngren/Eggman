from Cython.Build import cythonize
from setuptools import Extension, setup

ext = Extension(
    "eggman.cy_eggman",
    ["src/eggman/cy_eggman.pyx"],
    libraries=["gsl", "gslcblas"],
)

directives = {
    'language_level': "3",
    'embedsignature': True,
    'cdivision': True,
    'initializedcheck': False,
    'boundscheck': False,
    'binding': True,
}

setup(ext_modules=cythonize([ext], compiler_directives=directives, annotate=True))
