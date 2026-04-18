import os

from Cython.Build import cythonize
from setuptools import Extension, setup

ext = Extension(
    "eggman.cy_eggman",
    ["src/eggman/cy_eggman.pyx"],
    libraries=["gsl", "gslcblas"],
    # libraries=["cspice", "gsl", "gslcblas"],
    # library_dirs=[os.getcwd() + "/cspice/lib/"],
    # include_dirs=[os.getcwd() + "/cspice/include/"],
)

directives = {
    'language_level': "3",
    'embedsignature': True,
    'cdivision': True,
    'initializedcheck': False,
    'boundscheck': False,
    'binding': True,
}

# Hack to detect if coverage data is being generated
if os.path.exists("./htmlcov"):
    ext.define_macros = [("CYTHON_TRACE_NOGIL", "1")]
    directives['linetrace'] = True

setup(ext_modules=cythonize([ext], compiler_directives=directives, annotate=True))
