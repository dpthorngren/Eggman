from setuptools import setup, Extension
from Cython.Build import cythonize
import os


ext = Extension(
    "eggman", ["src/eggman.pyx"],
    libraries=[":cspice.a", "gsl"],
    library_dirs=[os.getcwd()+"/cspice/lib/"],
    include_dirs=[os.getcwd()+"/cspice/include/"],
)

setup(
    name="eggman",
    version='0.2',
    description="Code for calculating the geometry and transit depths of piecewise ellipsoidal objects.",
    author="Daniel Thorngren",
    ext_modules=cythonize([ext], compiler_directives={'embedsignature': True, 'language_level': "3"}),
    install_requires=['cython']
)
