from setuptools import setup, Extension
from Cython.Build import cythonize
import os


ext = Extension(
    "geom", ["src/geom.pyx"],
    libraries=[":cspice.a", "gsl"],
    library_dirs=[os.getcwd()+"/cspice/lib/"],
    include_dirs=[os.getcwd()+"/cspice/include/"],
)

setup(
    name="geom",
    version='0.1',
    description="A module for computing the geometry of transiting ellipsoidal planets.",
    author="Daniel Thorngren",
    ext_modules=cythonize([ext], compiler_directives={'embedsignature': True, 'language_level': "3"}),
    install_requires=['cython']
)
