from Cython.Build import cythonize
from setuptools import Extension, setup

ext = Extension(
    name="eggman.cy_eggman",
    sources=["src/eggman/cy_eggman.pyx"],
    libraries=["gsl"],
    include_dirs=["src/include/", "src/cpp_source/"],
    extra_compile_args=['-O3'],
)

directives = {
    'language_level': "3",
    'embedsignature': False,
    'cdivision': True,
    'initializedcheck': False,
    'boundscheck': False,
    'binding': True,
}

setup(ext_modules=cythonize([ext], compiler_directives=directives, annotate=True))
