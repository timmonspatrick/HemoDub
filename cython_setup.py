from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


##Usage:
##python3 cython_setup.py build_ext --inplace

setup(
    ext_modules=(cythonize("residue_distribution.pyx"))
    )
