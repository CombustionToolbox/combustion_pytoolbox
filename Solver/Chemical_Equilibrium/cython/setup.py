import os
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

dirname = os.path.dirname(__file__)
filename = os.path.join(dirname, 'GibbsMinimizationCython.pyx')
setup(ext_modules = cythonize(filename), include_dirs=[numpy.get_include()])