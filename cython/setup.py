#!/usr/bin/env python
from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("rnppy",
                sources=["/home/alblc/dev/RnP/cython/rnppy.pyx"],
                language="c++",
                include_dirs=["/home/alblc/dev/RnP/c++/",numpy.get_include()],
                libraries=["rnp"],
                extra_link_args=["-L/home/alblc/dev/RnP/build/bin/lib/"],
                )],
    
)
