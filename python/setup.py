#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
if sys.hexversion < 0x03000000: # uniform unicode handling for both Python 2.x and 3.x
    def u(x):
        return x.decode('utf-8')
    def textfileopen(filename):
        return open(filename, mode='r')
else:
    def u(x):
        return x
    def textfileopen(filename):
        return open(filename, mode='r', encoding='utf_8')
u('''
  fastcluster: Fast hierarchical clustering routines for R and Python

  Copyright © 2011 Daniel Müllner
  <http://math.stanford.edu/~muellner>
''')
#import distutils.debug
#distutils.debug.DEBUG = 'yes'
from numpy.distutils.core import setup, Extension

with textfileopen('fastcluster.py') as f:
    for line in f:
        if line.find('__version_info__ =')==0:
            version = '.'.join(line.split("'")[1:-1:2])
            break

print('Version: ' + version)

setup(name='fastcluster', \
      version=version, \
      py_modules=['fastcluster'], \
      description='Fast hierarchical clustering routines for R and Python.', \
      long_description="""
This library provides Python functions for hierarchical clustering. It generates hierarchical clusters from distance matrices or from vector data.

Part of this module is intended to replace the functions

    linkage, single, complete, average, weighted, centroid, median, ward

in the module scipy.cluster.hierarchy with the same functionality but much faster algorithms. Moreover, the function 'linkage_vector' provides memory-efficient clustering for vector data.

The interface is very similar to MATLAB's Statistics Toolbox API to make code easier to port from MATLAB to Python/Numpy. The core implementation of this library is in C++ for efficiency.
""",
      ext_modules=[Extension('_fastcluster',
                             ['../src/fastcluster_python.cpp'],
                  # Feel free to uncomment the line below if you use the GCC.
                  # This switches to more aggressive optimization and turns
                  # more warning switches on. No warning should appear in
                  # the compilation process.
                  #
                  # Also, the author's Python distribution generates debug
                  # symbols by default. This can be turned off, resulting a in
                  # much smaller compiled library.
                  #
                  #extra_compile_args=['-O3', '-Wall', '-Wconversion', '-Wsign-conversion', '-g0', '-march=native', '-mtune=native', '-fno-math-errno'],
                  # (no -pedantic -Wextra, -ansi)
     )],
      keywords=['dendrogram', 'linkage', 'cluster', 'agglomerative', 'hierarchical', 'hierarchy', 'ward'],
      author=u("Daniel Müllner"),
      author_email="fastcluster@math.stanford.edu",
      license="BSD <http://opensource.org/licenses/BSD-2-Clause>",
      classifiers = ["Topic :: Scientific/Engineering :: Information Analysis",
                     "Topic :: Scientific/Engineering :: Artificial Intelligence",
                     "Topic :: Scientific/Engineering :: Bio-Informatics",
                     "Topic :: Scientific/Engineering :: Mathematics",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 2",
                     "Programming Language :: Python :: 3",
                     "Programming Language :: C++",
                     "Operating System :: OS Independent",
                     "License :: OSI Approved :: BSD License",
                     "Intended Audience :: Science/Research",
                     "Development Status :: 5 - Production/Stable"],
      url = 'http://math.stanford.edu/~muellner',
      )
