from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np
import os
import cython_gsl

machine=os.environ['MACHINE']
print 'compiling on', machine


if machine=='mylaptop':
    mpi4py_include='/Library/Frameworks/Python.framework/Versions/Current/lib/python2.7/site-packages/mpi4py/include' 
elif machine=='gwln':
    mpi4py_include='/home/wangxin/software/python/lib/python2.7/site-packages/mpi4py/include'
elif machine=='HHPC':
    mpi4py_include='/home/wangxin/software/python/python/2.7.3/lib/python2.7/site-packages/mpi4py/include'
    


ext_modules = [ 

Extension('power', ['power.pyx'], 
          include_dirs=[np.get_include()])
 ]

setup( name = 'cython wrapper', cmdclass = {'build_ext': build_ext},
       ext_modules = ext_modules )
