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
Extension('cubicspline', ['cubicspline.pyx'], 
          include_dirs=[np.get_include()]), 

#Extension('power', ['power.pyx'], 
#          include_dirs=[np.get_include()]), 

Extension('ndintegration', ['ndintegration.pyx'], 
          include_dirs=[np.get_include()], 
          libraries=['cuba'] ), 

Extension('mymath', ['mymath.pyx'], 
          include_dirs=[np.get_include()]), 

Extension('mpiutil', ['mpiutil.pyx'], 
           include_dirs=[np.get_include(), 
	   mpi4py_include]), 

Extension("sf", ["sf.pyx"], 
          include_dirs= [np.get_include(), cython_gsl.get_cython_include_dir(),
	  '../../clib', ], 
          libraries=cython_gsl.get_libraries()+['mysf'],
          library_dirs=[cython_gsl.get_library_dir(), '../../clib/lib'] ),

Extension("alm", ["alm.pyx"], 
          include_dirs= [np.get_include(), '../../clib', ], 
          libraries=['mysf'],
          library_dirs=[ '../../clib/lib'] ),

Extension("cgalio", ["cgalio.pyx", "../c/cgalio.c"],
          libraries=cython_gsl.get_libraries()+["cuba"],
          include_dirs=[np.get_include(), cython_gsl.get_cython_include_dir(), "../c/"],
          library_dirs=[cython_gsl.get_library_dir(), './', "../c/"]    ) 

 ]

setup( name = 'cython wrapper', cmdclass = {'build_ext': build_ext},
       ext_modules = ext_modules )
