import numpy as np
cimport numpy as np
from cython_gsl cimport *
from libc.math cimport floor, sqrt
from libc.stdlib cimport malloc, free
cimport clib.alm as alm
cimport genscript.cyth.mymath as ma


cpdef getlm(int lmax, int i):
    cdef int l, m, *lm

    lm=ma.imat1(2)

    alm.alm_getlm(lmax, i, lm)
    l, m = lm[0], lm[1]

    free(lm)

    return l, m

    
cpdef getidx(int lmax, int l, int m):
    return alm.alm_getidx(lmax, l, m)
    
cpdef getsize(int lmax, int mmax ):
    return alm.alm_getsize(lmax, mmax)
    
cpdef getlmax(int s, int mmax):
    return alm.alm_getlmax(s, mmax)
