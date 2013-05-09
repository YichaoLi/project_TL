import numpy as np
cimport numpy as np
from cython_gsl cimport *
from libc.math cimport floor, sqrt

cimport clib.mysf as mysf


cpdef wigner3j_m0(int l1, int l2, int l3):
    ''' wigner 3j symbol with m1=m2=m3=0
    '''
    return mysf.ThreeJ_l1l2l3(l1,l2,l3)

cpdef wigner3j(int l1, int l2, int l3, int m1, int m2, int m3):
    return gsl_sf_coupling_3j(<int>2*l1,<int>2*l2,<int>2*l3,<int>2*m1,<int>2*m2,<int>2*m3)


cpdef Gaunt(int l1, int l2, int l3, int m1, int m2, int m3):
    cdef double n

    n=sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/4/np.pi)

    #return n*wigner3j(l1,l2,l3,0,0,0)*wigner3j(l1,l2,l3,m1,m2,m3)
    return n*wigner3j_m0(l1,l2,l3)*wigner3j(l1,l2,l3,m1,m2,m3)


cpdef Gaunt_v1(int l1, int l2, int l3, int m1, int m2, int m3):
    cdef double n
    n=sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/4/np.pi)

    return n*wigner3j(l1,l2,l3,0,0,0)*wigner3j(l1,l2,l3,m1,m2,m3)
