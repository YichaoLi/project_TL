""" API of myinterpolate.c """
from cython_gsl cimport *


cdef extern from "../../myinterpolate.h":
    ctypedef struct Interpar:
        pass
        #gsl_interp_accel *acc
        #gsl_spline *spl
    

    int myinterp_init( Interpar *f, double *x, double *y, int n)
    int gslinterp_init(Interpar *f, double *x, double *y, int n)

    void * gslinterp_void_init(double *x, double *y, int n)

    int myinterp_free( Interpar *f)
    double myinterp(Interpar *f, double x)
    
