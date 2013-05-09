import numpy as np
cimport numpy as np

cimport genscript.cyth.cubicspline as cs
from libc.stdlib cimport malloc, free
from libc.math cimport sqrt, floor





cdef int *imat1(int x):
    cdef:
        int i
        int *vec

    vec = <int *>malloc(sizeof(int) * x)

    return vec 








cdef double *dmat1(int x):
    cdef:
        int i
        double *vec

    vec = <double *>malloc(sizeof(double) * x)

    """
    for i in range(x):
        vec[i]=0
    """
  
    return vec 




cdef double **dmat2(int x, int y):
    cdef:
        int i, j
        double **mat


    mat = <double **>malloc(sizeof(double *) * x)
  
    for i in range(x):
        mat[i] = <double *> malloc(sizeof(double) * y)
    
    """
    for i in range(x):
        for j in range(y):
            mat[i][j]=0
    """

    return mat 


cdef double ***dmat3(int x, int y, int z):
    cdef:
        int i, j, k
        double ***mat

    mat=<double ***>malloc(sizeof(double **)*x)
  
    for i in range(x):
        mat[i]=<double **>malloc(sizeof(double *)*y)
        for j in range(y):
            mat[i][j]=<double *> malloc(sizeof(double)*z)

    """
    for i in range(x):
        for j in range(y):
            for k in range(z):
                mat[i][j][k] = 0
    """

    return mat 






cdef double ****dmat4(int x1, int x2, int x3, int x4):
    cdef:
        int i, j, k, l
        double ****mat

    mat=<double ****>malloc(sizeof(double ***)*x1)
  
    for i in range(x1):
        mat[i]=<double ***>malloc(sizeof(double **)*x2)
        for j in range(x2):
            mat[i][j]=<double **> malloc(sizeof(double *)*x3)
            for k in range(x3):
                mat[i][j][k]=<double *>malloc(sizeof(double)*x4)


    """
    for i in range(x1):
        for j in range(x2):
            for k in range(x3):
                for l in range(x4):
                    mat[i][j][k][l] = 0
    """
    return mat 







cdef double dmin(double x, double y):
    cdef double z
    
    if x-y<=0:
        z=x
    else:
        z=y
    
    return z


cdef double dmax(double x, double y):
    cdef double z
    
    if x-y>=0:
        z = x
    else:
        z = y
    
    return z



cdef int imin(int x, int y):
    return <int> dmin( <double>x, <double>y )
    


cdef int imax(int x, int y):
    return <int> dmax( <double>x, <double>y )




'''
cdef double max_seq(double *x, int xn): 
    cdef:
        int i
        double xm=x[0]
    
    for(i=0; i<xn; i++):
      if( x[i]- xm>0 )
        xm = x[i]
    
    return xm
'''
