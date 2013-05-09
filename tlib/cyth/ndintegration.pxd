"""
>>> API of n-dimensional Monte-Carlo integration 
"""

#cdef cuhre(void (* integrand), ndim, ncomp, userdata, erel, eabs)
cdef cuhre(void (* integrand), ndim, ncomp, void *userdata, erel, eabs)

cdef vegas(void (* integrand), ndim, ncomp, userdata, erel, eabs)


cdef double integrand_scale(int ndim, double *x, double *lower, double *upper)
