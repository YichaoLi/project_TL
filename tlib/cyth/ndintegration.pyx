import numpy as np
cimport numpy as np
cimport genscript.cyth.cuba as cuba




#cdef cuhre(void (* integrand), ndim, ncomp, userdata, erel, eabs):
cdef cuhre(void (* integrand), ndim, ncomp, void *userdata, erel, eabs):
    """
    >>> wrapper of cuhre caller from `cuba`
    """
    cdef:
        int nd, nc, verbose, nregions, neval, fail
        double epsrel, epsabs, integral[10], error[10], prob[10]

    nd, nc = ndim, ncomp 
    epsrel, epsabs = erel, eabs

    verbose=2
    key=0
    mineval, maxeval= 50, 50000000

    cuba.Cuhre(nd, nc, <cuba.integrand_t>integrand, <void *>userdata,
              epsrel, epsabs, verbose, mineval, maxeval, key,
              &nregions, &neval, &fail, integral, error, prob);

    res=[]
    for i in range(nc):
        res.append(integral[i])


    return np.array(res)



cdef vegas(void (* integrand), ndim, ncomp, userdata, erel, eabs):

    return




cdef double integrand_scale(int ndim, double *x, double *lower, double *upper):
    cdef double ran, jacobian, scaled
    jacobian=1
    for i in range(ndim):
        ran=upper[i]-lower[i]
        jacobian*=ran
        scaled=lower[i]+x[i]*ran
        x[i]=scaled

    return jacobian


