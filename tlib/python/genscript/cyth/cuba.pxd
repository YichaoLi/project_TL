# API definition, converted from cuba.h 
cdef extern from "cuba.h":
    ctypedef int (*integrand_t)(int *ndim, double *x,
      int *ncomp, double *f, void *userdata)
    

#    ctypedef int (*integrand_t)(int *ndim, double x[],
#      int *ncomp, double f[], void *userdata)
    
#    ctypedef void (*peakfinder_t)(int *ndim, double b[],
#      int *n, double x[])
    
    # -- Vegas -- #
    void Vegas(int ndim, int ncomp, integrand_t integrand, void *userdata,
      double epsrel, double epsabs, int flags, int seed,
      int mineval, int maxeval, int nstart, int nincrease, int nbatch,
      int gridno, char *statefile, int *neval, int *fail,
      double integral[], double error[], double prob[])

    # -- Cuhre -- #
    void Cuhre(int ndim, int ncomp, integrand_t integrand, void *userdata,
      double epsrel, double epsabs, int flags, int mineval, int maxeval,
      int key, int *nregions, int *neval, int *fail,
      double integral[], double error[], double prob[])
    

