cdef extern from "alm.h":

    int alm_getlm(int lmax, int i, int *lm)
    
    int alm_getidx(int lmax, int l, int m)
    
    int alm_getsize(int lmax, int mmax )
    
    int alm_getlmax(int s, int mmax)
