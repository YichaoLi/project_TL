#ifndef _H_ALM_
#define _H_ALM_


#include <stdio.h>
#include <math.h>

int alm_getlm(int lmax, int i, int *lm);

int alm_getidx(int lmax, int l, int m);


int alm_getsize(int lmax, int mmax );


int alm_getlmax(int s, int mmax);

#endif
