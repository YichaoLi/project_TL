#include <stdio.h>
#include <math.h>



int alm_getlm(int lmax, int i, int *lm) {

  /*     ---  Get the l and m from index and lmax. ---
      Parameters
      ----------
      lmax : int
        The maximum l defining the alm layout
      i : int 
        The index for which to compute the l and m.
  */

  int l, m;

  m = (int) (ceil(((2*lmax+1)-sqrt(pow((2*lmax+1),2)-8*(i-lmax)))/2)) ;
  l = i-m*(2*lmax+1-m)/2;

  lm[0]=l;  lm[1]=m;

  return 0;
  }


int alm_getidx(int lmax, int l, int m) {
  /* Returns index corresponding to (l,m) in an array describing alm up to lmax.
  
      lmax : int
        The maximum l, defines the alm layout
      l : int
        The l for which to get the index
      m : int
        The m for which to get the index
      
      Returns
      -------
      idx : int
        The index corresponding to (l,m)
  */

  return m*(2*lmax+1-m)/2+l;
  }



int alm_getsize(int lmax, int mmax ) {
    /* Returns the size of the array needed to store alm up to *lmax* and *mmax*

        Parameters
        ----------
        lmax : int
          The maximum l, defines the alm layout
        mmax : int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        size : int
          The size of the array needed to store alm up to lmax, mmax.
   */

  if ((mmax < 0) ||( mmax > lmax) )
     mmax = lmax;

  return mmax * (2 * lmax + 1 - mmax) / 2 + lmax + 1;
  }



int alm_getlmax(int s, int mmax) {
    /* Returns the lmax corresponding to a given array size.
        
        Parameters
        ----------
        s : int
          Size of the array
        mmax : None or int, optional
          The maximum m, defines the alm layout. Default: lmax.

        Returns
        -------
        lmax : int
          The maximum l of the array, or -1 if it is not a valid size.
    */
  double x;

  if (mmax >= 0)
      x = (2 * s + pow(mmax,2) - mmax - 2) / (2 * mmax + 2);
  else
      x = (-3 + sqrt(1 + 8 * s)) / 2;

  if (x != floor(x))
      return -1;
  else
      return (int) x;
  }
