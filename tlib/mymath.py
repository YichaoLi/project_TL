import math
import numpy as np
import pylab as pl



''' ->> array/list related routines <<- '''
def nearest_sval(arr, val):
    # ->> find nearest value to `val` in `arr` 
    idx =  (np.abs(arr-val)).argmin()
    return arr[idx]

def nearest(arr, val):
    if isinstance(val, list)|isinstance(val, np.ndarray):
        return np.array([ nearest_sval(arr, val[i]) for i in range(len(val)) ])
    else:
        return nearest_sval(arr, val)


def nearest_idx_sval(arr, val):
    # ->> find nearest index to `val` in `arr` 
    return (np.abs(arr-val)).argmin()

def nearest_idx(arr, val):
    if isinstance(val, list)|isinstance(val, np.ndarray):
        return np.array([ nearest_idx_sval(arr, val[i]) for i in range(len(val)) ])
    else:
        return nearest_idx_sval(arr, val)


def inrange(arr, val):
    # ->> judge whether `val` in range of `arr`, order of `arr` is unimportant
    min, max = arr.min(), arr.max()

    if (val-min>0)&(val-max<0):
        return True
    else:
        return False



''' ->> mathematical functions <<- '''

def function_abs(func, x, *args):

    try: 
        return func(x)
    except:
        if isinstance(x, list)|isinstance(x, np.ndarray):
	    try: 
	        f= [triangle_ele(i)  for i in x]
                return np.array(f)
	    except:
	        f=[triangle(i) for i in x]
		return np.array(f)
        else:
            raise Exception('error.')



def function2d_abs(func, x, y, *args):

    try: 
        return func(x,y)
    except:
        if isinstance(x, list)|isinstance(x, np.ndarray):
	    # only check x, assume y is the same
	    try: 
	        f= [func(x[i], y[i], *args)  for i in range(len(x))]
                return np.array(f)
	    except:
	        f= [function2d_abs(func, x[i], y[i], *args)  for i in range(len(x))]
		return np.array(f)
        else:
            raise Exception('error.')


def Heaviside_step(x):
    """ Use the right-continuous convention that H(0)= 1 """
    return (np.sign(x)+1)/2.

def boxfunc(x, xa, xb):
    """ Return the box function from xa to xb, assume xa<=xb """
    if isinstance(x, np.ndarray):
        b=Heaviside_step(x-xa)-Heaviside_step(x-xb)
	b[np.nonzero((x==xa)|(x==xb))]=1
	return b

    else:
        if (x==xa)|(x==xb):
            return 1
        else:
            return Heaviside_step(x-xa)-Heaviside_step(x-xb)


def triangle_ele(x):
    if abs(x)<1:
        f=1-abs(x)
    else:
        f=0
    return f

def triangle(x):

    try: 
        return triangle_ele(x)
    except:
        if isinstance(x, list)|isinstance(x, np.ndarray):
	    try: 
	        f= [triangle_ele(i)  for i in x]
                return np.array(f)
	    except:
	        f=[triangle(i) for i in x]
		return np.array(f)
        else:
            raise Exception('error.')




""" From utilis.py of Mark
"""

def trapz(y, x=None, ax=-1, method=pl.add.reduce):
    """trapz(y,x=None,ax=-1) integrates y along the given dimension of
    the data array using the trapezoidal rule.

    you can call it with method=pl.add.accumulate to yield partial sums
    """
    y = pl.asarray(y)
    if x is None:
        d = 1.0
    else:
        d = pl.diff(x,axis=ax)
    y = pl.asarray(y)
    nd = len(y.shape)
    s1 = nd*[slice(None)]
    s2 = nd*[slice(None)]
    s1[ax] = slice(1,None)
    s2[ax] = slice(None,-1)

    ans = method(0.5* d * (y[s1]+y[s2]),ax)
    return ans


def trapezoidRule(f,dx):
    """
    numerical integral of a tabulated function f using trapezoid rule

    the function is assumed to be sampled at uniform intervals separated
    by dx between the limits of the integral

    4.1.11 Num.Rec., first ed.
    """
    return (pl.sum(f) - 0.5*f[0]-0.5*f[-1])*dx

def trapezoidArray(f,dx):
    """
    same as trapezoidRule() but return the integral as fn of z
    (mostly in-place operations for optimized memory usage)
    """
    i = pl.add.accumulate(f)
    pl.add(i,-0.5*f,  i)
    pl.add(i,-0.5*f[0],  i)
    pl.multiply(i,dx, i)
    return i



def eig_sort(eigval, eigvec, axis=-1):
    i = list(np.ogrid[[slice(x) for x in eigval.shape]])
    i[axis] = eigval.argsort(axis)

    #eigvec_v=np.rollaxis(eigvec.view(), -1, -2)[i,:]
    eigvec_v=np.rollaxis( np.array([ eigvec[...,k,:][i] for k in \
              range(eigvec.shape[-2] )] ), 0, -1)
    #print eigvec_v.shape

    return eigval[i], eigvec_v 


def sort(d1, d2, axis=0, combine=True):
    i = list(np.ogrid[[slice(x) for x in d1.shape]])
    i[axis] = d1.argsort(axis)

    if combine==True:
        return np.array([d1[i], d2[i]])

    else:
        return d1[i], d2[i]

