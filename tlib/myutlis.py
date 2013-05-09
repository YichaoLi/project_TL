import os
import numpy as np
import pylab as pl
import scipy.misc as sc
import ConfigParser
import json
from collections import Sequence
from itertools import chain, count

import genclass as gc
import datanalysis.f77io as f77io


""" ------ >>> string related <<< ----- """
def add_quotation(string):
    # add quotation mark, useful in def_variable in parameters.py
    if isinstance(string, str):
        return '\''+string+'\''
    else:
        return string

def str2bool(v):
    if( (v.lower() in ("yes", "true", "t") ) ==True):
        return True 
    else:
        raise Exception

def convert_strdict(dic):
    for key in dic:
        try: 
	    res=int(dic[key])
	except: 
	    try: 
	        res=float(dic[key])
	    except: 
	        try: 
		    res=str2bool(dic[key])
		except: 
		    try:
                        res=json.loads(dic[key])
		    except:
		        res=dic[key]
	dic[key]=res
    return dic


""" ------ >>> array/list related <<< ----- """
def index_select(index):
    """
    >>> return a list of array index, each element represent the
        index of the element.
    """
    return zip(*index)


def rollarr(x, shift, ax='all'):
    """ shift the ndarray,  
    """
    if ax=='all':
        if not isinstance(shift, int):
            raise Excpetion('`shift` should be `int` when axis=`all`')
	    
	ax=np.arange(len(x.shape))
	s=[shift]*len(x.shape)

        if shift==0:
	    return x

    else:
        #if (isinstance(shift, list))|(isinstance(shift, np.ndarray)):
	#    raise Exception('shift should be a `list` or `ndarray`')

	if (len(x.shape))!=len(ax):
	    raise Exception

        s=shift


    y=x
    for i in range(len(ax)):
       y=np.roll(y, s[i], axis=ax[i])


    return y





def depth(seq):

    if isinstance(seq, np.ndarray)|isinstance(seq, list)|isinstance(seq, tuple):
        return len(np.array(seq).shape)

    '''
    if isinstance(seq, np.ndarray):
        return len(seq.shape)
    elif isinstance(seq, list)|isinstance(seq, tuple):
        pass
    else:
        print 'depth error: NOT `ndarray` or `list` or `tuple`'
        return -1

    seq = iter(seq)
    try:
        for level in count():
            seq = chain([next(seq)], seq)
            seq = chain.from_iterable(s for s in seq if isinstance(s, Sequence))
    except StopIteration:
        return level

    '''



""" ------ >>> ploting reated <<< ----- """
def boundary(data, fact=[1, 1]):
    return [fact[0]*np.amin(data), fact[1]*np.amax(data)]

# get the default plot line style
def get_ltype( ind=16, colorfirst=False, cset=['k','b','y','r'], lset=['-','--','-.',':']):
    if(colorfirst):
        ltype = np.ravel(np.array([ [cset[i]+lset[j] for i in range(len(cset))] \
                                for j in range(len(lset)) ]) )
    else:
        if(ind<=len(cset)*len(lset)):
            ltype = np.ravel(np.array([ [cset[i]+lset[j] for j in range(len(lset))] \
                                for i in range(len(cset)) ]) )
        else:
            rnd = ind/(len(cset)*len(lset)) + 1
            print "rnd={0}, ind={1}, len={2}".format(rnd, ind,  len(cset)*len(lset) )
    
            lt = np.ravel(np.array([ [cset[i]+lset[j] for j in range(len(lset))] \
                                    for i in range(len(cset)) ]) )
	    ltype = np.array(list(lt)*rnd )



    return ltype



""" ------ >>> I/O Modules <<< ----- """
# if not exist, create director dir
def createdir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)
        print 'created directory {0}'.format(dir)


# reading 
def read_npz(fname):
    dat=np.load(fname)
    data=[dat[dat.files[i]]  for i in range(len(dat.files))]
    return data

"""
def load_nparray(fname):
    data=
    return data
"""


def read_datacube(fname, data_form, grid=128):
    if(data_form=='npy'):
        data=np.load(fname) 

    elif(data_form=='npz'):
        data=read_npz(fname)

    elif(data_form=='binary'):
        data=f77io.getdatacube(ngrid=grid, prefix='', file=fname)
    else:
        raise Exception('read_datacube error.')

    return data



""" ------ >>> mathematical modules <<< ----- """
def taylorexp(x, b, nstart=0, array=False):
    if(array):
        f=np.zeros(x.shape)
    else:
        f=0

    ln=len(b)
    for n in range(ln):
        f+=b[n]*np.power(x, n+nstart)/sc.factorial(n+nstart)

    return f




""" ------ >>> other modules <<< ----- """
def getfname(prefix, suffix, n):
    
    fn = []
    for i in range(n):
        fn.append( prefix + str(i) + suffix )
    
    return fn
