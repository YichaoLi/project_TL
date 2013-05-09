import numpy as np
import pylab as pl
import genscript.myutlis as mut




def condition_matrix_onetype(condmat, typelist):
    ''' Select certain types from `typelist` of `condmat`.
    '''
    cond=np.array([False]*np.prod(condmat.shape)).reshape(condmat.shape)

    for i in range(len(typelist)):
        cond |= (condmat==typelist[i])

    return cond



def condition_matrix(condmat, typelist):

    d=mut.depth(typelist)

    if d==1:
        return condition_matrix_onetype(condmat, typelist)

    elif d==2:
        cond=[condition_matrix_onetype(condmat,typelist[i]) for i in range(len(typelist))]
	return np.array(cond)

    else:
        raise Exception('condition matrix error')

