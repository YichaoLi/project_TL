""" Routines for reading files
"""
import numpy as np
import array as arr



def readgrid_numpy_miguel(filename, ngrid, dtype):

    if dtype=='byte':
        d=np.fromfile(filename, np.int8)
        header=256
    elif dtype=='float':
        d=np.fromfile(filename, np.float32)
        header=256*8/32
    elif dtype=='double':
        d=np.fromfile(filename, np.float)
        header=256*8/64
        
    return d[:header], d[header:].reshape(ngrid, ngrid, ngrid)



def readgrid(filename,ngrid=512, dtype='float'):
    """ ?read the grid data, remember the `256` byte dummy header, 
    """

    if dtype=='double':
        array_type= 'd'
	numpy_type= np.float

    elif dtype=='float':
        array_type= 'f'
	numpy_type= np.float32

    elif dtype=='byte':
        #return readgrid_numpy_miguel(filename, ngrid, dtype)
        head, d=readgrid_numpy_miguel(filename, ngrid, dtype)
	return d

    else:
        raise Exception('unsupported datatype')

    F = open(filename,'rb')
    header=F.read(256)
    den = arr.array(array_type)
    den.fromfile(F,ngrid*ngrid*ngrid)
    F.close()


    den = np.array(den).reshape((ngrid,ngrid,ngrid)).astype(numpy_type)

    return den



def rgrid(filename, ngrid=512, dtype='float', comp=1):
    """ ?read the binary grid data without any header
    """

    if dtype=='double':
        array_type= 'd'
	numpy_type= np.float

    elif dtype=='float':
        array_type= 'f'
	numpy_type= np.float32

    else:
        raise Exception('unsupported datatype')

    F = open(filename,'rb')
    den = arr.array(array_type)

    den.fromfile(F,ngrid*ngrid*ngrid*comp)
    F.close()

    if comp==1:
        #den.fromfile(F,ngrid*ngrid*ngrid)
        #F.close()
        den = np.array(den).reshape((ngrid,ngrid,ngrid)).astype(numpy_type)
    else:
        den = np.array(den).reshape((ngrid,ngrid,ngrid,comp)).astype(numpy_type)


    return den




def get_flines(fname):
    with open(fname, 'r') as f:
        for i, l in enumerate(f):
            pass
    return i + 1




def read_stream(fname):
    dat=[]
    for lines in open(fname, 'r'):
        dat+=map(float, lines.split()) 

    return np.array(dat)



def read_miguel_cgal(fname):
    ''' integer: indicate # of particles
        then all doubles
    '''
    
    f=open(fname,'rb')
    n=f.read(2)
    print n
    #array_type= 'd'

    '''
    dat=arr.array(array_type)
    dat.fromfile(f, n)
    F.close()

    d = np.array(dat).reshape((ngrid,ngrid,ngrid)).astype(numpy_type)
    return d
    '''

