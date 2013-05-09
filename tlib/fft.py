import numpy as np
import scipy.fftpack as sfft


def kfftn(dk_shape, ndim, boxsize):
    ''' return frequency grid '''
    k=[sfft.fftfreq(dk_shape[i], 2.*np.pi/boxsize ) for i in \
      range(ndim-1) ] + [sfft.rfftfreq(dk_shape[-1],2.*np.pi/boxsize) ]
    return k

def krfftn(dk_shape, ndim, boxsize):
    ''' return the k grid used for rfftn, where the last axis is done 
        using rfft '''
    k=[sfft.fftfreq(dk_shape[i], 2.*np.pi/boxsize ) for i in \
      range(ndim-1) ] + [sfft.rfftfreq(dk_shape[-1],2.*np.pi/boxsize) ]
    return k



def gradient(d, boxsize, fft_axes=None):
    ''' Get the gradient of data, return the field in real space 
        `ndim`: could be the same as ndim of `d`, i.e. d is scalar field, 
	`fft_axes`:  tupple/list of axes of FFT spatial component
    '''

    if fft_axes==None:
        ndim==d.ndim
	ftype=0
	vcomp_idx=[]
    else:
        ndim=len(fft_axes)
	vcomp_idx=list(set(range(d.ndim))-set(fft_axes))
        ftype=d.ndim-ndim


    if ftype==0:
        #->> d is scalar
        print 'Taking gradient of scalar field'
        dd=d.view()
    
    elif ftype==1:
        #->> d is not vector, use the first `ndim` axis as spatial coord
        print 'Taking gradient of vector field'
	dd=np.rollaxis(d.view(), vcomp_idx[0], 0)
	fft_axes=range(1, d.ndim)
    
    else:
        raise Excpetion('ftype error')
 
 
    # -> rFFT, Fourier transform of each component
    dk=np.fft.rfftn(dd, axes=fft_axes) 
    print 'FFT dk shape:', dk.shape

    if ftype>1:
        dks=dk[0].shape
    else:
        dks=dk.shape


    # -> get k vector
    k_=krfftn(dks, ndim, boxsize)
    kf=np.frompyfunc(lambda x: k_[x], 1, 1)
    # -> k index grid
    kidx=np.indices(dks)

    # -> get the gradient in Fourier space
    dk_grd=np.zeros(dks)


    dgrd=[] #np.zeros([ndim]*ftype+list(d.shape) )
    for i in range(len(dk)):

        for j in range(ndim):
            dk_grd=kf(kidx[j])

        dgrd.append(np.irfftn(dk_grd))


    if ftype==0:
        return
    elif ftype==1:
        return np.array(dgrd)



def hessian(d, boxsize):
    # ->> get the Hessian matrix of scalar `d`
    g=gradient(d, boxsize, fft_axes=None)

    #g_axes= 
    h=gradient(g, boxsize, fft_axes=[])

    return h


def smoothgauss(d, sigma):
    return



class FFTdata(object):
    ''' FFT class, useful for combined operations in Fourier space '''

    def __init__(self, d, boxsize):
        self.data=d
	self.boxsize=float(boxsize)
	self.ngrid=d.shape[0]
        self.kmin=2.*np.pi/float(boxsize)


    def rfft(self, d=None, ngrid=None, boxsize=None, real=True):
        if d==None:
            d=self.data
	if ngrid==None:
	    ngrid=self.ngrid
	if boxsize==None:
	    boxsize=self.boxsize

        if real:
            self.dk=np.fft.rfftn(d)
	    self.k=[sfft.fftfreq(self.dk.shape[i], 2.*np.pi/boxsize ) for i in \
            range(d.ndim-1) ] + [sfft.rfftfreq(self.dk.shape[-1],2.*np.pi/boxsize) ]
	else:
	    raise Exception()

	return #self.dk


    def irfft(self, d=None):
        if d==None:
	    d=self.data

        return np.fft.irfftn(d)


    def smooth_gaussian(self, sigma):
        ''' ->> smooth radius sigma <<- '''
        self.sigmak = 1./sigma

        return


    def gradient(self):
        return


    def hessian(self):
        return
