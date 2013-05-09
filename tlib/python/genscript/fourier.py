"""
>>> Continuous Fourier Transform routines
"""
import numpy as np
import scipy.integrate as integrate
from scipy import interpolate


""" 
    >>>    CFT re-interpolate from FFT   <<<
"""
class FourierTransform(object):
    """ class for calculating continuous Fourier transformation 
    """

    def __init__(self, dimension, xlim=10, dx=0.1, data=None,
                 function=None, fargs=None, empty=False):
        """
	>>> `input`
               funtion: input function to generate data
	       fargs:   arguments of input function
	"""

	if (data==None)&(function==None):
	    raise Exception('`data` and `function` cannot be `None` at the same time')

        if dimension>2:
	    raise Exception('dimension is larger than 2')

	self.dimension=dimension

        if empty!=True:
            if function!=None:
                self.data_fromfunction(dimension, function, fargs, xlim, dx)
	    else:
	        self.import_data(dimension, data)

            """ get the frequency `k` of the data """
	    self.getfreq(dimension)



    """ ----- FFT methods: prepare data ----- """
    def data_fromfunction(self, dimension, function, fargs, xlim, dx):
        """ prepare the data from given function with arguments `fargs` """

	self.dx=dx
	self.xlim=xlim
	self.xlim_tot=2*xlim

        """ 1 dimensional """
        if dimension==1:
            self.x=np.arange(-self.xlim_tot, self.xlim_tot, self.dx)

            n1,n2=np.nonzero((self.x>=-self.xlim)&(self.x<-self.xlim+dx))[0], \
                   np.nonzero((self.x>=self.xlim-self.dx)&(self.x<self.xlim))[0]

	    self.data=np.zeros(len(self.x))
	    self.data[n1:n2]=function(self.x[n1:n2], fargs)

        """ 2 dimensional """
	if dimension==2:
	    if(len(dx)!=2)|(len(xlim)!=2):
	        raise Exception('dx or xlim should be 2-dimensional')

            self.x=np.array([np.arange(-self.xlim_tot[i], self.xlim_tot[i], \
	                    self.dx[i]) for i in range(2) ])
            self.xcoord, self.ycoord=np.meshgrid(self.x[0],self.x[1])

            xn1,xn2=np.nonzero((self.x[0]>=-self.xlim[0])&(self.x[0]<-self.xlim[0]+self.dx[0]))[0],\
	             np.nonzero((self.x[0]>=self.xlim[0]-self.dx[0])&(self.x[0]<self.xlim[0]))[0]
            yn1,yn2=np.nonzero((self.x[1]>=-self.xlim[1])&(self.x[1]<-self.xlim[1]+self.dx[1]))[0],\
	             np.nonzero((self.x[1]>=self.xlim[1]-self.dx[1])&(self.x[1]<self.xlim[1]))[0]

            print xn1, xn2, yn1, yn2

            self.data=np.zeros(self.xcoord.shape)
            self.data[xn1:xn2,yn1:yn2]=function(self.xcoord[xn1:xn2,yn1:yn2], \
	                                        self.ycoord[xn1:xn2, yn1:yn2], fargs)




    def import_data(dimension, data):
        """ prepare the data from existing array `data` """
        pass

	
    def getfreq(self, dimension):
        """ get the array of frequency """
        if dimension==1:
            self.k = np.fft.fftshift(np.fft.fftfreq(self.x.size, d=self.dx) )

        elif dimension==2:
	    pass



    def cft1d(self):
        """ get the 1D FFT """
        return np.fft.fftshift(np.fft.fft(self.data))*self.dx

    def icft1d(self):
        """ get the 1D FFT """
        return np.fft.fftshift(np.fft.ifft(self.data))*self.dx


    def cft2d(self):
        return np.fft.fftshift(np.fft.fft2(self.data))*self.dx[0]*self.dx[1]

    def icft2d(self):
        return np.fft.fftshift(np.fft.ifft2(self.data))*self.dx[0]*self.dx[1]








""" --------------------------------------------------"""

def cft_old(data, x, k, lim, with2pi=True):
    """ return the continuous Fourier transformation through numerical integration
        input: 
	      data: 2D array, 0th row is x-axis and 1st row the y-axis
	      lim: the limit of data
	      k:    wavenumber k
	output: 
    """
    print lim
    datrep=interpolate.splrep(x, data, s=0) 

    if with2pi:
        fact=2.*np.pi
    else:
        fact=1.

    datft_real, datft_imag=[], []
    for ki in k:
        intg_real=lambda t: interpolate.splev(t,datrep,der=0)*np.cos(fact*t*ki)
        intg_imag=lambda t: interpolate.splev(t,datrep,der=0)*np.sin(fact*t*ki)

        res_real, err_real = integrate.quad(intg_real,lim[0],lim[1],limit=int(1e2))
        res_imag, err_imag = integrate.quad(intg_imag,lim[0],lim[1],limit=int(1e2))

        print ki, res_real, res_imag

	datft_real.append(res_real)
	datft_imag.append(res_imag)

    return datft_real, datft_imag



def cft2d_old(data, x, k, lim ):
    pass

