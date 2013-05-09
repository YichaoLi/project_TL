import numpy as np
cimport genscript.cyth.cubicspline as cs




cdef class PowerSpectrum(object):


    def __init__(self, cp, data=None, data_y=None, window=None, R=None):

        """
	>>> Initialize cosmological parameters
	"""
        #self.cp=cp
        #self.window=window

        #(self.h, self.om0, self.omB, self.omL, self.omR, self.n) = (
        #self.cp.h, self.cp.omec+self.cp.omeb, self.cp.omeb,
        #self.cp.omex, 4.17e-5/(self.cp.h)**2, self.cp.ns[0])

        (self.h, self.om0, self.omB, self.omL, self.omR, self.n) = (
        cp.h, cp.omec+cp.omeb, cp.omeb, cp.omex, 4.17e-5/(cp.h)**2, cp.ns[0])


        self.omC = self.om0-self.omB
        self.delH = 1.94e-5*self.om0**(-0.785 - 0.05*np.log(self.om0))*np.exp(-0.95*(self.n - 1)\
                    - 0.169*(self.n - 1)**2)

        #self.t  = <double>self.bbks
        self.D0 = self.D1_unorm(0)

        ##   h^-1 normalization for sigma8=0.9
        self.bondEfsNorm = 216622.0
        ##   bare normalizatin by Pan (apparently correct for h_1)   
        self.bondEfsNorm = 122976.0


        if data!=None:
            """ Initialize Powerspectrum interpolator """

            if isinstance(data, str):
                dat = np.array(list(zip(* np.genfromtxt(data) )) )  
                k, pk=dat[0,:], dat[1,:]

                #print dat

            elif isinstance(data, np.ndarray):
                if len(data.shape)==2:
                   k, pk=data[0,:], data[1,:]
                elif len(data.shape)==1:
                   k, pk=data, data_y
                else:
                    raise Exception()

            else:
                raise Exception()

            self.power=cs.Interpolater(k, pk)

            if window!=None:
                self.power_smooth=cs.Interpolater(k, pk*window(k, R))

        

    def __call__(self, k, smooth=False):
        """ return initialized power spectrum from `data` """
        return self.evaluate(k, smooth)

   
    """ interpolate power spectrum """
    def evaluate(self, k, smooth):
        """ return the evaluation value of powerspectrum, surpport ndarray of `k`
	    from Interpolater class """

        if smooth:
            return self.power_smooth(k)
        else:
            return self.power(k)




    """ >>> Transfer functions <<< """
    def bbks(self,k):
        """
        BBKS transfer function
        """
        q = k/self.om0/self.h**2*np.exp(self.omB + np.sqrt(2*self.h)*self.omB/self.om0)
        return np.log(1.0 + 2.34*q)/(2.34 *q) *(1 + 3.89*q + (16.1*q)**2 + (5.46*q)**3 + (6.71*q)**4)**(-0.25)


    def bondEfs(self,k):
        """
        alternative transfer function fit by Bond end Efstathiou for
        gif simulations which apparently use it

        this is normalized differently from the bbks transfer function
        and can only be used at z=0
        """
        Alpha=6.4
        b=3.0
        c=1.7
        nu=1.13
        gamma=self.om0*self.h
        q=k/gamma

        ##   h^-1 normalization for sigma8=0.9
        self.bondEfsNorm = 216622.0
        ##  bare normalizatin by Pan (apparently correct for h_1)   
        #self.bondEfsNorm = 122976.0

        return self.bondEfsNorm/(1.0+(Alpha*q +(b*q)**1.5 + (c*q)**2)**nu)**(2.0/nu)


    """ >>> Growth function
    """

    def om0z(self,z):
        """
        Omega as a function of redshift
        """
        return self.om0*(1.0 + z)**3/(self.omL+self.omR*(1.0 + z)**2+self.om0*(1.0 + z)**3)


    def omLz(self,z):
        """
        Lambda as a function of redshift
        """
        return self.omL/(self.omL + self.omR*(1.0 + z)**2 + self.om0*(1.0 + z)**3)


    def D1_unorm(self,z):
        """
        Growth function
        """
        return 5*self.om0z(z)/(1.0 + z)/2.0/(self.om0z(z)**(4./7.0) - self.omLz(z) + \
               (1.0 + self.om0z(z)/2.0)*(1.0 + self.omLz(z)/70.0))


    def D1(self, z):
        """ Renormalized growth function
        """
        return self.D1_unorm(z)/self.D0
