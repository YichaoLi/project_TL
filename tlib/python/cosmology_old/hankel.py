import pylab as M
import numpy as N
import bessint as BI

class Hankel3D:
    def __init__(self,nMax):
        """

        """
        self.bessint = BI.BesselIntegrals(0.5,nMax)

    def transform1(self,f,x,n,h,pk2xi=1):
        """

        """
        bi = 1.0/x**3*self.bessint.besselInt(lambda z:z**1.5*f(z/x),n,h)
        pf = 1.0/(2.0*M.pi)**1.5
        if pk2xi==0:
            pf =1.0/pf
        return pf*bi

    def transform(self,f,x,n,h,pk2xi=1):
        """

        """
        if N.isscalar(x):
            return self.transform1(f,x,n,h,pk2xi)
        else:
            return M.array(map(lambda z:self.transform1(f,z,n,h,pk2xi),x))
