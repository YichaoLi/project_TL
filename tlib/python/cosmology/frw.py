"""             Friedman-Robertson-Walker Cosmology
    Utilities to compute various distance measures in the FRW cosmology
    Many of them are directly copied from Mark's cosmopy lib
"""
import pylab as pl
import numpy as np

import cosmology.cosparameter as cospar
import genscript.mymath as math
import genscript.units as units


class const: pass
const.c = 299792.458 # speed of light [km/s]


class Distance(object):
    '''
    Cosmological distances
    
    Fast version for arrays using linear interpolation of
    pre-calculated integral values for 1/E(z) and 1+w(z)

    all distances are in *bare* Mpc, for h^1Mpc units

    d -> d*h

    @@??? cp.w can be callable? How?
    '''
    def __init__(self,cosmoparam=None,zmax=20,_dz=1e-5):
        if cosmoparam==None:  cp = cospar.CosmoParams()
        else:                 cp = cosmoparam
        self.H0 = cp.H0
        self.h = self.H0/100.0
        self.m0 = cp.omec+ cp.omeb
        self.q0 = cp.omex
        self.k0 = 1 - self.m0 - self.q0    
        # hubble distance
        self.cH = const.c / self.H0 
        # w is a function of z
        w = cp.w
        self.w = w
        _z = pl.arange(0.0,zmax+_dz/2,_dz)
        if not callable(w):
            if w==-1:
                self._X = 0
            else:
                #X = lambda z: (1+w)*z/(1+z)
                #self._X = pl.array(map(X,_z))
                self._X = (1+w)*pl.log(1+_z)
            self.wfunc = lambda z: w(z)
        else:
            # trapezoid integral of 1+w(z) -> X
            o = pl.array(map(lambda z: (1+w(z))/(1+z), _z))
            self._X = math.trapezoidArray(o,_dz)
            self.wfunc = w
        # trapezoid integral 1/E(z)
        r = self.reE(_z)
        _iReE = math.trapezoidArray(r,_dz)
        # save
        self._z = _z
        self._dz = _dz
        self._rc = _iReE * self.cH
        
    def __str__(self):
        return 'Distance: %g %g %g %g %s' % (self.H0, self.m0, self.k0,
                                             self.q0, str(self.w))        

    def z2a(self, z):
        return 1./(1.+z)

    def a2z(self, a):
        return 1./a-1

    def E(self,z):
        "H(z) = H(0) E(z)"
        onez = 1+z
        onez2 = onez*onez
        onez3 = onez*onez2
        arg = self.m0*onez3 + self.k0*onez2 + self.q0*pl.exp(3*self._X)
        #arg = self.m0*onez3 + self.k0*onez2 + self.q0

        return pl.sqrt(arg)
    
    def H(self,z, unit_normalize=False):
        if unit_normalize==False:
	    fact=1.
	else:
	    fact=units.mega_parsec/1.e3

        return self.E(z)*self.H0*fact

    def mH(self, z, unit_normalize=False):
        return self.H(z, unit_normalize=unit_normalize)*self.z2a(z)

    def reE(self,z):
        return 1/self.E(z)
    
    def rc(self,z):
        """
        Comoving line of sight distance
        
        interpolate (linear for now)
        """
        i = pl.searchsorted(self._z,z)-1
        f = (z-self._z[i]) / self._dz
        d = self._rc[i] + f*(self._rc[i+1]-self._rc[i])        
        return pl.where(z==0,0,d)

    def rtc(self,z):
        "comoving transverse a.k.a. proper motion distance"
        Dh = self.cH
        Dc = self.rc(z)
        sqOmR = pl.sqrt(abs(self.k0))
        if self.k0 > 0:
            Dm = Dh / sqOmR * pl.sinh(sqOmR*Dc/Dh)
        elif self.k0 < 0:
            Dm = Dh / sqOmR * pl.sin (sqOmR*Dc/Dh);
        else:
            Dm = Dc
        return Dm

    def da(self,z):
        d = self.rtc(z)
        pl.divide(d,1+z, d)
        return d

    def dl(self,z):
        d = self.rtc(z)
        pl.multiply(d,1+z, d)
        return d

    def dm(self,z):
        return 5* pl.log10(self.dl(z)/(10e-6))
    
    def diffRc(self,z):
        return self.cH / self.E(z)

    def diffRtc(self,z):
        Dh = self.cH
        Dc = self.rc(z)
        sqOmR = pl.sqrt(abs(self.k0))
        diff = Dh/self.E(z)
        if self.k0 > 0:
            diff = diff * pl.cosh(sqOmR*Dc/Dh)
        elif self.k0 < 0:
            diff = diff * pl.cos (sqOmR*Dc/Dh);
        return diff

    def vc(self,z):
        Dm = self.rtc(z)
        return self.cH * Dm*Dm * self.reE(z)


    #aliases
    comovingLineOfSight = rc
    proper = rc
    diffProper = diffRc
    dac = rtc
    comovingTransverse = rtc
    diffComovingTransverse = diffRtc
    properMotion = rtc
    diffProperMotion = diffRtc
    angularDiameter=da
    luminosity=dl
    modulus=dm
    comovingVolumeElement=vc


