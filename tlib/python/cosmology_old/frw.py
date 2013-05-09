"""Friedman-Robertson-Walker Cosmology

Utilities to compute various distance measures in the FRW cosmology

Current revision:
    ID:         $Id: frw.py 59 2008-02-28 20:36:04Z neyrinck $
    Date:       $Date: 2008-02-28 10:36:04 -1000 (Thu, 28 Feb 2008) $
    Revision:   $Revision: 59 $
"""
REVISION = '$Revision: 59 $'

import unittest
import pylab as M
import numpy as N

from param import CosmoParams
import utils


class const: pass
const.c = 299792.458 # speed of light [km/s]


class Distance:
    '''
    Cosmological distances
    
    Fast version for arrays using linear interpolation of
    pre-calculated integral values for 1/E(z) and 1+w(z)

    all distances are in *bare* Mpc, for h^1Mpc units

    d -> d*h
    '''
    def __init__(self,cosmoparam=None,zmax=20,_dz=1e-5):
        if cosmoparam==None:  cp = CosmoParams()
        else:                 cp = cosmoparam
        self.H0 = cp.Hubble
        self.h = self.H0/100.0
        self.m0 = cp.Omega_CDM + cp.Omega_Baryon
        self.q0 = cp.Omega_Lambda
        self.k0 = 1 - self.m0 - self.q0    
        # hubble distance
        self.cH = const.c / self.H0 
        # w is a function of z
        w = cp.w
        self.w = w
        _z = M.arange(0.0,zmax+_dz/2,_dz)
        if not callable(w):
            if w==-1:
                self._X = 0
            else:
                #X = lambda z: (1+w)*z/(1+z)
                #self._X = M.array(map(X,_z))
                self._X = (1+w)*M.log(1+_z)
            self.wfunc = lambda z: w(z)
        else:
            # trapezoid integral of 1+w(z) -> X
            o = M.array(map(lambda z: (1+w(z))/(1+z), _z))
            self._X = utils.trapezoidArray(o,_dz)
            self.wfunc = w
        # trapezoid integral 1/E(z)
        r = self.reE(_z)
        _iReE = utils.trapezoidArray(r,_dz)
        # save
        self._z = _z
        self._dz = _dz
        self._rc = _iReE * self.cH
        
    def __str__(self):
        return 'Distance: %g %g %g %g %s' % (self.H0, self.m0, self.k0,
                                             self.q0, str(self.w))        
    def E(self,z):
        "H(z) = H(0) E(z)"
        onez = 1+z
        onez2 = onez*onez
        onez3 = onez*onez2
        arg = self.m0*onez3 + self.k0*onez2 + self.q0*M.exp(3*self._X)
        return M.sqrt(arg)
    
    def H(self,z):
        return self.E(z)*self.H0

    def reE(self,z):
        return 1/self.E(z)
    
    def rc(self,z):
        """
        Comoving line of sight distance
        
        interpolate (linear for now)
        """
        i = M.searchsorted(self._z,z)-1
        f = (z-self._z[i]) / self._dz
        d = self._rc[i] + f*(self._rc[i+1]-self._rc[i])        
        return M.where(z==0,0,d)

    def rtc(self,z):
        "comoving transverse a.k.a. proper motion distance"
        Dh = self.cH
        Dc = self.rc(z)
        sqOmR = M.sqrt(abs(self.k0))
        if self.k0 > 0:
            Dm = Dh / sqOmR * M.sinh(sqOmR*Dc/Dh)
        elif self.k0 < 0:
            Dm = Dh / sqOmR * M.sin (sqOmR*Dc/Dh);
        else:
            Dm = Dc
        return Dm

    def da(self,z):
        d = self.rtc(z)
        M.divide(d,1+z, d)
        return d

    def dl(self,z):
        d = self.rtc(z)
        M.multiply(d,1+z, d)
        return d

    def dm(self,z):
        return 5* M.log10(self.dl(z)/(10e-6))
    
    def diffRc(self,z):
        return self.cH / self.E(z)

    def diffRtc(self,z):
        Dh = self.cH
        Dc = self.rc(z)
        sqOmR = M.sqrt(abs(self.k0))
        diff = Dh/self.E(z)
        if self.k0 > 0:
            diff = diff * M.cosh(sqOmR*Dc/Dh)
        elif self.k0 < 0:
            diff = diff * M.cos (sqOmR*Dc/Dh);
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

class DistanceTest(unittest.TestCase):
    def runTest(self):
        z = M.arange(0,5,0.1)
        #cp = CosmoParams( hubble=70, omega_cdm=0.25,
        #                  omega_baryon=0.05, w=lambda z: -1)
        cp = CosmoParams( hubble=70, omega_cdm=0.25,
                          omega_baryon=0.05, w=-1.0)
        dist = Distance(cp)
        for s in z:
            d = dist.comovingLineOfSight(s)
            print s, d

        
if __name__=='__main__':
    unittest.main()
