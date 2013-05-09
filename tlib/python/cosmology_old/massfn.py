""">>> massfn.py <<<

This module holds stuff to calculate mass functions based on given cosmologies.

Current revision:
    ID:         $Id: param.py 17 2005-10-31 18:18 neyrinck $
    Date:       $Date: Oct 31 18:22:52 HST 2005 $
    Revision:   $Revision: 25 $

(C) 2005 The CosmoPy Team (see Copyright for details)
"""

import pylab as M
import numpy as N
import utils
import pt
import halo

import os

def msun2chimp(m_msun,c):
    """
    Converts a mass in solar masses to chimps.
    A chimp is the matter contained in a Cubic H-Inv MegaParsec,

    Arguably, omega_cdm could be changed to omega_matter.
    """
    #omega_matter = c.cp.omega_baryon + c.cp.omega_cdm
    return m_msun*(c.cp.hubble/100.)/(2.7755e11*c.cp.omega_cdm)

def msunh2chimp(m_msun,c):
    """
    Converts a mass in solar masses/h to chimps.
    A chimp is the matter contained in a Cubic H-Inv MegaParsec,

    Arguably, omega_cdm could be changed to omega_matter.
    """
    #omega_matter = c.cp.omega_baryon + c.cp.omega_cdm
    return m_msun/(2.7755e11*c.cp.omega_cdm)

def chimpRadius(m):
    """
    Gets the radius in Mpc/h of a sphere enclosing mass m, where m is in chimps
    """
    return (3.*m/(4.*M.pi))**(1./3.)

def getMassFunction(h,c):
    """
    Get n(m,z) from a halo model instance for which nu(m) has already been calculated,
    and a Camb instance.
    """
    
    nuprime2 = h.p.st_little_a * h.nu**2

    nufnu = 2.*(1.+ 1./nuprime2**h.p.stq)*M.sqrt(nuprime2/(2.*M.pi))* \
                M.exp(-nuprime2/2.) # hold off on normalization

    dlognu = h.m*0.

    for i in range(len(h.nu)):
        dlognu[i] = 0.5*M.log(h.nu_pad[i+2]/h.nu_pad[i])

    nmz_unnorm = (dlognu/h.dlogm)*nufnu/h.m**2

    w = N.where(nmz_unnorm < 1.7e308)[0]
    lw = len(w)
    if lw < len(nmz_unnorm):
        print "Warning! the mass function's blowing up!"

    h.nmz = nmz_unnorm*1.
    totaln = halo.generalIntOverMassFn(1,1,1.,h,whichp='mm')
    if h.p.st_big_a == 0.:
        h.nmz /= totaln
    else:
        h.nmz *= h.p.st_big_a

    print 'Normalization const (integrated):',1./totaln
    # if this isn't close to what you expect (~0.322 for Sheth-Tormen, 0.5 for Press-Schechter),
    # you need to expand the mass integration range, the mass bins per dex, or extrapolate c.pk.
    if h.p.st_big_a != 0.:
        print 'Used:',h.p.st_big_a

def poissonize(h,tol):
    """
    Monte-Carlo sample the mass function, using a Poisson distribution in each bin.
    h = halo model instance
    tol = threshold beyond which we don't bother with the Monte Carlism.  If there's a high number of
          haloes in a bin (>~1000?), Poissonizing will make little difference
    """
    
    poissonnmz = h.nmz*0.
    arraytopoissonize = h.nmz*(h.volume*h.m*h.dlogm)
    #pylab.clf()
    #pylab.loglog(h.m,1.+arraytopoissonize)
    #print arraytopoissonize
    tolm2 = tol**(-2)
    for j in range(len(h.nmz)):
        if (arraytopoissonize[j] < tolm2) and (arraytopoissonize[j]>1e-45):
            arraytopoissonize[j] = N.random.poisson(arraytopoissonize[j])
            #poissonnmz[j] = jimmy
            #print arraytopoissonize[j],poissonnmz[j]
        elif (arraytopoissonize[j] <= 1e-45):
            arraytopoissonize[j] = 0.

    #pylab.loglog(h.m,1.+poissonnmz)
    ans = arraytopoissonize/(h.volume*h.m*h.dlogm)

    #pylab.savefig('poissonize.ps')

    return(ans)
