""">>> proj.py <<<

This module holds routines to project the power spectrum on the
sky.

Initially the following approximations will be used:
1) small angle (flat sky) approximation
2) flat spacetime
3) narrow enough redshift distribution, such that the change
of the growth factor can be neglected over the redshift interval


Current revision:
    ID:         $Id: proj.py 60 2008-02-28 20:46:31Z neyrinck $
    Date:       $Date: 2008-02-28 10:46:31 -1000 (Thu, 28 Feb 2008) $
    Revision:   $Revision: 60 $

(C) 2005 The CosmoPy Team (see Copyright for details)
"""
import unittest
import pylab as M
import frw
import pt
import utils

def calcdNdz2mass(z,z0=0.043,beta=1.825,lam=1.574):
    """
    unnormalized dN/dz for the 2Mass catalog,

    fit from Afshordi, Loh and Strauss (2004)
    default is for 12 < K_20 < 12.5

    other values:

    for 12.5 < K_20 < 13
    z0=0.054,beta=1.800, lambda=1.600
    
    for 13   < K_20 < 13.5
    z0=0.067,beta=1.765, lambda=1.636

    for 13.5 < K_20 < 14
    z0=0.084,beta=1.723, lambda=1.684 
    """
    return z**(beta*lam-1.0)*M.exp(-(z/z0)**beta)


def projectCl(lvec,P,D,dNdz,z,growthFac=None):
    """
    project C_l's given a power spectrum class P (Camb or
    BBKS) and Distance class D together

    arguments:
    lvec: vector of l values
    P: p.pk p.k contains the power spectrum, e.g. pt.Camb instance
    D: frw.Distance instance
    dNdz,z, growthFac: vectors suitable for trapezoid z-integration

    presently it crashes if z=0.0 is included, start from a small z value
    """
    lvec = M.asarray(lvec)
    dNdz2 = M.asarray(dNdz)**2
    z = M.asarray(z)
    da1 = 1./D.rtc(z)/D.h #comoving Da in h^1Mpc

    dNdz2vc = dNdz2/D.vc(z)/D.h**3 # comovin volume in (h^-1Mpc)^3
    #`use growth factor if given
    if growthFac:
        dNdz2vc = dNdz2vc*(growthFac**2)
    lk = M.log(P.k)
    pk = P.pk
    
##     return M.asarray([utils.trapz(utils.splineResample(pk,lk,
##                      M.log(l*da1))*dNdz2vc,z) for l in lvec])
    return M.asarray([utils.trapz(utils.interpolateLin(pk,lk,
                     M.log(l*da1))*dNdz2vc,z) for l in lvec])


def plotCl(ob=0.05):
    from pylab import loglog,show
    #z = M.arange(0.0005,0.5,0.0005)
    #dNdz = calcdNdz2mass(z)
    tamas = utils.readColumns("simit3n.dndz")
    z = tamas[0]
    dNdz = tamas[1]
    dNdz = dNdz/utils.trapz(dNdz,z)
    c = pt.Camb( hubble=70, omega_cdm=0.3-ob,
              omega_baryon=ob, w=-1.0,transfer_kmax=500,
                 use_physical='F',do_nonlinear=1)
    c.run()
    dist = frw.Distance(c.cp,zmax=1)
    l = M.arange(10.0,1000,1.0)
    cl = projectCl(l,c,dist,dNdz,z)
    loglog(l,l**2*cl)
    show()

class TestCl(unittest.TestCase):
    """
    test C_l projection
    """
    from pylab import loglog, show

    def runTest(self):
#        from pylab import loglog,show
        z = M.arange(0.0005,0.5,0.0005)
        dNdz = calcdNdz2mass(z)
        dNdz = dNdz/utils.trapz(dNdz,z)
        c = pt.Camb( hubble=70, omega_cdm=0.25,
                  omega_baryon=0.05, w=-1.0,transfer_kmax=500,
                     use_physical='F',do_nonlinear=1)
        c.run()
        dist = frw.Distance(c.cp,zmax=1)
        l = M.arange(10.0,1000,10.0)
        cl = projectCl(l,c,dist,dNdz,z)
#        loglog(l,l**2*cl)
#        show()
        print cl
        return True
        

if __name__=='__main__':
    unittest.main()
