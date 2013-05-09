import pylab as M
import numpy as N
import pt
import halo
import info
import copy

def hod():
    """
    Displays various galaxy power spectra w/ different HOD's.
    See param.py for explanations of hod params.

    camb needs to be in ../CAMB by default, but that can be changed in pt.py
    """

    c = pt.Camb(hubble = 70., ombh2 = 0.05*(0.7)**2, omch2 = 0.25*(0.7)**2) #may want to go out to transfer_kmax=100 for high accuracy
    c.run()
    pt.normalizePk(c,0.8) #sigma_8
    c.kextend(-10,60) #needed so that sigma(r) integral converges for a wide range of r

    #Sheth-Tormen
    h = halo.HaloModel(c,st_big_a = 0., st_little_a=0.707, stq = 0.3, k = 10**M.arange(-2,2.,0.1),massdivsperdex=5)

    h.pmm = halo.getHaloPknl(c,h)
    
    h.p.whichp='gg' # by default, h.p.whichp = 'mm', returning the matter power spectrum

    # by default, h.p.k_mmin_msun = 1e11, h.p,k_betas = 1
    mmins = 10.**M.arange(10.,12.01,0.2)
    for i in range(len(mmins)):
        h.p.k_mmin_msun = mmins[i]
        h.refreshHOD(c)
        h.pgg = halo.getHaloPknl(c,h)
        M.loglog(h.k,h.pgg,'y')

    M.loglog(h.k,h.pk,'k')
    M.loglog(h.k,h.pmm,'b')
    M.show()

    betas = M.arange(0.5,1.5,0.1)
    h.p.k_mmin_msun = 1e11
    for i in range(len(betas)):
        h.p.k_betas = betas[i]
        h.refreshHOD(c)
        h.pgg = halo.getHaloPknl(c,h)
        M.loglog(h.k,h.pgg,'y')

    M.loglog(h.k,h.pk,'k')
    M.loglog(h.k,h.pmm,'b')
    M.show()

def redshift():
    """
    Evolution with redshift of matter power spectrum
    """
    zs = M.arange(0.,5.,2.)

    for z in zs:
        print z
        c = pt.Camb(hubble = 70., ombh2 = 0.05*(0.7)**2, omch2 = 0.25*(0.7)**2,transfer_redshift = [z])
        c.run()
        ps = pt.PowerSpectrum(c.cp)
        c.kextend(-10,60) #To ensure accurate sigma(r) -- if it doesn't, a warning will ensue
        pt.normalizePk(c,0.8*ps.d1(z)/ps.d1(0.)) #sigma_8 at redshift z
        
        #Sheth-Tormen
        h = halo.HaloModel(c,st_big_a = 0., st_little_a=0.707, stq = 0.3, k = 10**M.arange(-2,2.01,0.2),massdivsperdex=5)

        h.pmm = halo.getHaloPknl(c,h)
        M.loglog(h.k, h.pmm, label='z='+str(z))
        M.loglog(h.k, h.pk,'k:',label='linear')

        cp_halofit = c.cp
        cp_halofit['do_nonlinear'] = 1 # Halofit (Smith et al) fit
        chalofit = pt.Camb(cambParam=cp_halofit)
        chalofit.run()
        wheretoplot = N.where(chalofit.k > 1e-2)[0]
        M.loglog(chalofit.k[wheretoplot[::10]],chalofit.pk[wheretoplot[::10]],'--',label='halofit')

    M.legend()
    M.show()

def getInfoCurve():
    """
    Various functions to calculate example parameter error bars as in
    Neyrinck & Szapudi 2007, MNRAS 375, L51
    """

    c = pt.Camb(hubble = 70., ombh2 = 0.05*(0.7)**2, omch2 = 0.25*(0.7)**2)
    c.run()
    c.kextend(-10,60) # necessary to make sigma(m) integral converge well.
    pt.normalizePk(c,0.8) #sigma_8

    outputdir = 'example/'
    #Sheth-Tormen
    h = halo.HaloModel(c,st_big_a = 0., st_little_a=0.707, stq = 0.3, k = 10.**M.arange(-2,1.01,0.25),massdivsperdex=5)
    #For final calculations, should use more massdivsperdex, e.g. 20 (maybe 10 is ok)
    #also, k is really coarse, as you'll see if you run this.

    # get covariance matrix from halo-model trispectrum (saves it in the 'prefix' directory)
    # it also automatically runs halo.getHaloPknl
    halo.getHaloCov(outputdir,c,h)

    # power spectrum at h.k (range of k at which halo model quantities are evaluated)
    M.loglog(h.k,h.pnl)
    M.show()

    # get derivs wrt ln A, tilt
    h.dloga = halo.getdlogPnldCosmoParam(c,h,'scalar_amp',linlog='log')
    h.dtilt = halo.getdlogPnldCosmoParam(c,h,'scalar_spectral_index',linlog='lin')
    M.loglog(h.k,h.dloga**2,label='ln A')
    M.loglog(h.k,h.dtilt**2,label='tilt')
    M.legend()
    M.show()
    M.save(outputdir+'dlogpnldloga.dat',M.transpose([h.k,h.dloga]),fmt='%6.5e')
    M.save(outputdir+'dlogpnldtilt.dat',M.transpose([h.k,h.dtilt]),fmt='%6.5e')

    # get parameter covariance matrix (just a function of k, since there's only one variable)
    k, covmat = info.getParamCovMat(outputdir,dlogfilenames=['dlogpnldloga.dat','dlogpnldtilt.dat'])

    # plot the unmarginalized error bars in ln A and the tilt,
    # if the matter power spectrum is analyzed from k= k[0] to k.

    M.loglog(k, M.sqrt(covmat[0,0,:]),label='ln A')
    M.loglog(k, M.sqrt(covmat[1,1,:]),label='tilt')

    M.legend()
    M.show()
