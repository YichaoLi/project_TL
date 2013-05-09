""">>> Halo.py <<<

This module holds routines to calculate things using the halo model.

Current revision:
    ID:         $ID: $
    Date:       $Date: 2008-12-14 16:07:39 -1000 (Sun, 14 Dec 2008) $
    Revision:   $Revision: 63 $

(C) 2007 The CosmoPy Team (see Copyright for details)
"""
import pylab as M
import numpy as N
import utils
import massfn
import pt
import param
import os
import copy

class HaloModel(utils.SpecialFunctions):
    """
    class to initialize halo model
    """
    def __init__(self, c, haloModelParam=None,**kws):
        """
        kws is a set of arguments to change (see params in param.py)
        e.g. HaloModel(stA=0.) etc
        """
        if haloModelParam==None:
            self.p = param.HaloModelParams(**kws)
        else:
            self.p = haloModelParam
        self.run(c)

    def run(self,c): # c is a camb object
        self.z = c.cp.transfer_redshift[0]
        self.dlog10m = 1./self.p.massdivsperdex
        self.dlogm = self.dlog10m * M.log(10.)
        if self.p.lomass == 0.:
            self.p.lomass = (4.*M.pi/3.)*c.cp.omega_cdm*(2.*M.pi/c.k[-1])**3
        self.m_pad = self.p.lomass*\
                M.exp(M.arange(-self.dlogm,M.log(self.p.himass/self.p.lomass)+self.dlogm*2, self.dlogm))
        self.m = self.m_pad[1:-2]
        if self.p.integrateto == 0.:
            self.integratetoindex = len(self.m)
        else:
            self.integratetoindex = N.where(self.m < self.p.integrateto)[0][-1]

        ps = pt.PowerSpectrum(c.cp)
        if self.p.delta == 0.:
            self.delta = ps.deltac(self.z)
        else:
            self.delta = self.p.delta
        print 'delta = ',self.delta
        self.deltaz0 = ps.deltac(0.) #the ~1.67 factor, without redshift correction
        sigma_pad = M.sqrt(pt.sigma2fromPk(c,massfn.chimpRadius(self.m_pad))) * ps.d1(0.)/ps.d1(self.z)

        # calculate sigma with ps normalized to z=0
        self.nu_pad = self.delta/sigma_pad
        self.sigma = sigma_pad[1:-2]
        self.nu = self.nu_pad[1:-2]
        self.b = getHaloBiasCoef(self)
        massfn.getMassFunction(self,c) # sets h.nmz

        bint0 = generalIntOverMassFn(0,1,1.,self,whichp='mm')
        bint1 = generalIntOverMassFn(1,1,1.,self,whichp='mm')
        bint2 = generalIntOverMassFn(2,1,1.,self,whichp='mm')
        #test 1st,2nd-order bias.  bint0,1 should be ~1, bint2 should be ~0
        if (M.fabs(M.log(bint0)) > 0.01) + (M.fabs(M.log(bint1)) > 0.01) + (M.fabs(bint2) > 0.01):
            print 'Warning! (bint0,1,2) = ',bint0,bint1,bint2
            print 'Should be about (1,1,0).  If worried, increase range of mass integration,'
            print 'mass bins per dex, or extrapolate c.pk (which affects sigma(m))'
        self.mstar = mStar(self.m,self.nu)
        self.logcbar = getLogcbar(self)
        #self.k = c.k[self.p.startki::self.p.strideki]

        self.k = self.p.k*1.
        
        if self.p.use_sea_2h: #Use nonlinear power spectrum for 2h term
            # not recently tested
            cp_cnl = copy.copy(c.cp)
            cp_cnl['do_nonlinear'] = 1
            cnl = pt.Camb(cambParam = cp_cnl)
            cnl.run()
            self.pk = M.exp(utils.splineIntLinExt(M.log(cnl.pk),M.log(cnl.k),M.log(self.k)))

        else: #Use linear power spectrum for 2h term
            self.pk = c.pkInterp(self.k)

        if self.p.dofthp == 1:
            self.fthp = deltah(self)

        # Convert Kravtsov HOD quantities to internal units (chimps, etc.)
        self.k_mmin = massfn.msun2chimp(self.p.k_mmin_msun,c)
        self.k_m1 = self.k_mmin * self.p.k_m1overmmin
        # for the non-Kravtsov HOD, for some reason, we use units of msun/h
        self.mcut = massfn.msunh2chimp(self.p.mcut_msunh,c)
        self.m11 = massfn.msunh2chimp(self.p.m11_msunh,c)
        self.m13 = massfn.msunh2chimp(self.p.m13_msunh,c)
        self.logsqrtfact = M.log(M.sqrt(self.m13/self.m11))

        self.ngalbar = generalIntOverMassFn(0,1, 1.,self,whichp='gg')

    def refreshHOD(self,c):
        """
        Do this if you change a HOD parameter
        """
        self.k_mmin = massfn.msun2chimp(self.p.k_mmin_msun,c)
        self.k_m1 = self.k_mmin * self.p.k_m1overmmin
        self.mcut = massfn.msunh2chimp(self.p.mcut_msunh,c)
        self.m11 = massfn.msunh2chimp(self.p.m11_msunh,c)
        self.m13 = massfn.msunh2chimp(self.p.m13_msunh,c)
        self.logsqrtfact = M.log(M.sqrt(self.m13/self.m11))

        self.ngalbar = generalIntOverMassFn(0,1, 1.,self,whichp = 'gg')
        
def hodMoment(n,h):
    """
    Returns the nth moment of the Halo Occupation Distribution
    <Ngal(Ngal-1) ... (Ngal-n)|m>.
    Notation initially as in Scoccimarro et al (2001)

    The only really tested function is the other function hodNsMoment, the Satellite HOD of K04.
    """
    if h.p.hodkrav == 1:
        ans = 1.*M.zeros(len(h.m))
        h.wk_mmin = N.where(h.m > h.k_mmin)[0][0]
        h.wk_m1 = N.where(h.m > h.k_m1)[0][0]
        n_s = (h.m[h.wk_mmin:]/h.k_m1)**h.p.k_betas
        ans[h.wk_mmin:] = (n_s**(n-1))*(n_s + n)

        return ans
        
    elif h.p.hodkrav == 2:
        print "Slightly more complicated/accurate Krav. not implemented yet."

    else:
        """
        Binomial HOD not fully implemented
        """

        alpha2 = 1.*M.zeros(len(h.m))
        h.w11 = N.where(h.m > h.m11)[0][0]
        h.wcut = N.where(h.m >= h.mcut)[0][0]
        h.w13 = N.where(h.m >= h.m13)[0][0]
        # don't need to where it each time

        alpha2[h.w11:h.w13] = M.log(M.sqrt(h.m[h.w11:h.w13]/h.m11))/h.logsqrtfact
        alpha2[h.w13:] += 1.
    
        momnom = 1.*M.ones(len(h.m))
        # Moment w/o m, i.e. <Ngal(Ngal-1)...(Ngal-n-1)>
        # binomial approx. of Scoccimarro et al (2001)
        for j in range(1,n):
            momnom *= j*alpha2 - (j - 1.)

        ngalm = 1.*M.zeros(len(h.m))
        ngalm[h.w11:h.wcut] += h.p.n0
        ngalm[h.wcut:] = h.p.n0*(h.m[h.wcut:]/h.mcut)**h.p.alphaexp

        #M.loglog(h.m[h.w11:],ans[h.w11:],'k')
        #M.loglog(h.m[h.w11:],momnom[h.w11:]*ngalm[h.w11:]**n,'b')
        #M.show()

        return ans
        #return momnom * ngalm**n

def hodNsMoment(n,h):
    """
    Returns the nth moment of the Satellite HOD
    <Ns(Ngal-1) ... (Ngal-n)|m>.
    Notation as in Kravtsov et al (2004)
    """
    if h.p.hodkrav == 1:
        ans = 1.*M.zeros(len(h.m))
        h.wk_mmin = N.where(h.m > h.k_mmin)[0][0]
        h.wk_m1 = N.where(h.m > h.k_m1)[0][0]
        n_s = (h.m[h.wk_mmin:]/h.k_m1)**h.p.k_betas
        ans[h.wk_mmin:] = n_s**n

    else:
        print "Satellite HOD only for Kravtsov-style satellite HOD."

    return ans


def getHaloBiasCoef(h):
    """
    Halo bias parameters (Mo et al (1997), Cooray & Sheth, sect. 3.4)
    """
    a1 = 1.
    a2 = -17./21.
    a3 = 341./567.
    a4 = -55805./130977.

    delta2 = h.delta*h.delta
    sigma2 = h.sigma*h.sigma
    nu2 = delta2/sigma2

    # notation as in Cooray & Sheth, sect. 3.4
    # ST a -> CS q
    # ST q -> CS p
    # ST A -> CS A
    cs_q = h.p.st_little_a
    cs_p = h.p.stq

    epsilon1 = (cs_q*nu2 - 1.)/h.deltaz0
    epsilon2 = cs_q*nu2*(cs_q*nu2-3.)/h.deltaz0**2
    e1 = 2.*cs_p/h.deltaz0/(1.+(cs_q*nu2)**cs_p)
    e2 = e1 * ((1.+2.*cs_p)/h.deltaz0 + 2.*epsilon1)
    # use h.deltaz0, not h.delta here (bug fix from Sebastien Heinis)

    b = M.reshape(1.*M.ones(3*len(h.sigma)),(3,len(h.sigma)))
    # already took care of b[0,:]
    b[1,:] = 1. + epsilon1 + e1
    b[2,:] = 2.*(1.+a2)*(epsilon1+e1) + epsilon2 + e2

    # This is what's in Mo et all (1997); this is just for the Press-Schechter case
    #b[1,:] = (nu2-1.)/h.delta + 1.
    #b[2,:] = 2.*(1.+a2)*(nu2-1)/h.delta + (nu2-3.)/sigma2
    #b[3,:] =  6.*(a2+a3)*(nu2-1)/h.delta + 3.*(1.+2.*a2)*(nu2-3.)/sigma2 + \
    #        (nu2**2-6.*nu2+3.)/(h.delta*sigma2)
    # b[4,:] = 24.*(a3+a4)*(nu2-1.)/h.delta + \
    #           12.*(a2**2 + 2.*(a2+a3))*(nu2-3.)/sigma2 + \
    #           4.*(1.+3.*a2)*(nu2**2-6*nu2+3)/(h.delta*sigma2) + \
    #           (nu2**2-10.*nu2+15.)/sigma2
    
    return b

def mStar(m,nu):
    """
    Returns M* based on an array of nu(M)'s.
    M* is defined to be the mass at which nu(M) = 1.
    Used for concentration distribution.
    """
    closest = N.where(nu < 1.)[0][-1] #nu increases with M

    logmstar = M.log(m[closest]) + M.log(m[closest+1]/m[closest])/M.log(nu[closest+1]/nu[closest])*\
               M.fabs(M.log(nu[closest]))
    return M.exp(logmstar)

def haloDensProf(r,rs,h):
    """
    returns unnormalized halo density profile
    rs: scale radius
    r: radius at which it's evaluated
    """
    return 1./((r/rs)**h.p.dpalpha * (1. + r/rs)**h.p.dpbetanfw)

def deltah(h):
    """
    Gets delta_h (notation as in Cooray & Hu 2001),
    the Fourier-transformed halo profile
    """

    if h.p.nfw == 'y':
        # Use cisi
        ans = M.outer(h.k*0.,h.m*0.)
        cbar = M.exp(h.logcbar)
        rvir = (3./(4.*M.pi)*h.m/200.)**(1./3.)
        rs = rvir/cbar
        lf = M.log(1+cbar)-cbar/(1+cbar)
        for ki in range(len(h.k)):
            ci1c, si1c = utils.cisiarr((1+cbar)*h.k[ki]*rs)
            ci0, si0 = utils.cisiarr(h.k[ki]*rs)
            ans[ki,:] = (M.sin(h.k[ki]*rs)*(si1c-si0) - \
                         M.sin(cbar*h.k[ki]*rs)/((1+cbar)*h.k[ki]*rs) + \
                         M.cos(h.k[ki]*rs)*(ci1c-ci0))/lf
            
        else:
            g = utils.HGQ(5) # 5-pt Hermite-Gauss quadrature class
            ans = 1.*M.zeros((len(h.k),len(h.m),5))
            
            for mi in range(len(h.m)):
                logc = h.logcbar[mi] + g.abscissas*h.p.sigmalogc #array of c's
                c = M.exp(logc)
                rvir = ((3./4.)*h.m[mi]/(200.*M.pi))**(1./3.)
                rs = rvir/c
                lf = M.log(1+c)-c/(1+c)
                for ki in range(len(h.k)):
                    ci1c, si1c = utils.cisiarr((1+c)*h.k[ki]*rs)
                    ci0, si0 = utils.cisiarr(h.k[ki]*rs)
                    ans[ki,mi,:] = (M.sin(h.k[ki]*rs)*(si1c-si0) - \
                                    M.sin(c*h.k[ki]*rs)/((1+c)*h.k[ki]*rs) + \
                                    M.cos(h.k[ki]*rs)*(ci1c-ci0))/lf
                    
    return ans

def getLogcbar(h):
    """
    Returns log(cbar), where cbar is the mean concentration parameter
    e.g. Bullock et al: cbar = 9/(1+z) (m/mstar)^-0.13
    """
    return M.log(h.p.cbarcoef/(1.+h.z)*(h.m/h.mstar)**h.p.cbarslope)

def generalIntOverMassFn(beta, mExponent, extra_integrand, h, whichp = 0):
    """       beta
    Computes I (k) (as in Cooray & Hu (2001)), integrating over the halo mass function.
              mu

    beta: order of bias coefficient (b_beta in the integrand)
    mExponent: additional (over the power of m from the logarithmic integration)
                mass exponent (put mu here)
    extra_integrand: e.g. halo density profiles
    h: halo model instance

    e.g. i11[k] = integralOverMassFn(1,1, g.int(h.fthp[k,:,:]),h)
     where the two : indices in h.fthp are for mass and concentration parameter;
     g.int integrates over concentration parameter, so g.int(h.fthp[k,:,:]) is an array
     for each mass.
    """
    if whichp == 0:
        whichp == h.p.whichp
    
    if whichp == 'gg':
        integrand = h.b[beta,:]*h.m*hodMoment(mExponent,h)*h.nmz * extra_integrand

    elif whichp == 'mm':
        integrand = h.b[beta,:]*h.m*h.m**(mExponent)*h.nmz * extra_integrand
        
    ans = utils.simpsonRule(integrand, h.dlogm)
    return ans

def intOverMassFn(beta, mExponent, karrayunnumed, h, whichp = 0, plot = 0, show = 0,colstr = ''):
    """
    Get an integral of, e.g., halo density profiles over the mass function.
    More automatic than generalIntOverMassFn.
    
              beta
    Computes I (k) (as in Cooray & Hu (2001)), integrating over the halo mass function.
              mu
    beta: order of bias coefficient (b_beta in the integrand)
    mExponent: additional (over the power of m from the logarithmic integration)
                mass exponent (put mu here).
    karrayunnumed: Array of k indices to evaluate Fourier-transformed halo density profile
    h: halo model instance
    whichp: for now, 'gg' -> galaxy statistics.
                     'mm' -> matter stats.   Default: h.p.whichp

    """
    if whichp == 0:
        whichp = h.p.whichp

    karray = M.array(karrayunnumed)

    if len(karray) != mExponent:
        print "Warning! Different number of k indices than mExponent!"
    
    g = utils.HGQ(5) # 5-pt Hermite-Gauss quadrature class

    uM = h.fthp[0,:,:]*0. + 1.
    g_integrand = h.fthp[0,:,:]*0. + 1.
            
    if whichp == 'gg':
        """
        Right now, we've only implemented the Kravtsov-style 'satellite' halo occupation distribution,
        i.e. the HOD excluding the central galaxy.  For example, in Smith, Watts & Sheth 2006, sect. 5.4, the
        halo integrands simplify in the case of a satellite HOD in eq. 72 (one-halo P_gg) to
        <Ns(Ns-1)> |u(k)|^2 + 2 <Ns> u(k).  This second term includes the contribution of the central galaxy.
        In general, using SWS's notation, the n-galaxy convolution kernel (e.g. eq. 76) is

        W^ng(k_0, ..., k_n) = <Ns(Ns-1)...(Ns-n)> u(k_0)...u(k_n) + Sum_i <Ns...(Ns-(n-1))> u(k_0)...u(k_n)/u(k_i).
        """                                                   
        
        hodNsMomentM = M.outer(hodNsMoment(mExponent,h),1.*M.ones(5))
        hodNsMomentMminus1 = M.outer(hodNsMoment(mExponent-1,h),1.*M.ones(5))
        
        # Leading term in mExponent; fthp(k)^mExponent
        for ki in karray:
            uM *= h.fthp[ki,:,:]

        g_integrand = hodNsMomentM*uM*1.
        
        #mExponent * fthp^(mExponent-1) term
        for ki in karray:
            g_integrand += hodNsMomentMminus1 * uM/h.fthp[ki,:,:]
            
        extra_integrand = g.int(g_integrand)
                

    elif whichp == 'mm':
        for ki in karray:
            g_integrand *= h.fthp[ki,:,:]
        
        extra_integrand = h.m**(mExponent)* g.int(g_integrand)

    integrand = h.b[beta,:]*h.m*h.nmz * extra_integrand

    if plot == 1:
        M.loglog(h.m[::10],integrand[::10],colstr)
        #M.semilogx(h.m,integrand)

    if show == 1:
        M.show()
        
    ans = utils.simpsonRule(integrand[:h.integratetoindex], h.dlogm)
    if whichp == 'gg':
        ans /= h.ngalbar**mExponent
        
    return ans

def getHaloPknl(c,h):
    """
    Returns halo model power spectrum (if h.p.whichp == 'gg', galaxy; if h.p.whichp == 'mm', matter)
    """
    g = utils.HGQ(5) # 5-pt Hermite-Gauss quadrature class

    i11 = 0.*h.k
    i02 = 0.*h.k
    extra2h = 0.*h.k
    for k1 in range(len(h.k)):
        i11[k1] = intOverMassFn(1,1, [k1],h)
        i02[k1] = intOverMassFn(0,2, [k1,k1],h)

    #print 'i''s:',i11[0],i02[0]
    #if h.largescalenormal == 1:
    #    i11 /= i11[0]

    p1h = i02
    #one-halo term
    p2h = h.pk*i11**2
    #two-halo term
    #df = h.k**3/(2.*M.pi**2)
    #df = 1.

    #M.loglog(h.k,p1h*df,'r',linewidth=0.5)
    #M.loglog(h.k,p2h*df,'y',linewidth=0.5)
    #M.loglog(h.k,h.pk*df,'r',linewidth=0.5)
    #M.loglog(h.k,(p1h+p2h)*df,'k',linewidth=0.25)
   
    return p1h + p2h

def getdlogPnldCosmoParam(c,h,paramstring,epsilon = 2.**(-6),linlog='log',sig8 = 0.,usepnl=True, usehalofit=False):
    """
    Returns dlog(Pnl(k))/d(parameter)

    linlog: if 'lin' (or anything but 'log'), do dlogP/dparam; if 'log', do dlogP/dlogparam

    epsilon = interval in param or log param

    Warning; this keeps the scalar amplitude A constant, not sigma_8.
    To keep sigma_8 constant, set sig8 to be nonzero.
    To get linear power spectrum dlog's, set usepnl = False
    Should use Smith et al. HALOFIT power spectrum if usepnl=False, usehalofit=True (though this hasn't been tested).
    """

    cp = copy.copy(c.cp)

    if (usehalofit):
        cp['do_nonlinear'] = 1

    if linlog == 'log':
        multfact = M.exp(epsilon)
        divfact = M.exp(-epsilon)
    else:
        multfact = 1. + epsilon
        divfact = 1. - epsilon

    if (paramstring == 'barf'):
        ommh2 = cp['ombh2']+cp['omch2']
        cp['ombh2'] *= multfact
        cp['omch2'] = ommh2 - cp['ombh2']
    elif isinstance(cp[paramstring],float):
        cp[paramstring] *= multfact
    else:
        for i in range(len(cp[paramstring])):
            cp[paramstring][i] *= multfact

    cplus = pt.Camb(cambParam = cp)
    cplus.run()

    if (paramstring == 'barf'):
        cp['ombh2'] *= divfact/multfact #since we're starting from the cplus params
        cp['omch2'] = ommh2 - cp['ombh2']
    elif isinstance(cp[paramstring],float):
        cp[paramstring] *= multfact
    else:
        for i in range(len(cp[paramstring])):
            cp[paramstring][i] *= divfact/multfact #since we're starting from the cplus params

    cminus = pt.Camb(cambParam = cp)
    cminus.run()

    cplus.kextend(-10,60)
    cminus.kextend(-10,60)

    if (sig8 > 0.):
        pt.normalizePk(cplus,sig8)
        pt.normalizePk(cminus,sig8)

    if (usepnl):
        hplus = HaloModel(cplus,haloModelParam = h.p)
        hplus.pknl = getHaloPknl(cplus,hplus)
        hminus = HaloModel(cminus,haloModelParam = h.p)
        hminus.pknl = getHaloPknl(cminus,hminus)
    else:
        hplus = copy.deepcopy(h)
        hplus.pknl = cplus.pkInterp(h.k)
        hminus = copy.deepcopy(h)
        hminus.pknl = cminus.pkInterp(h.k)

    dlog = M.log(hplus.pknl/hminus.pknl)/(2.*epsilon)    

    #M.loglog(h.k,hplus.pknl)
    #M.loglog(h.k,hminus.pknl)
    #M.show()

    return dlog

def getHaloTrispec(c,h, startki = 0, endki = 0, strideki = 1, adder = 0., highprecisionthresh = 0.03):
    """
    Calculates the halo model trispectrum T(k1,-k1,k2,-k2),
    as in Cooray & Hu (2001)

    adder: something that will be added later, e.g. gaussian
    part of the covariance.  In there for precision checks

    Should work for galaxies, too, but shot noise isn't implemented.  Also, for some HOD parameters,
    may expose PT trispectrum on scales where it may not work
    """

    #p = c.cp

    tribi = pt.TriBiSpectrumFromCamb(c)
    g = utils.HGQ(5)

    # Could change a bunch of these to scalars

    i11 = 0.*h.k
    i21 = 0.*h.k
    i02 = 0.*h.k
    pk1plusk3_perp = 0.*h.k
    t3hB_perp = 0.*h.k
    t4hT_perp = 0.*h.k
    dsq4hnoT = 0.*h.k
    dsq4hshouldbe = 0.*h.k
    dsq1h = 0.*h.k
    dsq2h = 0.*h.k
    dsq2h31 = 0.*h.k
    dsq3h = 0.*h.k
    dsq4h = 0.*h.k
    qsq = 0.*h.k
    i04 = M.outer(0.*h.k,0.*h.k)
    i12 = M.outer(0.*h.k,0.*h.k)
    i13_112 = M.outer(0.*h.k,0.*h.k)
    i13_122 = M.outer(0.*h.k,0.*h.k)
    i22 = M.outer(0.*h.k,0.*h.k)
    i114 = M.outer(0.*h.k,0.*h.k)
    i1114 = M.outer(0.*h.k,0.*h.k)
    
    
    pk1plusk3 = M.outer(0.*h.k,0.*h.k)
    t2h31 = M.outer(0.*h.k,0.*h.k)
    t3hnoB = M.outer(0.*h.k,0.*h.k)
    t3hB = M.outer(0.*h.k,0.*h.k)
    t4hnoT = M.outer(0.*h.k,0.*h.k)
    t4hT = M.outer(0.*h.k,0.*h.k)
    b10 = M.outer(0.*h.k,0.*h.k)
    t10 = M.outer(0.*h.k,0.*h.k)

    for k1 in range(len(h.k)):
        i11[k1] = intOverMassFn(1,1, [k1],h)
        i21[k1] = intOverMassFn(2,1, [k1],h)
        i02[k1] = intOverMassFn(0,2, [k1,k1],h)

    if endki == 0:
        endki = len(h.k)

    for k1 in range(startki, endki, strideki):
        for k3 in range(k1,len(h.k)):
        #for k3 in range(k1,endki,strideki):
            i04[k1,k3] = intOverMassFn(0,4, [k1,k1,k3,k3], h)
            i13_112[k1,k3] = intOverMassFn(1,3, [k1,k1,k3], h)
            i13_122[k1,k3] = i13_112[k3,k1]
            i12[k1,k3] = intOverMassFn(1,2, [k1,k3], h)
            i22[k1,k3] = intOverMassFn(2,2, [k1,k3], h)

            t2h31[k1,k3] = 2.*(h.pk[k1]*i13_122[k1,k3]*i11[k1] + \
                             h.pk[k3]*i13_112[k1,k3]*i11[k3])

            t3hnoB[k1,k3] = (i11[k1]*h.pk[k1])**2 * i22[k3,k3] + \
                            (i11[k3]*h.pk[k3])**2 * i22[k1,k1] + \
                            4.*(i11[k1]*h.pk[k1])*(i11[k3]*h.pk[k3])*i22[k1,k3]

            t4hnoT[k1,k3] = 2.*i11[k1]*i11[k3]*h.pk[k1]*h.pk[k3] *\
                          (i21[k1]*i11[k3]*h.pk[k3] + \
                             i21[k3]*i11[k1]*h.pk[k1])

            # First Romberg-integrate explicitly-angular-averaged things to low precision
            pk1plusk3[k1,k3] = utils.openRomb(\
                lambda cth:c.pkInterp(M.sqrt(h.k[k1]**2 + h.k[k3]**2 + \
                                 2.*h.k[k1]*h.k[k3]*cth)), -1.,1.,eps=0.3,k=3)/2.

            b10[k1,k3] = utils.openRomb(lambda cth:tribi.b\
                                    (h.k[k1],h.k[k3], cth, c),-1.,1.,eps=0.3,k=3)/2.
            t3hB[k1,k3] = 4. * b10[k1,k3] * i12[k1,k3]*i11[k1]*i11[k3]


            #if k1 == k3:
            #t10[k1,k3] = 32.*h.pk[k1]**2*utils.openRomb(lambda cth: (3.+10*cth)**2*c.pkInterp(h.pk[k1]*M.sqrt(2.*(1-cth))),-1.,1.,eps = 0.3,k=2)/2. - 11./378.*h.pk[k1]**3
            # could change to this if we wanted to; quicker, but less uniform
            
            t10[k1,k3] = utils.openRomb(lambda cth:tribi.tk1mk1k2mk2_array\
                                    (h.k[k1],h.k[k3], cth, c,0),-1.,1.,eps=0.3,k=3)/2.

            t4hT[k1,k3] = t10[k1,k3] * i11[k1]**2 * i11[k3]**2

            tentativetotal = M.fabs(i04[k1,k3]+2*pk1plusk3[k1,k3]*i12[k1,k3]+t2h31[k1,k3]+\
                             t3hnoB[k1,k3]+t3hB[k1,k3]+t4hnoT[k1,k3]+t4hT[k1,k3] +\
                             adder[k1,k3])

            if (adder[k1,k3] != 0.):
                print 'adder = ',adder[k1,k3]

            #calculate Romberg-integrated things to high precision, if they are >1/2 of total
            if M.fabs(2*pk1plusk3[k1,k3]*i12[k1,k3]) > highprecisionthresh*tentativetotal:
                print 't2h22: ',pk1plusk3[k1,k3],
                pk1plusk3[k1,k3] = utils.openRomb(
                    lambda cth:c.pkInterp(M.sqrt(h.k[k1]**2 + h.k[k3]**2 + \
                                                 2.*h.k[k1]*h.k[k3]*cth)), -1.,1.,eps=0.03,k=7,jmax=18)/2.
                print pk1plusk3[k1,k3]
                
            if M.fabs(t3hB[k1,k3]) > highprecisionthresh*tentativetotal:
                print 't3hB: ',b10[k1,k3],
                b10[k1,k3] = utils.openRomb(lambda cth:tribi.b\
                                        (h.k[k1],h.k[k3], cth, c),-1.,1.,eps=0.01,k=5,jmax=30)/2.
                print b10[k1,k3]
                t3hB[k1,k3] = 4. * b10[k1,k3] * i12[k1,k3]*i11[k1]*i11[k3]
                   
            if M.fabs(t4hT[k1,k3]) > highprecisionthresh*tentativetotal:
                print 't4hT:', t10[k1,k3],
                    
                t10[k1,k3] = utils.openRomb(lambda cth:tribi.tk1mk1k2mk2_array\
                                        (h.k[k1],h.k[k3], cth, c,0),-1.,1.,eps=0.01,k=5)/2.
                print t10[k1,k3]
                t4hT[k1,k3] = t10[k1,k3] * i11[k1]**2 * i11[k3]**2

            nrm = 2.*h.pk[k1]*h.pk[k3]

            #output some stuff at each entry in the CovMat
            print k1,k3,i04[k1,k3]/nrm, (2.*pk1plusk3[k1,k3]*i12[k1,k3]+t2h31[k1,k3])/nrm, \
                  (t3hnoB[k1,k3]+t3hB[k1,k3])/nrm, t4hT[k1,k3]/nrm, t4hnoT[k1,k3]/nrm, (t4hT[k1,k3]+t4hnoT[k1,k3])/nrm,\
                  (i04[k1,k3]+ 2.*pk1plusk3[k1,k3]*i12[k1,k3]+t2h31[k1,k3]+ \
                   t3hnoB[k1,k3]+t3hB[k1,k3]+t4hT[k1,k3]+t4hnoT[k1,k3])/nrm
            
        pk1plusk3_perp[k1] = c.pkInterp(M.sqrt(2.)*h.k[k1])

        t3hB_perp[k1] = 4.*tribi.b(h.k[k1],h.k[k1],0.,c) *\
                           i12[k1,k1]*i11[k1]**2
        squaretri = tribi.tk1mk1k2mk2(h.k[k1],h.k[k1], 0., c,0)
        t4hT_perp[k1] = i11[k1]**4 * squaretri
        qsq[k1] = squaretri/(4.*h.pk[k1]**2 * (2.*pk1plusk3_perp[k1] +\
                                               h.pk[k1]))
        s = pk1plusk3_perp[k1]/h.pk[k1]
        dsq4hshouldbe[k1] = 0.085*(4.*h.pk[k1]**2 * \
                                   (2.*pk1plusk3_perp[k1] + h.pk[k1]))
        
        dsq1h[k1] = i04[k1,k1]
        dsq2h[k1] = 2.*pk1plusk3_perp[k1]*i12[k1,k1]**2. + t2h31[k1,k1]
        dsq2h31[k1] = t2h31[k1,k1]
        dsq3h[k1] = t3hnoB[k1,k1] + t3hB_perp[k1]
        dsq4hnoT[k1] = 4.*(i11[k1]*h.pk[k1])**3*i21[k1]
        dsq4h[k1] = t4hnoT[k1,k1] + t4hT_perp[k1]

    dsq = dsq1h + dsq2h + dsq3h + dsq4h

    df = h.k**3/(2.*M.pi**2)
    ot = 1./3.

    # These are debugging files; they output the square-configuration reduced trispectrum.
    #M.save(h.prefix+'dsq1h.dat',M.transpose([h.k,dsq1h]))
    #M.save(h.prefix+'dsq2h.dat',M.transpose([h.k,dsq2h]))
    #M.save(h.prefix+'dsq2h31.dat',M.transpose([h.k,dsq2h31]))
    #M.save(h.prefix+'dsq3h.dat',M.transpose([h.k,dsq3h]))
    #M.save(h.prefix+'dsq4h.dat',M.transpose([h.k,dsq4h]))
    rat = M.fabs(dsq4hnoT/t4hT_perp)

    t1h = i04
    t2h22 = 2.*pk1plusk3*i12**2
    for k1 in range(len(h.k)):
        for k3 in range(k1+1,len(h.k)):
            t10[k3,k1] = t10[k1,k3]
            t1h[k3,k1] = t1h[k1,k3]
            t2h22[k3,k1] = t2h22[k1,k3]
            t2h31[k3,k1] = t2h31[k1,k3]
            t3hnoB[k3,k1] = t3hnoB[k1,k3]
            t3hB[k3,k1] = t3hB[k1,k3]
            t4hnoT[k3,k1] = t4hnoT[k1,k3]
            t4hT[k3,k1] = t4hT[k1,k3]
    
    t2h = t2h22 + t2h31
    t3h = t3hnoB + t3hB
    t4h = t4hnoT + t4hT

    ans = t1h+t2h+t3h+t4h
    
    if h.p.outputalltterms == 0:
        return ans
    elif h.p.outputalltterms == 1:
        return ans,t10,t1h,t2h,t3h,t4h
    elif h.p.outputalltterms == 2:
        return ans,t10,t1h,t2h22,t2h31,t3hB,t3hnoB,t4hT,t4hnoT
    else:
        return

def getHaloCov(prefix,c,h):
    """
    Output halo model covariance matrix, correlation matrix into the directory 'prefix'
    """
    os.system('mkdir '+prefix)

    h.pnl = getHaloPknl(c,h)

    M.save(prefix+'pnl.dat',M.transpose([h.k,h.pnl]), fmt = '%18.16e')

    h.prefix = prefix
    #h.dlogPnldlogA = getdlogPnldlogA(c,h)
    #M.save(prefix+'dlogpnldloga.dat',M.transpose([h.k,h.dlogPnldlogA]),fmt='%6.5e')
    vk = h.k*0.
    vk[0] = (h.k[0]*h.k[1])**1.5 - (h.k[0]**3/h.k[1])**1.5
    for k1 in M.arange(1,len(h.k)-1):
        vk[k1] = (h.k[k1]*h.k[k1+1])**1.5 - (h.k[k1]*h.k[k1-1])**1.5
    vk[-1] = (h.k[-1]**3/h.k[-2])**1.5 - (h.k[-1]*h.k[-2])**1.5
    vk *= 4.*M.pi/3.

    gausspart = M.outer(h.k*0.,h.k*0.)
    for k1 in M.arange(len(h.k)):
        gausspart[k1,k1] = (2.*M.pi)**3 * 2.*(h.pnl[k1]**2)/vk[k1]

    if h.p.outputalltterms == 0:
        t = getHaloTrispec(c,h, adder=gausspart)
    elif h.p.outputalltterms == 1:
        t,t10,t1h,t2h,t3h,t4h = getHaloTrispec(c,h, adder=gausspart)
    elif h.p.outputalltterms == 2:
        t,t10,t1h,t2h22,t2h31,t3hB,t3hnoB,t4hT,t4hnoT = getHaloTrispec(c,h,adder=gausspart)
    
    covar = t*1.

    cocg = h.k*0.
    for k1 in M.arange(len(h.k)):
        cocg[k1] = M.sqrt(covar[k1,k1]/gausspart[k1,k1])

    covar += gausspart
    M.save(prefix+'covar.dat',covar, fmt = '%18.16e')
    M.save(prefix+'gausspart.dat',gausspart, fmt = '%18.16e')
    #t10->pt.dat is the perturbation theory trispectrum by itself.  However,
    #it might not be calculated at high precision in the nonlinear
    #regime, since where other things dominate on small scales, the
    #first pass at calculating it is done with low precision.

    #
    if h.p.outputalltterms == 1:
        M.save(prefix+'pt.dat',t10, fmt = '%18.16e')
        M.save(prefix+'t1h.dat',t1h, fmt = '%18.16e')
        M.save(prefix+'t2h.dat',t2h, fmt = '%18.16e')
        M.save(prefix+'t3h.dat',t3h, fmt = '%18.16e')
        M.save(prefix+'t4h.dat',t4h, fmt = '%18.16e')
    if h.p.outputalltterms == 2:
        M.save(prefix+'pt.dat',t10, fmt = '%18.16e')
        M.save(prefix+'t1h.dat',t1h, fmt = '%18.16e')
        M.save(prefix+'t2h22.dat',t2h22, fmt = '%18.16e')
        M.save(prefix+'t2h31.dat',t2h31, fmt = '%18.16e')
        M.save(prefix+'t3hB.dat',t3hB, fmt = '%18.16e')
        M.save(prefix+'t3hnoB.dat',t3hnoB, fmt = '%18.16e')
        M.save(prefix+'t4hT.dat',t4hT, fmt = '%18.16e')
        M.save(prefix+'t4hnoT.dat',t4hnoT, fmt = '%18.16e')
        
    correl = 0.*covar

    tnorm = t
    for i in M.arange(len(h.k)):
        for j in M.arange(len(h.k)):
            correl[i,j] = covar[i,j]/M.sqrt(covar[i,i]*covar[j,j])

    M.save(prefix+'nbins.dat',M.array([len(h.k)]), fmt = '%d')
    M.save(prefix+'correl.dat',correl, fmt = '%4.3f')

def degraderesolution(prefix,factor,dlogstring):
    covar = M.load(prefix+'covar.dat')
    pnl = M.load(prefix+'pnl.dat')
    dlog = M.load(prefix+dlogstring)[:,1]
    k = pnl[:,0]*1.
    p = pnl[:,1]*1.
    gausspart = M.load(prefix+'gausspart.dat')
    nbins = len(k)

    nongausspart = covar - gausspart

    nongausspartnew = nongausspart[:nbins-factor:factor,:nbins-factor:factor]*0.
    knew = k[:nbins-factor:factor]*0.
    pnew = p[:nbins-factor:factor]*0.
    gausspartnew = gausspart[:nbins-factor:factor,:nbins-factor:factor]*0.
    nbinsnew = len(knew)
    dlognew = dlog[:nbins-factor:factor]*0.

    for i1 in range(0,nbins-factor,factor):
        i1new = i1/factor
        print i1,i1+factor-1,nbins
        print i1new,nbinsnew
        weights = k[i1:i1+factor-1]**3
        sumweights = M.sum(weights)
        pnew[i1new] = M.sum(p[i1:i1+factor-1]*weights)/sumweights
        knew[i1new] = M.sum(k[i1:i1+factor-1]*weights)/sumweights
        dlognew[i1new] = M.sum(dlog[i1:i1+factor-1]*weights)/sumweights

    sqrtkfact = M.sqrt(k[1]/k[0])
        
    for i1 in range(0,nbins-factor,factor):
        i1new = i1/factor
        for i2 in range(0,nbins-factor,factor):
            i2new = i2/factor
                                                                       
            weights2 = M.outer(k[i1:i1+factor-1]**3,k[i2:i2+factor-1]**3)
            sumweights2 = M.sum(M.sum(weights2))
            nongausspartnew[i1new,i2new] = M.sum(M.sum(nongausspart[i1:i1+factor-1,i2:i2+factor-1]*weights2))/sumweights2

            if i1new == i2new:
                vk = (4.*M.pi/3.)*((k[i1+factor-1]*sqrtkfact)**3 - (k[i1]/sqrtkfact)**3)
                gausspartnew[i1new,i2new] = (2.*M.pi)**3 * 2.*(pnew[i1new]**2)/vk
                                                                       
    covarnew = gausspartnew + nongausspartnew

    prefixnew = prefix+'degrade'+str(factor)+'/'
    os.system('mkdir '+prefixnew)
    M.save(prefixnew+'pnl.dat',M.transpose([knew,pnew]), fmt = '%18.16e')
    M.save(prefixnew+'covar.dat',covarnew, fmt = '%18.16e')
    M.save(prefixnew+'gausspart.dat',gausspartnew, fmt = '%18.16e')
    M.save(prefixnew+dlogstring,M.transpose([knew,dlognew]), fmt = '%18.16e')
    M.save(prefix+'nbins.dat',M.array([nbinsnew],shape=(1,1,)), fmt = '%d')
