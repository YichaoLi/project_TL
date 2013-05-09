import pylab as M
import numpy as N
import f77io
import sys
import myload
#import scipy.stats as SS
import scipy.special as SSp
import power
import cosmology.pt
import fontsize
#import cepstrum
import itertools as I

def getds():
    d63 = f77io.getdatacube()
    d45 = f77io.getdatacube(file='densmesh45.f77')
    d0 = f77io.getdatacube(file='densmesh0.f77')

    return d0,d45,d63

def dists(dlist):
    for d in dlist:
        M.plot(d.flatten(),N.log(d.flatten()),',')

    M.show()

def decreaseres(d,f=2):
    print d.shape, len(d.shape)

    if (f == 1):
        return(d)

    if (len(d.shape) == 3):
        ng = len(d[:,0,0])
        newng = ng/f
        newd = N.zeros(newng**3).reshape((newng,newng,newng))
        for i in range(newng):
            for j in range(newng):
                for k in range(newng):
                    newd[i,j,k] = N.mean(d[f*i:f*(i+1),f*j:f*(j+1),f*k:f*(k+1)])

    elif (len(d.shape) == 2):
        ng = len(d[:,0])
        newng = ng/f
        newd = N.zeros(newng**2).reshape((newng,newng))
        for i in range(newng):
            for j in range(newng):
                newd[i,j] = N.mean(d[f*i:f*(i+1),f*j:f*(j+1)])

    return(newd)
                

def findshift(deni,denf_g,sti=0,endi=0):
    ng = len(deni[:,0,0])
    if endi ==0:
        endi = ng
    deni_copy = 1.*deni
    deni_copy = N.roll(deni_copy,sti-1,axis=0)
    for i in N.arange(sti,endi):
        deni_copy = N.roll(deni_copy,1,axis=0)
        for j in range(ng):
            deni_copy = N.roll(deni_copy,1,axis=1)
            for k in range(ng):
                deni_copy = N.roll(deni_copy,1,axis=2)
                print N.sum((deni_copy-denf_g).flatten()**2), i,j,k
                sys.stdout.flush()

def showshift(di,dfg):
    shifts = M.load('findshift.64.head')[:,1:].astype(int)
    #shifts *= 4
    for i in range(10):
        print shifts[i,:]
        M.subplot(221)
        M.pcolor(N.roll(N.roll(N.roll(di,shifts[i,0],axis=0),
                               shifts[i,1],axis=1),shifts[i,2],axis=2)[:,:,0])
        M.axis('tight')
        M.subplot(222)
        M.pcolor(dfg[:,:,0])
        M.axis('tight')
        M.subplot(223)
        M.pcolor(di[:,:,0]-dfg[:,:,0])
        M.axis('tight')
        M.show()


def makedensgaussanal(denf,uselog=False,avgrepeats=True, just0=False, tiegauss=False, sigmagauss = 0.,spatialsigma=8.):
    denshape = denf.shape

    denff = denf.flatten()
    if (tiegauss): 
        # use actual density first to order, then to break ties use 
        # density smoothed in 8-Mpc radius
        denny = N.empty(N.prod(denf.shape),dtype='float32,float32')
        denny['f0'] = denf.flatten()
        denny['f1'] = smoothgauss(denf,sigma=spatialsigma).flatten()
        o_f = N.argsort(denny,order=('f0','f1'))

    else:
        o_f = N.argsort(denff)

    gaussf = 0.*denff.astype(N.float)
    lenny = len(gaussf)


    if (sigmagauss == 0.):
        if (uselog):
            sigmagauss = N.sqrt(N.var(N.log(denff)))
        else:
            sigmagauss = N.sqrt(N.var(denff))

    step = 1./lenny

    gaussf[o_f] = N.sqrt(2.)*sigmagauss*SSp.erfinv(2.*N.arange(0.5*step,1,step)-1.)

    # average together repeated elements
    if (avgrepeats):
        i = 0
        while (i < lenny-1):
            if (denff[o_f[i]] == denff[o_f[i+1]]):
                numdups = N.searchsorted(denff[o_f[i:]],denff[o_f[i]],
                                         side='right')
                if (numdups > 1):
                    gaussf[o_f[i:i+numdups]] = N.zeros(numdups) + \
                        N.mean(gaussf[o_f[i:i+numdups]])
            else: 
                numdups = 1

            i += numdups
            if (just0):
                break

    gaussf = gaussf.reshape(denshape)

    return gaussf

def makedensgaussanal2(denf,uselog=False,avgrepeats=True, just0=False, tiegauss=False, sigmagauss = 0.,spatialsigma=8.):
    denshape = denf.shape

    denff = denf.flatten()
    mind = min(denff)
    denff = denff - mind + 0.1

    gaussf = 0.*denff.astype(N.float)
    lenny = len(gaussf)

    o_f = N.argsort(denff)

    if (sigmagauss == 0.):
        if (uselog):
            sigmagauss = N.sqrt(N.var(N.log(denff)))
        else:
            sigmagauss = N.sqrt(N.var(denff))

    step = 1./lenny

    gaussf[o_f] = N.sqrt(2.)*sigmagauss*SSp.erfinv(2.*N.arange(0.5*step,1,step)-1.)

    # copied with little comprehension from 
    # http://www.weask.us/entry/find-length-sequences-identical-values-numpy-array
    if (avgrepeats):
        place = 0
        print len(denff[o_f])
        for a, b in I.groupby(denff[o_f]):
            sumb=sum(b)
            numdup = sumb/a
            #print a,sumb,numdup,place,denff[o_f[place:place+numdup]]
            if (numdup > 1):
                meanie = N.mean(gaussf[o_f[place:place+numdup]])
                #print numdup,meanie,gaussf[o_f[place:place+numdup]],denff[o_f[place:place+numdup]]
                gaussf[o_f[place:place+numdup]] = N.repeat(meanie,numdup)
            place += numdup

    gaussf = gaussf.reshape(denshape)

    return gaussf

def histogramGauss2d(d1,d2,bins=256,sigmagauss = 1.,plotty=False,uniform=False):
    

    o1 = N.argsort(d1)
    o2 = N.argsort(d2)

    d1g = 0.*d1
    d2g = 0.*d2
    lenny = len(d1)

    step = 1./lenny

    if (uniform):
        d1g[o1] = N.arange(0.5*step,1,step).astype(N.float32)
        d2g[o2] = N.arange(0.5*step,1,step).astype(N.float32)
    else:
        d1g[o1] = (N.sqrt(2.)*sigmagauss*SSp.erfinv(2.*N.arange(0.5*step,1,step)-1.)).astype(N.float32)
        d2g[o2] = N.sqrt(2.)*sigmagauss*SSp.erfinv(2.*N.arange(0.5*step,1,step)-1.).astype(N.float32)

    print d1g[o1[-5:]],d2g[o2[-5:]]
    print d1[o1[-5:]],d2[o2[-5:]]

    hist2d, xedgesg, yedgesg = N.histogram2d(d1g,d2g,bins=bins,normed=True)

    where_xedgesg = d1g[o1].searchsorted(xedgesg)
    where_yedgesg = d2g[o2].searchsorted(yedgesg)
    
    xcenters = 0.*xedgesg[:-1]
    ycenters = 0.*yedgesg[:-1]
    xbinsize = 0.*xedgesg[:-1]
    ybinsize = 0.*yedgesg[:-1]
    xedges = d1[o1[where_xedgesg]]
    yedges = d2[o2[where_yedgesg]]
    
    print xedgesg[-5:]
    print where_xedgesg[-5:]
    print xedges[-5:]
    print hist2d[-5:,-5:]

    for i in range(len(xcenters)):
        xcenters[i] = N.mean(d1[o1[where_xedgesg[i]:where_xedgesg[i+1]]])
        ycenters[i] = N.mean(d2[o2[where_yedgesg[i]:where_yedgesg[i+1]]])
        xbinsize[i] = xedges[i+1]-xedges[i]
        ybinsize[i] = yedges[i+1]-yedges[i]

    if (plotty):
        M.subplot(221)
        M.pcolor(xedgesg,yedgesg,hist2d)
        M.colorbar()
        M.subplot(222)
        M.pcolor(xedgesg,yedgesg,hist2d*N.outer(xbinsize,ybinsize))
        M.colorbar()
        M.subplot(223)
        M.pcolor(xedges,yedges,hist2d)

    return hist2d,xedgesg,yedgesg,xedges,yedges,xcenters,ycenters

def makedensgauss(deni,denf):
    ng = len(denf[:,0,0])
    o_f = N.argsort(denf.flatten())
    o_i = N.argsort(deni.flatten())

    gaussf = 0.*denf.flatten()

    gaussf[o_f] = deni.flatten()[o_i]
    gaussf = gaussf.reshape((ng,ng,ng))    

    return gaussf

def var(dlist):
    for d in dlist:
        vard = N.var(d.flatten())
        varlogd = N.var(N.log(d.flatten()))
        print vard,varlogd,vard/varlogd

def showf(dlist,r=2,c=2,make2d=False):
    M.clf()
    
    pn = 0
    print make2d
    if make2d:
        for d in dlist:
            pn += 1
            M.subplot(r,c,pn)
            M.pcolor(d[:,:,0])
            M.axis('tight'); M.colorbar()

    else:
        for d in dlist:
            pn += 1
            M.subplot(r,c,pn)
            M.pcolor(d)
            M.axis('tight'); M.colorbar()

def show3(dlist,r=2,c=2,greyscale=False,output=False,samerange=True):

    """ For paper figure, dlist should be (d63-1, d0-1, log(d63), d63ga) """

#distrib.show3((d63[:128,:128,0]-1,d0[:128,:128,0]-1,N.log(d63[:128,:128,0]),d63ga[:128,:128,0]),greyscale=True)

    M.clf()

    fig = M.figure(figsize=(6.4, 6.4), dpi=100) 
    axesarr=N.array([[0.01,0.51,0.4,0.4],
                     [0.51,0.51,0.4,0.4],
                     [0.01,0.01,0.4,0.4],
                     [0.51,0.01,0.4,0.4]])

    print axesarr
    colorbax = 1.*axesarr
    print colorbax
    colorbax[:,2] = 0.*colorbax[:,2] + 0.03
    colorbax[:,0] += 0.4

    print colorbax

    if greyscale:
        colorscheme='binary'
    else:
        colorscheme='jet'

    # d63, d0, log d63, d63g
    titlearr=[r'$\delta$',r'$\delta_{\rm initial}$',r'$\log(1+\delta)$',r'$\delta_{\rm Gauss}$']

    if (dlist[1] != None):
        min23 = min(min(dlist[2].flatten()),min(dlist[3].flatten()))
        max23 = max(max(dlist[2].flatten()),max(dlist[3].flatten()))

        max0 = max(dlist[1].flatten())
        min0 = min(dlist[1].flatten())

        initfact = min(max23/max0,min23/min0)
        print min23,max23, initfact

    sc = 0
    for d in dlist:
        if (d != None):
            M.axes(axesarr[sc])
            M.title(titlearr[sc],fontsize=23)
            if (sc > 1):
                print titlearr[sc]
                if (samerange):
                    M.pcolor(d,cmap=M.get_cmap(colorscheme),vmin = min23,vmax=max23)
                else:
                    M.pcolor(d,cmap=M.get_cmap(colorscheme))
            elif (sc == 1):
            #print min(d.flatten()*initfact),max(d.flatten()*initfact)
                if (samerange):
                    M.pcolor(d*initfact,cmap=M.get_cmap(colorscheme),vmin = min23,vmax=max23)
                else:
                    M.pcolor(d,cmap=M.get_cmap(colorscheme))

            else:
                M.pcolor(d,cmap=M.get_cmap(colorscheme))

#        if (sc == 1):
#            M.colorbar(ticks=[-0.1,-0.05,0,0.05,0.1])
#        else:

            M.axis('tight')
            M.axis('equal')
            M.axis('tight')
            M.xticks([])
            M.yticks([])

            cax = M.axes(colorbax[sc])
            M.colorbar(cax=cax)

        sc += 1

    #M.savefig('showdens.eps',dpi=8)
    #M.gcf().set_size_inches((6.4,6.4))
    #M.gcf().set_size_inches((15.,12.))
    if (output):
        if greyscale:
            M.savefig('showdens_grey.png',dpi=100)
            M.savefig('showdens_grey.pdf')
        else:
            fig.savefig('showdens.png',dpi=100)
            M.savefig('showdens.pdf')

    #M.show()

def testeps(d):
    M.clf()
    M.pcolor(d)
    M.axis('tight')
    M.colorbar()
    M.gcf().set_size_inches((7.5,6.))
    M.savefig('test.png',dpi=240)


def showpplog(prefix,color,fact=1.,xis=1.,xislog=1.,sumto=10,camb=0,cellsize=0):
    p = M.load(prefix+'/pm.pnl.dat')
    plog = M.load(prefix+'/plogm.pnl.dat')

    # all times xis
    sump = N.sum(M.load(prefix+'c11')[:sumto,1]*p[:sumto,2])
    c21xis = N.sum(M.load(prefix+'c21')[:sumto,1]*p[:sumto,2])/sump
    c22xis = N.sum(M.load(prefix+'c22')[:sumto,1]*p[:sumto,2])/sump
    c31xis = N.sum(M.load(prefix+'c31')[:sumto,1]*p[:sumto,2])/sump
    
    #bias = N.sum(p[:sumto,1]*p[:sumto,2])/N.sum(plog[:sumto,1]*plog[:sumto,2])
    bias = N.sum((p[:sumto,1]/plog[:sumto,1]) * p[:sumto,2])/N.sum(p[:sumto,2])
    biaserror = N.std(p[:sumto,1]/plog[:sumto,1])
    simpleapprox = 1./(1.-0.44*xis)

    c21 = camb.c21(cellsize)
    c22 = camb.c22(cellsize)
    c31 = camb.c31(cellsize)
    s3 = camb.s3(cellsize)
    #print cellsize,c21, c21/xis
    approx = 1./(1+xis*(2.-c21))
    approx2 = 1./(1+xis*(2-c21)+xis**2*(7-2*s3-4*c21 + 2.*c31/3. + c22/4.))
    #print bias,simpleapprox,approx,approx2
    print bias,biaserror

    M.loglog([cellsize],[simpleapprox-1],'yo')
    M.loglog([cellsize],[approx-1],'rp')
    M.loglog([cellsize],[approx2-1.],'bh')
    M.loglog([cellsize],[fact-1.],'gD')
    M.loglog([cellsize],[bias-1],'k.')
    #M.loglog([cellsize],[approx2],'r.')

def rejuice_hv():
    meanlogrho = N.exp(-N.array([-8.5677281E-03, -3.7504282E-02, -0.1114288, 
                           -0.2148134, -0.4872689, -1.163443]))
    xis = N.array([0.017202929, 0.077699677, 0.2589330, 0.606343639, 2.14528, 6.77099])

    mlr_mill = N.exp(-N.array([-0.2402180, -0.49094526, -0.8136205, -1.161325]))
    vd_mill = N.array([0.706927, 2.51522, 9.590195, 35.7166])
    vld_mill = N.array([0.4493474,0.827824,1.197856,1.49259])

    #showpplog('hv/r512','c',fact=meanlogrho[4],vardelta =vardelta[4])
    #showpplog('hv/r1024','y',fact=meanlogrho[5],vardelta =vardelta[5])
    c_hv = pt.Camb(hubble=70.,ombh2=0.0196,omch2=0.3*0.7**2)
    c_hv.run()
    sig82_hv = pt.normalizePk(c_hv,0.9)
    c_mill = pt.Camb(hubble=73., ombh2 = (0.25-0.045)*0.7**2, omch2 = 0.25*0.7**2)
    c_mill.run()
    sig82_mill = pt.normalizePk(c_mill,0.9)
    print sig82_hv, sig82_mill

    M.clf()

    M.xlabel('cellsize (h/Mpc)',fontsize=20)
    M.ylabel(r'$P_\delta(k)/P_{\log (1+\delta)}(k)-1$',fontsize=20)
    M.axes([0.15,0.15,0.8,0.8])

    bhv32=showpplog('hv/r32/','b',fact=meanlogrho[0],xis =xis[0],camb=c_hv,
              cellsize=3000/32.)
    bhv64=showpplog('hv/r64/','g',fact=meanlogrho[1],xis =xis[1],camb=c_hv,
              cellsize=3000/64.)
    bhv128=showpplog('hv/r128/','r',fact=meanlogrho[2],xis =xis[2],camb=c_hv,
              cellsize=3000/128.)
    bhv256=showpplog('hv/r256/','m',fact=meanlogrho[3],xis =xis[3],camb=c_hv,
              cellsize=3000/256.)
  
    bm32=showpplog('mill/s63r8/','b--',fact=mlr_mill[0],xis=vd_mill[0],camb=c_mill,
              cellsize=500./32.)
    bm64=showpplog('mill/s63r4/','g--',fact=mlr_mill[1],xis =vd_mill[1],camb=c_mill,
              cellsize=500./64.)
    bm128=showpplog('mill/s63r2/','r--',fact=mlr_mill[2],xis =vd_mill[2],camb=c_mill,
              cellsize=500./128.)
    bm256=showpplog('mill/s63/','m--',fact=mlr_mill[3],xis =vd_mill[3],camb=c_mill,
              cellsize=500./256.)

    M.show()

def rejuice_hv():
    meanlogrho = N.exp(-N.array([-8.5677281E-03, -3.7504282E-02, -0.1114288, 
                           -0.2148134, -0.4872689, -1.163443]))
    xis = N.array([0.017202929, 0.077699677, 0.2589330, 0.606343639, 2.14528, 6.77099])

    mlr_mill = N.exp(-N.array([-0.2402180, -0.49094526, -0.8136205, -1.161325]))
    vd_mill = N.array([0.706927, 2.51522, 9.590195, 35.7166])
    vld_mill = N.array([0.4493474,0.827824,1.197856,1.49259])

    #showpplog('hv/r512','c',fact=meanlogrho[4],vardelta =vardelta[4])
    #showpplog('hv/r1024','y',fact=meanlogrho[5],vardelta =vardelta[5])

    c_hv = pt.Camb(hubble=70.,ombh2=0.0196,omch2=0.3*0.7**2)
    c_hv.run()
    sig82_hv = pt.normalizePk(c_hv,0.9)
    c_mill = pt.Camb(hubble=73., ombh2 = (0.25-0.045)*0.7**2, omch2 = 0.25*0.7**2)
    c_mill.run()
    sig82_mill = pt.normalizePk(c_mill,0.9)
    print sig82_hv, sig82_mill

    bhv32=showpplog('hv/r32/','b',fact=meanlogrho[0],xis =xis[0],camb=c_hv,
              cellsize=3000/32.)
    bhv64=showpplog('hv/r64/','g',fact=meanlogrho[1],xis =xis[1],camb=c_hv,
              cellsize=3000/64.)
    bhv128=showpplog('hv/r128/','r',fact=meanlogrho[2],xis =xis[2],camb=c_hv,
              cellsize=3000/128.)
    bhv256=showpplog('hv/r256/','m',fact=meanlogrho[3],xis =xis[3],camb=c_hv,
              cellsize=3000/256.)
  
    bm32=showpplog('mill/s63r8/','b--',fact=mlr_mill[0],xis=vd_mill[0],camb=c_mill,
              cellsize=500./32.)
    bm64=showpplog('mill/s63r4/','g--',fact=mlr_mill[1],xis =vd_mill[1],camb=c_mill,
              cellsize=500./64.)
    bm128=showpplog('mill/s63r2/','r--',fact=mlr_mill[2],xis =vd_mill[2],camb=c_mill,
              cellsize=500./128.)
    bm256=showpplog('mill/s63/','m--',fact=mlr_mill[3],xis =vd_mill[3],camb=c_mill,
              cellsize=500./256.)

    M.show()


def rejuice(d63,d63_2,d63_4,d63_8):
    #pinit = M.load('mill/s63/pm.pnl.dat')

    p1 = M.load('mill/s63/pm.pnl.dat')
    plog1 = M.load('mill/s63/plogm.pnl.dat')
    p2 = M.load('mill/s63r2/pm.pnl.dat')
    plog2 = M.load('mill/s63r2/plogm.pnl.dat')
    p4 = M.load('mill/s63r4/pm.pnl.dat')
    plog4 = M.load('mill/s63r4/plogm.pnl.dat')
    p8 = M.load('mill/s63r8/pm.pnl.dat')
    plog8 = M.load('mill/s63r8/plogm.pnl.dat')

    f63= N.exp(-N.mean(N.log(d63.flatten())))
    f63_2= N.exp(-N.mean(N.log(d63_2.flatten())))
    f63_4= N.exp(-N.mean(N.log(d63_4.flatten())))
    f63_8= N.exp(-N.mean(N.log(d63_8.flatten())))

    #M.loglog(p1[:,0],p1[:,1]/(plog1[:,1]*f63),'b--')
    #M.loglog(p2[:,0],p2[:,1]/(plog2[:,1]*f63_2),'g--')
    #M.loglog(p4[:,0],p4[:,1]/(plog4[:,1]*f63_4),'r--')
    #M.loglog(p8[:,0],p8[:,1]/(plog8[:,1]*f63_8),'y--')

    #xis = N.mean(d63.flatten()**2)
    #xis_2 = N.mean(d63_2.flatten()**2)
    #xis_4 = N.mean(d63_4.flatten()**2)
    #xis_8 = N.mean(d63_8.flatten()**2)

    xis = (1.+ 0.5*N.sqrt(N.var(d63.flatten())))
    xis_2 = (1.+0.5*N.sqrt(N.var(d63_2.flatten())))
    xis_4 = (1.+0.5*N.sqrt(N.var(d63_4.flatten())))
    xis_8 = 1.+0.5*N.sqrt(N.var(d63_8.flatten()))
    
    print 'exps:',f63,f63_2,f63_4,f63_8
    print 'xis:',xis, xis_2,xis_4,xis_8

    M.loglog(plog1[:,0],p1[:,1]/(plog1[:,1]*f63)*(1.+2.*xis**2),'b')
    M.loglog(plog2[:,0],p2[:,1]/(plog2[:,1]*f63_2)*(1.+2.*xis_2**2),'g')
    M.loglog(plog4[:,0],p4[:,1]/(plog4[:,1]*f63_4)*(1.+2.*xis_4**2),'r')
    M.loglog(plog8[:,0],p8[:,1]/(plog8[:,1]*f63_8)*(1.+2.*xis_8**2),'y')

    M.loglog(plog1[:,0],p1[:,1]/(plog1[:,1]*xis),'b')
    M.loglog(plog2[:,0],p2[:,1]/(plog2[:,1]*xis_2),'g')
    M.loglog(plog4[:,0],p4[:,1]/(plog4[:,1]*xis_4),'r')
    M.loglog(plog8[:,0],p8[:,1]/(plog8[:,1]*xis_8),'y')


    M.xlabel(r'$k\ [\rm{Mpc}/h]$',fontsize=20)
    M.ylabel(r'$P_\delta(k)/P_{\log (1+\delta)}(k)$',fontsize=20)

    bias1 = N.sum(p1[:5,1]*p1[:5,2])/N.sum(plog1[:5,1]*plog1[:5,2])
    bias2 = N.sum(p2[:5,1]*p2[:5,2])/N.sum(plog2[:5,1]*plog2[:5,2])
    bias4 = N.sum(p4[:5,1]*p4[:5,2])/N.sum(plog4[:5,1]*plog4[:5,2])
    bias8 = N.sum(p8[:5,1]*p8[:5,2])/N.sum(plog8[:5,1]*plog8[:5,2])

    print bias1,bias2,bias4,bias8#, N.log(bias1),N.log(bias2),N.log(bias4)       
    M.show()

def histo(denlist):

    for d in denlist:
        h,bin_edges = N.histogram(N.log10(d.flatten()),normed=True,range=(-2,2),bins=50)
        M.semilogx(10.**(bin_edges),h)

    M.show()



def logskew(d,floors = [],col='',divvypk=[]):
    maxd = max(d.flatten())
    if len(floors) == 0:
        floors = 1./maxd * 2.**N.arange(-5,5.01,1.)
    print floors

    for fl in floors:
        dcopy = 1.*d
        wltf = N.where(dcopy < fl)
        dcopy[wltf] = 0.*dcopy[wltf] + fl

        logrho = N.log(dcopy)
        var  = N.var(decreaseres(logrho,f=16).flatten())
        #skew = SS.skew(logrho.flatten())
        #print fl,'log:',var,'lin:',N.var(dcopy.flatten())
        print fl,'log:',var#,SS.skew(dcopy.flatten())

        k,pk = power.pk(logrho)

        if (len(divvypk) > 0):
            M.loglog(k,pk/divvypk/var,col)

    if (len(divvypk) == 0):
        return pk,pk/var
    else:
        return

def shotnoise(d1,d4): # for hv
    M.clf()
    #logskew(d2,floors=N.arange(-2,2.1,2),col='r')

    realpk,realpkovar = logskew(d1,floors=[1e-10],col='k--')
    logskew(d4,floors=0.2*2**N.arange(-2,2.1,2),col='',divvypk=realpkovar)

    M.show()

def gaussy(d):
    colors = ['b','g','r','c','m','y','b--','g--','r--','c--','m--','y--']
    multi = 2.**N.arange(0,4.1,1)

    M.clf()

    for i in range(len(multi)):
        dinflate = (d-1.)*multi[i]
        k,p = power.pk(dinflate)
        w0 = N.where(dinflate <=-0.9999)
        dinflate[w0] = dinflate[w0]*0. -0.9999
        k,plog = power.pk(N.log(dinflate+1.))
        print p[0]/plog[0], p[-1]/plog[-1]

        M.loglog(k,p/plog,colors[i])

    M.show()

def expy(d):
    colors = ['b','g','r','c','m','y','k','b--','g--','r--','c--','m--','y--']
    multi = 2.**N.arange(2,3.1,1.)
    #multi = [77.]
    for i in range(len(multi)):
        dinflate = N.exp((d)*multi[i])-1.
        k,p = power.pk(dinflate)

        var = N.var((d)*multi[i])
        lognocor = ((N.exp(var)-1)*N.exp(var)/var)
        lognovar = (N.exp(var)-1)*N.exp(var)
        print 1+2*var,'lognocor:',lognocor
        #varlog = N.var(N.log(dinflate+1.).flatten())
        k,plog = power.pk(N.log(dinflate+1))
        print p[0]/plog[0], p[-1]/plog[-1]

        #M.subplot(121)
        M.loglog(k,p/plog/lognocor,colors[i])
        #M.subplot(122)

        #M.loglog(k,plog,colors[i])
        #M.loglog([k[0],k[-1]],[lognovar,lognovar],colors[i])
        #M.loglog(k,1./plog**mul,colors[i])
        

    
    #preal = M.load('mill/s63/pm.pnl.dat')
    #plogreal = M.load('mill/s63/plogm.pnl.dat')    
    #M.loglog(preal[:,0],preal[:,1]/plogreal[:,1],'b')

    #M.xlabel(r'$k\ [\rm{Mpc}/h]$',fontsize=20)
    #M.ylabel(r'$P_\delta(k)/P_{\log (1+\delta)}(k)$',fontsize=20)

    M.show()


def vars(listd):

    M.clf()

    for d in listd:
        print N.var(d-1.)

def momentous(listd):

    M.clf()
    for d in listd:

        momentnum = N.arange(2,10)
        moment = 0.*momentnum
        partialsum = 0.*moment
        for i in range(len(momentnum)):
            moment[i] = N.mean((d.flatten()-1)**momentnum[i])

        print moment
        M.subplot(211)
        M.semilogy(momentnum,moment,'k.')
        for i in N.arange(len(momentnum)-1,-1,-1):
            print moment[i],
            moment[i] /= moment[0]**(momentnum[i])
            print moment[i]
            partialsum[i:] +=moment[i]/float(momentnum[i])*(-1)**(momentnum[i]+1)
        print moment
        M.semilogy(momentnum,moment,'b.')
        wpos = N.where(partialsum > 0.)
        wneg = N.where(partialsum < 0.)

        print moment
        print partialsum
        M.subplot(212)
        M.semilogy(momentnum[wpos],partialsum[wpos],'k.')
        M.semilogy(momentnum[wneg],-partialsum[wneg],'r.')

        M.show()

def ngp(p, boxsize=2.,ngrid=256):

    pint = N.floor(p*ngrid/boxsize).astype(N.int)
    print 'shape:',p.shape
    dim = p.shape[1]
    print 'dim=',dim
    if dim == 2:
        grid = N.zeros(ngrid**2).reshape(ngrid,ngrid)
        for i in range(len(p[:,0])):
            grid[pint[i,0],pint[i,1]] += 1.

    elif dim == 3:
        grid = N.zeros(ngrid**3).reshape(ngrid,ngrid,ngrid)
        for i in range(len(p[:,0])):
            grid[pint[i,0],pint[i,1],pint[i,2]] += 1.
    else:
        print 'Croptersplank!'

    gridmean = N.mean(grid.flatten())
    print 'gridmean = ',gridmean
    grid /= N.sum(grid.flatten())
    #grid /= gridmean

    return grid

def histmom(dlist,redshifts,nbins=50):
    """ d=rho, not delta """
    #h,b = N.histogram(d.flatten()-1.,bins=200,normed=True,range=(-1.,5.))
    #b += (b[1]-b[0])/2.
    #M.plot(b,h)
    #M.show()
    colors = ['g','k','b']
    thicks = [1,1,1]

    pn = 0
    M.clf()
    M.gcf().set_size_inches((6,2.4))


    axes12 = [[0.14,0.255,0.85,0.83-0.25],[0.14,0.83,0.85,0.145]]
    for d in dlist:
        print 'gaddytime?'
        M.axes(axes12[0])
        std = N.sqrt(N.var(d.flatten()))
        h,b = N.histogram((d.flatten()-1.)/std,bins=nbins,range=(-5,5),normed=True)
        b += (b[1]-b[0])/2.
        M.plot(b,h,'b',linewidth=1.5,label=r'$\delta$')

        print max(h)

        logd = N.log(d.flatten())
        meanlog = N.mean(logd)
        stdlog = N.sqrt(N.var(logd))
        hlog,blog = N.histogram((N.log(d.flatten())-meanlog)/stdlog,bins=nbins,range=(-5,5),normed=True)
        blog += (blog[1]-blog[0])/2.
        M.plot(blog,hlog,'g--',linewidth=2,label=r'$\log(1+\delta)$')


        M.plot(blog,M.exp(-blog**2/2.)/N.sqrt(2.*M.pi),'k:',linewidth=1.5,label=r'${\rm Gaussian}$')
        M.axis([-4,4,0,0.5])


        M.axes(axes12[1])
        M.plot(b,h,'b',linewidth=1.5)
        M.xticks([])
        M.yticks([1,2,3,4],('','','','4'))
        fontsize.mySetFontSize(ylabel=17)


        M.axis([-4,4,0.5,max(h)])

        print 'z=',str(redshifts[pn])
        #print 'Log skewness, kurtosis:',SS.skew(logd),SS.kurtosis(logd)
        #print 'skewness, kurtosis:',SS.skew(d.flatten()),SS.kurtosis(d.flatten())
        pn += 1
        
    

    M.axes(axes12[0])
    M.xlabel(r'$(x-\bar{x}>)/\sigma(x)$',fontsize=25)
    M.ylabel(r'$P[(x-\bar{x})/\sigma(x)]$',fontsize=25)
    M.legend(pad=0.02,labelspacing=0)
    fontsize.mySetFontSize(xlabel=17,ylabel=17,legend=17)

    #M.savefig('lognormality.eps')
    #M.savefig('lognormality.png')
    M.show()

def remechoes(d, nsigma=1.):
    ans = 1.*d - N.mean(d.flatten())
    sig = N.std(d.flatten())
    return ans * (N.abs(ans) > nsigma*sig)


def abovethresh(d,thresh=0.):
    return d*(d > 0.)

def smoothgauss(d,sigma=8.,boxsize=1200.):

#    sigmak = 2.*N.pi/sigma
    sigmak = 1./sigma
    kmin = 2.*N.pi/boxsize

    dk = N.fft.rfftn(d)
    s = d.shape
    sk = dk.shape

    if (len(s) == 3):
        a = N.fromfunction(lambda x,y,z:x, sk)
        a[N.where(a > s[0]/2)] -= s[0]
        b = N.fromfunction(lambda x,y,z:y, sk)
        b[N.where(b > s[1]/2)] -= s[1]
        c = N.fromfunction(lambda x,y,z:z, sk)
        c[N.where(c > s[2]/2)] -= s[2]

        k2 = kmin**2*(a**2+b**2+c**2)

    elif (len(s) == 2):
        b = N.fromfunction(lambda y,z:y, sk)
        b[N.where(b > s[0]/2)] -= s[0]
        c = N.fromfunction(lambda y,z:z, sk)
        #c[N.where(c > s[1]/2)] -= s[1]

        k2 = kmin**2*(b**2+c**2)
    elif (len(s) == 1):
        c = N.fromfunction(lambda z:z, sk)
        #c[N.where(c > s[0]/2)] -= s[0]
        
        k2 = kmin**2*c**2

    else:
        raise Exception('smoothgauss error: dimension>3')
        
#    gaussian = 1./(N.sqrt(2.*N.pi)*sigmak)*N.exp(-k2/sigmak**2/2.)
#    gaussian = 1./(N.sqrt(2.*N.pi))*N.exp(-k2/sigmak**2/2.)
    gaussian = N.exp(-k2/sigmak**2/2.)

    print dk.shape, gaussian.shape
    print N.sum(gaussian.flatten())


    return N.fft.irfftn(dk*gaussian)
    
def potential(d,boxsize=500.):

    kmin = 2.*N.pi/boxsize

    dk = N.fft.rfftn(d)
    s = d.shape
    sk = dk.shape
    if (len(s) == 3):
        a = N.fromfunction(lambda x,y,z:x, sk)
        a[N.where(a > s[0]/2)] -= s[0]
        b = N.fromfunction(lambda x,y,z:y, sk)
        b[N.where(b > s[1]/2)] -= s[1]
        c = N.fromfunction(lambda x,y,z:z, sk)
        c[N.where(c > s[2]/2)] -= s[2]

        k2 = kmin**2*(a**2+b**2+c**2)

    elif (len(s) == 2):
        b = N.fromfunction(lambda y,z:y, sk)
        b[N.where(b > s[0]/2)] -= s[0]
        c = N.fromfunction(lambda y,z:z, sk)
        c[N.where(c > s[1]/2)] -= s[1]

        k2 = kmin**2*(b**2+c**2)
        k2[0,0] = 1.

    elif (len(s) == 1):
        c = N.fromfunction(lambda z:z, sk)
        c[N.where(c > s[0]/2)] -= s[0]
        
        k2 = kmin**2*c**2
        
    return N.fft.irfftn(dk/k2**1.5)
    
