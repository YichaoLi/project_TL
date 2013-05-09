"""
   Langrangian correlaion funtion

   (c) 2012 Istvan Szapudi
"""

import numpy as N
import scipy, scipy.spatial
import f77io

dataDir = '/data/nb64/'
posFName = dataDir+'z0.pos'
volFName= dataDir+'z0.vol'

def readData(pf = posFName, vf = volFName, nMax = 1000, scale=200):
    """
    read nMax data points and the corresponding volumes
    normalize volumes

    assume that volume and positions scale the same way
    """
    pos = f77io.vobozin(pf)
    vol = f77io.readvol(vf)
    v = sum(vol)
#    print len(vol), len(pos)
    ntot = len(vol) #not checked if same, but assumed
    if nMax > 0:
        inds = N.random.permutation(range(len(vol)))[:nMax]
        pos = N.take(pos,inds,axis=0)
        vol = N.take(vol,inds)
        
    
    vol = vol*v/ntot
    pos = pos*scale
    print sum(vol), scale**3
    return pos, vol

def bruteForceLCF(pos,vol,rMax=20,bNum=20):
    """
    Langrangian estimator corresponding to DD/RR-1 and (DD-2DR)/RR+1
    brute force algorithm w/o kdtree
    
    20x faster than _p
    """

    dd = N.zeros(bNum)
    dv = N.zeros(bNum)
    vv = N.zeros(bNum)

    i = 0
    for p in pos:
        rs = N.sqrt(N.sum((pos[i+1:]-p)**2,axis=1))
        bs = (N.floor(rs*bNum/rMax)).astype(int)
#       2x slower than python loop
#        for b in range(bNum):
#            dd[b] = N.sum(N.where(bs==b,1,0))
#            dv[b] = N.sum(N.where(bs==b,vol[i+1:]+vol[i],0))
#            vv[b] = N.sum(N.where(bs==b,vol[i+1:]*vol[i],0))
        j = 0
        for b in bs:
            if b < bNum:
                dd[b] = dd[b]+1.0
                vv[b] = vv[b]+vol[i]*vol[i+1+j]
                dv[b] = dv[b]+vol[i]+vol[i+1+j] # dv+vd
            j = j+1
        i = i+1
        print i

    return dd,dv,vv,dd/vv-1.0, (dd-dv)/vv+1.0


def bruteForceLCF_p(pos,vol,rMax=20,bNum=20):
    """
    Langrangian estimator corresponding to DD/RR-1 and (DD-2DR)/RR+1
    brute force algorithm w/o kdtree

    """

    dd = N.zeros(bNum)
    dv = N.zeros(bNum)
    vv = N.zeros(bNum)

    np = len(vol)
    for i in range(np):
        for j in range(i+1,np):
            r = N.sqrt(N.sum((pos[i]-pos[j])**2))
            b = int(N.floor(r*bNum/rMax))
            if b < bNum:
                dd[b] = dd[b]+1.0
                vv[b] = vv[b]+vol[i]*vol[j]
                dv[b] = dv[b]+vol[i]+vol[j] # dv+vd
#        print i

    return dd,dv,vv,dd/vv-1.0, (dd-dv)/vv+1.0

def calcSimpleLCF(pos,vol,rMax=20,bNum=20):
    """
    Langrangian estimator corresponding to DD/RR-1 and (DD-2DR)/RR+1

    final result needs to be multiplied with

    dr 4pi/3(L/l)**3 /r**2

    L = true sim. size (Mpc h-1)
    l = sim. size in units that Langrangian density is 1
    dr bin width
    r bin size (Mpc h^-1)

    """
    kdt = scipy.spatial.KDTree(pos)
    pairs = kdt.query_pairs(rMax)
    print len(pairs)
#   work around kdtree crashes on nangaku
#    n = len(vol)
#    pairs = [(i,j) for i in range(n) for j in range(n)]


    dd = N.zeros(bNum)
    dv = N.zeros(bNum)
    vv = N.zeros(bNum)
    dvi = N.zeros(bNum)
    vivi = N.zeros(bNum)

    for (i,j) in pairs:
         r = N.sqrt(N.sum((pos[i]-pos[j])**2))
         b = int(N.floor(r*bNum/rMax))
         if b < bNum:
             dd[b] = dd[b]+1.0
             vv[b] = vv[b]+vol[i]*vol[j]
             dv[b] = dv[b]+vol[i]+vol[j] # dv+vd
             vivi[b] = vivi[b]+1/(vol[i]*vol[j])
             dvi[b] = dvi[b]+1/vol[i]+1/vol[j] # dv+vd
             
    return dd,dv,vv,dd/vv-1.0, (dd-dv)/vv+1.0,dvi,vivi,vivi/dd-1,(vivi-dvi)/dd+1.0

def loadCF(fRoot='nb64_brute'):
    dd = N.loadtxt(fRoot+'.dd')
    dv = N.loadtxt(fRoot+'.dv')
    vv = N.loadtxt(fRoot+'.vv')
    cf = N.loadtxt(fRoot+'.cf')
    cfls = N.loadtxt(fRoot+'.cfls')
    return dd,dv,vv,cf,cfls


def full(nMax=-1):
    outFile='nb64_brute'
    pos,vol=readData(nMax=nMax)
    dd,dv,vv,cf,cfls=bruteForceLCF(pos,vol,rMax=120,bNum=120)
#    dd,dv,vv,cf,cfls=calcSimpleLCF(pos,vol,rMax=20,bNum=20)
    N.savetxt(outFile+'.dd',dd)
    N.savetxt(outFile+'.dv',dv)
    N.savetxt(outFile+'.vv',vv)
    N.savetxt(outFile+'.cf',cf)
    N.savetxt(outFile+'.cfls',cfls)

if __name__=='__main__':
    full()
