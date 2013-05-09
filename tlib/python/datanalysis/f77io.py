import numpy as N
import array as A
import struct as S
import os, sys
import pylab as M

def getdatacube(ngrid = 256, prefix='', file='', galtrans=False,ndim=3):

    F = open(prefix+file,'rb')
    f77dummy = S.unpack('i',F.read(4))[0]
    ng3 = S.unpack('i',F.read(4))[0]
    f77dummy = S.unpack('i',F.read(4))[0]
    f77dummy = S.unpack('i',F.read(4))[0]

    print 'ng3=',ng3

    den = A.array('f')
    den.fromfile(F,ng3)
    if ndim == 3:
        ngrid = int(round(ng3**(1./3.)))
        den = N.array(den).reshape((ngrid,ngrid,ngrid))#.astype(N.float32)
    elif ndim == 2:
        ngrid = int(round(ng3**(1./2.)))
        den = N.array(den).reshape((ngrid,ngrid))#.astype(N.float32)
    
    F.close()

    if(galtrans):
        den = den.transpose(1,2,0)

    return den

def putdatacube(den, prefix, file, dtype='float32'):

    lenden = len(den.flatten())
    F = open(prefix+file,'w')
    f77dummy = N.array([4L,lenden,4L,lenden*4]).astype('int32')
    #f77dummy = N.array([4L]).astype('int32')
    print f77dummy
    f77dummy.tofile(F)
    F.flush()
    #os.system('od -i '+prefix+file)
    (((den.flatten()).astype(dtype)).squeeze()).tofile(F)
    #os.system('od -f '+prefix+file+' | head')
    f77dummy = N.array([4L*lenden]).astype('int32')
    f77dummy.tofile(F)
    print f77dummy
    F.close()

def readfield(file,typestring,arraydim=0,bytes=4,dummybytes=4,arrayshape=(),dodum=False,save=True,byteswap=False):
    if dodum:
        f77dummy = S.unpack('i',file.read(4))[0]
        if dummybytes == 8:
            f77dummy = S.unpack('i',file.read(4))[0]
        print f77dummy,
    if (arraydim > 0):
        ans = A.array(typestring)
        ans.fromfile(file,arraydim)
        if (byteswap):
            ans.byteswap()
        if (len(arrayshape) == 0):
            if (typestring == 'f'):
                ans = N.array(ans,dtype=N.float32)
            else:
                ans = N.array(ans)
        else:
            if (typestring == 'f'):
                ans = N.array(ans,dtype=N.float32).reshape(arrayshape)
            else:
                N.array(ans,dtype=N.float32).reshape(arrayshape)

    else:
        ans = S.unpack(typestring,file.read(bytes))[0]
    if dodum:
        f77dummy = S.unpack('i',file.read(4))[0]
        if dummybytes == 8:
            f77dummy = S.unpack('i',file.read(4))[0]
        print f77dummy
    if save:
        return ans

def readpos(file='MMF200M.GAD'):

    F = open(file,'rb')

    npart = readfield(F,'i',arraydim=1,bytes=4,dodum=False)[0]
    print npart
    pos = N.zeros((npart,3))
    pos[:,0] = readfield(F,'f',arraydim=npart,dodum=False)
    pos[:,1] = readfield(F,'f',arraydim=npart,dodum=False)
    pos[:,2] = readfield(F,'f',arraydim=npart,dodum=False)
    return pos

def read2lptcrappos(file):

    F = open(file,'rb')

    npart = readfield(F,'i',arraydim=1,bytes=4,dodum=True,dummybytes=8)[0]
    print npart
    pos = N.zeros((npart,3))
    pos[:,0] = readfield(F,'f',arraydim=npart,dodum=True,dummybytes=8)
    pos[:,1] = readfield(F,'f',arraydim=npart,dodum=True,dummybytes=8)
    pos[:,2] = readfield(F,'f',arraydim=npart,dodum=True,dummybytes=8)
    return pos


def readgrafic(prefix, file, read_vel=False,read_id=False,debug=False):
    F = open(prefix+file,'rb')

    header = readfield(F,'i',arraydim=11,bytes=4,dodum=True)
    ng = header[0]
    cube = N.zeros((ng,ng,ng))
    for i in range(ng):
        cube[:,:,0]=readfield(F,'f',arraydim=ng**2,bytes=4,dodum=True).reshape(ng,ng)
    F.close()
    return cube


def readgadget(prefix, file, read_vel=False,read_id=False,debug=False):

    F = open(prefix+file,'rb')

    f77dummy = S.unpack('i',F.read(4))[0]
    if (debug):
        print f77dummy
    npart = readfield(F,'i',arraydim=6,bytes=4,dodum=False)
    massarr = readfield(F,'d',arraydim=6,bytes=8,dodum=False)
    time = S.unpack('d',F.read(8))[0]
    redshift = S.unpack('d',F.read(8))[0]
    flag_sfr = S.unpack('i',F.read(4))[0]
    flag_feedback = S.unpack('i',F.read(4))[0]
    npartTotal = readfield(F,'i',arraydim=6,bytes=4,dodum=False)

    bytesleft=f77dummy - 6*4 - 6*8 - 8 - 8 - 2*4-6*4
    la = readfield(F,'i',arraydim=bytesleft/4,dodum=False,bytes=4)
    f77dummy = S.unpack('i',F.read(4))[0]

    if (debug):
        print npart
        print massarr
        print time,redshift,flag_sfr,flag_feedback
        print npartTotal

    pos = readfield(F,'f',arraydim=3*npart[1],dodum=True,arrayshape=(npart[1],3))
    if (read_vel):
        pos = readfield(F,'f',arraydim=3*npart[1],dodum=True,arrayshape=(npart[1],3))
    else:
        readfield(F,'f',arraydim=3*npart[1],dodum=True,save=False)

    id = readfield(F,'i',arraydim=npart[1],dodum=True)
    print max(id), min(id)
    pos_reord = N.empty(shape=pos.shape,dtype=N.float32)
    pos_reord[id-1,:] = 1.*pos

    return pos_reord

def readgadgetmad(prefix, file, read_vel=False,debug=True,nfiles=1,byteswap=False,partlabel=0):
    #if nfiles==1:
    #    filenames=[file]
    #else:
    #    filenames=file+'.'+N.arange(nfiles).astype(string)
    pos_reord = None
    for i in N.arange(nfiles):
        print 'FILE#: ',i
        sys.stdout.flush()
        if (nfiles == 1):
            filename = file
        else:
            filename = file+'.'+str(i)

        F = open(prefix+filename,'rb')
        dummy = readfield(F,'i',arraydim=2,bytes=4,dodum=False,byteswap=byteswap)
        #f77dummy = S.unpack('>i',F.read(4))[0]
        if (debug):
            print 'dummy=',dummy
        npart = readfield(F,'i',arraydim=6,bytes=4,dodum=False,byteswap=byteswap)
        massarr = readfield(F,'d',arraydim=6,bytes=8,dodum=False,byteswap=byteswap)
        time = S.unpack('>d',F.read(8))[0]
        redshift = S.unpack('>d',F.read(8))[0]
        flag_sfr = S.unpack('>i',F.read(4))[0]
        flag_feedback = S.unpack('>i',F.read(4))[0]
        npartTotal = readfield(F,'i',arraydim=6,bytes=4,dodum=False,byteswap=byteswap)

        if (debug):
            print 'debug info:'
            print 'npart=',npart
            print 'massarr=',massarr
            print 'crap=',time,redshift,flag_sfr,flag_feedback
            print 'npartTotal=',npartTotal
            print N.sum(npart)

        if pos_reord == None:
            pos_reord = N.empty((npartTotal[partlabel],3),dtype=N.float32)
            print 'allocating pos_reord'

        sys.stdout.flush()

        bytesleft=256 - 6*4 - 6*8 - 8 - 8 - 2*4-6*4
        la = readfield(F,'i',arraydim=bytesleft/4,dodum=False,bytes=4,byteswap=byteswap)
        f77dummy1 = S.unpack('>i',F.read(4))[0]
        #f77dummy2 = S.unpack('>i',F.read(4))[0]
        #f77dummy3 = S.unpack('>i',F.read(4))[0]

        if (debug):
            print 'before pos:',la, f77dummy1#,f77dummy2,f77dummy3

        if (read_vel):
            readfield(F,'f',arraydim=3*npart[partlabel],dodum=True,save=False)
            pos = readfield(F,'f',arraydim=3*npart[partlabel],dodum=False,arrayshape=(npart[partlabel],3),byteswap=byteswap)
        else:
            pos = readfield(F,'f',arraydim=3*npart[0],dodum=False,arrayshape=(npart[0],3),byteswap=byteswap)

            if (debug):
                print 'pos read:',pos.shape
                print min(pos[:,0]),max(pos[:,0])
                print min(pos[:,1]),max(pos[:,1])
                print min(pos[:,2]),max(pos[:,2])
                print pos
            dum1=S.unpack('>i',F.read(4))[0]

            if (debug):
                vel=readfield(F,'f',arraydim=3*npart[0],dodum=True,byteswap=byteswap)
                print vel,min(vel),max(vel)
            else:
                readfield(F,'f',arraydim=3*npart[0],dodum=True,byteswap=byteswap,save=False)
            

        id = readfield(F,'i',arraydim=npart[0],dodum=True,byteswap=byteswap)
        print 'id:',id
        print 'min,maxid:',min(id),max(id)
        sys.stdout.flush()
        print pos_reord[id-1,:].shape, pos.shape
        pos_reord[id-1,:] = 1.*pos

    return pos_reord

def writegrid(den, prefix, file):

    F = open(prefix+file,'w')
    (((den.flatten()).astype('float32')).squeeze()).tofile(F)
    F.close()

def readgrid(shape=(286,286), prefix='', file=''):

    F = open(prefix+file,'r')
    den = A.array('f')
    den.fromfile(F,N.prod(shape))
    F.close()
    return(N.array(den).reshape(shape).astype(N.float32))

def vobozout(pos,filename,f77=False):
    F = open(filename,'w')

    
    np = len(pos[:,0])
    print 'np=',np
    if (f77):
        (4*N.array([1],dtype=N.int32)).tofile(F)
    N.array([np],dtype=N.int32).tofile(F)
    if (f77):
        (4*N.array([1,np])).tofile(F)
    (pos[:,0]).astype(N.float32).tofile(F)
    if (f77):
        (4*N.array([np,np])).tofile(F)
    (pos[:,1]).astype(N.float32).tofile(F)
    if (f77):
        (4*N.array([np,np])).tofile(F)
    (pos[:,2]).astype(N.float32).tofile(F)
    if (f77):
        (4*N.array([np])).tofile(F)

    F.close()

    #os.system('od -i '+filename+' | head -1')
    #os.system('od -f '+filename+' | head')

def vobozin(filename,dodum=False,dchar='f'):
    F = open(filename,'r')

    if (dodum):
        dummy = S.unpack('i',F.read(4))[0]
        dummy = S.unpack('i',F.read(4))[0]
    np = S.unpack('i',F.read(4))[0]
    if (dodum):
        dummy = S.unpack('i',F.read(4))[0]
        dummy = S.unpack('i',F.read(4))[0]

    print 'np=',np
    pos = N.empty((np,3),dtype=N.float32)
    pos[:,0] = readfield(F,dchar,arraydim=np,dodum=dodum,dummybytes=8)
    pos[:,1] = readfield(F,dchar,arraydim=np,dodum=dodum,dummybytes=8)
    pos[:,2] = readfield(F,dchar,arraydim=np,dodum=dodum,dummybytes=8)

    return pos

def vozisol_ascii(p,npreal=6619136,prefix='',root='lg2048_dm2',randseed=12345):
    np = len(p[:,0])

    #min0 = min(p[:,0])
    #min1 = min(p[:,1])
    #min2 = min(p[:,2])
    #max0 = max(p[:,0])
    #max1 = max(p[:,1])
    #max2 = max(p[:,2])
    #maxrange=max([max0-min0,max1-min1,max2-min2])

    #print min0,min1,min2,maxrange

    rs = N.random.RandomState(randseed)

    p1=p#/1e4
    #if (randseed != None):
    #    p1 += rs.rand(np,3)/2.**20
    #p1[:,0] = (p[:,0]-min0)/maxrange
    #p1[:,1] = (p[:,1]-min1)/maxrange
    #p1[:,2] = (p[:,2]-min2)/maxrange

    F=open(prefix+root+'.pos','w')
    F.write('%d %d\n'%(np, npreal))
    for i in range(np):
        F.write('%10.8f %10.8f %10.8f\n'%(p1[i,0],p1[i,1],p1[i,2]))

    F.close()

def makegridgad(katrin='/panfs/scratch2/vol9/heitmann/Grid2/',sim=0,pref2='',filepref='Gadget_1.0000.',nfiles=32,byteswap=False,ng=256,out=None):

    if out == None:
        out = 's%03d.%d.grid'%(sim,ng)

    pref1='M%03d'%sim
    path=katrin+pref1+'/'+pref2+'/'

    littleh = M.load('lanparam2.txt')[sim,3]
    boxsize=1.3e3*littleh
    print 'boxsize=',boxsize

    d=N.zeros((ng,ng,ng),dtype=N.int32)
    for i in range(nfiles):
    #for i in range(1):
        print i
        sys.stdout.flush()
        pos=readgadget(prefix='',file=path+filepref+str(i))
        h,edges = N.histogramdd(pos,bins=(ng,ng,ng),\
                           range=((0.,boxsize),(0.,boxsize),(0.,boxsize)))
        d += h
        mean = N.mean(d.astype(N.float64).flatten())
        print 'mean = ',mean

    d = d.astype('float32')/mean
        
    putdatacube(d,prefix='/scratch2/neyrinck/'+pref1+'/',file=out)
    return d

def readid(file,lonid=False):
    F = open(file,'r')

    np = S.unpack('i',F.read(4))[0]
    print np
    print 'que?'
    if lonid:
        id = N.empty(np,dtype=N.uint64)
        id = readfield(F,'Q',dodum=False,arraydim=np)
    else:
        id = N.empty(np,dtype=N.uint32)
        id = readfield(F,'i',dodum=False,arraydim=np)
    F.close()

    return id

def readvectors(ngrid = 128, prefix='',file='', galtrans=False,ndim=3):

    F = open(prefix+file,'rb')
    ng = S.unpack('i',F.read(4))[0]
    ng3 = ng**3
    print 'ng3=',ng3

    vec = A.array('f')
    vec.fromfile(F,ng3*3)
    F.close()

    return vec.reshape(ng,ng,ng,3)

def makeasciipos(p,filename="gaddy.pos",npreal=None):

    p=p.astype(N.float32)
    F=open(filename,'w')
    np = len(p[:,0])
    if npreal == None:
        npreal=1*np

    F.write(str(np)+' '+str(npreal)+'\n')

    for i in range(np):
        F.write(str(p[i,0])+' '+str(p[i,1])+' '+str(p[i,2])+'\n')

    F.close()
    
