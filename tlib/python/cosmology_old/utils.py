""">>> utils.py <<<
 
   this is a loose collection of  utilities, mainly IO

   Current revision:

   ID:         $Id: utils.py 62 2008-08-21 21:57:49Z neyrinck $
   Date:       $Date: 2008-08-21 11:57:49 -1000 (Thu, 21 Aug 2008) $
   Revision:   $Revision: 62 $

(C) 2008 The CosmoPy Team (see Copyright for details)
   
"""
import numpy as N
import pylab as M
import scipy.signal as SS
import scipy.interpolate as SI
import unittest

def myFloat(str):
    """
    return the float value, or
    float('NaN')
    """
    try:
        return float(str)
    except ValueError:
        return float('NaN')
    
def readRaws(fname,ignore=0):
    """
      read and return the raws of a file

    optionally ignore the first ingnore lines
    comments with # allowed
    """

    data = []
    i = 0
    for line in file(fname):
         if line[0]=='#' or i < ignore:
             i += 1
             continue
         row = [myFloat(val) for val in line.split()]
         data.append(row)
    return M.array(data)


def readColumns(fname,ignore=0):
    """
    read columns of a file into a (newly created) array

    optionally ignore the first ingnore lines
    comments with # allowed
    """
         
    data = M.transpose(readRaws(fname,ignore))
    return data

def splineIntLinExt(y,x,xNew,splinecoeff = 0.):
    """
    Use the scipy spline interpolation, but linearly extrapolate at the edges,
    since scipy.signal.cspline1d assumes periodic boundary conditions
    """    
    if len(x) < 4:
        return interpolateLin(y,x,xNew)

    if isinstance(xNew,float):
        wasfloat = 1
        xNew = M.array([xNew])
    else:
        wasfloat = 0

    whereSpline = N.where((xNew > x[0]) * (xNew < x[-1]))[0]
    whereLin = N.where((xNew <= x[0]) + (xNew >= x[-1]))[0]
    
    ans = xNew * 0.
    if len(whereSpline) > 0:
        if isinstance(splinecoeff,float): # not pre-calculated.
            splinecoeff = SS.cspline1d(y)
        ans[whereSpline] = SS.cspline1d_eval(splinecoeff, xNew[whereSpline], dx=x[1]-x[0], x0 = x[0])

    if len(whereLin) > 0:
        ans[whereLin] = interpolateLin(y,x,xNew[whereLin])

    if wasfloat:
        return ans[0]
    else:
        return ans

def splineIntLinExt_LogLog(y,x,xNew):
    """
    Log-log polation.  Spline coefficients in log.
    """
    return M.exp(splineIntLinExt(M.log(y),M.log(x),M.log(xNew)))

def interpolateLin(y,x,xNew):
    """
    linear interpolation of y[x] onto y[xNew]
    Linearly extrapolates if outside range
    """
    xInd = M.clip(M.searchsorted(x,xNew)-1,0,len(x)-2)
    xFract = (xNew-x[xInd])/(x[xInd+1]-x[xInd])
    return y[xInd]+xFract*(y[xInd+1]-y[xInd])

def interpolateLinLog(y,x,xNew):
    """
    linear interpolation in LOG space of y[x] onto y[xNew]
    Linearly extrapolates if outside range
    """
    logx = M.log(x)
    logy = M.log(y)
    logxNew = M.log(xNew)
    
    logxInd = M.clip(M.searchsorted(logx,logxNew)-1,0,len(logx)-2)
    logxFract = (logxNew-logx[logxInd])/(logx[logxInd+1]-logx[logxInd])

    return M.exp(logy[logxInd]+logxFract*(logy[logxInd+1]-logy[logxInd]))

def trapz(y, x=None, ax=-1, method=M.add.reduce):
    """trapz(y,x=None,ax=-1) integrates y along the given dimension of
    the data array using the trapezoidal rule.

    you can call it with method=M.add.accumulate to yield partial sums
    """
    y = M.asarray(y)
    if x is None:
        d = 1.0
    else:
        d = M.diff(x,axis=ax)
    y = M.asarray(y)
    nd = len(y.shape)
    s1 = nd*[slice(None)]
    s2 = nd*[slice(None)]
    s1[ax] = slice(1,None)
    s2[ax] = slice(None,-1)

    ans = method(0.5* d * (y[s1]+y[s2]),ax)
    return ans


def trapezoidRule(f,dx):
    """
    numerical integral of a tabulated function f using trapezoid rule

    the function is assumed to be sampled at uniform intervals separated
    by dx between the limits of the integral

    4.1.11 Num.Rec., first ed.
    """
    return (M.sum(f) - 0.5*f[0]-0.5*f[-1])*dx

def trapezoidArray(f,dx):
    """
    same as trapezoidRule() but return the integral as fn of z
    (mostly in-place operations for optimized memory usage)
    """
    i = M.add.accumulate(f)
    M.add(i,-0.5*f,  i)
    M.add(i,-0.5*f[0],  i)
    M.multiply(i,dx, i)
    return i

def simpsonRule(f,dx):
    """
    numerical integral of a tabulated function f using Simpson's rule

    the function is assumed to be sampled at uniform intervals separated
    by dx between the limits of the integral

    4.1.14 Num.Rec., first ed.
    """
    return (M.sum(f) +(-31.0*f[0]+11.0*f[1]-5.0*f[2]+f[3]
                          -31.0*f[-1]+11.0*f[-2]-5.0*f[-3]+f[-4])/48.0)*dx

def openRomb(integrand, a, b,eps=1e-6,jmax=14,k=5):
    """
    Returns the integral on the _open_interval_ (a,b).
    Integration is performed by Romberg's method of order 2k,
    where, e.g., k=2 is Simpson's rule.
    """
    jmaxp=jmax+1

    s = 0.*M.zeros(jmaxp)
    h = 0.*M.zeros(jmaxp+1)
    ss = 0.
    dss = 0.

    h[0]=1.0
    for j in range(0,jmax):
        s[j]=tripleInt(integrand,a,b,s[j],j)
        if j >= k:
            ss,dss = interpPoly(h[j-k:j],s[j-k:j],k,0.0)
            if M.fabs(dss) <= eps*M.fabs(ss):
                return ss
        s[j+1]=s[j]
        h[j+1]=h[j]/9.

    print 'Non-convergence in openRomb'
    return ss

def interpPoly(xin,yin,n,x):
    """
    Polynomial interpolation of degree n.  yin = f(xin) is the interpolated function;
    the output is y = f(x) and an estimate of the error, dy.
    """

    nmax = 20

    c = 0.*M.zeros(nmax)
    d = 0.*M.zeros(nmax)
    ns=1
    dif=M.fabs(x-xin[1])
    for i in range(0,n):
        dift=M.fabs(x-xin[i])
        if dift < dif:
            ns=i
            dif=dift
        c[i]=yin[i]
        d[i]=yin[i]

    y=yin[ns]
    ns=ns-1
    for m in range(1,n):
        for i in range(0,n-m):
            ho=xin[i]-x
            hp=xin[i+m]-x
            w=c[i+1]-d[i]
            den=ho-hp
            if den == 0.:
                print 'interpPoly failed!'
                return y,dy
            den=w/den
            d[i]=hp*den
            c[i]=ho*den

        if (2*ns < n-m):
            dy=c[ns+1]
        else:
            dy=d[ns]
            ns=ns-1

        y=y+dy

    return y,dy

def adaptIntPlot(f,a,b):
    """
    Adaptive (doubling partition) integration.
    Minimizes function evaluations at the expense of some tedious
    array minipulations.
    """
    maxiter = 20
    miniter = 5
    tolerance = 0.1
    maxnx = 2**maxiter
    minnx = 2**miniter
    x = 0.*M.zeros(maxnx)
    
    dx = (b-a)/2.#**minsteps
    nx = 2
    x[0] = a
    x[1] = a+dx
    integral = M.sum(f(x[1:2]))*dx # 1 so we don't include the first endpt
    dx /= 2.
    newintegral = integral/2. + M.sum(f(x[:nx]+dx))*dx
    for i in range(nx-1,-1,-1):
        x[2*i] = x[i]
        x[2*i+1] = x[i] + dx
    nx *= 2
    keepgoing = 1
    while keepgoing == 1:
        integral = newintegral
        dx /= 2.
        eff = f(x[:nx]+dx)
        M.plot(x[:nx]+dx,f(x[:nx]+dx))
        newintegral = integral/2. + M.sum(eff)*dx#M.sum(f(x[:nx]+dx))*dx
        print newintegral*nx/(nx-1)
        for i in range(nx-1,-1,-1):
            x[2*i] = x[i]
            x[2*i+1] = x[i] + dx
        nx *= 2
        keepgoing = 0
        if integral*newintegral > 0.:
            if ((M.fabs(M.log(integral*(nx/2)/(nx/2-1)/(newintegral*nx/(nx-1)))) >\
                 tolerance) and (nx < maxnx/2)) or (nx < minnx):
                keepgoing = 1
        elif integral*newintegral == 0.:
            print "Hmmm, we have a zero integral here.  Assuming convergence."
        else:
            keepgoing = 1
    
    M.show()
    print nx,
    if nx == maxnx/2:
        print 'No convergence in utils.adaptInt!'
    return newintegral*nx/(nx-1)

def adaptInt(f,a,b,tolerance):
    """
    Adaptive (doubling partition) integration.
    Minimizes function evaluations at the expense of some tedious
    array minipulations.
    """
    maxiter = 25
    miniter = 5
    maxnx = 2**maxiter
    minnx = 2**miniter
    x = 0.*M.zeros(maxnx)
    
    dx = (b-a)/2.#**minsteps
    nx = 2
    x[0] = a
    x[1] = a+dx
    integral = M.sum(f(x[1:2]))*dx # 1 so we don't include the first endpt
    dx /= 2.
    newintegral = integral/2. + M.sum(f(x[:nx]+dx))*dx
    for i in range(nx-1,-1,-1):
        x[2*i] = x[i]
        x[2*i+1] = x[i] + dx
    nx *= 2
    keepgoing = 1
    while keepgoing == 1:
        integral = newintegral
        dx /= 2.
        eff = f(x[:nx]+dx)
        newintegral = integral/2. + M.sum(eff)*dx#M.sum(f(x[:nx]+dx))*dx
        for i in range(nx-1,-1,-1):
            x[2*i] = x[i]
            x[2*i+1] = x[i] + dx
        nx *= 2
        keepgoing = 0
        if integral*newintegral > 0.:
            if ((M.fabs(M.log(integral*(nx/2)/(nx/2-1)/(newintegral*nx/(nx-1)))) >\
                 tolerance) and (nx < maxnx/2)) or (nx < minnx):
                keepgoing = 1
        elif integral*newintegral == 0.:
            print "Hmmm, we have a zero integral here.  Assuming convergence."
        else:
            keepgoing = 1
    
    print nx,
    if nx == maxnx/2:
        print 'No convergence in utils.adaptInt!'
    return newintegral*nx/(nx-1)

def tripleInt(funct,a,b,s,n):
    """
    This computes an integral of funct by tripling the number
    of function evaluations at each step.
    """
    
    if n == 0:
        s=(b-a)*M.sum(funct(M.array([0.5*(a+b)])))
    else:
        it=3**(n-1)
        tnm=it
        dell3 = (b-a)/tnm
        dell=dell3/3.
        dell2 = dell*2.
        x=a+0.5*dell
        array = a+0.5*dell + dell3*M.arange(0.,it)
        summy = M.sum(funct(array) + funct(array+dell2))

        s=(s+(b-a)*summy/tnm)/3.

    return s

class LogInterp:
    """
    
    """

    def __init__(self,x,f):
        """

        """
        self.xa = x
        self.fa = f
        n = x.shape[0]
        self.n = n
        self.xMin = x[0]
        self.xMax = x[n-1]
        lnx = M.log(x)
        lnf = M.log(f)
        self.fSpl = SI.splrep(lnx,lnf,s=0)
        self.nAtxMin = (lnf[1]-lnf[0])/(lnx[1]-lnx[0])
        self.fAtxMin = f[0]
        self.nAtxMax = (lnf[n-1]-lnf[n-2])/(lnx[n-1]-lnx[n-2])
        self.fAtxMax = f[n-1]

    def f1(self,x):
        return self.fAtxMin*(x/self.xMin)**self.nAtxMin        

    def f2(self,x):
        if len(x)<1:
            return
        lnx = M.log(x)
        lnf = SI.splev(lnx,self.fSpl,der=0)
        return M.exp(lnf)
    
    def f3(self,x):
        return self.fAtxMax*(x/self.xMax)**self.nAtxMax

    def f(self,x):
        if N.isscalar(x):
            x = M.array([x])
        return N.piecewise(x,[x<=self.xMin,x>=self.xMax],[self.f1,self.f3,self.f2])



class LinInterp:
    def __init__(self,x,f):
        """
        
        """
        self.xa = x
        self.fa = f
        n = x.shape[0]
        self.n = n
        self.xMin = x[0]
        self.xMax = x[n-1]
        self.fSpl = SI.splrep(x,f,s=0)
        self.nAtxMin = (M.log(f[1])-M.log(f[0]))/(M.log(x[1])-M.log(x[0]))
        self.fAtxMin = f[0]
        self.nAtxMax = (M.log(abs(f[n-1]))-M.log(abs(f[n-2])))/(M.log(x[n-1])-M.log(x[n-2]))
        self.fAtxMax = f[n-1]

    def f1(self,x):
        return self.fAtxMin*(x/self.xMin)**self.nAtxMin        

    def f2(self,x):
        if len(x)<1:
            return
        return SI.splev(x,self.fSpl,der=0)
    
    def f3(self,x):
        return self.fAtxMax*(x/self.xMax)**self.nAtxMax

    def f(self,x):
        if N.isscalar(x):
            x = M.array([x])
        return N.piecewise(x,[x<=self.xMin,x>=self.xMax],[self.f1,self.f3,self.f2])


class integralTests(unittest.TestCase):
    """
    unit tests for the simple integral formulas
    """
    def testTrapezoid(self):
        x = M.arange(0,1+0.000005,0.00001)
        y = x*x
        self.failIf(trapezoidRule(y,0.00001) != 0.33333333335000059)

    def testSimpson(self):
        x = M.arange(0,1+0.000005,0.00001)
        y = x*x
        self.failIf(simpsonRule(y,0.00001) != 0.33333333333333387)

def cisiarr(xarr):
    """
    Cosine and sine integrals Ci(xarr) and Si(xarr) for arrays.
    Uses a Numerical Recipes algorithm.
    """

    EPS=6e-8
    EULER=.57721566
    MAXIT=101
    PIBY2=1.5707963
    FPMIN=1.e-30
    TMIN=2.

    ciarr = xarr*0.
    siarr = xarr*0.
    for el in range(len(xarr)):
        x = xarr[el]
        t=M.fabs(x)
        if t == 0.:
            si=0.
        else:
            if t > TMIN:
                b=1.+t*1j
                c=1./FPMIN
                d=1./b
                h=d
                for i in range(2,MAXIT):
                    a=-(i-1)**2
                    b += 2.
                    d=1./(a*d+b)
                    c=b+a/c
                    dell=c*d
                    h *= dell
                    if (M.fabs(dell.real - 1.) + M.fabs(dell.imag)) < EPS:
                        i = 0
                        break
                    
                if i > 0:
                    print 'Continued fraction method failed in cisiarr'
                        
                h *= (M.cos(t)-M.sin(t)*1j)
                ci=-h.real
                si=PIBY2+h.imag
            else:
                if t < M.sqrt(FPMIN):
                    sumc=0.
                    sums=t
                else:
                    summy=0.
                    sums=0.
                    sumc=0.
                    sign=1.
                    fact=1.
                    odd=1
                    for k in range(1,MAXIT):
                        fact *= t/k
                        term=fact/k
                        summy += sign*term
                        err=term/M.fabs(summy)
                        if odd:
                            sign=-sign
                            sums=summy
                            summy=sumc
                        else:
                            sumc=summy
                            summy=sums

                        if err < EPS:
                            k = 0
                            break
                        
                        odd=1-odd

                    if k > 0:
                        print "Too many iterations in cisiarr!"

                si=sums
                ci=sumc+M.log(t)+EULER

            if x < 0.:
                si=-si

        ciarr[el] = ci
        siarr[el] = si

    return ciarr,siarr


class SpecialFunctions:
    """
    mainly exact bessels

    """

    def j0(self,x):
        """
        spherical bessel l=0
        """
        return M.choose(M.greater(x,0.001),(1.0-x**2/6.0,M.sin(x)/x))

    def j1(self,x):
        """
        spherical bessel l=1
        """
        return M.choose(M.greater(x,0.001),(x/3.0,(M.sin(x)/x - M.cos(x))/x))

    def j2(self,x):
        """
        spherical bessel l=2
        """
        return M.choose(M.greater(x,0.001),(x**2/15.0,((M.cos(x) * -3. / x -
                     M.sin(x) * (1. - 3. /x**2))/x)))
        return 

class SpecialFunctionTests(unittest.TestCase):
    """
    unit tests for the SpecialFunction class
    """
    def testjs(self):
        sf = SpecialFunctions()
        
        self.failIf(sf.j0(0.1) != 0.99833416646828155 or
                    sf.j0(0.00001) != 0.99999999998333333  )

        self.failIf(sf.j1(0.1) != 0.033300011902557269 or
                    sf.j1(0.00001) != 3.3333333333333337e-06 )

        self.failIf(sf.j2(0.1) != 0.00066619060838490896 or
                    sf.j2(0.00001) != 6.6666666666666679e-12  )

        self.failIf(M.fabs(sf.j0(M.array([0.1,0.00001]))-
                    M.array([ 0.99833417,  1.        ])).max() > 10e-8)

        self.failIf(M.fabs(sf.j1(M.array([0.1,0.00001]))-
                    M.array([ 3.33000119e-02,3.33333333e-06])).max() > 10e-8)

        self.failIf(M.fabs(sf.j2(M.array([0.1,0.00001]))-
                    M.array([6.66190608e-04,6.66666667e-12 ])).max() > 10e-8)


def meanVar(data,outMean=None):
    """
    each data[i] is a vector, return the average vector and the variance
    vector
 
    each vector is assumed to be the same length of course
 
    if a mean is available from outside, set outMean which is a vector
    of len(data[i]); then the mean returned is a residual
    """

    mean = []
    var = []
    n = len(data)
    print n,len(data[0])
    try:
        for i in range(n):
            data[i] -= outMean
    except TypeError:
        pass
    for i in range(len(data[0])):
        m = 0.0
        v = 0.0
        for j in range(n):
            m += data[j][i]
            v += data[j][i]*data[j][i]
        if outMean == None:
            v -= m*m/float(n)
        m /= float(n)
        v /= (n-1)
        mean.append(m)
        var.append(v)
    return (M.array(mean),M.array(var))

class HGQ:
    """
    Hermite-Gauss Quadrature (Gaussian weight function)
    accepts number of abscissas
    """
    def __init__(self, nabscissas):
        if nabscissas == 1:
            self.abscissas = M.array([0.])
            self.weights = M.array([1.])
        elif nabscissas == 3:
            self.abscissas = M.array([-0.5, 0., 0.5])*M.sqrt(6.)
            self.weights = M.array([1./6., 2./3., 1./6.])
        elif nabscissas == 5:
            self.abscissas = M.array([-2.02018,-0.958572,0.,0.958572,2.02018])
            self.weights = M.array([0.0199532,0.393619,0.945309,0.393619,0.0199532])
            self.weights /= M.sum(self.weights)
        else:
            print 'Invalid number of abscissas (1,3, or 5).'

    def int(self, integrand):
        """
        Integrates over second argument of an array
        with Gaussian Quadrature weights
        """
        return M.sum(M.outer(1.*M.ones(integrand.shape[0]),self.weights) * \
                   integrand,1)

def main():
    unittest.main()
if __name__=='__main__' :
    main()
