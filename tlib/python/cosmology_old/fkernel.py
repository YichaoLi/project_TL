"""
   Perturbation theory kernels for test
"""

import pylab as M
import cex

class DivergentPTKernelError(cex.Error):
    """
    exception for when the kernel is divergent due to a too
    small vector in alpha or beta
    """
    pass


class PTKernels:
    """
    recursive definition of perturbation theory kernels 
    for tests (possibly slow)
    """

    def __init__(self,eps = 0.0001):
	"""
	set parameters
	"""
	self.eps = eps # tolerance for small vectors

    def alpha(self,k1,k2):
	"""
	alpha coupling function

	k1,k2 are numarray vectors
	"""
	if M.sum(k1*k1) < self.eps or M.sum(k2*k2) < self.eps:
	    raise DivergentPTKernelError
	return M.sum((k1+k2)*k1)/M.sum(k1*k1)
	

    def beta(self,k1,k2):
	"""
	alpha coupling function

	k1,k2 are numarray vectors
	"""
	if M.sum(k1*k1) < self.eps or M.sum(k2*k2) < self.eps:
	    raise DivergentPTKernelError
	return (M.sum((k1+k2)*(k1+k2)))*M.sum((k1*k2))/M.sum(k1*k1)/M.sum(k2*k2)/2.
	

    def FG(self,kList):
	"""
	recursive function for the F and G kernels
	
	kList is a List of numarray vectors
	"""

	if len(kList) == 1:
	    return (1.0,1.0)

	sgg = 0.0
	sgf = 0.0
	n = len(kList)
	for m in range(1,n):
	    k1 = M.sum(kList[:m])
	    k2 = M.sum(kList[m:])
	    k1Sqr = M.sum(k1*k1)
	    k2Sqr = M.sum(k2*k2)
	    print k1,k2,k1Sqr,k2Sqr
	    if k1Sqr > self.eps and k2Sqr > self.eps:
		locF1,locG1 =  self.FG(kList[:m])
		locF2,locG2 =  self.FG(kList[m:])
		sgg += locG1*locG2*self.beta(k1,k2)
		sgf += locG1*locF2*self.alpha(k1,k2)

	return (((2*n+1.)*sgf+2*sgg)/(2*n+3.)/(n-1.),
		(3*sgf+2*n*sgg)/(2*n+3.)/(n-1.))

