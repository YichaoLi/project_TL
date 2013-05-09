import numpy as np
import genscript.cyth.cubicspline as cs



class CMBPowerSpectrum(object):
    def __init__(self, cp=None, data=None, data_y=None):

        if cp!=None: #
	    self.cp=cp


        if data!=None:
            """ Initialize Powerspectrum interpolator """
            if isinstance(data, str):
                dat = np.array(zip(* np.genfromtxt(data) ))   
                l, c=dat[0,:], dat[1,:]


            elif isinstance(data, np.ndarray):
                if len(data.shape)==2:
                   l, c=data[0,:], data[1,:]
                elif len(data.shape)==1:
                   l, c=data, data_y
                else:
                    raise Exception()

            else:
                raise Exception()

            self.cl=cs.Interpolater(l, c)




    def __call__(self, l):
        """ return initialized power spectrum from `data` """
        return self.evaluate(l)

   
    """ interpolate power spectrum """
    def evaluate(self, l):
        return self.cl(l)







