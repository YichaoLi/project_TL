import sys, os
import numpy as np


from genscript.extendclass import *
import genscript.progcontrol as pc







''' ------ >>> run the script <<< ------ '''
pcontrol_dict={
    'bao_fisher': False, 
    'cmb_fisher': False 
    }
param_dict={
    'telpara_fname': 'parameter.cfg',
    'telpara_sec': 'Cylinder',
    'fishpara_fname': 'parameter.cfg',
    'fishpara_sec':   'Fisher_Matrix'
    }
cosmoparam_dict={
    'power_spectrum_fname': None
    }

if __name__=='__main__':

    # ->> Initialization <<- #
    sec='Fisher_Matrix_General'
    init_dict = myDict(pcontrol_dict)+myDict(param_dict)+myDict(cosmoparam_dict)
    p=pc.prog_init(section=sec, **init_dict)


    # -> 
    print p.paramdict



    ''' ->> Finalize <<- '''
    p.finalize()



