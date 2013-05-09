import numpy as np

from genscript.extendclass import *
import genscript.progcontrol as pc
import genscript.parameter as par







''' ->> !!!the multiple inheritance of class `Parameters` has some problem, fix it.
'''

class First(par.Parameters):
    def __init__():
        pass


class Second(par.Parameters):
    def __init__():
        pass



class User(First, Second):
    def __init__():
        super(User, self).__init__()



if __name__=='__main__':
    p=pc.prog_init()



    p.finalize()
