import genscript.parameter as par
import genscript.myutlis as mut

"""
>>> NOTE: One cannot automatically generate C variables with names from given string.
    NOT SURE what to do yet.
                          >>> !!!  NOT Working ... !!!  <<<
"""

def typeconvert(vlist):
    typelist=[]

    for i in range(len(vlist)):
        if isinstance(vlist[i], int):
            typelist.append('int')

        elif isinstance(vlist[i], float):
            typelist.append('double')

        elif isinstance(vlist[i], str):
            typelist.append('char *')

        elif isinstance(vlist[i], list):
            if isinstance(vlist[i][0], int ):
                typelist.append('int *')
            elif isinstance(vlist[i][0], float):
                typelist.append('double *')
        else:
            raise Exception()

    return typelist



cdef class cythonParameter:
    def __init__(self, vardict):
        """ define the cython version of variables, given value dictionary and types.
        """

        keys_list, values_list = [], []
        for i in range(len(vardict)):
            keys_list.append(vardict.keys()[i])
            values_list.append(vardict.values()[i])

        types_list=typeconvert(values_list) 

 
        for i in range(len(vardict)):
            if not isinstance(values_list[i], list):
                exec('cdef '+types_list[i]+' self.'+keys_list[i]+'='+str(mut.add_quotation(vardict.values()[i])))
            else:
                exec('cdef '+types_list[i]+' self.'+keys_list[i])
                exec('self.'+keys_list[i]+'='+str(values_list[i]))

    

