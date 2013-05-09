""" >>>  Parameters Classes <<< """
import sys
import time
import abc
import ConfigParser
import numpy as np
import pylab as pl

import myutlis as mut


class ParamDict(dict):
    """ 
        >>> Base class of parameter dictionary <<< 
    """
    def __init__(self, valuedict, **pardict):
        super(ParamDict,self).__init__(valuedict)
        for i in pardict:
            if i in self: self[i] = pardict[i]
            else: raise Exception('ParamDict failed')

    def __getattr__(self,name):
        return self[name]

    def __str__(self):
        r = 'parameters in ParamDict:\n'
        for i in self:
            r += "\t -> %s = %s\n" % (i, str(self[i]))
        return r[:-1]



class Parameters(object):
    """ 
    >>> Abstract Base Class for cosmological parameters, survey parameters etc. <<<
    >>> If the paramfname & section != None, then update from given 
	parameter file.  If def_var==True, then define all variables from 
	the dictionary.  
    >>> NOTE that when update from dictionary--'pardict' get a higher
        priority than 'parameter file' ***              
    >>> Get the update from parameter file if exists. 
        However, NOTE that the 'defulat_values' & 'names_dict' SHOULD always
	be the superset of the updated dictionary from file.  
    >>> paradict will always be a 'varible_dict', it is updated to 
        varible_dict.
    """

    __metaclass__ = abc.ABCMeta  # Enforce Abstract class

    def __init__(self, values_dict, names_dict=None, paramfname=None, section=None, 
                 def_var=True, **pardict):
	"""
	values_dict: dictionary of variable values, keys are long `human-readble` 
	             names.  
	names_dict:  dictionary of variable names, values are variable names 
	             used in code.
	def_var:     if True, then define all variables from dictionary.
	paramfname:  the name of parameter file where all parameters could be updated
	pardict:     specify the parameter dictionary need to be updated.
	"""

	self.values_dict=values_dict.copy()
	if(names_dict!=None):
	    self.names_dict=names_dict.copy()
	else:
	    self.names_dict=names_dict  # just get None

        """ Update from parameter files, only existing default keys
	    NO NEW keys will be added.  
	"""
        self.values_dict=self.update_values_dict_fromfile(self.values_dict, paramfname, section)

        """ convert `values_dict` + `names_dict` to varibles_dict """
	vardict=self.get_varibledict(self.values_dict, self.names_dict)

        """ update """
        vardict=self.update_dict_values(vardict, **pardict)

	""" define the parameter dictionary """
	self.paramdict=ParamDict(vardict)

        """ >>> convert dictionary into varibles <<< """
        if(def_var):
            self.def_varibles(vardict)



    def get_varibledict(self, values_dict, names_dict):
        """ 
	    >>> Return the varialbe dictionary from values_dict & names_dict
	"""
        if(names_dict!=None):
            """ >>> get the variable dictionary, i.e. (variable_name, values) pairs <<< """
            vardict=dict((v, k) for k, v in names_dict.items() )
	    for i in range(len(vardict.values())):
	        vardict[vardict.keys()[i]]=values_dict[vardict.values()[i]]
	else:
            vardict=values_dict

        return vardict


        """ >>>  Upadate  methods  <<<"""
    def add_new_everything(self, update_values, update_names, paramfname=None, section=None, 
                   def_var=True, **pardict):
        """   This is useful since one can define a subclass without knowing
	        lots of details of the super class.  """
        # update from files, again, NO NEW keys
        update_values=self.update_values_dict_fromfile(update_values, paramfname, section)

        # update values & names here, so NEW keys CAN be added. 
	# Useful for subclass
	self.values_dict.update(update_values)
	if(self.names_dict!=None):
	    self.names_dict.update(update_names)

        #update paramdict
	varupdate=self.get_varibledict(update_values, update_names)
	self.paramdict.update(varupdate)

        # further update from arguments
	# self.paramdict.update(**pardict)  # also add new keys
        self.paramdict=self.update_dict_values(self.paramdict, **pardict)

        if(def_var):
            self.def_varibles(varupdate)



    def update_values_dict_fromfile(self, values_dict, paramfname, section):
        """ update values_dict from parameters file.
	"""
        if((paramfname!=None) & (section!=None) ):
	    fp=open(paramfname)
            parser=ConfigParser.ConfigParser()
            parser.readfp(fp)

            """ >>>  update ONLY the DEFAULT KEYS from parameter file <<< """
	    update=mut.convert_strdict(dict(parser.items(section)) )
            newvalues_dict=self.update_dict_values(values_dict, **update)
            #self.values_dict.update(**update)

	    fp.close()
	else:
            newvalues_dict=values_dict

	return newvalues_dict


    def update_dict_values(self, d, **update):
        """ Updating the values dictionary, ONLY for default EXISTing keys. """
        for key in update:
            if key in d:
                d[key]=update[key]
	return d



    def update_params(self, **update):
        ''' Update both paramdict as well as varialbes
	'''
        self.paramdict.update(**update)
        self.def_varibles(self.paramdict)
	self.get_derived()

    def def_varibles(self, values_dict):
	for i in range(len(values_dict)):
            exec('self.'+values_dict.keys()[i]+'='+str(mut.add_quotation(values_dict.values()[i])))



    """ >>> Abstract properties <<< """
    @abc.abstractmethod
    def get_derived(self):
        """ Initialize other parameters """
        pass



