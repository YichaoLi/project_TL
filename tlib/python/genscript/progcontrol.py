""" >>> program control <<<
"""
import os, sys, inspect
import time
import ConfigParser
import numpy as np
import parameter as par
import genclass as gc

import cosmology.cosparameter as cospar
import cosmology.frw 
import cosmology.power 

import mpiutil as mpi

defaultProgramControlDict={
    'debug': True,
    'debug_level': 10,
    'mpi': False
    }



class prog(par.Parameters):
    """ >>> program control class
    """

    def __init__(self, paramfname, section, varname=None, **pardict):

        # start clocking
	self.start_time=time.time()

	# debug informations
	#section='PROGRAM_CONTROL'

        # update and include extra value dictionary from `pardict`
	value_dict=defaultProgramControlDict
	value_dict.update(pardict)

        super(prog,self).__init__(value_dict, names_dict=None, 
                     paramfname=paramfname, section=section, def_var=True)


        self.paramfname=paramfname
	self.section=section

        # ini file reader
	if paramfname!=None:
            self.param=ConfigParser.ConfigParser()
            self.param.read(paramfname)
	    self.get_workspace(varname=varname)

        # get derived 
	self.get_derived()



    def get_workspace(self, varname=None):

        if varname==None:
            self.parafile=gc.emptydata(self.param, section=self.section)
	elif isinstance(varname, str):
	    exec('self.'+varname+'=gc.emptydata(self.param, section=self.section)')
	else:
	    raise Exception('program initialization error')


    def get_derived(self):
        pass


    def getattribute(self, obj, level=5):
        if((self.debug==True) & (self.debug_level>=level)):
	    print '-'*70
	    print 'All attribute:\n{0}'.format(obj.__dict__.keys() )
	    print 'Attribute values:\n{0}'.format(obj.__dict__)
	    print '-'*70



    def finalize(self):
	mpi.barrier()  # no effect without MPI
        if mpi.rank0:
	    t=time.time()-self.start_time

	    h=int(t/3600)
	    m=int((t-h*3600)/60)
	    s=(t-h*3600-m*60)

            print '\nWe are done, and execution time = {0}h/{1}min/{2}sec'.format(h,m,s)





def prog_init(paramfile=None, section=None, varname=None, init_cosmology=True, **dict):
    ''' ->> Program initialization routines <<- '''
    if len(sys.argv)==3:
        paramfile=sys.argv[1]
        section=sys.argv[2]

    elif len(sys.argv)==2:
        if section==None:
            raise Exception('section name of parameter file `{0}` is needed'.\
	                 format(paramfile))
        paramfile=sys.argv[1]

    elif len(sys.argv)!=1:
        raise Exception('script input error')


    print '>> Importing parameters from section `{0}` of file \n -- `{1}`\n'.format(\
          section, paramfile)

    if mpi.size>1:
        if mpi.rank0:
            print '------ >>> Starting MPI <<< ------'
        print "MPI process %i of %i." % (mpi.rank, mpi.size)

    p=prog(paramfile, section, varname=varname, **dict)

    if (mpi.size>1):
        p.mpi = True
	p.paramdict.update(mpi=True)

    if init_cosmology==True:
        p.cp=cospar.CosmoParams()
        p.dist=frw.Distance(cosmoparam=p.cp)
	#p.power=power.PowerSpectrum(p.cp)

        try:
            pk= np.genfromtxt(p.power_spectrum_fname)
            p.pk=power.PowerSpectrum(p.cp, data=p.power_spectrum_fname)
	except:
            p.pk=power.PowerSpectrum(p.cp)

        #p.pk=power.PowerSpectrum(p.cp)

    if (p.debug==True)&(mpi.rank0):
        print p.paramdict

        if varname==None:
	    varname='parafile'

	try:
            exec('p.getattribute(p.'+varname +')')
	except:
	    pass

    add_root_path()


    return p




def add_root_path():

    folder = os.path.realpath(os.path.abspath( os.getcwd() ) )
    if folder not in sys.path:
        sys.path.insert(0, folder)





def current_path():
    # realpath() with make your script run, even if you symlink it :)
    cmd_folder = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile( inspect.currentframe() ))[0]))
    if cmd_folder not in sys.path:
        sys.path.insert(0, cmd_folder)
    
    # use this if you want to include modules from a subforder
    cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0],"subfolder")))
    if cmd_subfolder not in sys.path:
        sys.path.insert(0, cmd_subfolder)
    
    
    # Info:
    # cmd_folder = os.path.dirname(os.path.abspath(__file__)) # DO NOT USE __file__ !!!
    # __file__ fails if script is called in different ways on Windows
    # __file__ fails if someone does os.chdir() before
    # sys.argv[0] also fails because it doesn't not always contains the path
