""" >>>  General Classes  <<< """
import sys
import time
import ConfigParser
import numpy as np
import pylab as pl


"""  --- >>> mydata class <<< ---  """
# emptydata only get the necessary information of 'folder', 'file_number', 'fname' etc.
class emptydata(object):
    def __init__(self, config, section='DEFAULT'):
	self.section=section
	self.getfolder(config)
	self.getfn(config)
	self.getfname_set(config)
        try: 
	    self.get_dataformat(config)
	except: 
	    print 'no format information of the data is imported.'
	    pass

        self.get_fn_output(config)
        self.get_output_set(config)

	print 'reading from: section={0}; folder={1}'.format(self.section, self.folder)

    def getfolder(self, config, field='folder'):
	self.folder = config.get(self.section, field)

    def getfn(self, config, field = 'file_number'):
	self.fn = config.getint(self.section, field)

    def getfname_set(self, config, prefix='fname_'):
	self.fname = [ config.get(self.section, prefix+str(i+1)) for i in \
		       range(self.fn)]
	self.fnamefull = [self.folder+self.fname[i] for i in range(len(self.fname))]


    def get_dataformat(self, config, prefix='data_format'):
	self.dform_support=set(['gadget', 'npy', 'npz', 'binary'])
        datf=config.get(self.section, prefix)
	self.dataformat=[x.strip() for x in datf.split(',')]
	if(len(self.dataformat)!=self.fn):
	    raise Exception('emptydata: # of dataformat inconsistent.')

        if( (set(self.dataformat)<=self.dform_support) != True):
	    raise Exception('emptydata: unsupported dataformat.')


    def get_fn_output(self, config, field = 'output_file_number'):
	self.fn_output = config.getint(self.section, field)

    def get_output_set(self, config, prefix= 'output_'):
	self.output = [ config.get(self.section, prefix+str(i+1)) for i in \
		       range(self.fn_output)]

	self.outputfull = [self.folder+self.output[i] for i in range(self.fn_output)]
	self.output_split=[self.outputfull[i].split('/') for i in range(self.fn_output)]



class data(emptydata):

    def get_data(self, fname, dtype='col'):
	if(dtype=='col'):
	    data = np.array(list(zip(* np.genfromtxt( fname ) )) )  
	elif(dtype=='row'):
	    data = np.genfromtxt(fname) 
	else:
	    raise Exception
	return data


    def get_dataset(self, fnameset, dtype='col', arrtype='numpy'):
	dat = [self.get_data(fnameset[i], dtype=dtype) for i in range(self.fn)]

	if(arrtype=='numpy'):
	    data = np.array(dat)
	    self.arrtype = arrtype

	return data


class mydata():
    """ ??? maybe  absorbed into 'data' and 'emptydata' ??? """
    def __init__(self, fname):
	self.fname=fname

    def getdata(self, dtype='col', arrtype='list'):
        nfile = len(self.fname)
	if(dtype=='col'):
	    # colum data 
	    self.data = [ np.array(list(zip(* np.genfromtxt(self.fname[i]) )) ) for 
                          i in range(nfile) ] 
            self.ndat = [ len(self.data[i]) for i in range(nfile) ]
	elif(dtype=='row'):
	    self.data = [ np.genfromtxt(self.fname[i]) for i in range(nfile) ] 
            self.ndat = [ len(self.data[i]) for i in range(nfile) ]
	else:
	    raise Exception

	if(arrtype=='numpy'):
	    self.data = np.array(self.data)

	self.arrtype = arrtype

    # get data from only one file
    #def getdata_one(self, dtype='col', arrtype='list'):



#--- >>> this class stores information of the workspace where the script run <<< ---
class workspace():
    """ --- >>> ??? still useful ??? <<< --- """
    def __init__(self, config, section='DEFAULT'):

        self.folder_in  = config.get(section, 'input_folder')
        self.folder_out = config.get(section, 'output_folder')
        print 'input_folder:', self.folder_in, 'output_folder:', self.folder_out

        self.infbundle = config.getboolean(section, 'input_bundle_file')
        self.infn      = config.getint(section, 'input_file_number')
        self.infn_bund = config.getint(section, 'bundle_file_number')
        self.outfn      = config.getint(section, 'output_file_number')
	
	if(self.infbundle == True):
	    self.prefix_in = config.get(section, 'input_file_prefix')
            self.suffix_in = config.get(section, 'input_file_suffix')

	    try:
		self.bindex_min = config.getint(section, 'bundle_index_min')
	    except ValueError:
		print 'Can NOT import min bundle index from parameter file.'
	    """else: 
		self.bindex_min = 0 """
	    self.bindex_max = self.bindex_min + self.infn_bund

	    print self.bindex_min, self.bindex_max

	    self.bfname_in = [ self.folder_in +self.prefix_in + str(i) + \
	                       self.suffix_in for i in range(self.bindex_min, self.bindex_max)]

	self.fname_in = [self.folder_in + config.get(section, 'input_file_name_'\
	                  +str(i+1) ) for i in range(self.infn)]
		

	self.fname_out = [self.folder_out+config.get(section, \
	                 'output_file_name_'+str(i+1) ) for i in range(self.outfn)]

	self.section = config.get(section, 'main_section')

	"""
	print self.fname_in
	print self.bfname_in
	print self.fname_out
	print self.section 
	"""
        print 'initialization is done.'




"""  --- >>>  Data cube class, e.g. 3D density field <<< ---  """
class datcube(object):
    def __init__(self, data, boxsize):
        self.data = data
	self.boxsize = boxsize
        self.ngrid = np.array(self.data.shape)
	self.getgrid2d()

    def getgrid2d(self):
        if(self.ngrid[0]==self.ngrid[1] & self.ngrid[0]==self.ngrid[2]):
	    b1, b2, grd = -self.boxsize/2., self.boxsize/2., self.ngrid[0]*1j   
	self.grid2d = np.mgrid[b1:b2:grd, b1:b2:grd]




class parameterfile(object):
    """ --- >>> general class of parameters files, ??? Possible Gernalization:  
                list of section?? then 'allparams' may have inconsistent items....??  <<< --- """

    def __init__(self, fname, section='DEFAULT'):

        if(isinstance(fname, 'str')==False):
	    raise Exception('parameter file class error, fname is not a string')

	if(isinstance(section, 'list')):
            raise Exception('list of section are NOT supported yet.')
	elif(isinstance(section, 'str')):
            pass
	else:
	    raise Exception('parameter file class error, section is neither string nor list.')

        self.fname=fname
	self.section = section

        #--- >>> initialize the parser <<< ---
        self.parser=ConfigParser.ConfigParser()
        self.parser.read(fname)

	# >>> get all the items <<< #
	self.allparams = dict(self.parser.items(section))
	self.allparams_keys = self.allparams.keys()
	self.allparams_values = self.allparams.values()



class myprog(object):
    """ --- >>> Program Control <<< --- """
    def __init__(self, argv, section='DEFAULT'):
        if len(sys.argv)!=2:
            print 'input parameter file needed.'
            sys.exit(0)
        
        self.param=ConfigParser.ConfigParser()
        self.param.read(sys.argv[1])
	self.get_progcontrol()
	self.section = section

	self.start_time=time.time()



    def get_progcontrol(self, section='PROGRAM_CONTROL'):
        self.debug=self.param.getboolean(section, 'debug')
	self.debug_level=self.param.getint(section, 'debug_level')
	if(self.debug==True):
	    print 'myprog: debug={0}, debug level={1}'.format(self.debug, 
	           self.debug_level)

    def getattribute(self, obj, level=5):
        if((self.debug==True) & (self.debug_level>=level)):
	    print '-'*70
	    print 'All attribute:\n{0}'.format(obj.__dict__.keys() )
	    print 'Attribute values:\n{0}'.format(obj.__dict__)
	    print '-'*70

    def finalize(self):
        print '\nWe are done, and execution time = {0} sec.'.format(time.time()-self.start_time)
