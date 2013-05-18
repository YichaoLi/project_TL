#! /usr/bin/env python 
'''
    This is an example of the module. We have a framework for calling and 
    executing the modules. That is using the manager.py or mpimanager.py 
    with a XXX.pipe. If you want to work within this framework, your code
    need follow the structure defined here. 
'''

import scipy as np
import numpy as np

from kiyopy import parse_ini


'''
    We use a dict manage the parameters. The dict should be definded here
    and the value can be give in your XXX.pipe. To distinguish parameter
    dicts of different module, we use a prefix. 
'''
prefix = 'test_'
params_init = {
    'param_1' : 1.,
    'param_2' : 2.,
}

class TestExample(object):
    def __init__(self, parameter_file_or_dict=None, feedback=1):
        '''
            The initialize function here. 
            parameter_file_or_dict : This can be your XXX.pipe 
            feedback : how much information givin by the code.
        '''
        # call parse_ini.parse() to initialize the parameters. 
        # all parameters with prefix 'test_' will be read out,
        # and save into a dict self.params.
        self.params = parse_ini.parse(parameter_file_or_dict,
                                      params_init,
                                      prefix=prefix,
                                      feedback=feedback)
        self.feedback = feedback

    def execute(self, nprocesses=1):
        '''
            Do what you want to do. The manager will automatically call
            this method to start the program. 

        '''

        print self.params['param_1']
        print self.params['param_2']

        self.test()

    def test(self):
        print self.params['param_1'] + self.params['param_2']
