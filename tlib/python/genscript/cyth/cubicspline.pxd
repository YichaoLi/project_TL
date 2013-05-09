"""
"""
import numpy as np
cimport numpy as np


dbltype = np.int
ctypedef np.float64_t dbltype_t

cdef extern from "math.h":
    double exp(double x)
    double log(double x)


cdef class Interpolater(object):

    cdef np.ndarray __data_
    cdef np.ndarray __y2_


    cdef dbltype_t value_cdef(self, dbltype_t x)


