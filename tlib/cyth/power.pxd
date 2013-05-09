""" 
"""
cimport genscript.cyth.cubicspline as cs


cdef class PowerSpectrum(object):

    cdef double h, om0, omB, omL, omR, n
    cdef double omC, delH, D0, bondEfsNorm

    cdef public cs.Interpolater power
    cdef public cs.Interpolater power_smooth

