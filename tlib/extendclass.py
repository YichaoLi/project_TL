""" >>> Extended class
"""
import sys
import numpy as np



class myDict(dict):
    """ my extended class of dictionary
    """
    def __add__(self, other):
        copy = self.copy()
        copy.update(other)
        return copy

    def __radd__(self, other):
        copy = other.copy()
        copy.update(self)
        return copy






