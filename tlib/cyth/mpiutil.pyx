import numpy as np
cimport numpy as np


from mpi4py cimport MPI as mpi
from mpi4py cimport mpi_c 

_comm = mpi.COMM_WORLD
world = _comm

rank = _comm.Get_rank()
size = _comm.Get_size()

_rank = rank if rank else 0
_size = size if size else 1

rank0 = True if _rank == 0 else False


cpdef np.ndarray[np.int_t, ndim=1] partition_list_alternate(\
              np.ndarray[np.int_t, ndim=1]full_list, i, n):
    '''Partition a list into `n` pieces. Return the `i`th partition.'''
    return full_list[i::n]



"""
#cdef np.ndarray[np.int_t, ndim=1] partition_list_alternate(full_list, int i, \
#                   n, ndarray=True):
def partition_list_alternate(full_list, i, n, ndarray=True):
    '''Partition a list into `n` pieces. Return the `i`th partition.'''
    if ndarray:
        return np.array(full_list[i::n])
    else:
        return full_list[i::n]


def partition_list_mpi(full_list, ndarray=True):
    '''Return the partition of a list specific to the current MPI process.'''
    return partition_list_alternate(full_list, _rank, _size, ndarray=ndarray)


def mpirange(*args):
    '''An MPI aware version of `range`, each process gets its own sub section.
    '''
    full_list = range(*args)
    
    #if alternate:
    return partition_list_alternate(full_list, _rank, _size, ndarray=False)
    #else:
    #    return np.array_split(full_list, _size)[rank_]
"""

cpdef barrier():
    if _size > 1:
        _comm.Barrier()
    
