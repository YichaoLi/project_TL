''' >> Partially borrowed from Richard's code <<- '''
import warnings
import numpy as np

_rank = 0
_size = 1
_comm = None
world = None

rank0 = True

## Try to setup MPI and get the comm, rank and size.
## If not they should end up as rank=0, size=1.
try:
    from mpi4py import MPI

    _comm = MPI.COMM_WORLD
    world = _comm
    
    rank = _comm.Get_rank()
    size = _comm.Get_size()

    _rank = rank if rank else 0
    _size = size if size else 1

    #if rank:
    #    print "MPI process %i of %i." % (_rank, _size)

    rank0 = True if _rank == 0 else False
    
except ImportError:
    warnings.warn("Warning: mpi4py not installed.")


def partition_list_alternate(full_list, i, n, ndarray=True):
    """Partition a list into `n` pieces. Return the `i`th partition."""
    if ndarray:
        return np.array(full_list[i::n])
    else:
        return full_list[i::n]


def partition_list_mpi(full_list, ndarray=True):
    """Return the partition of a list specific to the current MPI process."""
    return partition_list_alternate(full_list, _rank, _size, ndarray=ndarray)


def mpirange(*args):
    """An MPI aware version of `range`, each process gets its own sub section.
    """
    full_list = range(*args)
    
    #if alternate:
    return partition_list_alternate(full_list, _rank, _size, ndarray=False)
    #else:
    #    return np.array_split(full_list, _size)[rank_]




def barrier():
    if _size > 1:
        _comm.Barrier()
    



def gather(d, n, root=0):
    # gather distributed data

    if _size>1:
        dat_=_comm.gather(d, root=root) 

        if rank0:
            dat=np.zeros(n)
            drange=range(n)

	    for i in range(_size):
	        dat[drange[i::_size]] = dat_[i]
	else:
	    dat=None

	dat=_comm.bcast(dat, root=0)

        return dat

    else:
        return d



def gather_unify(d, root=0):

    if _size>1:

        _d = _comm.gather(d, root=0)

        if rank0:
            data = _d[0]
            for i in range(1, _size):
                data += _d[i]

        else:
            data=None
    
        data=_comm.bcast(data, root=0)

	return data
    
    else:
        return _d

