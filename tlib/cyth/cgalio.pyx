import numpy as np
cimport numpy as np
cimport cgal as cg




def read_position(fname, npart):
    cdef np.ndarray[np.float32_t, ndim=1, mode = 'c'] x, y, z

    x=np.zeros(npart, dtype=np.float32) 
    y=np.zeros(npart, dtype=np.float32)
    z=np.zeros(npart, dtype=np.float32)

    cg.read_particle_positions(<char *>fname, <int>npart, <float *>x.data, \
                               <float *>y.data, <float *>z.data)


    return x, y, z



def read_velocity(fname, npart):
    cdef np.ndarray[np.float_t, ndim=1, mode = 'c'] d
    d=np.zeros(npart, dtype=np.float)

    cg.read_particle_velocity(<char *>fname, <int>npart, <double *>d.data)

    return d



def write_position(fname, xx, yy, zz):
    cdef:
        np.ndarray[np.float32_t, ndim=1, mode = 'c'] x, y, z
        int npart

    npart = <int>xx.shape[0]
    x=np.ascontiguousarray(xx, dtype=np.float32)
    y=np.ascontiguousarray(yy, dtype=np.float32)
    z=np.ascontiguousarray(zz, dtype=np.float32)

    cg.write_particle_positions(<char *>fname, <int>npart, <float *>x.data, \
                                <float *>y.data, <float *>z.data)
    return




def write_velocity(fname, dd):
    cdef:
        np.ndarray[np.float64_t, ndim=1, mode = 'c'] d
        int npart

    npart = <int>dd.shape[0]
    d=np.ascontiguousarray(dd, dtype=np.float64)

    cg.write_particle_velocity(<char *>fname, <int>npart, <double *>d.data)
    return


def read_header(fname):
    cdef:
        int npart 
        float boxsize

    cg.read_particle_header(<char *>fname, &npart, &boxsize)

    return npart, boxsize
