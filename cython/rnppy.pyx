cimport numpy as np
import numpy as np


np.import_array()

cdef extern from "rnp.h":
    int r6pDoubleLin(double * X, double * u, int direction, double r0, int maxpow, double * Rout, double * wout, double * Cout, double * tout)

def r6ppy(np.ndarray[double, ndim=2, mode="c"] X not None,
            np.ndarray[double, ndim=2, mode="c"] u not None
            ):
    cdef np.ndarray[double, ndim=2, mode="fortran"] R = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] C = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] w = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] t = np.empty((3,64),order='fortran')


    nres = r6pDoubleLin(<double*> np.PyArray_DATA(X), <double*> np.PyArray_DATA(u), 0,0,1, <double*> np.PyArray_DATA(R), <double*> np.PyArray_DATA(w), <double*> np.PyArray_DATA(C), <double*> np.PyArray_DATA(t))

    return R[:,:nres], C[:,:nres], w[:,:nres], t[:,:nres]