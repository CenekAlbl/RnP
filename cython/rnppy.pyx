cimport numpy as np
import numpy as np


np.import_array()

cdef extern from "rnp.h":
    int r6pSingleLin(double * X, double * u, int direction, double r0, int maxpow, double * Rout, double * wout, double * Cout, double * tout)

cdef extern from "rnp.h":
    int r6pDoubleLin(double * X, double * u, double r0, double * C, double * t, double *v, double * w, int direction)

def r6pSingleLinpy(np.ndarray[double, ndim=2, mode="c"] X not None,
            np.ndarray[double, ndim=2, mode="c"] u not None
            ):
    cdef np.ndarray[double, ndim=2, mode="fortran"] v = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] C = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] w = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] t = np.empty((3,64),order='fortran')

    nres = r6pSingleLin(<double*> np.PyArray_DATA(X), <double*> np.PyArray_DATA(u), 0,0,1, <double*> np.PyArray_DATA(v), <double*> np.PyArray_DATA(w), <double*> np.PyArray_DATA(C), <double*> np.PyArray_DATA(t))

    return v[:,:nres], C[:,:nres], w[:,:nres], t[:,:nres]

def r6pDoubleLinpy(np.ndarray[double, ndim=2, mode="c"] X not None,
            np.ndarray[double, ndim=2, mode="c"] u not None
            ):
    cdef np.ndarray[double, ndim=2, mode="fortran"] v = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] C = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] w = np.empty((3,64),order='fortran')
    cdef np.ndarray[double, ndim=2, mode="fortran"] t = np.empty((3,64),order='fortran')

    nres = r6pDoubleLin(<double*> np.PyArray_DATA(X), <double*> np.PyArray_DATA(u), 0, <double*> np.PyArray_DATA(C), <double*> np.PyArray_DATA(t), <double*> np.PyArray_DATA(v), <double*> np.PyArray_DATA(w), 0)

    return v[:,:nres], C[:,:nres], w[:,:nres], t[:,:nres]