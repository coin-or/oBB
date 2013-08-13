import numpy as np
cimport numpy as np
from cpython cimport PyObject, Py_INCREF

cdef extern from "IQuadProg++.hh":
    cdef cppclass DoubleVector:
        void initialize(double a,unsigned int n)
        void setArray(double* data, unsigned int n)
        PyObject * getArray()
        #DoubleVector *new_DoubleVector "new DoubleVector" ()
    cdef cppclass DoubleMatrix:
        PyObject * getArray()
        void setArray(double* data, unsigned int nrows, unsigned int ncols)
        unsigned int nrows()
        unsigned int ncols()

    double solve_quadprog(DoubleMatrix* G, DoubleVector* g0,
            DoubleMatrix* CE, DoubleVector* ce0,
            DoubleMatrix* CI, DoubleVector* ci0,
            DoubleVector* x)


cdef class PyVector:
    cdef DoubleVector* CppSelf
    #cdef setCppSelf(self, Vector[double]* s)


cdef class PyMatrix:
    cdef DoubleMatrix* CppSelf
