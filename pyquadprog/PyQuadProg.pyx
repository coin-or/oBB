import numpy as np
cimport PyQuadProg


class PyQuadProg:
    '''
    Please see README.rst to find out what these arguments are.
    '''
    def __init__(self, G, g0, CE, ce0, CI, ci0):
        self.G = PyMatrix(G)
        self.CE = PyMatrix(CE)
        self.CI = PyMatrix(CI)
        self.g0 = PyVector(g0)
        self.ce0 = PyVector(ce0)
        self.ci0 = PyVector(ci0.astype(np.double))
        self.x = PyVector(np.empty(1, dtype=np.double))
        cy_solve_quadprog(self.G, self.g0, self.CE, self.ce0,
                                 self.CI, self.ci0, self.x)


cdef class PyVector:
    def __cinit__(self, np.ndarray[np.double_t, ndim=1]data):
        self.CppSelf = new DoubleVector()
        self.setArray(data)

    def __dealloc__(self):
        del self.CppSelf

    def initialize(self, a, n):
        self.CppSelf.initialize(a, n)

    def getArray(self):
        return <object>self.CppSelf.getArray()

    def setArray(self, np.ndarray[np.double_t, ndim=1]data):
        Py_INCREF(data)
        self.CppSelf.setArray(<double*>data.data, data.shape[0])

    def __repr__(self):
        return repr(self.getArray())

    def __str__(self):
        return str(self.getArray())


cdef class PyMatrix:
    def __cinit__(self, np.ndarray[np.double_t, ndim=2]data):
        self.CppSelf = new DoubleMatrix()
        self.setArray(data)

    def __dealloc__(self):
        del self.CppSelf

    def getArray(self):
        mat = <object>self.CppSelf.getArray()
        return mat.reshape((self.nRows, self.nCols))

    def setArray(self, np.ndarray[np.double_t, ndim=2]data):
        Py_INCREF(data)
        self.CppSelf.setArray(<double*>data.data, data.shape[0], data.shape[1])

    property nRows:
        def __get__(self):
            return self.CppSelf.nrows()

    property nCols:
        def __get__(self):
            return self.CppSelf.ncols()

    def __str__(self):
        return str(self.getArray())


def cy_solve_quadprog(PyMatrix G, PyVector g0,
            PyMatrix CE, PyVector ce0,
            PyMatrix CI, PyVector ci0,
            PyVector x):
    return solve_quadprog(G.CppSelf, g0.CppSelf,
            CE.CppSelf, ce0.CppSelf,
            CI.CppSelf, ci0.CppSelf,
            x.CppSelf)



