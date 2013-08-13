#include "Array.hh"
#include "QuadProg++.hh"
#include "Python.h"
#include <numpy/arrayobject.h>

using namespace QuadProgPP;


class DoubleVector : public Vector<double>{
    public:
    DoubleVector(){v = NULL; _import_array();}
    ~DoubleVector(){if (v) delete[] v; v = NULL; n = 0;}
    void initialize(const double a, const unsigned int n){
        this->n = n;
        v = new double[n];
        for (unsigned int i = 0; i < n; i++)
            v[i] = a;
    }

    void setArray(double *data, unsigned int len){
        if (v) delete[] v;
        v = data;
        n = len;
    }

    PyObject* getArray(){
        npy_intp dims = this->n;
        PyObject *Arr = PyArray_SimpleNewFromData(1, &dims,
                PyArray_DOUBLE, this->v);
    return Arr;}
};


class DoubleMatrix : public Matrix<double>{
    public:
    DoubleMatrix(){v = NULL; _import_array();}
    ~DoubleMatrix(){if (v) delete[] v; v = NULL; n = 0; m = 0;}
    PyObject* getArray(){
        npy_intp dims = this->n * this->m;
        PyObject *Arr = PyArray_SimpleNewFromData(1, &dims,
                PyArray_DOUBLE, this->v[0]);
    return Arr;
    }

    void setArray(double *data, unsigned int nrows, unsigned int ncols){
        if (v) delete[] v;
        v = new double*[nrows * ncols];
        v[0] = data;
        for (int i = 1 ; i <= nrows ; i++)
            v[i] = v[i-1] + ncols;
        n = nrows;
        m = ncols;
    }

};


double solve_quadprog(DoubleMatrix* G, DoubleVector* g0,
            const DoubleMatrix* CE, const DoubleVector* ce0,
            const DoubleMatrix* CI, const DoubleVector* ci0,
            DoubleVector* x){
    return solve_quadprog(*G, *g0, *CE, *ce0, *CI, *ci0, *x);
}

