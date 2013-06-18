#include <Python.h>
#include <numpy/arrayobject.h>
#include "untrendy.h"

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

static PyObject
*untrendy_find_discontinuities(PyObject *self, PyObject *args)
{
    double dt, Q, thresh;
    PyObject *t_obj, *chi_obj;

    if (!PyArg_ParseTuple(args, "OOddd", &t_obj, &chi_obj, &dt, &Q, &thresh))
        return NULL;

    PyArrayObject *t_array = (PyArrayObject*)PyArray_FROM_OTF(t_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyArrayObject *chi_array = (PyArrayObject*)PyArray_FROM_OTF(chi_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    if (t_array == NULL || chi_array == NULL) {
        Py_XDECREF(t_array);
        Py_XDECREF(chi_array);
        return NULL;
    }

    int n = (int)PyArray_DIM(t_array, 0);
    if ((int)PyArray_DIM(chi_array, 0) != n) {
        PyErr_SetString(PyExc_ValueError, "Dimension mismatch");
        Py_DECREF(t_array);
        Py_DECREF(chi_array);
        return NULL;
    }

    double *t = (double*)PyArray_DATA(t_array),
           *chi = (double*)PyArray_DATA(chi_array);

    int ind = find_discontinuities(n, t, chi, dt, Q, thresh);

    Py_DECREF(t_array);
    Py_DECREF(chi_array);

    PyObject *ret = Py_BuildValue("i", ind);
    return ret;
}

static PyMethodDef untrendy_methods[] = {
    {"find_discontinuities",
     (PyCFunction)untrendy_find_discontinuities,
     METH_VARARGS,
     "Find discontinuities in a time series."},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int untrendy_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int untrendy_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "_untrendy",
        NULL,
        sizeof(struct module_state),
        untrendy_methods,
        NULL,
        untrendy_traverse,
        untrendy_clear,
        NULL
};

#define INITERROR return NULL

PyObject *PyInit__untrendy(void)

#else
#define INITERROR return

void init_untrendy(void)

#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("_untrendy", untrendy_methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("_untrendy.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    import_array();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
