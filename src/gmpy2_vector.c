/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_vector.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
 *                                                                         *
 * This file is part of GMPY2.                                             *
 *                                                                         *
 * GMPY2 is free software: you can redistribute it and/or modify it under  *
 * the terms of the GNU Lesser General Public License as published by the  *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2 is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    *
 * License for more details.                                               *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with GMPY2; if not, see <http://www.gnu.org/licenses/>    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* The following code was a test case for creating a gmpy2 function that would
 * apply an MPFR function to all the elements of a list. It was slightly faster
 * than map(gmpy2.list, <<list>>) but not enough to justify the effort. The
 * code is left for possible future use.
 */

PyDoc_STRVAR(GMPy_doc_function_vector,
"vector(iterable) -> list\n\n"
"Template for applying a function to an iterable.");

PyDoc_STRVAR(GMPy_doc_context_vector,
"vector(iterable) -> list\n\n"
"Template for applying a function to an iterable.");

static PyObject *
GMPy_Context_Vector(PyObject *self, PyObject *other)
{
    PyObject *result, *tempres;
    Py_ssize_t i, seq_length;

    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(other = PySequence_List(other))) {
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }

    /* other contains a new list containing all the values from the
     * iterable. Create a list to store the results.
     */

    seq_length = PyList_GET_SIZE(other);
    if (!(result = PyList_New(seq_length))) {
        Py_DECREF(other);
        return NULL;
    }

    /* Iterate through the list. */

    for (i=0; i < seq_length; i++) {
        if (!(tempres = GMPy_Number_Sin(PyList_GET_ITEM(other, i), context))) {
            Py_DECREF(other);
            Py_DECREF(result);
            TYPE_ERROR("all items in iterable must be numbers");
            return NULL;
        }

        if (PyList_SetItem(result, i, tempres) < 0) {
            Py_DECREF(other);
            Py_DECREF(result);
            return NULL;
        }
    }

    Py_DECREF(other);

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_vector2,
"vector2(iterable, iterable) -> list\n\n"
"Template for applying a function to a pair of iterables.");

PyDoc_STRVAR(GMPy_doc_context_vector2,
"vector2(iterable) -> list\n\n"
"Template for applying a function to a pair of iterables.");

static PyObject *
GMPy_Context_Vector2(PyObject *self, PyObject *args)
{
    PyObject *arg1, *arg2, *result;
    MPFR_Object *tempres;

    Py_ssize_t i, seq_length;

    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("vector2() requires 2 arguments");
        return NULL;
    }

    if (!(arg1 = PySequence_List(PyTuple_GET_ITEM(args, 0)))) {
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }

    if (!(arg2 = PySequence_List(PyTuple_GET_ITEM(args, 1)))) {
        Py_DECREF(arg1);
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }


    /* other contains a new list containing all the values from the
     * iterable. Create a list to store the results.
     */

    if (PyList_GET_SIZE(arg1) != PyList_GET_SIZE(arg2)) {
        Py_DECREF(arg1);
        Py_DECREF(arg2);
        TYPE_ERROR("arguments must be the same length");
        return NULL;
    }

    seq_length = PyList_GET_SIZE(arg1);
    if (!(result = PyList_New(seq_length))) {
        Py_DECREF(arg1);
        Py_DECREF(arg2);
        return NULL;
    }

    /* Iterate through the list. */

    for (i=0; i < seq_length; i++) {
        //~ if (!(tempres = GMPy_Number_Mul(PyList_GET_ITEM(arg1, i),
                                        //~ PyList_GET_ITEM(arg2, i),
                                        //~ context))) {
            //~ Py_DECREF(arg1);
            //~ Py_DECREF(arg2);
            //~ Py_DECREF(result);
            //~ TYPE_ERROR("all items in iterable must be numbers");
            //~ return NULL;
        //~ }

        tempres = GMPy_MPFR_New(0, context);
        tempres->rc = mpfr_mul(tempres->f,
                               MPFR(PyList_GET_ITEM(arg1, i)),
                               MPFR(PyList_GET_ITEM(arg2, i)),
                               GET_MPFR_ROUND(context));

        if (PyList_SetItem(result, i, (PyObject*)tempres) < 0) {
            Py_DECREF(arg1);
            Py_DECREF(arg2);
            Py_DECREF(result);
            return NULL;
        }
    }

    Py_DECREF(arg1);
    Py_DECREF(arg2);

    return (PyObject*)result;
}
