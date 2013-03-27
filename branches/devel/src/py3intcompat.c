/* This file contains source code copied from Python's source code.

   "Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009,
           2010 Python Software Foundation; All Rights Reserved"

   For the full text of the license, see http://www.python.org/psf/license/.

   =======================================

   This file provides support to older versions of Python for several new C
   API functions introduced in Python 2.7, Python 3.0, or Python 3.2.

   PyLong_AsLongAndOverflow accepts a PyLong object and either returns a C
   long or sets an overflow flag. A PyInt is also accepted on Python 2.x so
   the macro PyIntOrLong_Check() is defined in gmpy.c to accept either type
   of integer. Actually, PyLong_AsLongAndOverflow will try to convert almost
   any object to an integer if it is not already a PyInt or PyLong. If you
   don't verify the object is either a PyInt or PyLong first, you MUST check
   for an error in the return value.

   PyLong_AsLongLongAndOverflow behaves similarly except it returns a C
   long long.

   Example
   -------
   #if PY_MAJOR_VERSION == 3
   #define PyIntOrLong_Check(op) PyLong_Check(op)
   #else
   #define PyIntOrLong_Check(op) (PyInt_Check(op) || PyLong_Check(op))
   #endif

   #include "py3intcompat.c"

   ......

   long temp;
   int overflow;

   ......

   if(PyIntOrLong_Check(obj)) {
       temp = PyLong_AsLongAndOverflow(obj, &overflow);
       if(overflow) {
           obj is a PyLong that won't fit into a C long.
       } else {
           Process a PyInt or a PyLong that will fit into a C long.
       }
   }

   ......

   # Don't check the object type first.
   temp = PyLong_AsLongAndOverflow(obj, &overflow);
   if (temp == -1 && PyErr_Occurred()) {
       raise an error message
   } else if (overflow) {
       obj is a PyLong that won't fit into a C long.
   } else {
       Process a PyInt or a PyLong that will fit into a C long.
   }


*/

#include "longintrepr.h"

/* Get a C long int from a long int object.
   Returns -1 and sets an error condition if overflow occurs. */

#if (PY_VERSION_HEX < 0x02070000)

#ifndef PyLong_SHIFT
#define PyLong_SHIFT SHIFT
#endif

#ifndef PyLong_MASK
#define PyLong_MASK MASK
#endif

#define PY_ABS_LONG_MIN     (0-(unsigned long)LONG_MIN)

typedef short sdigit;

static long
PyLong_AsLongAndOverflow(PyObject *vv, int *overflow)
{
    /* This version by Tim Peters */
    register PyLongObject *v;
    unsigned long x, prev;
    long res;
    Py_ssize_t i;
    int sign;
    int do_decref = 0; /* if nb_int was called */

    *overflow = 0;
    if (vv == NULL) {
        PyErr_BadInternalCall();
        return -1;
    }

    if(PyInt_Check(vv))
        return PyInt_AsLong(vv);

    if (!PyLong_Check(vv)) {
        PyNumberMethods *nb;
        nb = vv->ob_type->tp_as_number;
        if (nb == NULL || nb->nb_int == NULL) {
            PyErr_SetString(PyExc_TypeError,
                    "an integer is required");
            return -1;
        }
        vv = (*nb->nb_int) (vv);
        if (vv == NULL)
            return -1;
        do_decref = 1;
        if(PyInt_Check(vv)) {
            res = PyInt_AsLong(vv);
            goto exit;
        }
        if (!PyLong_Check(vv)) {
            Py_DECREF(vv);
            PyErr_SetString(PyExc_TypeError,
                    "nb_int should return int object");
            return -1;
        }
    }

    res = -1;
    v = (PyLongObject *)vv;
    i = Py_SIZE(v);

    switch (i) {
    case -1:
        res = -(sdigit)v->ob_digit[0];
        break;
    case 0:
        res = 0;
        break;
    case 1:
        res = v->ob_digit[0];
        break;
    default:
        sign = 1;
        x = 0;
        if (i < 0) {
            sign = -1;
            i = -(i);
        }
        while (--i >= 0) {
            prev = x;
            x = (x << PyLong_SHIFT) + v->ob_digit[i];
            if ((x >> PyLong_SHIFT) != prev) {
                *overflow = Py_SIZE(v) > 0 ? 1 : -1;
                goto exit;
            }
        }
        /* Haven't lost any bits, but casting to long requires extra care
         * (see comment above).
             */
        if (x <= (unsigned long)LONG_MAX) {
            res = (long)x * sign;
        }
        else if (sign < 0 && x == PY_ABS_LONG_MIN) {
            res = LONG_MIN;
        }
        else {
            *overflow = sign;
            /* res is already set to -1 */
        }
    }
 exit:
    if (do_decref) {
        Py_DECREF(vv);
    }
    return res;
}
#endif   /* if (PY_VERSION_HEX < 0x02070000) */

#ifdef _WIN64
#if HAVE_LONG_LONG
#if (PY_VERSION_HEX < 0x02070000)
#define PY_ABS_LLONG_MIN    (0-(unsigned PY_LONG_LONG)PY_LLONG_MIN)
static PY_LONG_LONG
PyLong_AsLongLongAndOverflow(PyObject *vv, int *overflow)
{
    /* This version by Tim Peters */
    register PyLongObject *v;
    unsigned PY_LONG_LONG x, prev;
    PY_LONG_LONG res;
    Py_ssize_t i;
    int sign;
    int do_decref = 0; /* if nb_int was called */

    *overflow = 0;
    if (vv == NULL) {
        PyErr_BadInternalCall();
        return -1;
    }

    if(PyInt_Check(vv))
        return PyInt_AsLong(vv);

    if (!PyLong_Check(vv)) {
        PyNumberMethods *nb;
        nb = vv->ob_type->tp_as_number;
        if (nb == NULL || nb->nb_int == NULL) {
            PyErr_SetString(PyExc_TypeError,
                    "an integer is required");
            return -1;
        }
        vv = (*nb->nb_int) (vv);
        if (vv == NULL)
            return -1;
        do_decref = 1;
        if(PyInt_Check(vv)) {
            res = PyInt_AsLong(vv);
            goto exit;
        }
        if (!PyLong_Check(vv)) {
            Py_DECREF(vv);
            PyErr_SetString(PyExc_TypeError,
                    "nb_int should return int object");
            return -1;
        }
    }

    res = -1;
    v = (PyLongObject *)vv;
    i = Py_SIZE(v);

    switch (i) {
    case -1:
        res = -(sdigit)v->ob_digit[0];
        break;
    case 0:
        res = 0;
        break;
    case 1:
        res = v->ob_digit[0];
        break;
    default:
        sign = 1;
        x = 0;
        if (i < 0) {
            sign = -1;
            i = -(i);
        }
        while (--i >= 0) {
            prev = x;
            x = (x << PyLong_SHIFT) + v->ob_digit[i];
            if ((x >> PyLong_SHIFT) != prev) {
                *overflow = sign;
                goto exit;
            }
        }
        /* Haven't lost any bits, but casting to long requires extra
         * care (see comment above).
         */
        if (x <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
            res = (PY_LONG_LONG)x * sign;
        }
        else if (sign < 0 && x == PY_ABS_LLONG_MIN) {
            res = PY_LLONG_MIN;
        }
        else {
            *overflow = sign;
            /* res is already set to -1 */
        }
    }
 exit:
    if (do_decref) {
        Py_DECREF(vv);
    }
    return res;
}
#endif   /* if (PY_VERSION_HEX < 0x02070000) */

#if ((PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION < 2))
#define PY_ABS_LLONG_MIN    (0-(unsigned PY_LONG_LONG)PY_LLONG_MIN)
static PY_LONG_LONG
PyLong_AsLongLongAndOverflow(PyObject *vv, int *overflow)
{
    /* This version by Tim Peters */
    register PyLongObject *v;
    unsigned PY_LONG_LONG x, prev;
    PY_LONG_LONG res;
    Py_ssize_t i;
    int sign;
    int do_decref = 0; /* if nb_int was called */

    *overflow = 0;
    if (vv == NULL) {
        PyErr_BadInternalCall();
        return -1;
    }

    if (!PyLong_Check(vv)) {
        PyNumberMethods *nb;
        nb = vv->ob_type->tp_as_number;
        if (nb == NULL || nb->nb_int == NULL) {
            PyErr_SetString(PyExc_TypeError,
                    "an integer is required");
            return -1;
        }
        vv = (*nb->nb_int) (vv);
        if (vv == NULL)
            return -1;
        do_decref = 1;
        if (!PyLong_Check(vv)) {
            Py_DECREF(vv);
            PyErr_SetString(PyExc_TypeError,
                    "nb_int should return int object");
            return -1;
        }
    }

    res = -1;
    v = (PyLongObject *)vv;
    i = Py_SIZE(v);

    switch (i) {
    case -1:
        res = -(sdigit)v->ob_digit[0];
        break;
    case 0:
        res = 0;
        break;
    case 1:
        res = v->ob_digit[0];
        break;
    default:
        sign = 1;
        x = 0;
        if (i < 0) {
            sign = -1;
            i = -(i);
        }
        while (--i >= 0) {
            prev = x;
            x = (x << PyLong_SHIFT) + v->ob_digit[i];
            if ((x >> PyLong_SHIFT) != prev) {
                *overflow = sign;
                goto exit;
            }
        }
        /* Haven't lost any bits, but casting to long requires extra
         * care (see comment above).
         */
        if (x <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
            res = (PY_LONG_LONG)x * sign;
        }
        else if (sign < 0 && x == PY_ABS_LLONG_MIN) {
            res = PY_LLONG_MIN;
        }
        else {
            *overflow = sign;
            /* res is already set to -1 */
        }
    }
 exit:
    if (do_decref) {
        Py_DECREF(vv);
    }
    return res;
}
#endif   /* if (PY_VERSION_HEX < 0x03020000) */

#endif   /* if HAVE_LONG_LONG */
#endif   /* ifdef _WIN64 */
