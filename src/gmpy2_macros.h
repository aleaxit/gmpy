/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_macros.h                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018 Case Van Horsen                        *
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

/* This file contains a collection of macros that can be used to reduce
 * repetitive code. As new macros are written to support of refactoring,
 * they should be placed here.
 */

/* NAME is used as part of the GMPy function name. It usually uses an upper-
 *      case first character.
 * FUNC is the component of the actual name used by MPFR and MPC.
 *
 * GMPY_MPFR_MPC_UNIOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, context)
 *     GMPy_Complex_NAME(x, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_MPC_UNIOP_EX(NAME, FUNC) creates the following functions:
 *     GMPy_MPFR_NAME(x, context)
 *     GMPy_Real_NAME(x, context)
 *     GMPy_MPC_NAME(x, context)
 *     GMPy_Complex_NAME(x, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_MPC_UNIOP_TEMPLATE(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_Real_NAME & GMPy_Complex_NAME exist
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_MPC_TRIOP_EX(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, y, Z, context)
 *     GMPy_Complex_NAME(x, y, Z, context)
 *     GMPy_Number_NAME(x, y, Z, context)
 *     - assumes GMPy_Integer_NAME & GMPy_Rational_NAME also exist
 *     GMPy_Context_NAME(self, args)
 *     - called with METH_VARARGS
 *
 * GMPY_MPFR_MPC_UNIOP_TEMPLATE(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_Real_NAME & GMPy_Complex_NAME exist
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_UNIOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_UNIOP_TEMPLATE(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_Real_NAME exists
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 *
 * GMPY_MPFR_BINOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, y, context)
 *     GMPy_Number_NAME(x, y, context)
 *     GMPy_Context_NAME(self, args)
 *     - called with METH_VARARGS
 *
 */

#define GMPY_MPFR_MPC_UNIOP(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \

#define GMPY_MPFR_MPC_UNIOP_EX(NAME, FUNC) \
static PyObject * \
_GMPy_MPFR_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result; \
    CHECK_CONTEXT(context); \
    if (!(result = GMPy_MPFR_New(0, context))) { \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, MPFR(x), GET_MPFR_ROUND(context)); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *tempx; \
    PyObject *result; \
    CHECK_CONTEXT(context); \
    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context))) { \
        return NULL; \
    } \
    result = _GMPy_MPFR_##NAME((PyObject*)tempx, context); \
    Py_DECREF((PyObject*)tempx); \
    return (PyObject*)result; \
} \
static PyObject * \
_GMPy_MPC_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPC_Object *result; \
    CHECK_CONTEXT(context); \
    if (!(result = GMPy_MPC_New(0, 0, context))) { \
        return NULL; \
    } \
    result->rc = mpc_##FUNC(result->c, MPC(x), GET_MPC_ROUND(context)); \
    _GMPy_MPC_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Complex_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPC_Object *tempx; \
    PyObject *result; \
    CHECK_CONTEXT(context); \
    if (!(tempx = GMPy_MPC_From_Complex(x, 1, 1, context))) { \
        return NULL; \
    } \
    result = _GMPy_MPC_##NAME((PyObject*)tempx, context); \
    Py_DECREF((PyObject*)tempx); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (MPFR_Check(x)) \
        return _GMPy_MPFR_##NAME(x, context); \
    if (MPC_Check(x)) \
        return _GMPy_MPC_##NAME(x, context); \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    if (IS_COMPLEX(x)) \
        return GMPy_Complex_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_MPC_UNIOP_TEMPLATE(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    if (IS_COMPLEX(x)) \
        return GMPy_Complex_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_MPC_UNIOP_TEMPLATE_EX(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (MPFR_Check(x)) \
        return _GMPy_MPFR_##NAME(x, context); \
    if (MPC_Check(x)) \
        return _GMPy_MPC_##NAME(x, context); \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    if (IS_COMPLEX(x)) \
        return GMPy_Complex_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_MPC_TRIOP_TEMPLATE(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context) \
{ \
    if (MPZ_Check(x) && MPZ_Check(y) && MPZ_Check(z)) \
        return _GMPy_MPZ_##NAME(x, y, z, context); \
    if (MPQ_Check(x) && MPQ_Check(y) && MPQ_Check(z)) \
        return _GMPy_MPQ_##NAME(x, y, z, context); \
    if (MPFR_Check(x) && MPFR_Check(y) && MPFR_Check(z)) \
        return _GMPy_MPFR_##NAME(x, y, z, context); \
    if (MPC_Check(x) && MPC_Check(y) && MPC_Check(z)) \
        return _GMPy_MPC_##NAME(x, y, z, context); \
    if (IS_INTEGER(x) && IS_INTEGER(y) && IS_INTEGER(z)) \
        return GMPy_Integer_##NAME(x, y, z, context); \
    if (IS_RATIONAL(x) && IS_RATIONAL(y) && IS_RATIONAL(z)) \
        return GMPy_Rational_##NAME(x, y, z, context); \
    if (IS_REAL(x) && IS_REAL(y) && IS_REAL(z)) \
        return GMPy_Real_##NAME(x, y, z, context); \
    if (IS_COMPLEX(x) && IS_COMPLEX(y) && IS_COMPLEX(z)) \
        return GMPy_Complex_##NAME(x, y, z, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 3) { \
        TYPE_ERROR(#FUNC"() requires 3 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), \
                              PyTuple_GET_ITEM(args, 2), context); \
}

#define GMPY_MPFR_UNIOP(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result, *tempx; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_UNIOP_NOROUND(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result, *tempx; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_MPFR_Method_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT(context); \
    return GMPy_Number_##NAME(self, context); \
}\
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_UNIOP_NOROUND_NOMETHOD(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result, *tempx; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_UNIOP_EX(NAME, FUNC) \
static PyObject * \
_GMPy_MPFR_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result; \
    CHECK_CONTEXT(context); \
    if (!(result = GMPy_MPFR_New(0, context))) { \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, MPFR(x), GET_MPFR_ROUND(context)); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    PyObject *result, *tempx; \
    CHECK_CONTEXT(context); \
    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) { \
        return NULL; \
    } \
    result = _GMPy_MPFR_##NAME(tempx, context); \
    Py_DECREF(tempx); \
    return result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (MPFR_Check(x)) \
        return _GMPy_MPFR_##NAME(x, context); \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_UNIOP_TEMPLATE(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_UNIOP_TEMPLATE_EX(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (MPFR_Check(x)) \
        return _GMPy_MPFR_##NAME(x, context); \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *other) \
{ \
    CTXT_Object *context = NULL; \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(other, context); \
}

#define GMPY_MPFR_BINOP(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL, *tempy = NULL; \
    CHECK_CONTEXT(context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    tempy = GMPy_MPFR_From_Real(y, 1, context); \
    result = GMPy_MPFR_New(0, context); \
    if (!result || !tempx || !tempy) { \
        Py_XDECREF((PyObject*)tempx); \
        Py_XDECREF((PyObject*)tempy); \
        Py_XDECREF((PyObject*)result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    Py_DECREF((PyObject*)tempy); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    if (IS_REAL(x) && IS_REAL(y)) \
        return GMPy_Real_##NAME(x, y, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

/* Macro to support functions that require ('mpfr', 'int').
 * More precisely, the first argument must pass IS_REAL() and the second
 * argument must pass IS_INTEGER(). */

#define GMPY_MPFR_BINOP_REAL_ULONG(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    unsigned long n; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    n = c_ulong_From_Integer(y); \
    if (!result || !tempx || (n == (unsigned long)(-1) && PyErr_Occurred())) { \
        Py_XDECREF((PyObject*)tempx); \
        Py_XDECREF((PyObject*)result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, n, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    if (IS_REAL(x) && PyIntOrLong_Check(y)) \
        return GMPy_Real_##NAME(x, y, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

/* Macro to support functions that require ('mpfr', 'int').
 * More precisely, the first argument must pass IS_REAL() and the second
 * argument must pass IS_INTEGER(). The calling sequence passes n first
 * to the MPFR library.*/

#define GMPY_MPFR_BINOP_REAL_LONG(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    long n; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    n = c_long_From_Integer(y); \
    if (!result || !tempx || (n == -1 && PyErr_Occurred())) { \
        Py_XDECREF((PyObject*)tempx); \
        Py_XDECREF((PyObject*)result); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, n, tempx->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    if (IS_REAL(x) && PyIntOrLong_Check(y)) \
        return GMPy_Real_##NAME(x, y, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

#define GMPY_MPFR_BINOP_TEMPLATE(NAME, FUNC) \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    if (IS_REAL(x) && IS_REAL(y)) \
        return GMPy_Real_##NAME(x, y, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
} \

#define GMPY_MPFR_BINOP_EX(NAME, FUNC) \
static PyObject * \
_GMPy_MPFR_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    MPFR_Object *result; \
    CHECK_CONTEXT(context); \
    if (!(result = GMPy_MPFR_New(0, context))) { \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context)); \
    _GMPy_MPFR_Cleanup(&result, context); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    PyObject *result, *tempx, *tempy; \
    CHECK_CONTEXT(context); \
    tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context); \
    tempy = (PyObject*)GMPy_MPFR_From_Real(y, 1, context); \
    if (!tempx || !tempy) { \
        Py_XDECREF(tempx); \
        Py_XDECREF(tempy); \
        return NULL; \
    } \
    result = _GMPy_MPFR_##NAME(tempx, tempy, context); \
    Py_DECREF(tempx); \
    Py_DECREF(tempy); \
    return result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, PyObject *y, CTXT_Object *context) \
{ \
    if (MPFR_Check(x) && MPFR_Check(y)) \
        return _GMPy_MPFR_##NAME(x, y, context); \
    if (IS_REAL(x) && IS_REAL(y)) \
        return GMPy_Real_##NAME(x, y, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(#FUNC"() requires 2 arguments"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context); \
}
