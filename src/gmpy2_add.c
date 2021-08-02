/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_add.c                                                             *
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

/* This file implements the + operator, gmpy2.add(), and context.add().
 */

/* Add two Integer objects (see gmpy2_convert.h). If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted
 * into an mpz, a type error is returned.
 */

static PyObject *
GMPy_Integer_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype, 
                         CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPZANY(xtype)) {
        if (IS_TYPE_MPZANY(ytype)) {
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_add(result->z, MPZ(x), MPZ(y));
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            return (PyObject*)result;
        }

        if (IS_TYPE_PyInteger(ytype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(y, &error);

            if (!error) {
                if (temp >= 0) {
                    mpz_add_ui(result->z, MPZ(x), temp);
                }
                else {
                    mpz_sub_ui(result->z, MPZ(x), -temp);
                }
            }
            else {
                mpz_set_PyIntOrLong(result->z, y);
                GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
                mpz_add(result->z, MPZ(x), result->z);
                GMPY_MAYBE_END_ALLOW_THREADS(context);
            }
            return (PyObject*)result;
        }
    }

    if (IS_TYPE_MPZANY(ytype)) {
        if (IS_TYPE_PyInteger(xtype)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(x, &error);

            if (!error) {
                if (temp >= 0) {
                    mpz_add_ui(result->z, MPZ(y), temp);
                }
                else {
                    mpz_sub_ui(result->z, MPZ(y), -temp);
                }
            }
            else {
                mpz_set_PyIntOrLong(result->z, x);
                GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
                mpz_add(result->z, result->z, MPZ(y));
                GMPY_MAYBE_END_ALLOW_THREADS(context);
            }
            return (PyObject*)result;
        }
    }

    if (IS_TYPE_INTEGER(xtype) && (IS_TYPE_INTEGER(ytype))) {
        MPZ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_IntegerWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPZ_From_IntegerWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_add(result->z, tempx->z, tempy->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    TYPE_ERROR("add() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Rational_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                          CTXT_Object *context)
{
    MPQ_Object *result = NULL;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPQ(xtype) && IS_TYPE_MPQ(ytype)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_add(result->q, MPQ(x), MPQ(y));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        return (PyObject*)result;
    }

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype)) {
        MPQ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_RationalWithType(x, xtype, context)) ||
            !(tempy = GMPy_MPQ_From_RationalWithType(y, ytype, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpq_add(result->q, tempx->q, tempy->q);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    TYPE_ERROR("add() argument type not supported");
    return NULL;
}

/* Addition can be performed by the equivalent of mpfr.__add__ or by
 * gmpy2.add().
 *
 *   GMPy_Real_Add(x, y, context) returns x+y using the provided context. If
 *   provided context is NULL, then the current context is used. If an error
 *   occurs, NULL is returned and an exception is set. If either x or y can't
 *   be converted to an mpfr, then Py_NotImplemented is returned.
 *    GMPy_Real_Add() will not try to promote the result to a different type
 *   (i.e. mpc).
 */

/* Attempt to add two numbers and return an mpfr. The code path is optimized by
 * checking for mpfr objects first. Returns Py_NotImplemented if both objects
 * are not valid reals.  */

static PyObject *
GMPy_Real_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                      CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPFR(xtype) && IS_TYPE_MPFR(ytype)) {
        mpfr_clear_flags();
        result->rc = mpfr_add(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }
    
    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype)) {
        MPFR_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_RealWithType(x, xtype, 1, context)) ||
            !(tempy = GMPy_MPFR_From_RealWithType(y, ytype, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpfr_clear_flags();
        result->rc = mpfr_add(result->f, MPFR(tempx), MPFR(tempy), GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPFR_Cleanup(&result, context);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    TYPE_ERROR("add() argument type not supported");
    return NULL;
}

/* GMPy_Complex_Add(x, y, context) returns x+y using the provided context. If
 * context is NULL, then the current context is used. If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted to
 * an mpc, then Py_NotImplemented is returned. */


static PyObject *
GMPy_Complex_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype,
                         CTXT_Object *context)
{
    MPC_Object *result = NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (IS_TYPE_MPC(xtype) && IS_TYPE_MPC(ytype)) {
        result->rc = mpc_add(result->c, MPC(x), MPC(y), GET_MPC_ROUND(context));
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype)) {
        MPC_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPC_From_ComplexWithType(x, xtype, 1, 1, context)) ||
            !(tempy = GMPy_MPC_From_ComplexWithType(y, ytype, 1, 1, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        result->rc = mpc_add(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    TYPE_ERROR("add() argument type not supported");
    return NULL;
}

/* Implement all the slot methods here. */

static PyObject *
GMPy_Number_Add_Slot(PyObject *x, PyObject *y)
{
    CTXT_Object *context = NULL;
    CHECK_CONTEXT(context);

    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);
    
    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_AddWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_AddWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_AddWithType(x, xtype, y, ytype, context);
        
    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_AddWithType(x, xtype, y, ytype, context);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement context.add() and gmpy2.add(). */

PyDoc_STRVAR(GMPy_doc_function_add,
"add(x, y) -> number\n\n"
"Return x + y.");

static PyObject *
GMPy_Number_Add(PyObject *x, PyObject *y, CTXT_Object *context)
{
    int xtype = GMPy_ObjectType(x);
    int ytype = GMPy_ObjectType(y);
    
    if (IS_TYPE_INTEGER(xtype) && IS_TYPE_INTEGER(ytype))
        return GMPy_Integer_AddWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_RATIONAL(xtype) && IS_TYPE_RATIONAL(ytype))
        return GMPy_Rational_AddWithType(x, xtype, y, ytype, context);

    if (IS_TYPE_REAL(xtype) && IS_TYPE_REAL(ytype))
        return GMPy_Real_AddWithType(x, xtype, y, ytype, context);
        
    if (IS_TYPE_COMPLEX(xtype) && IS_TYPE_COMPLEX(ytype))
        return GMPy_Complex_AddWithType(x, xtype, y, ytype, context);

    TYPE_ERROR("add() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_add,
"context.add(x, y) -> number\n\n"
"Return x + y.");

static PyObject *
GMPy_Context_Add(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("add() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Add(PyTuple_GET_ITEM(args, 0),
                           PyTuple_GET_ITEM(args, 1),
                           context);
}

