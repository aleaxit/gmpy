/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_add.c                                                             *
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

/* This file implements the + operator, gmpy2.add(), and context.add().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL
 * value for context implies the function should use the currently active
 * context.
 *
 *   GMPy_Number_Add(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Add_Slot
 *   GMPy_MPQ_Add_Slot
 *   GMPy_MPFR_Add_Slot
 *   GMPy_MPC_Add_Slot
 *
 *   GMPy_Integer_Add(Integer, Integer, context|NULL)
 *   GMPy_Rational_Add(Rational, Rational, context|NULL)
 *   GMPy_Real_Add(Real, Real, context|NULL)
 *   GMPy_Complex_Add(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Add(context, args)
 *
 */

/* Add two Integer objects (see gmpy2_convert.h). If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted
 * into an mpz, Py_NotImplemented is returned.
 */

static PyObject *
GMPy_Integer_Add(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPZ_Check(x)) {
        if (PyIntOrLong_Check(y)) {
            int error;
            native_si temp = GMPy_Integer_AsNative_siAndError(y, &error);

            if (!error) {
                if (temp >= 0) {
                    mpz_add_ui(result->z, MPZ(x), temp);
                }
                else {
                    mpz_sub_ui(result->z, MPZ(x), -temp);
                }
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, y);
                mpz_add(result->z, MPZ(x), global.tempz);
            }
            return (PyObject*)result;
        }

        if (MPZ_Check(y)) {
            mpz_add(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
    }

    if (MPZ_Check(y)) {
        if (PyIntOrLong_Check(x)) {
            int error;
            native_si temp = GMPy_Integer_AsNative_siAndError(x, &error);

            if (!error) {
                if (temp >= 0) {
                    mpz_add_ui(result->z, MPZ(y), temp);
                }
                else {
                    mpz_sub_ui(result->z, MPZ(y), -temp);
                }
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, x);
                mpz_add(result->z, MPZ(y), global.tempz);
            }
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        MPZ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPZ_From_Integer(x, context)) ||
            !(tempy = GMPy_MPZ_From_Integer(y, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpz_add(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Integer_Add().");
    Py_DECREF((PyObject*)result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Implement __add__ for MPZ_Object. On entry, one of the two arguments must
 * be an MPZ_Object. If the other object is an Integer, add and return an
 * MPZ_Object. If the other object isn't an MPZ_Object, call the appropriate
 * function. If no appropriate function can be found, return NotImplemented.
 */

static PyObject *
GMPy_MPZ_Add_Slot(PyObject *x, PyObject *y)
{
    if (MPZ_Check(x) && MPZ_Check(y)) {
        MPZ_Object *result = NULL;

        if ((result = GMPy_MPZ_New(NULL))) {
            mpz_add(result->z, MPZ(x), MPZ(y));
        }
        return (PyObject*)result;
    }

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Add(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Add(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Add two Rational objects (see gmpy2_convert.h). Returns None and
 * raises TypeError if both objects are not valid rationals. GMPy_Rational_Add
 * is intended to be called from GMPy_Number_Add(). */

static PyObject *
GMPy_Rational_Add(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPQ_Object *result = NULL;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPQ_Check(x) && MPQ_Check(y)) {
        mpq_add(result->q, MPQ(x), MPQ(y));
        return (PyObject*)result;
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        MPQ_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPQ_From_Number(x, context)) ||
            !(tempy = GMPy_MPQ_From_Number(y, context))) {
            /* LCOV_EXCL_START */
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
            /* LCOV_EXCL_STOP */
        }

        mpq_add(result->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Rational_Add().");
    Py_DECREF((PyObject*)result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Implement __add__ for Pympq. On entry, one of the two arguments must
 * be a Pympq. If the other object is a Rational, add and return a Pympq.
 * If the other object isn't a Pympq, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_MPQ_Add_Slot(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Add(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
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
 *
 *   GMPy_mpfr_add_fast(x, y) is the entry point for mpfr.__add__.
 */

/* Attempt to add two numbers and return an mpfr. The code path is optimized by
 * checking for mpfr objects first. Returns Py_NotImplemented if both objects
 * are not valid reals.  */

static PyObject *
GMPy_Real_Add(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPFR_Check(x)) {
        if (MPFR_Check(y)) {
            mpfr_clear_flags();

            result->rc = mpfr_add(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
            goto done;
        }

        if (PyIntOrLong_Check(y)) {
            int error;
            long temp = GMPy_Integer_AsLongAndError(y, &error);

            if (!error) {
                mpfr_clear_flags();

                result->rc = mpfr_add_si(result->f, MPFR(x), temp, GET_MPFR_ROUND(context));
                goto done;
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, y);
                mpfr_clear_flags();

                result->rc = mpfr_add_z(result->f, MPFR(x), global.tempz, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(y)) {
            mpfr_clear_flags();

            result->rc = mpfr_add_z(result->f, MPFR(x), MPZ(y), GET_MPFR_ROUND(context));
            goto done;
        }

        if (IS_RATIONAL(y)) {
            MPQ_Object *tempy = NULL;

            if (!(tempy = GMPy_MPQ_From_Number(y, context))) {
                /* LCOV_EXCL_START */
                Py_DECREF((PyObject*)result);
                return NULL;
                /* LCOV_EXCL_STOP */
            }

            mpfr_clear_flags();

            result->rc = mpfr_add_q(result->f, MPFR(x), tempy->q, GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempy);
            goto done;
        }

        if (PyFloat_Check(y)) {
            mpfr_clear_flags();

            result->rc = mpfr_add_d(result->f, MPFR(x), PyFloat_AS_DOUBLE(y), GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (MPFR_Check(y)) {
        if (PyIntOrLong_Check(x)) {
            int error;
            long temp = GMPy_Integer_AsLongAndError(x, &error);

            if (!error) {
                mpfr_clear_flags();

                result->rc = mpfr_add_si(result->f, MPFR(y), temp, GET_MPFR_ROUND(context));
                goto done;
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, x);
                mpfr_clear_flags();

                result->rc = mpfr_add_z(result->f, MPFR(y), global.tempz, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(x)) {
            mpfr_clear_flags();

            result->rc = mpfr_add_z(result->f, MPFR(y), MPZ(x), GET_MPFR_ROUND(context));
            goto done;
        }

        if (IS_RATIONAL(x)) {
            MPQ_Object *tempx = NULL;

            if (!(tempx = GMPy_MPQ_From_Number(x, context))) {
                /* LCOV_EXCL_START */
                Py_DECREF((PyObject*)result);
                return NULL;
                /* LCOV_EXCL_STOP */
            }

            mpfr_clear_flags();

            result->rc = mpfr_add_q(result->f, MPFR(y), tempx->q, GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempx);
            goto done;
        }

        if (PyFloat_Check(x)) {
            mpfr_clear_flags();

            result->rc = mpfr_add_d(result->f, MPFR(y), PyFloat_AS_DOUBLE(x), GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (IS_REAL(x) && IS_REAL(y)) {
        MPFR_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)) ||
            !(tempy = GMPy_MPFR_From_Real(y, 1, context))) {
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
        goto done;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    SYSTEM_ERROR("Internal error in GMPy_Real_Add().");
    return NULL;
    /* LCOV_EXCL_STOP */

  done:
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

/* Implement __add__ for Pympfr. On entry, one of the two arguments must
 * be a Pympfr. If the other object is a Real, add and return a Pympfr.
 * If the other object isn't a Pympfr, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_MPFR_Add_Slot(PyObject *x, PyObject *y)
{
    if ((MPFR_Check(x)) && MPFR_Check(y)) {
        MPFR_Object *result;
        CTXT_Object *context = NULL;

        CHECK_CONTEXT(context);

        if ((result = GMPy_MPFR_New(0, context))) {
            mpfr_clear_flags();

            result->rc = mpfr_add(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
            _GMPy_MPFR_Cleanup(&result, context);
        }
        return (PyObject*)result;
    }

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* GMPy_Complex_Add(x, y, context) returns x+y using the provided context. If
 * context is NULL, then the current context is used. If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted to
 * an mpc, then Py_NotImplemented is returned. */

static PyObject *
GMPy_Complex_Add(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPC_Check(x) && MPC_Check(y)) {

        result->rc = mpc_add(result->c, MPC(x), MPC(y), GET_MPC_ROUND(context));

        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    if (IS_COMPLEX(x) && IS_COMPLEX(y)) {
        MPC_Object *tempx = NULL, *tempy = NULL;

        if (!(tempx = GMPy_MPC_From_Complex(x, 1, 1, context)) ||
            !(tempy = GMPy_MPC_From_Complex(y, 1, 1, context))) {
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

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    SYSTEM_ERROR("Internal error in GMPy_Complex_Add().");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* GMPy_MPC_Add_Slot() is called by mpc.__add__. It just gets a borrowed reference
 * to the current context and call Pympc_Add_Complex(). Since mpc is the last
 * step of the numeric ladder, the NotImplemented return value from
 * Pympc_Add_Complex() is correct and is just passed on. */

static PyObject *
GMPy_MPC_Add_Slot(PyObject *x, PyObject *y)
{
    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Number_Add(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Add(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Add(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, context);

    TYPE_ERROR("add() argument type not supported");
    return NULL;
}

/* Implement context.add() and gmpy2.add(). */

PyDoc_STRVAR(GMPy_doc_function_add,
"add(x, y) -> number\n\n"
"Return x + y.");

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

