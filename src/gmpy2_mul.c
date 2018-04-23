/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mul.c                                                             *
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

/* This file implements the * operator, gmpy2.mul() and context.mul().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_Mul(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Mul_Slot
 *   GMPy_MPQ_Mul_Slot
 *   GMPy_MPFR_Mul_Slot
 *   GMPy_MPC_Mul_Slot
 *
 *   GMPy_Integer_Mul(Integer, Integer, context|NULL)
 *   GMPy_Rational_Mul(Rational, Rational, context|NULL)
 *   GMPy_Real_Mul(Real, Real, context|NULL)
 *   GMPy_Complex_Mul(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Mul(context, args)
 *
 */

/* Multiply two Integer objects (see gmpy2_convert.c). If an error occurs,
 * NULL is returned and an exception is set. If either x or y can't be
 * converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
GMPy_Integer_Mul(PyObject *x, PyObject *y, CTXT_Object *context)
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
                mpz_mul_si(result->z, MPZ(x), temp);
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, y);
                mpz_mul(result->z, MPZ(x), global.tempz);
            }
            return (PyObject*)result;
        }

        if (MPZ_Check(y)) {
            mpz_mul(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
    }

    if (MPZ_Check(y)) {
        if (PyIntOrLong_Check(x)) {
            int error;
            native_si temp = GMPy_Integer_AsNative_siAndError(x, &error);

            if (!error) {
                mpz_mul_si(result->z, MPZ(y), temp);
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, x);
                mpz_mul(result->z, MPZ(y), global.tempz);
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

        mpz_mul(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Integer_Mul().");
    Py_DECREF((PyObject*)result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Implement __mul__ for MPZ_Object. On entry, one of the two arguments must
 * be an MPZ_Object. If the other object is an Integer, add and return an
 * MPZ_Object. If the other object isn't an MPZ_Object, call the appropriate
 * function. If no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_MPZ_Mul_Slot(PyObject *x, PyObject *y)
{
    if (MPZ_Check(x) && MPZ_Check(y)) {
        MPZ_Object *result = NULL;

        if ((result = GMPy_MPZ_New(NULL))) {
            mpz_mul(result->z, MPZ(x), MPZ(y));
        }
        return (PyObject*)result;
    }

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Mul(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mul(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Multiply two Rational objects (see convert.c/IS_RATIONAL). Returns None and
 * raises TypeError if both objects are not valid rationals. GMPy_Rational_Mul
 * is intended to be called from GMPy_Number_Mul. */

static PyObject *
GMPy_Rational_Mul(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPQ_Object *result = NULL;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPQ_Check(x) && MPQ_Check(y)) {
        mpq_mul(result->q, MPQ(x), MPQ(y));
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

        mpq_mul(result->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Rational_Mul().");
    Py_DECREF((PyObject*)result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* Implement __mul__ for Pympq. On entry, one of the two arguments must
 * be a Pympq. If the other object is a Rational, add and return a Pympq.
 * If the other object isn't a Pympq, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_MPQ_Mul_Slot(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mul(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Attempt to multiply two numbers and return an mpfr. The code path is
 * optimized by checking for mpfr objects first. Returns Py_NotImplemented if
 * both objects are not valid reals.  */

static PyObject *
GMPy_Real_Mul(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPFR_Check(x) && MPFR_Check(y)) {
        mpfr_clear_flags();

        result->rc = mpfr_mul(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
        goto done;
    }

    if (MPFR_Check(x)) {
        if (PyIntOrLong_Check(y)) {
            int error;
            long temp = GMPy_Integer_AsLongAndError(y, &error);

            if (!error) {
                mpfr_clear_flags();

                result->rc = mpfr_mul_si(result->f, MPFR(x), temp, GET_MPFR_ROUND(context));
                goto done;
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, y);
                mpfr_clear_flags();

                result->rc = mpfr_mul_z(result->f, MPFR(x), global.tempz, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(y)) {
            mpfr_clear_flags();

            result->rc = mpfr_mul_z(result->f, MPFR(x), MPZ(y), GET_MPFR_ROUND(context));
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

            result->rc = mpfr_mul_q(result->f, MPFR(x), tempy->q, GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempy);
            goto done;
        }

        if (PyFloat_Check(y)) {
            mpfr_clear_flags();

            result->rc = mpfr_mul_d(result->f, MPFR(x), PyFloat_AS_DOUBLE(y), GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (MPFR_Check(y)) {
        if (PyIntOrLong_Check(x)) {
            int error;
            long temp = GMPy_Integer_AsLongAndError(x, &error);

            if (!error) {
                mpfr_clear_flags();

                result->rc = mpfr_mul_si(result->f, MPFR(y), temp, GET_MPFR_ROUND(context));
                goto done;
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, x);
                mpfr_clear_flags();

                result->rc = mpfr_mul_z(result->f, MPFR(y), global.tempz, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(x)) {
            mpfr_clear_flags();

            result->rc = mpfr_mul_z(result->f, MPFR(y), MPZ(x), GET_MPFR_ROUND(context));
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

            result->rc = mpfr_mul_q(result->f, MPFR(y), tempx->q, GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempx);
            goto done;
        }

        if (PyFloat_Check(x)) {
            mpfr_clear_flags();

            result->rc = mpfr_mul_d(result->f, MPFR(y), PyFloat_AS_DOUBLE(x), GET_MPFR_ROUND(context));
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

        result->rc = mpfr_mul(result->f, MPFR(tempx), MPFR(tempy), GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    SYSTEM_ERROR("Internal error in GMPy_Real_Mul().");
    return NULL;
    /* LCOV_EXCL_STOP */

  done:
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

/* Implement __mul__ for Pympfr. On entry, one of the two arguments must
 * be a Pympfr. If the other object is a Real, add and return a Pympfr.
 * If the other object isn't a Pympfr, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_MPFR_Mul_Slot(PyObject *x, PyObject *y)
{
    if (MPFR_Check(x) && MPFR_Check(y)) {
        MPFR_Object *result;
        CTXT_Object *context = NULL;

        CHECK_CONTEXT(context);

        if ((result = GMPy_MPFR_New(0, context))) {
            mpfr_clear_flags();

            result->rc = mpfr_mul(result->f, MPFR(x), MPFR(y), GET_MPFR_ROUND(context));
            _GMPy_MPFR_Cleanup(&result, context);
        }
        return (PyObject*)result;
    }

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* GMPy_Complex_Mul(x, y, context) returns x*y using the provided context. If
 * an error occurs, NULL is returned and an exception is set. If either x or
 * y can't be converted to an mpc, then Py_NotImplemented is returned. */

static PyObject *
GMPy_Complex_Mul(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (MPC_Check(x) && MPC_Check(y)) {

        result->rc = mpc_mul(result->c, MPC(x), MPC(y), GET_MPC_ROUND(context));

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

        result->rc = mpc_mul(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);

        _GMPy_MPC_Cleanup(&result, context);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    Py_DECREF((PyObject*)result);
    SYSTEM_ERROR("Internal error in GMPy_Complex_Mul().");
    return NULL;
    /* LCOV_EXCL_STOP */
}

/* GMPy_MPZ_Mul_Slot() is called by mpc.__mul__. It just gets a borrowed reference
 * to the current context and call Pympc_Mul_Complex(). Since mpc is the last
 * step of the numeric ladder, the NotImplemented return value from
 * Pympc_Add_Complex() is correct and is just passed on. */

static PyObject *
GMPy_MPC_Mul_Slot(PyObject *x, PyObject *y)
{
    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Number_Mul(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Mul(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mul(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, context);

    TYPE_ERROR("mul() argument type not supported");
    return NULL;
}

/* Implement context.mul() and gmpy2.mul(). */

PyDoc_STRVAR(GMPy_doc_function_mul,
"mul(x, y) -> number\n\n"
"Return x * y.");

PyDoc_STRVAR(GMPy_doc_context_mul,
"context.mul(x, y) -> number\n\n"
"Return x * y.");

static PyObject *
GMPy_Context_Mul(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mul() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Mul(PyTuple_GET_ITEM(args, 0),
                           PyTuple_GET_ITEM(args, 1),
                           context);
}

