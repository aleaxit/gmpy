/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_floordiv.c                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Case Van Horsen      *
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

/* This file implements the // operator, gmpy2.floor_div() and
 * context.floor_div().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_FloorDiv(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_FloorDiv_Slot
 *   GMPy_MPQ_FloorDiv_Slot
 *   GMPy_MPFR_FloorDiv_Slot
 *   GMPy_MPC_FloorDiv_Slot
 *
 *   GMPy_Integer_FloorDiv(Integer, Integer, context|NULL)
 *   GMPy_Rational_FloorDiv(Rational, Rational, context|NULL)
 *   GMPy_Real_FloorDiv(Real, Real, context|NULL)
 *   GMPy_Complex_FloorDiv(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_FloorDiv(context, args)
 *
 */

static PyObject *
GMPy_Integer_FloorDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPZ_Object *result;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_fdiv_q(result->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_q_ui(result->z, MPZ(x), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else {
                mpz_cdiv_q_ui(result->z, MPZ(x), -temp_si);
                mpz_neg(result->z, result->z);
            }
            return (PyObject*)result;
        }

        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpz_fdiv_q(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
    }

    if (CHECK_MPZANY(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (PyIntOrLong_Check(x)) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, x);
            mpz_fdiv_q(result->z, tempz, MPZ(y));
            mpz_cloc(tempz);
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        MPZ_Object *tempx, *tempy;

        tempx = GMPy_MPZ_From_Integer(x, context);
        tempy = GMPy_MPZ_From_Integer(y, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpz_fdiv_q(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement floor division for MPZ_Object. On entry, one of the two arguments
 * must be an MPZ_Object. If the other object is an Integer, return an
 * MPZ_Object. If the other object isn't an MPZ_Object, call the appropriate
 * function. If no appropriate function can be found, return NotImplemented.
 */

static PyObject *
GMPy_MPZ_FloorDiv_Slot(PyObject *x, PyObject *y)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_FloorDiv(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_FloorDiv(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_FloorDiv(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_FloorDiv(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Rational_FloorDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPZ_Object *result;
    MPQ_Object *tempq;

    CHECK_CONTEXT_SET_EXPONENT(context);

    result = GMPy_MPZ_New(context);
    tempq = GMPy_MPQ_New(context);
    if (!result || !tempq) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempq);
        return NULL;
    }

    if (MPQ_Check(x) && MPQ_Check(y)) {
        if (mpq_sgn(MPQ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }
        mpq_div(tempq->q, MPQ(x), MPQ(y));
        mpz_fdiv_q(result->z, mpq_numref(tempq->q), mpq_denref(tempq->q));
        Py_DECREF((PyObject*)tempq);
        return (PyObject*)result;
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        MPQ_Object *tempx, *tempy;

        tempx = GMPy_MPQ_From_Number(x, context);
        tempy = GMPy_MPQ_From_Number(y, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            goto error;
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            goto error;
        }

        mpq_div(tempq->q, tempx->q, tempy->q);
        mpz_fdiv_q(result->z, mpq_numref(tempx->q), mpq_denref(tempy->q));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        Py_DECREF((PyObject*)tempq);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_DECREF((PyObject*)result);
    Py_DECREF((PyObject*)tempq);
    return NULL;
}

static PyObject *
GMPy_MPQ_FloorDiv_Slot(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_FloorDiv(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_FloorDiv(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_FloorDiv(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Attempt floor division of two numbers and return an mpfr. The code path is
 * optimized by checking for mpfr objects first. Returns Py_NotImplemented if
 * both objects are not valid reals.  */

static PyObject *
GMPy_Real_FloorDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    /* This only processes mpfr if the exponent is still in-bounds. Need
     * to handle the rare case at the end. */

    if (MPFR_CheckAndExp(x) && MPFR_CheckAndExp(y)) {
        mpfr_clear_flags();
        result->rc = mpfr_div(result->f, MPFR(x), MPFR(y),
                              GET_MPFR_ROUND(context));
        result->rc = mpfr_floor(result->f, result->f);
        goto done;
    }

    if (MPFR_CheckAndExp(x)) {
        if (PyIntOrLong_Check(y)) {
            mpz_t tempz;
            mpir_si temp_si;
            int overflow;

            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpfr_clear_flags();
                result->rc = mpfr_div_z(result->f, MPFR(x),
                                        tempz, GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                result->rc = mpfr_floor(result->f, result->f);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_div_si(result->f, MPFR(x),
                                         temp_si, GET_MPFR_ROUND(context));
                result->rc = mpfr_floor(result->f, result->f);
                goto done;
            }
        }

        if (CHECK_MPZANY(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_div_z(result->f, MPFR(x),
                                    MPZ(y), GET_MPFR_ROUND(context));
            result->rc = mpfr_floor(result->f, result->f);
            goto done;
        }

        if (IS_RATIONAL(y)) {
            MPQ_Object *tempy;

            if (!(tempy = GMPy_MPQ_From_Number(y, context))) {
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_div_q(result->f, MPFR(x), tempy->q,
                                    GET_MPFR_ROUND(context));
            result->rc = mpfr_floor(result->f, result->f);
            Py_DECREF((PyObject*)tempy);
            goto done;
        }

        if (PyFloat_Check(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_div_d(result->f, MPFR(x),
                                    PyFloat_AS_DOUBLE(y), GET_MPFR_ROUND(context));
            result->rc = mpfr_floor(result->f, result->f);
            goto done;
        }
    }

    if (MPFR_CheckAndExp(y)) {
        if (PyIntOrLong_Check(x)) {
            mpir_si temp_si;
            int overflow;

            temp_si = PyLong_AsSIAndOverflow(x, &overflow);
            if (!overflow) {
                mpfr_clear_flags();
                result->rc = mpfr_si_div(result->f, temp_si, MPFR(y),
                                         GET_MPFR_ROUND(context));
                result->rc = mpfr_floor(result->f, result->f);
                goto done;
            }
        }

        /* Since mpfr_z_div does not exist, this combination is handled at the
         * end by converting x to an mpfr. Ditto for rational.*/

        if (PyFloat_Check(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_d_div(result->f, PyFloat_AS_DOUBLE(x),
                                    MPFR(y), GET_MPFR_ROUND(context));
            result->rc = mpfr_floor(result->f, result->f);
            goto done;
        }
    }

    /* In addition to handling PyFloat + PyFloat, the rare case when the
     * exponent bounds have been changed is handled here. See
     * Pympfr_From_Real() for details. */

    if (IS_REAL(x) && IS_REAL(y)) {
        MPFR_Object *tempx, *tempy;

        tempx = GMPy_MPFR_From_Real(x, 1, context);
        tempy = GMPy_MPFR_From_Real(y, 1, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpfr_clear_flags();
        result->rc = mpfr_div(result->f, MPFR(tempx), MPFR(tempy),
                              GET_MPFR_ROUND(context));
        result->rc = mpfr_floor(result->f, result->f);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    MPFR_CLEANUP_2(result, context, "division");
    return (PyObject*)result;
}

static PyObject *
GMPy_MPFR_FloorDiv_Slot(PyObject *x, PyObject *y)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_FloorDiv(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_FloorDiv(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Complex_FloorDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    TYPE_ERROR("can't take floor of complex number");
    return NULL;
}

static PyObject *
GMPy_MPC_FloorDiv_Slot(PyObject *x, PyObject *y)
{
    return GMPy_Complex_FloorDiv(x, y, NULL);
}

PyDoc_STRVAR(GMPy_doc_floordiv,
"floor_div(x, y) -> number\n\n"
"Return x // y; uses floor division.");

static PyObject *
GMPy_Number_FloorDiv(PyObject *x, PyObject *y, CTXT_Object *context)
{
    LOAD_CONTEXT_SET_EXPONENT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_FloorDiv(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_FloorDiv(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_FloorDiv(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_FloorDiv(x, y, context);

    TYPE_ERROR("floor_div() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_floordiv,
"context.floor_div(x, y) -> number\n\n"
"Return x // y; uses floor division.");

static PyObject *
GMPy_Context_FloorDiv(PyObject *self, PyObject *args)
{
    PyObject *result;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("floor_div() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        /* If we are passed a read-only context, make a copy of it before
         * proceeding. Remember to decref context when we're done. */

        if (((CTXT_Object*)self)->ctx.readonly) {
            context = (CTXT_Object*)GMPy_CTXT_Copy(self, NULL);
            if (!context)
                return NULL;
        }
        else {
            context = (CTXT_Object*)self;
            Py_INCREF((PyObject*)context);
        }
    }
    else {
        CHECK_CONTEXT_SET_EXPONENT(context);
        Py_INCREF((PyObject*)context);
    }

    result = GMPy_Number_FloorDiv(PyTuple_GET_ITEM(args, 0),
                                  PyTuple_GET_ITEM(args, 1),
                                  context);
    Py_DECREF((PyObject*)context);
    return result;
}

