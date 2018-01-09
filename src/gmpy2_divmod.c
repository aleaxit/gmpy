/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_divmod.c                                                          *
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

/* This file implements __divmod__ and context.divmod().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_DivMod(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_DivMod_Slot
 *   GMPy_MPQ_DivMod_Slot
 *   GMPy_MPFR_DivMod_Slot
 *   GMPy_MPC_DivMod_Slot
 *
 *   GMPy_Integer_DivMod(Integer, Integer, context|NULL)
 *   GMPy_Rational_DivMod(Rational, Rational, context|NULL)
 *   GMPy_Real_DivMod(Real, Real, context|NULL)
 *   GMPy_Complex_DivMod(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_DivMod(context, args)
 *
 */

static PyObject *
GMPy_Integer_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    PyObject *result = NULL;
    MPZ_Object *tempx = NULL, *tempy = NULL, *rem = NULL, *quo = NULL;

    if (!(result = PyTuple_New(2)) ||
        !(rem = GMPy_MPZ_New(context)) ||
        !(quo = GMPy_MPZ_New(context))) {

        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            int error;
            native_si temp = GMPy_Integer_AsNative_siAndError(y, &error);

            if (error) {
                mpz_set_PyIntOrLong(global.tempz, y);
                mpz_fdiv_qr(quo->z, rem->z, MPZ(x), global.tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_qr_ui(quo->z, rem->z, MPZ(x), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("division or modulo by zero");
                goto error;
            }
            else {
                mpz_cdiv_qr_ui(quo->z, rem->z, MPZ(x), -temp);
                mpz_neg(quo->z, quo->z);
            }
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }

        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                goto error;
            }
            mpz_fdiv_qr(quo->z, rem->z, MPZ(x), MPZ(y));
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }
    }

    if (CHECK_MPZANY(y) && PyIntOrLong_Check(x)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, x);
            mpz_fdiv_qr(quo->z, rem->z, global.tempz, MPZ(y));
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {

        if (!(tempx = GMPy_MPZ_From_Integer(x, context)) ||
            !(tempy = GMPy_MPZ_From_Integer(y, context))) {

            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }
        mpz_fdiv_qr(quo->z, rem->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Integer_DivMod().");
    /* LCOV_EXCL_STOP */
  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)rem);
    Py_XDECREF((PyObject*)quo);
    Py_XDECREF(result);
    return NULL;
}

static PyObject *
GMPy_MPZ_DivMod_Slot(PyObject *x, PyObject *y)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_DivMod(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_DivMod(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_DivMod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_DivMod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Rational_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPQ_Object *tempx = NULL, *tempy = NULL, *rem = NULL;
    MPZ_Object *quo = NULL;
    PyObject *result = NULL;

    if (!(result = PyTuple_New(2)) ||
        !(rem = GMPy_MPQ_New(context)) ||
        !(quo = GMPy_MPZ_New(context))) {

        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {

        if (!(tempx = GMPy_MPQ_From_Number(x, context)) ||
            !(tempy = GMPy_MPQ_From_Number(y, context))) {

            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }

        mpq_div(rem->q, tempx->q, tempy->q);
        mpz_fdiv_q(quo->z, mpq_numref(rem->q), mpq_denref(rem->q));
        /* Need to calculate x - quo * y. */
        mpq_set_z(rem->q, quo->z);
        mpq_mul(rem->q, rem->q, tempy->q);
        mpq_sub(rem->q, tempx->q, rem->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Rational_DivMod().");
  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)rem);
    Py_XDECREF((PyObject*)quo);
    Py_XDECREF(result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_MPQ_DivMod_Slot(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_DivMod(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_DivMod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_DivMod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Real_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL, *quo = NULL, *rem = NULL, *temp;
    PyObject *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = PyTuple_New(2)) ||
        !(rem = GMPy_MPFR_New(0, context)) ||
        !(quo = GMPy_MPFR_New(0, context))) {

        /* LCOV_EXCL_START */
        goto error;
        /* LCOV_EXCL_STOP */
    }

    if (IS_REAL(x) && IS_REAL(y)) {

        if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)) ||
            !(tempy = GMPy_MPFR_From_Real(y, 1, context))) {

            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        if (mpfr_zero_p(tempy->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("divmod() division by zero");
                goto error;
            }
            mpfr_set_nan(quo->f);
            mpfr_set_nan(rem->f);
            goto okay;
        }

        if (mpfr_nan_p(tempx->f) || mpfr_nan_p(tempy->f) || mpfr_inf_p(tempx->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("divmod() invalid operation");
                goto error;
            }
            mpfr_set_nan(quo->f);
            mpfr_set_nan(rem->f);
            goto okay;
        }

        if (mpfr_inf_p(tempy->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("divmod() invalid operation");
                goto error;
            }
            if (mpfr_zero_p(tempx->f)) {
                mpfr_set_zero(quo->f, mpfr_sgn(tempy->f));
                mpfr_set_zero(rem->f, mpfr_sgn(tempy->f));
            }
            else if ((mpfr_signbit(tempx->f)) != (mpfr_signbit(tempy->f))) {
                mpfr_set_si(quo->f, -1, MPFR_RNDN);
                mpfr_set_inf(rem->f, mpfr_sgn(tempy->f));
            }
            else {
                mpfr_set_si(quo->f, 0, MPFR_RNDN);
                rem->rc = mpfr_set(rem->f, tempx->f, MPFR_RNDN);
            }
            goto okay;
        }

        if (!(temp = GMPy_MPFR_New(0, context))) {
            /* LCOV_EXCL_START */
            goto error;
            /* LCOV_EXCL_STOP */
        }
        mpfr_fmod(rem->f, tempx->f, tempy->f, MPFR_RNDN);
        mpfr_sub(temp->f, tempx->f, rem->f, MPFR_RNDN);
        mpfr_div(quo->f, temp->f, tempy->f, MPFR_RNDN);

        if (!mpfr_zero_p(rem->f)) {
            if ((mpfr_sgn(tempy->f) < 0) != (mpfr_sgn(rem->f) < 0)) {
                mpfr_add(rem->f, rem->f, tempy->f, MPFR_RNDN);
                mpfr_sub_ui(quo->f, quo->f, 1, MPFR_RNDN);
            }
        }
        else {
            mpfr_copysign(rem->f, rem->f, tempy->f, MPFR_RNDN);
        }

        if (!mpfr_zero_p(quo->f)) {
            mpfr_round(quo->f, quo->f);
        }
        else {
            mpfr_setsign(quo->f, quo->f, mpfr_sgn(tempx->f) * mpfr_sgn(tempy->f) - 1, MPFR_RNDN);
        }
        Py_DECREF((PyObject*)temp);

        GMPY_MPFR_CHECK_RANGE(quo, context);
        GMPY_MPFR_CHECK_RANGE(rem, context);
        GMPY_MPFR_SUBNORMALIZE(quo, context);
        GMPY_MPFR_SUBNORMALIZE(rem, context);
      okay:
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return (PyObject*)result;
    }

    /* LCOV_EXCL_START */
    SYSTEM_ERROR("Internal error in GMPy_Real_DivMod_1().");
  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)rem);
    Py_XDECREF((PyObject*)quo);
    Py_XDECREF(result);
    return NULL;
    /* LCOV_EXCL_STOP */
}

static PyObject *
GMPy_MPFR_DivMod_Slot(PyObject *x, PyObject *y)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_DivMod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_DivMod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Complex_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    TYPE_ERROR("can't take floor or mod of complex number.");
    return NULL;
}

static PyObject *
GMPy_MPC_DivMod_Slot(PyObject *x, PyObject *y)
{
    return GMPy_Complex_DivMod(x, y, NULL);
}

static PyObject *
GMPy_Number_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_DivMod(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_DivMod(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_DivMod(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_DivMod(x, y, context);

    TYPE_ERROR("divmod() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_divmod,
"context.div_mod(x, y) -> (quotient, remainder)\n\n"
"Return div_mod(x, y); uses alternate spelling to avoid naming conflicts.\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"mpfr arguments to context.div_mod().");

static PyObject *
GMPy_Context_DivMod(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("div_mod() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_DivMod(PyTuple_GET_ITEM(args, 0),
                              PyTuple_GET_ITEM(args, 1),
                              context);
}

