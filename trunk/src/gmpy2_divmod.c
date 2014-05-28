/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_divmod.c                                                          *
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
    PyObject *result;
    MPZ_Object *tempx, *tempy, *rem, *quo;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    CHECK_CONTEXT(context);

    result = PyTuple_New(2);
    rem = GMPy_MPZ_New(context);
    quo = GMPy_MPZ_New(context);
    if (!result || !rem || !quo) {
        Py_XDECREF((PyObject*)rem);
        Py_XDECREF((PyObject*)quo);
        Py_XDECREF(result);
        return NULL;
    }

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_fdiv_qr(quo->z, rem->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_qr_ui(quo->z, rem->z, MPZ(x), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rem);
                Py_DECREF((PyObject*)quo);
                Py_DECREF(result);
                return NULL;
            }
            else {
                mpz_cdiv_qr_ui(quo->z, rem->z, MPZ(x), -temp_si);
                mpz_neg(quo->z, quo->z);
            }
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }
        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rem);
                Py_DECREF((PyObject*)quo);
                Py_DECREF(result);
                return NULL;
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
            Py_DECREF((PyObject*)rem);
            Py_DECREF((PyObject*)quo);
            Py_DECREF(result);
            return NULL;
        }
        mpz_inoc(tempz);
        mpz_set_PyIntOrLong(tempz, x);
        mpz_fdiv_qr(quo->z, rem->z, tempz, MPZ(y));
        mpz_cloc(tempz);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return (PyObject*)result;
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        tempx = GMPy_MPZ_From_Integer(x, context);
        tempy = GMPy_MPZ_From_Integer(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("could not convert Integer to mpz");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)rem);
            Py_DECREF((PyObject*)quo);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)rem);
            Py_DECREF((PyObject*)quo);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(quo->z, rem->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
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
    MPQ_Object *tempx, *tempy, *rem;
    MPZ_Object *quo;
    PyObject *result;

    CHECK_CONTEXT(context);

    result = PyTuple_New(2);
    rem = GMPy_MPQ_New(context);
    quo = GMPy_MPZ_New(context);
    if (!result || !rem || !quo) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)rem);
        Py_XDECREF((PyObject*)quo);
        return NULL;
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        tempx = GMPy_MPQ_From_Number(x, context);
        tempy = GMPy_MPQ_From_Number(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("could not convert Rational to mpq");
            goto error;
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

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)rem);
    Py_DECREF((PyObject*)quo);
    Py_DECREF(result);
    return NULL;
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

/* Compute the quotient and remainder of two mpfr numbers. Match Python's
 * behavior for handling signs and special values. Returns Py_NotImplemented
 * if both objects are not valid reals.
 */

static PyObject *
GMPy_Real_DivMod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL, *quo, *rem;
    PyObject *result;

    CHECK_CONTEXT(context);

    result = PyTuple_New(2);
    rem = GMPy_MPFR_New(0, context);
    quo = GMPy_MPFR_New(0, context);
    if (!result || !rem || !quo) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)quo);
        Py_XDECREF((PyObject*)rem);
        return NULL;
    }

    if (IS_REAL(x) && IS_REAL(y)) {
        tempx = GMPy_MPFR_From_Real(x, 1, context);
        tempy = GMPy_MPFR_From_Real(y, 1, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("could not convert Real to mpfr");
            goto error;
        }
        if (mpfr_zero_p(tempy->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("divmod() division by zero");
                goto error;
            }
        }

        if (mpfr_nan_p(tempx->f) || mpfr_nan_p(tempy->f) || mpfr_inf_p(tempx->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("divmod() invalid operation");
                goto error;
            }
            else {
                mpfr_set_nan(quo->f);
                mpfr_set_nan(rem->f);
            }
        }
        else if (mpfr_inf_p(tempy->f)) {
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
                rem->rc = mpfr_set(rem->f, tempx->f, GET_MPFR_ROUND(context));
            }
        }
        else {
            MPFR_Object *temp;

            if (!(temp = GMPy_MPFR_New(0, context))) {
                goto error;
            }
            mpfr_fmod(rem->f, tempx->f, tempy->f, GET_MPFR_ROUND(context));
            mpfr_sub(temp->f, tempx->f, rem->f, GET_MPFR_ROUND(context));
            mpfr_div(quo->f, temp->f, tempy->f, GET_MPFR_ROUND(context));

            if (!mpfr_zero_p(rem->f)) {
                if ((mpfr_sgn(tempy->f) < 0) != (mpfr_sgn(rem->f) < 0)) {
                    mpfr_add(rem->f, rem->f, tempy->f, GET_MPFR_ROUND(context));
                    mpfr_sub_ui(quo->f, quo->f, 1, GET_MPFR_ROUND(context));
                }
            }
            else {
                mpfr_copysign(rem->f, rem->f, tempy->f, GET_MPFR_ROUND(context));
            }

            if (!mpfr_zero_p(quo->f)) {
                mpfr_round(quo->f, quo->f);
            }
            else {
                mpfr_setsign(quo->f, quo->f,
                             mpfr_sgn(tempx->f) * mpfr_sgn(tempy->f) - 1,
                             GET_MPFR_ROUND(context));
            }
            Py_DECREF((PyObject*)temp);
        }
        GMPY_MPFR_CHECK_RANGE(quo, context);
        GMPY_MPFR_CHECK_RANGE(rem, context);
        GMPY_MPFR_SUBNORMALIZE(quo, context);
        GMPY_MPFR_SUBNORMALIZE(rem, context);


        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)rem);
    Py_DECREF((PyObject*)quo);
    Py_DECREF(result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)rem);
    Py_DECREF((PyObject*)quo);
    Py_DECREF(result);
    return NULL;
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

PyDoc_STRVAR(GMPy_doc_divmod,
"div_mod(x, y) -> (quotient, remainder)\n\n"
"Return divmod(x, y); uses alternate spelling to avoid conflicts.\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"mpfr arguments to div_mod().");

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
"context.divmod(x, y) -> (quotient, remainder)\n\n"
"Return divmod(x, y).\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"mpfr arguments to context.divmod().");

static PyObject *
GMPy_Context_DivMod(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("divmod() requires 2 arguments.");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_DivMod(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1),
                              context);
}

