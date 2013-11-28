/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mod.c                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013 Case Van Horsen            *
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

/* This file implements __mod__, gmpy2.mod(), and context.mod().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_Mod(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_mpz_mod_fast
 *   GMPy_mpq_mod_fast
 *   GMPy_mpfr_mod_fast
 *   GMPy_mpc_mod_fast
 *
 *   GMPy_Integer_Mod(Integer, Integer, context|NULL)
 *   GMPy_Rational_Mod(Rational, Rational, context|NULL)
 *   GMPy_Real_Mod(Real, Real, context|NULL)
 *   GMPy_Complex_Mod(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Mod(context, args)
 *
 */

static PyObject *
GMPy_Integer_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPZ_Object *tempx, *tempy, *result;
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
                mpz_fdiv_r(result->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_r_ui(result->z, MPZ(x), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(result->z, MPZ(x), -temp_si);
            }
            return (PyObject*)result;
        }
        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpz_fdiv_r(result->z, MPZ(x), MPZ(y));
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
            mpz_fdiv_r(result->z, tempz, MPZ(y));
            mpz_cloc(tempz);
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        tempx = GMPy_MPZ_From_Integer_Temp(x, context);
        tempy = GMPy_MPZ_From_Integer_Temp(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Integer to mpz.");
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
        mpz_fdiv_r(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_mpz_mod_fast(PyObject *x, PyObject *y)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Mod(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mod(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Rational_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    mpz_t tempz;
    MPQ_Object *tempx, *tempy, *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        tempx = GMPy_MPQ_From_Number_Temp(x, context);
        tempy = GMPy_MPQ_From_Number_Temp(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Rational to mpq.");
            goto error;
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }

        mpz_inoc(tempz);
        mpq_div(result->q, tempx->q, tempy->q);
        mpz_fdiv_q(tempz, mpq_numref(result->q), mpq_denref(result->q));
        /* Need to calculate x - tempz * y. */
        mpq_set_z(result->q, tempz);
        mpq_mul(result->q, result->q, tempy->q);
        mpq_sub(result->q, tempx->q, result->q);
        mpz_cloc(tempz);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PyObject *
GMPy_mpq_mod_fast(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mod(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Compute the quotient and remainder of two mpfr numbers. Match Python's
 * behavior for handling signs and special values. Returns Py_NotImplemented
 * if both objects are not valid reals.
 */

static PyObject *
GMPy_Real_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx, *tempy, *temp, *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(0, context);
    if (!result || !temp) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)temp);
        return NULL;
    }

    if (IS_REAL(x) && IS_REAL(y)) {
        tempx = GMPy_MPFR_From_Real_Temp(x, 0, context);
        tempy = GMPy_MPFR_From_Real_Temp(y, 0, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            goto error;
        }
        if (mpfr_zero_p(tempy->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("'mpfr' division by zero in modulo");
                goto error;
            }
        }
        mpfr_clear_flags();
        if (mpfr_nan_p(tempx->f) || mpfr_nan_p(tempy->f) || mpfr_inf_p(tempx->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("'mpfr' invalid operation in modulo");
                goto error;
            }
            else {
                mpfr_set_nan(result->f);
            }
        }
        else if (mpfr_inf_p(tempy->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("'mpfr' invalid operation in modulo");
                goto error;
            }
            if (mpfr_signbit(tempy->f)) {
                mpfr_set_inf(result->f, -1);
            }
            else {
                result->rc = mpfr_set(result->f, tempx->f,
                                      GET_MPFR_ROUND(context));
            }
        }
        else {
            mpfr_div(temp->f, tempx->f, tempy->f, MPFR_RNDD);
            mpfr_floor(temp->f, temp->f);
            result->rc = mpfr_fms(result->f, temp->f, tempy->f, tempx->f,
                                  GET_MPFR_ROUND(context));
            mpfr_neg(result->f, result->f, GET_MPFR_ROUND(context));
        }
        Py_DECREF((PyObject*)temp);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        MPFR_CLEANUP_RESULT("mod");
        return (PyObject*)result;
    }

    Py_DECREF(result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)result);
    Py_DECREF((PyObject*)temp);
    return NULL;
}

static PyObject *
GMPy_mpfr_mod_fast(PyObject *x, PyObject *y)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Complex_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    TYPE_ERROR("can't take mod of complex number.");
    return NULL;
}

static PyObject *
GMPy_mpc_mod_fast(PyObject *x, PyObject *y)
{
    return GMPy_Complex_Mod(x, y, NULL);
}

PyDoc_STRVAR(GMPy_doc_mod,
"mod(x, y) -> number\n\n"
"Return mod(x, y).");

static PyObject *
GMPy_Number_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Mod(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mod(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mod(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mod(x, y, context);

    TYPE_ERROR("mod() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_mod,
"context.mod(x, y) -> number\n\n"
"Return mod(x, y).");

static PyObject *
GMPy_Context_Mod(PyObject *self, PyObject *args)
{
    PyObject *result;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mod() requires 2 arguments.");
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

    result = GMPy_Number_Mod(PyTuple_GET_ITEM(args, 0),
                                  PyTuple_GET_ITEM(args, 1),
                                  context);
    Py_DECREF((PyObject*)context);
    return result;
}

