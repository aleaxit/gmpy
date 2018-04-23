/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_pow.c                                                             *
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

/* This file implements the ** operator, Python's pow() function,
 * gmpy2.powmod(), and context.pow().
 *
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_Pow(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_MPANY_Pow_Slot
 *
 *   GMPy_Integer_Pow(Integer, Integer, Integer|Py_None, context|NULL)
 *   GMPy_Integer_PowMod(Integer, Integer, Integer|Py_None, context|NULL)
 *   GMPy_Rational_Pow(Rational, Rational, context|NULL)
 *   GMPy_Real_Pow(Real, Real, context|NULL)
 *   GMPy_Complex_Pow(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Pow(context, args)
 */


/* Pympz_Pow_Integer is called by GMPy_Number_Pow() after verifying that the
 * first two arguments are integers, but not necessarily mpz. The third
 * argument must either be an integer or Py_None. The context argument is not
 * currently used but may be used in the future.
 */

/* TODO: refactor to improve performance.
 *   - dont't convert the exponent to an MPZ if there is no modulus
 */

static PyObject *
GMPy_Integer_Pow(PyObject *b, PyObject *e, PyObject *m, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *tempb = NULL, *tempe = NULL, *tempm = NULL;
    int has_mod;

    /* Try to parse the modulus value first. */

    if (m == Py_None) {
        has_mod = 0;
    }
    else {
        has_mod = 1;
        if (!IS_INTEGER(m)) {
            Py_RETURN_NOTIMPLEMENTED;
        }
        else {
            if (!(tempm = GMPy_MPZ_From_Integer(m, context))) {
                goto err;
            }
        }
    }

    result = GMPy_MPZ_New(context);
    tempb = GMPy_MPZ_From_Integer(b, context);
    tempe = GMPy_MPZ_From_Integer(e, context);

    if (!tempb || !tempe || !result) {
        goto err;
    }

    if (!has_mod) {
        /* When no modulo is present, the exponent must fit in unsigned long.
         */
        unsigned long el;

        if (mpz_sgn(tempe->z) < 0) {
            VALUE_ERROR("pow() exponent cannot be negative");
            goto err;
        }

        /* Catch -1, 0, 1 getting raised to large exponents. */

        if (mpz_cmp_si(tempb->z, 0) == 0) {
            if (mpz_cmp_si(tempe->z, 0) == 0) {
                mpz_set_ui(result->z, 1);
            }
            else {
                mpz_set_ui(result->z, 0);
            }
            goto done;
        }

        if (mpz_cmp_si(tempb->z, 1) == 0) {
            mpz_set_ui(result->z, 1);
            goto done;
        }

        if (mpz_cmp_si(tempb->z, -1) == 0) {
            if (mpz_even_p(tempe->z)) {
                mpz_set_ui(result->z, 1);
            }
            else {
                mpz_set_si(result->z, -1);
            }
            goto done;
        }

        if (!mpz_fits_ulong_p(tempe->z)) {
            VALUE_ERROR("pow() outrageous exponent");
            goto err;
        }

        el = (unsigned long) mpz_get_ui(tempe->z);
        mpz_pow_ui(result->z, tempb->z, el);
        goto done;
    }
    else {
        /* Modulo is present. */
        int sign;
        mpz_t mm, base, exp;

        sign = mpz_sgn(tempm->z);
        if (sign == 0) {
            VALUE_ERROR("pow() 3rd argument cannot be 0");
            goto err;
        }

        mpz_init(mm);
        mpz_abs(mm, tempm->z);

        /* A negative exponent is allowed if inverse exists. */
        if (mpz_sgn(tempe->z) < 0) {
            mpz_init(base);
            mpz_init(exp);

            if (!mpz_invert(base, tempb->z, mm)) {
                VALUE_ERROR("pow() base not invertible");
                mpz_clear(base);
                mpz_clear(exp);
                mpz_clear(mm);
                goto err;
            }
            else {
                mpz_abs(exp, tempe->z);
            }

            mpz_powm(result->z, base, exp, mm);
            mpz_clear(base);
            mpz_clear(exp);
        }
        else {
            mpz_powm(result->z, tempb->z, tempe->z, mm);
        }
        mpz_clear(mm);

        /* Python uses a rather peculiar convention for negative modulos
         * If the modulo is negative, result should be in the interval
         * m < r <= 0 .
         */
        if ((sign < 0) && (mpz_sgn(MPZ(result)) > 0)) {
            mpz_add(result->z, result->z, tempm->z);
        }
    }

  done:
    Py_XDECREF((PyObject*)tempb);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempm);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)tempb);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempm);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PyObject *
GMPy_Rational_Pow(PyObject *base, PyObject *exp, PyObject *mod, CTXT_Object *context)
{
    MPQ_Object *tempbq = NULL, *resultq = NULL;
    MPZ_Object *tempez = NULL;
    int bsign;
    long tempexp;

    if (mod != Py_None) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    /* Only support mpq**int. Everything else gets converted to mpf. */
    if (IS_RATIONAL(base) && IS_INTEGER(exp)) {

        resultq = GMPy_MPQ_New(context);
        tempbq = GMPy_MPQ_From_Rational(base, context);
        tempez = GMPy_MPZ_From_Integer(exp, context);
        if (!resultq || !tempbq || !tempez) {
            Py_XDECREF((PyObject*)resultq);
            Py_XDECREF((PyObject*)tempbq);
            Py_XDECREF((PyObject*)tempez);
            return NULL;
        }

        if (!mpz_fits_slong_p(tempez->z)) {
            VALUE_ERROR("mpq.pow() outrageous exponent");
            Py_DECREF((PyObject*)resultq);
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return NULL;
        }
        tempexp = (long) mpz_get_si(tempez->z);

        if (tempexp == 0) {
            mpq_set_si(resultq->q, 1, 1);
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return (PyObject*)resultq;
        }

        bsign = mpq_sgn(tempbq->q);
        if (tempexp < 0) {
            if (bsign == 0) {
                ZERO_ERROR("pow() 0 base to negative exponent");
                Py_DECREF((PyObject*)resultq);
                Py_DECREF((PyObject*)tempbq);
                Py_DECREF((PyObject*)tempez);
                return NULL;
            }
            if (bsign < 0) {
                mpz_neg(mpq_numref(resultq->q), mpq_denref(tempbq->q));
            }
            else {
                mpz_set(mpq_numref(resultq->q), mpq_denref(tempbq->q));
            }
            mpz_abs(mpq_denref(resultq->q), mpq_numref(tempbq->q));
            tempexp = -tempexp;
        }
        else {
            mpq_set(resultq->q, tempbq->q);
        }

        if (tempexp > 1) {
            mpz_pow_ui(mpq_numref(resultq->q), mpq_numref(resultq->q), tempexp);
            mpz_pow_ui(mpq_denref(resultq->q), mpq_denref(resultq->q), tempexp);
        }
        Py_DECREF((PyObject*)tempbq);
        Py_DECREF((PyObject*)tempez);
        return (PyObject*)resultq;
    }
    else {
        return GMPy_Real_Pow(base, exp, Py_None, context);
    }
}

static PyObject *
GMPy_Real_Pow(PyObject *base, PyObject *exp, PyObject *mod, CTXT_Object *context)
{
    MPFR_Object *tempb = NULL, *tempe = NULL, *result = NULL;
    MPZ_Object *tempz = NULL;
    MPC_Object *mpc_result = NULL;

    if (mod != Py_None) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempb = GMPy_MPFR_From_Real(base, 1, context);
    if (!result || !tempb) {
        goto err;
    }

    mpfr_clear_flags();

    if (PyIntOrLong_Check(exp)) {
        int error;
        long temp_exp = GMPy_Integer_AsLongAndError(exp, &error);

        if (!error) {
            result->rc = mpfr_pow_si(result->f, tempb->f, temp_exp, GET_MPFR_ROUND(context));
        }
        else {
            mpz_t tempzz;
            mpz_init(tempzz);
            mpz_set_PyIntOrLong(tempzz, exp);
            result->rc = mpfr_pow_z(result->f, tempb->f, tempzz, GET_MPFR_ROUND(context));
            mpz_clear(tempzz);
        }
    }
    else if (IS_INTEGER(exp)) {
        if (!(tempz = GMPy_MPZ_From_Integer(exp, context))) {
            goto err;
        }
        result->rc = mpfr_pow_z(result->f, tempb->f, tempz->z, GET_MPFR_ROUND(context));
    }
    else {
        if (!(tempe = GMPy_MPFR_From_Real(exp, 1, context))) {
            goto err;
        }
        result->rc = mpfr_pow(result->f, tempb->f, tempe->f, GET_MPFR_ROUND(context));
    }

    /* If the result is NaN, check if a complex result works. */
    if (mpfr_nanflag_p() && context->ctx.allow_complex) {
        mpc_result = (MPC_Object*)GMPy_Complex_Pow(base, exp, Py_None, context);
        if (!mpc_result || MPC_IS_NAN_P(mpc_result)) {
            Py_XDECREF((PyObject*)mpc_result);
            context->ctx.invalid = 1;
            GMPY_INVALID("pow() invalid operation");
            goto err;
        }
        /* return a valid complex result */
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempz);
        Py_XDECREF((PyObject*)tempb);
        Py_XDECREF((PyObject*)result);
        return (PyObject*)mpc_result;
    }

    _GMPy_MPFR_Cleanup(&result, context);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return NULL;
}

static PyObject *
GMPy_Complex_Pow(PyObject *base, PyObject *exp, PyObject *mod, CTXT_Object *context)
{
    MPC_Object *tempb = NULL, *tempe = NULL, *result= NULL;
    MPFR_Object *tempf = NULL;
    MPZ_Object *tempz = NULL;

    if (mod != Py_None) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    CHECK_CONTEXT(context);

    result = GMPy_MPC_New(0, 0, context);
    tempb = GMPy_MPC_From_Complex(base, 1, 1, context);
    if (!result || !tempb) {
        goto err;
    }

    mpfr_clear_flags();

    if (PyIntOrLong_Check(exp)) {
        int error;
        long temp_exp = GMPy_Integer_AsLongAndError(exp, &error);

        if (!error) {
            result->rc = mpc_pow_si(result->c, tempb->c, temp_exp, GET_MPC_ROUND(context));
        }
        else {
            mpz_t tempzz;
            mpz_init(tempzz);
            mpz_set_PyIntOrLong(tempzz, exp);
            result->rc = mpc_pow_z(result->c, tempb->c, tempzz, GET_MPC_ROUND(context));
            mpz_clear(tempzz);
        }
    }
    else if (IS_INTEGER(exp)) {
        if (!(tempz = GMPy_MPZ_From_Integer(exp, context))) {
            goto err;
        }
        result->rc = mpc_pow_z(result->c, tempb->c, tempz->z, GET_MPC_ROUND(context));
    }
    else if (IS_REAL(exp)) {
        if (!(tempf = GMPy_MPFR_From_Real(exp, 1, context))) {
            goto err;
        }

        result->rc = mpc_pow_fr(result->c, tempb->c, tempf->f, GET_MPC_ROUND(context));
    }
    else {
        if (!(tempe = GMPy_MPC_From_Complex(exp, 1, 1, context))) {
            goto err;
        }

        result->rc = mpc_pow(result->c, tempb->c, tempe->c, GET_MPC_ROUND(context));
    }

    _GMPy_MPC_Cleanup(&result, context);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempf);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempz);
    Py_XDECREF((PyObject*)tempf);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_integer_powmod,
"powmod(x, y, m) -> mpz\n\n"
"Return (x**y) mod m. Same as the three argument version of Python's\n"
"built-in pow(), but converts all three arguments to mpz.");

static PyObject *
GMPy_Integer_PowMod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *m;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("powmod() requires 3 arguments.");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    m = PyTuple_GET_ITEM(args, 2);

    if (IS_INTEGER(x) && IS_INTEGER(y) && IS_INTEGER(m))
        return GMPy_Integer_Pow(x, y, m, NULL);

    TYPE_ERROR("powmod() argument types not supported");
    return NULL;
}

static PyObject *
GMPy_Number_Pow(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Pow(x, y, z, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Pow(x, y, z, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Pow(x, y, z, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Pow(x, y, z, context);

    TYPE_ERROR("pow() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_pow,
"context.pow(x, y) -> number\n\n"
"Return x ** y.");

static PyObject *
GMPy_Context_Pow(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("pow() requires 2 arguments.");
        return NULL;
    }
    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Pow(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1),
                           Py_None, context);
}

static PyObject *
GMPy_MPANY_Pow_Slot(PyObject *base, PyObject *exp, PyObject *mod)
{
    if (IS_INTEGER(base) && IS_INTEGER(exp))
        return GMPy_Integer_Pow(base, exp, mod, NULL);

    if (IS_RATIONAL(base) && IS_RATIONAL(exp))
        return GMPy_Rational_Pow(base, exp, mod, NULL);

    if (IS_REAL(base) && IS_REAL(exp))
        return GMPy_Real_Pow(base, exp, mod, NULL);

    if (IS_COMPLEX(base) && IS_COMPLEX(exp))
        return GMPy_Complex_Pow(base, exp, mod, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPFR_Pow_Slot(PyObject *base, PyObject *exp, PyObject *mod)
{
    if (MPFR_Check(base) && (PyIntOrLong_Check(exp)) && (mod == Py_None)) {
        MPFR_Object *result = NULL;
        int error = 0;
        long temp_exp;
        CTXT_Object *context = NULL;

        CHECK_CONTEXT(context);

        if ((result = GMPy_MPFR_New(0, context))) {
            mpfr_clear_flags();

            temp_exp = GMPy_Integer_AsLongAndError(exp, &error);

            if (!error) {
                result->rc = mpfr_pow_si(result->f, MPFR(base), temp_exp, GET_MPFR_ROUND(context));
            }
            else {
                mpz_t tempzz;
                mpz_init(tempzz);
                mpz_set_PyIntOrLong(tempzz, exp);
                result->rc = mpfr_pow_z(result->f, MPFR(base), tempzz, GET_MPFR_ROUND(context));
                mpz_clear(tempzz);
            }
            _GMPy_MPFR_Cleanup(&result, context);
        }
        return (PyObject*)result;
    }

    if (IS_REAL(base) && IS_REAL(exp))
        return GMPy_Real_Pow(base, exp, mod, NULL);

    if (IS_COMPLEX(base) && IS_COMPLEX(exp))
        return GMPy_Complex_Pow(base, exp, mod, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}


