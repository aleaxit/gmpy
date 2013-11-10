/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_pow.c                                                              *
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

/* This file implements the ** operator, Python's pow() function,
 * gmpy2.powmod(), and context.pow().
 *
 * Private API
 * ===========
 * The Python ** operator and the pow() function both call the nb_add slot
 * of a numeric type. This file implements the following private function:
 *
 *   GMPy_mpany_pow_fast; called via the nb_power slot of mp*
 *   GMPy_Integer_Pow(Integer, Integer, Integer|Py_None, context|NULL)
 *   GMPy_Rational_Pow(Rational, Rational, context|NULL)
 *   GMPy_Real_Pow(Real, Real, context|NULL)
 *   GMPy_Complex_Pow(Complex, Complex, context|NULL)
 *   GMPy_Context_Pow; called by gmpy2.pow() and context.pow()
 *
 * Public API
 * ==========
 * The following functions are availabe as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 * The first four functions check the type of the first argument and will set
 * an exception and return NULL if the check fails.
 *
 *   Not yet implemented!
 *
 *   GMPy_Number_Pow(Number, Number, context|NULL)
 *
 */


/* Pympz_Pow_Integer is called by GMPy_Number_Pow() after verifying that the
 * first two arguments are integers, but not necessarily mpz. The third
 * argument must either be an integer or Py_None. The context argument is not
 * currently used but may be used in the future.
 */

static PyObject *
GMPy_Integer_Pow(PyObject *b, PyObject *e, PyObject *m, GMPyContextObject *context)
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
            if (!(tempm = GMPy_MPZ_From_Integer_Temp(m))) {
                goto err;
            }
        }
    }

    result = GMPy_MPZ_New();
    tempb = GMPy_MPZ_From_Integer_Temp(b);
    tempe = GMPy_MPZ_From_Integer_Temp(e);

    if (!tempb || !tempe || !result) {
        goto err;
    }

    if (!has_mod) {
        /* When no modulo is present, the exponent must fit in mpir_ui and
         * the exponent must be positive.
         */
        mpir_ui el;

        if (mpz_sgn(tempe->z) < 0) {
            VALUE_ERROR("pow() exponent cannot be negative");
            goto err;
        }

        if (!mpz_fits_ui_p(tempe->z)) {
            VALUE_ERROR("pow() outrageous exponent");
            goto err;
        }

        el = mpz_get_ui(tempe->z);
        mpz_pow_ui(result->z, tempb->z, el);
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

        mpz_inoc(mm);
        mpz_abs(mm, tempm->z);

        /* A negative exponent is allowed if inverse exists. */
        if (mpz_sgn(tempe->z) < 0) {
            mpz_inoc(base);
            mpz_inoc(exp);

            if (!mpz_invert(base, tempb->z, mm)) {
                VALUE_ERROR("pow() base not invertible");
                mpz_cloc(base);
                mpz_cloc(exp);
                mpz_cloc(mm);
                goto err;
            }
            else {
                mpz_abs(exp, tempe->z);
            }

            mpz_powm(result->z, base, exp, mm);
            mpz_cloc(base);
            mpz_cloc(exp);
        }
        else {
            mpz_powm(result->z, tempb->z, tempe->z, mm);
        }
        mpz_cloc(mm);

        /* Python uses a rather peculiar convention for negative modulos
         * If the modulo is negative, result should be in the interval
         * m < r <= 0 .
         */
        if ((sign < 0) && (mpz_sgn(MPZ(result)) > 0)) {
            mpz_add(result->z, result->z, tempm->z);
        }
    }
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
GMPy_Rational_Pow(PyObject *base, PyObject *exp, GMPyContextObject *context)
{
    MPQ_Object *tempbq = NULL, *resultq = NULL;
    MPZ_Object *tempez = NULL;
    int esign, bsign;
    mpir_si tempexp;

    /* Only support mpq**int. Everything else gets converted to mpf. */
    if (IS_RATIONAL(base) && IS_INTEGER(exp)) {

        resultq = GMPy_MPQ_New();
        tempbq = Pympq_From_Rational(base);
        tempez = GMPy_MPZ_From_Integer_Temp(exp);
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

        esign = mpz_sgn(tempez->z);
        if (esign == 0) {
            mpq_set_si(resultq->q, 1, 1);
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return (PyObject*)resultq;
        }

        bsign = mpq_sgn(tempbq->q);
        if (esign < 0) {
            if (bsign == 0) {
                ZERO_ERROR("mpq.pow() 0 base to negative exponent");
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
            tempexp = -mpz_get_si(tempez->z);
        }
        else {
            mpq_set(resultq->q, tempbq->q);
            tempexp = mpz_get_si(tempez->z);
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
        return GMPy_Real_Pow(base, exp, context);
    }
}

static PyObject *
GMPy_Real_Pow(PyObject *base, PyObject *exp, GMPyContextObject *context)
{
    MPFR_Object *tempb = NULL, *tempe = NULL, *result = NULL;
    MPC_Object *mpc_result = NULL;

    if (!context)
        CURRENT_CONTEXT(context);

    SET_EXPONENT(context);

    result = (MPFR_Object*)Pympfr_new_context(context);
    tempb = GMPy_MPFR_From_Real_Temp(base, context);
    tempe = GMPy_MPFR_From_Real_Temp(exp, context);

    if (!result || !tempe || !tempb) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempb);
        return NULL;
    }

    if (mpfr_zero_p(tempb->f) && (mpfr_sgn(tempe->f) < 0)) {
        context->ctx.divzero = 1;
        if (context->ctx.traps & TRAP_DIVZERO) {
            GMPY_DIVZERO("zero cannot be raised to a negative power");
            goto err;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_pow(result->f, tempb->f, tempe->f,
                          GET_MPFR_ROUND(context));

    if (result && mpfr_nanflag_p() && context->ctx.allow_complex) {
        /* If the result is NaN, check if a complex result works. */

        mpc_result = (MPC_Object*)GMPy_Complex_Pow(base, exp, context);
        if (!mpc_result || MPC_IS_NAN_P(mpc_result)) {
            Py_XDECREF((PyObject*)mpc_result);
            context->ctx.invalid = 1;
            GMPY_INVALID("invalid operation in 'mpfr' pow()");
            goto err;
        }
        /* return a valid complex result */
        Py_DECREF((PyObject*)tempe);
        Py_DECREF((PyObject*)tempb);
        Py_DECREF((PyObject*)result);
        return (PyObject*)mpc_result;
    }

    MPFR_CLEANUP_2(result, context, "pow()");
    Py_DECREF((PyObject*)tempe);
    Py_DECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)result);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempb);
    return NULL;
}

static PyObject *
GMPy_Complex_Pow(PyObject *base, PyObject *exp, GMPyContextObject *context)
{
    MPC_Object *tempb = NULL, *tempe = NULL, *result= NULL;

    result = (MPC_Object*)Pympc_new_context(context);
    tempb = GMPy_MPC_From_Complex_Temp(base, context);
    tempe = GMPy_MPC_From_Complex_Temp(exp, context);

    if (!result || !tempe || !tempb) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempb);
        return NULL;
    }

    if (MPC_IS_ZERO_P(tempb) && MPC_IS_ZERO_P(tempe)) {
        mpc_set_ui(result->c, 1, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempe);
        Py_DECREF((PyObject*)tempb);
        return (PyObject*)result;
    }

    if (MPC_IS_ZERO_P(tempb) &&
        (!mpfr_zero_p(mpc_imagref(tempe->c)) ||
         mpfr_sgn(mpc_realref(tempe->c)) < 0)) {

        context->ctx.divzero = 1;
        if (context->ctx.traps & TRAP_DIVZERO) {
            GMPY_DIVZERO("zero cannot be raised to a negative or complex power");
            Py_DECREF((PyObject*)tempe);
            Py_DECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }

    result->rc = mpc_pow(result->c, tempb->c,
                         tempe->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempe);
    Py_DECREF((PyObject*)tempb);

    MPC_CLEANUP_2(result, context, "pow()");
    return (PyObject*)result;
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

PyDoc_STRVAR(GMPy_doc_context_pow,
"context.pow(x, y) -> number\n\n"
"Return x ** y.");

static PyObject *
GMPy_Context_Pow(PyObject *self, PyObject *args)
{
    Py_ssize_t argc;
    PyObject *arg0, *arg1;
    GMPyContextObject *context;

    argc = PyTuple_GET_SIZE(args);
    if (self && GMPyContext_Check(self)) {
        if (argc != 2) {
            TYPE_ERROR("context.pow() requires 2 arguments.");
            return NULL;
        }
        /* If we are passed a read-only context, make a copy of it before
         * proceeding. */

        if (((GMPyContextObject*)self)->ctx.readonly)
            context = (GMPyContextObject*)GMPyContext_context_copy(self, NULL);
        else
            context = (GMPyContextObject*)self;
    }
    else {
        /* pow() is only supported as a context method so this branch should
         * be impossible to reach.
         */
        SYSTEM_ERROR("pow() only supported as method of a context.");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);
    arg1 = PyTuple_GET_ITEM(args, 1);

    if (IS_INTEGER(arg0) && IS_INTEGER(arg1))
        return GMPy_Integer_Pow(arg0, arg1, Py_None, context);

    if (IS_RATIONAL(arg0) && IS_RATIONAL(arg1))
        return GMPy_Rational_Pow(arg0, arg1, context);

    if (IS_REAL(arg0) && IS_REAL(arg1))
        return GMPy_Real_Pow(arg0, arg1, context);

    if (IS_COMPLEX(arg0) && IS_COMPLEX(arg1))
        return GMPy_Complex_Pow(arg0, arg1, context);

    TYPE_ERROR("pow() argument types not supported");
    return NULL;
}

static PyObject *
GMPy_mpany_pow_fast(PyObject *base, PyObject *exp, PyObject *mod)
{
    if (IS_INTEGER(base) && IS_INTEGER(exp))
        return GMPy_Integer_Pow(base, exp, mod, NULL);

    if (IS_RATIONAL(base) && IS_RATIONAL(exp))
        return GMPy_Rational_Pow(base, exp, NULL);

    if (IS_REAL(base) && IS_REAL(exp))
        return GMPy_Real_Pow(base, exp, NULL);

    if (IS_COMPLEX(base) && IS_COMPLEX(exp))
        return GMPy_Complex_Pow(base, exp, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}


