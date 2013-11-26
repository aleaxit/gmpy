/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_abs.c                                                              *
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

/* This file implements abs(), gmpy2.abs(), and context.abs().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL
 * value for context implies the function should use the currently active
 * context.
 *
 *   GMPy_Number_Abs(Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_mpz_abs_fast; called by abs() via mpz.nb_absolute
 *   GMPy_mpq_abs_fast; called by abs() via mpq.nb_absolute
 *   GMPy_mpfr_abs_fast; called by abs() via mpfr.absolute
 *   GMPy_mpc_abs_fast; called by abs() via mpc.absolute
 *
 *   GMPy_Integer_Abs(Integer, context|NULL)
 *   GMPy_Rational_Abs(Rational, context|NULL)
 *   GMPy_Real_Abs(Real, context|NULL)
 *   GMPy_Complex_Abs(Complex, context|NULL)
 *
 *   GMPy_Context_Abs(context, args)
 */

/* Logical analysis
 *
 * If x already is an mpz, the sign is checked. If x is greater than or equal
 * to 0, the reference count is incremented and x is returned. If x is less
 * than 0, then a new object is created, set to the absolute value and the
 * new object is returned.
 *
 * Since we call GMPy_MPZ_From_Integer_New(), we could skip the check for mpz
 * argument but it is a reasonable optimization.
 */

static PyObject *
GMPy_Integer_Abs(PyObject *x, CTXT_Object *context)
{
    MPZ_Object *result;

    if (MPZ_Check(x)) {
        if (mpz_sgn(MPZ(x)) >= 0) {
            Py_INCREF(x);
            return x;
        }
        else {
            if ((result = GMPy_MPZ_New(context)))
                mpz_abs(result->z, MPZ(x));
            return (PyObject*)result;
        }
    }

    if ((result = GMPy_MPZ_From_Integer_New(x, context))) {
        mpz_abs(result->z, result->z);
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_mpz_abs_fast(MPZ_Object *x)
{
    return GMPy_Integer_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Rational_Abs(PyObject *x, CTXT_Object *context)
{
    MPQ_Object *result;

    if (MPQ_Check(x)) {
        if (mpz_sgn(mpq_numref(MPQ(x))) >= 0) {
            Py_INCREF(x);
            return x;
        }
        else {
            if ((result = GMPy_MPQ_New(context))) {
                mpq_set(result->q, MPQ(x));
                mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
            }
            return (PyObject*)result;
        }
    }

    if ((result = GMPy_MPQ_From_Number_New(x, context))) {
        mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_mpq_abs_fast(MPQ_Object *x)
{
    return GMPy_Rational_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Real_Abs(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPFR_From_Real_New(x, 0, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_abs(result->f, result->f, GET_MPFR_ROUND(context));
        MPFR_CLEANUP_2(result, context, "abs()");
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_mpfr_abs_fast(MPFR_Object *x)
{
    return GMPy_Real_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Complex_Abs(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;
    MPC_Object *tempx;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(tempx = GMPy_MPC_From_Complex_Temp(x, 0, 0, context)))  {
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpc_abs(result->f, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    MPFR_CLEANUP_2(result, context, "abs()");

    return (PyObject*)result;
}

static PyObject *
GMPy_mpc_abs_fast(MPC_Object *x)
{
    return GMPy_Complex_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Number_Abs(PyObject *x, CTXT_Object *context)
{
    if (IS_INTEGER(x))
        return GMPy_Integer_Abs(x, context);

    if (IS_RATIONAL_ONLY(x))
        return GMPy_Rational_Abs(x, context);

    if (IS_REAL_ONLY(x))
        return GMPy_Real_Abs(x, context);

    if (IS_COMPLEX_ONLY(x))
        return GMPy_Complex_Abs(x, context);

    TYPE_ERROR("abs() argument type not supported");
    return NULL;
}

/* Implement context.abs(). The following code assumes it used a as method of
 * a context. */

PyDoc_STRVAR(GMPy_doc_context_abs,
"context.abs(x) -> number\n\n"
"Return abs(x), the context is applied to the result.");

static PyObject *
GMPy_Context_Abs(PyObject *self, PyObject *args)
{
    assert(CTXT_Check(self));

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("context.abs() requires 1 argument.");
        return NULL;
    }

    return GMPy_Number_Abs(PyTuple_GET_ITEM(args, 0), (CTXT_Object*)self);
}

