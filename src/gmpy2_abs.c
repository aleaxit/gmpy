/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_abs.c                                                             *
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

/* This file implements __abs__, gmpy2.abs(), and context.abs().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. If the value
 * of context is NULL, then the function should use the currently active
 * context.
 *
 *   GMPy_Number_Abs(Number, context)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Abs_Slot
 *   GMPy_MPQ_Abs_Slot
 *   GMPy_MPFR_Abs_Slot
 *   GMPy_MPC_Abs_Slot
 *
 *   GMPy_Integer_Abs(Integer, context|NULL)
 *   GMPy_Rational_Abs(Rational, context|NULL)
 *   GMPy_Real_Abs(Real, context|NULL)
 *   GMPy_Complex_Abs(Complex, context|NULL)
 *
 *   GMPy_Context_Abs(context, obj)
 */

static PyObject *
GMPy_Integer_Abs(PyObject *x, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

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

    /* This is safe because result is not an incremented reference to an
     * existing value. Why?
     *   1) No values are interned like Python's integers.
     *   2) MPZ is already handled so GMPy_MPZ_From_Integer() can't return
     *      an incremented reference to an existing value (which it would do
     *      if passed an MPZ).
     */

    if ((result = GMPy_MPZ_From_Integer(x, context))) {
        mpz_abs(result->z, result->z);
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Abs_Slot(MPZ_Object *x)
{
    return GMPy_Integer_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Rational_Abs(PyObject *x, CTXT_Object *context)
{
    MPQ_Object *result = NULL;

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

    /* This is safe because result is not an incremented reference to an
     * existing value. MPQ is already handled so GMPy_MPQ_From_Rational()
     * can't return an incremented reference to an existing value (which it
     * would do if passed an MPQ).
     */

    if ((result = GMPy_MPQ_From_Rational(x, context))) {
        mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Abs_Slot(MPQ_Object *x)
{
    return GMPy_Rational_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Real_Abs(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)) ||
        !(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_abs(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPFR_Abs_Slot(MPFR_Object *x)
{
    return GMPy_Real_Abs((PyObject*)x, NULL);
}

static PyObject *
GMPy_Complex_Abs(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPC_Object *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPC_From_Complex(x, 1, 1, context)) ||
        !(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpc_abs(result->f, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;

}

static PyObject *
GMPy_MPC_Abs_Slot(MPC_Object *x)
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
GMPy_Context_Abs(PyObject *self, PyObject *other)
{
    return GMPy_Number_Abs(other, (CTXT_Object*)self);
}

