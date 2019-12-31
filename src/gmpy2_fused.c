/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_fused.c                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018, 2019 Case Van Horsen                  *
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

static PyObject *
_GMPy_MPZ_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_add(result->z, result->z, MPZ(z));
    return (PyObject*)result;
}

static PyObject *
GMPy_Integer_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_Integer(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_Integer(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_Integer(z, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMA(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPQ_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_add(result->q, result->q, MPQ(z));
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_Rational(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_Rational(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_Rational(z, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMA(tempx, tempy, tempz, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempz);
    return (result);
}

static PyObject *
_GMPy_MPFR_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fma(result->f, MPFR(x), MPFR(y), MPFR(z), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_Real(y, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_Real(z, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMA(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPC_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result->rc = mpc_fma(result->c, MPC(x), MPC(y), MPC(z), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPC_From_Complex(y, 1, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPC_From_Complex(z, 1, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPC_FMA(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fma,
"context.fma(x, y, z) -> number\n\n"
"Return correctly rounded result of (x * y) + z.");

PyDoc_STRVAR(GMPy_doc_function_fma,
"fma(x, y, z) -> number\n\n"
"Return correctly rounded result of (x * y) + z.");

GMPY_MPFR_MPC_TRIOP_TEMPLATE(FMA, fma);

static PyObject *
_GMPy_MPZ_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_sub(result->z, result->z, MPZ(z));
    return (PyObject*)result;
}

static PyObject *
GMPy_Integer_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_Integer(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_Integer(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_Integer(z, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPQ_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_sub(result->q, result->q, MPQ(z));
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_Rational(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_Rational(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_Rational(z, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPFR_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fms(result->f, MPFR(x), MPFR(y), MPFR(z), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_Real(y, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_Real(z, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

static PyObject *
_GMPy_MPC_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* TODO: We shouldn't temporarily mutate an mpc object. */

    mpc_neg(MPC(z), MPC(z), GET_MPC_ROUND(context));
    result->rc = mpc_fma(result->c, MPC(x), MPC(y), MPC(z), GET_MPC_ROUND(context));
    mpc_neg(MPC(z), MPC(z), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPC_From_Complex(y, 1, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPC_From_Complex(z, 1, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPC_FMS(tempx, tempy, tempz, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fms,
"context.fms(x, y, z) -> number\n\n"
"Return correctly rounded result of (x * y) - z.");

PyDoc_STRVAR(GMPy_doc_function_fms,
"fms(x, y, z) -> number\n\n"
"Return correctly rounded result of (x * y) - z.");

GMPY_MPFR_MPC_TRIOP_TEMPLATE(FMS, fms);

/* Add support for new fmma and fmms functions from MPFr 4. */\

#if MPFR_VERSION_MAJOR > 3

static PyObject *
_GMPy_MPZ_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPZ_New(context)) ||
        !(temp = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_mul(temp->z, MPZ(z), MPZ(t));
    mpz_add(result->z, result->z, temp->z);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_Integer_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_Integer(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_Integer(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_Integer(z, context)) ||
        !(tempt = (PyObject*)GMPy_MPZ_From_Integer(t, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMMA(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

static PyObject *
_GMPy_MPQ_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPQ_New(context)) ||
        !(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_mul(temp->q, MPQ(z), MPQ(t));
    mpq_add(result->q, result->q, temp->q);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_Rational(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_Rational(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_Rational(z, context)) ||
        !(tempt = (PyObject*)GMPy_MPQ_From_Rational(t, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMMA(tempx, tempy, tempz, tempt, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempz);
    Py_DECREF((PyObject*)tempt);
    return (result);
}

static PyObject *
_GMPy_MPFR_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fmma(result->f, MPFR(x), MPFR(y), MPFR(z), MPFR(t), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_FMMA(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_Real(y, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_Real(z, 1, context)) ||
        !(tempt = (PyObject*)GMPy_MPFR_From_Real(t, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMMA(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fmma,
"context.fmma(x, y, z, t) -> number\n\n"
"Return correctly rounded result of (x * y) + (z * t).");

PyDoc_STRVAR(GMPy_doc_function_fmma,
"fmma(x, y, z, t) -> number\n\n"
"Return correctly rounded result of (x * y) + (z + t).");

GMPY_MPFR_QUADOP_TEMPLATE(FMMA, fmma);

static PyObject *
_GMPy_MPZ_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPZ_New(context)) ||
        !(temp = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_mul(result->z, MPZ(x), MPZ(y));
    mpz_mul(temp->z, MPZ(z), MPZ(t));
    mpz_sub(result->z, result->z, temp->z);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_Integer_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPZ_From_Integer(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPZ_From_Integer(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPZ_From_Integer(z, context)) ||
        !(tempt = (PyObject*)GMPy_MPZ_From_Integer(t, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPZ_FMMS(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

static PyObject *
_GMPy_MPQ_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    if (!(result = GMPy_MPQ_New(context)) ||
        !(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF(temp);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_mul(result->q, MPQ(x), MPQ(y));
    mpq_mul(temp->q, MPQ(z), MPQ(t));
    mpq_sub(result->q, result->q, temp->q);
    Py_DECREF(temp);
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    if (!(tempx = (PyObject*)GMPy_MPQ_From_Rational(x, context)) ||
        !(tempy = (PyObject*)GMPy_MPQ_From_Rational(y, context)) ||
        !(tempz = (PyObject*)GMPy_MPQ_From_Rational(z, context)) ||
        !(tempt = (PyObject*)GMPy_MPQ_From_Rational(t, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPQ_FMMS(tempx, tempy, tempz, tempt, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempz);
    Py_DECREF((PyObject*)tempt);
    return (result);
}

static PyObject *
_GMPy_MPFR_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpfr_clear_flags();

    result->rc = mpfr_fmms(result->f, MPFR(x), MPFR(y), MPFR(z), MPFR(t), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_FMMS(PyObject *x, PyObject *y, PyObject *z, PyObject *t, CTXT_Object *context)
{
    PyObject *result, *tempx = NULL, *tempy = NULL, *tempz = NULL, *tempt = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context)) ||
        !(tempy = (PyObject*)GMPy_MPFR_From_Real(y, 1, context)) ||
        !(tempz = (PyObject*)GMPy_MPFR_From_Real(z, 1, context)) ||
        !(tempt = (PyObject*)GMPy_MPFR_From_Real(t, 1, context))) {
        /* LCOV_EXCL_START */
        Py_XDECREF(tempx);
        Py_XDECREF(tempy);
        Py_XDECREF(tempz);
        Py_XDECREF(tempt);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    result = _GMPy_MPFR_FMMS(tempx, tempy, tempz, tempt, context);
    Py_DECREF(tempx);
    Py_DECREF(tempy);
    Py_DECREF(tempz);
    Py_DECREF(tempt);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_fmms,
"context.fmms(x, y, z, t) -> number\n\n"
"Return correctly rounded result of (x * y) - (z * t).");

PyDoc_STRVAR(GMPy_doc_function_fmms,
"fmms(x, y, z, t) -> number\n\n"
"Return correctly rounded result of (x * y) - (z + t).");

GMPY_MPFR_QUADOP_TEMPLATE(FMMS, fmms);

#endif


