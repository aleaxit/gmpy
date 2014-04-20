/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_minus.c                                                           *
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

/* This file implements __neg__ and context.minus().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. If the value
 * of context is NULL, then the function should use the currently active
 * context.
 *
 *   GMPy_Number_Minus(Number, context)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Minus_Slot
 *   GMPy_MPQ_Minus_Slot
 *   GMPy_MPFR_Minus_Slot
 *   GMPy_MPC_Minus_Slot
 *
 *   GMPy_Integer_Minus(Integer, context|NULL)
 *   GMPy_Rational_Minus(Rational, context|NULL)
 *   GMPy_Real_Minus(Real, context|NULL)
 *   GMPy_Complex_Minus(Complex, context|NULL)
 *
 *   GMPy_Context_Minus(context, args)
 */

static PyObject *
GMPy_Integer_Minus(PyObject *x, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    result = GMPy_MPZ_New(context);
    tempx = GMPy_MPZ_From_Integer(x, context);
    if (!result || !tempx) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_neg(result->z, tempx->z);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Minus_Slot(MPZ_Object *x)
{
    return GMPy_Integer_Minus((PyObject*)x, NULL);
}

static PyObject *
GMPy_Rational_Minus(PyObject *x, CTXT_Object *context)
{
    MPQ_Object *result, *tempx;

    CHECK_CONTEXT(context);

    result = GMPy_MPQ_New(context);
    tempx = GMPy_MPQ_From_Rational(x, context);
    if (!result || !tempx) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpq_neg(result->q, tempx->q);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Minus_Slot(MPQ_Object *x)
{
    return GMPy_Rational_Minus((PyObject*)x, NULL);
}

static PyObject *
GMPy_Real_Minus(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPFR_From_Real(x, 1, context);
    result = GMPy_MPFR_New(0, context);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_neg(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPFR_CLEANUP(result, context, "minus()");
    return (PyObject*)result;
}

static PyObject *
GMPy_MPFR_Minus_Slot(MPFR_Object *x)
{
    return GMPy_Real_Minus((PyObject*)x, NULL);
}

static PyObject *
GMPy_Complex_Minus(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    result = GMPy_MPC_New(0, 0, context);
    if (!tempx || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpc_neg(result->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPC_CLEANUP(result, context, "minus()");
    return (PyObject*)result;
}

static PyObject *
GMPy_MPC_Minus_Slot(MPC_Object *x)
{
    return GMPy_Complex_Minus((PyObject*)x, NULL);
}

static PyObject *
GMPy_Number_Minus(PyObject *x, CTXT_Object *context)
{
    if (IS_INTEGER(x))
        return GMPy_Integer_Minus(x, context);

    if (IS_RATIONAL_ONLY(x))
        return GMPy_Rational_Minus(x, context);

    if (IS_REAL_ONLY(x))
        return GMPy_Real_Minus(x, context);

    if (IS_COMPLEX_ONLY(x))
        return GMPy_Complex_Minus(x, context);

    TYPE_ERROR("minus() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_context_minus,
"context.minus(x) -> number\n\n"
"Return -x. The context is applied to the result.");

static PyObject *
GMPy_Context_Minus(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("minus() requires 1 argument.");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Minus(PyTuple_GET_ITEM(args, 0), context);
}

