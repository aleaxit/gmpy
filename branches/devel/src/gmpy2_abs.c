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

/* This file implements abs() and context.abs().
 *
 * Private API
 * ===========
 * The Python function abs() calls the __abs__ method of a numeric type. This
 * file implements the following private functions:
 *
 *   GMPy_mpz_abs_fast; called by abs() via mpz.__abs__
 *   GMPy_mpq_abs_fast; called by abs() via mpq.__abs__
 *   GMPy_mpfr_abs_fast; called by abs() via mpfr.__abs__
 *   GMPy_mpc_abs_fast; called by abs() via mpc.__abs__
 *
 *   GMPY_Context_Abs; called by context.abs()
 *
 * Public API
 * ==========
 * The following functions are availabe as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 * The first four functions check the type of the first argument and will set
 * an exception and return NULL if the check fails.
 *
 *   GMPy_Abs_Integer(Integer, context|NULL)
 *   GMPy_Abs_Rational(Rational, context|NULL)
 *   GMPy_Abs_Real(Real, context|NULL)
 *   GMPy_Abs_Complex(Complex, context|NULL)
 *   GMPy_Abs_Number(Number, context|NULL)
 *
 */

static PyObject *
GMPy_mpz_abs_fast(PympzObject *self)
{
    PympzObject *result;

    if ((result = (PympzObject*)Pympz_new()))
        mpz_abs(result->z, self->z);

    return (PyObject*)result;
}

static PyObject *
GMPy_Abs_Integer(PyObject *self)
{
    PympzObject *result;

    if (IS_INTEGER(self)) {
        if (result = Pympz_From_Integer(self)) {
            mpz_abs(result->z, self->z);
        }
    }
    return (PyObject*)result;
}
static PyObject *
Pympq_abs(PympqObject *self)
{
    PympqObject *result;

    if ((result = (PympqObject*)Pympq_new())) {
        mpq_set(result->q, self->q);
        mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
    }

    return (PyObject*)result;
}


/* Implement context.abs(). The following code assumes it used a as method of
 * a context. */

PyDoc_STRVAR(doc_context_abs,
"context.abs(x) -> number\n\n"
"Return abs(x), the context is applied to the result.");

static PyObject *
Pympany_abs(PyObject *self, PyObject *args)
{
    Py_ssize_t argc;
    PyObject *arg0;
    GMPyContextObject *context;

    argc = PyTuple_GET_SIZE(args);
    if (self && GMPyContext_Check(self)) {
        if (argc != 1) {
            TYPE_ERROR("context.abs() requires 1 argument.");
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
        SYSTEM_ERROR("function is not supported");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    if (IS_INTEGER(arg0))
        return GMPy_Abs_Integer(arg0);

    if (IS_RATIONAL(arg0))
        return Pympq_abs(arg0);

    if (IS_REAL(arg0))
        return Pympfr_Abs_Real(arg0, context);

    if (IS_COMPLEX(arg0))
        return Pympc_Abs_Complex(arg0, context);

    TYPE_ERROR("abs() argument types not supported");
    return NULL;
}

/* Pympfr_Abs_Real is expected to be called by Pympany_abs (when used as a
 * context method) or Pympfr_abs_fast (when used as mpfr.__abs__). */

static PyObject *
Pympfr_Abs_Real(PyObject *x, GMPyContextObject *context)
{
    PympfrObject *result;

    if (!(result = (PympfrObject*)Pympfr_new(0)))
        return NULL;

    if (Pympfr_CheckAndExp(x)) {
        mpfr_clear_flags();
        result->rc = mpfr_abs(result->f, Pympfr_AS_MPFR(x), GET_MPFR_ROUND(context));
        MERGE_FLAGS;
        CHECK_FLAGS("abs()");
        goto done;
    }
    else if (IS_REAL(x)) {
        PympfrObject *tempx;

        tempx = Pympfr_From_Real_context(x, 0, context);
        if (!tempx) {
            SYSTEM_ERROR("Can not covert Real to 'mpfr'");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpfr_clear_flags();
        result->rc = mpfr_abs(result->f, tempx->f, GET_MPFR_ROUND(context));
        MERGE_FLAGS;
        CHECK_FLAGS("abs()");
        Py_DECREF((PyObject*)tempx);
        goto done;
    }
    else {
        TYPE_ERROR("abs() called with invalid type");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

  done:
    MPFR_CLEANUP_RESULT("abs()");
    return (PyObject*)result;
}

static PyObject *
Pympfr_abs_fast(PympfrObject *x)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);
    return Pympfr_Abs_Real((PyObject*)x, context);
}

static PyObject *
Pympc_abs(PyObject *self)
{
    PympfrObject *result = 0;
    PympcObject *tempx = 0;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    result = (PympfrObject*)Pympfr_new(0);
    tempx = Pympc_From_Complex(self, 0, 0);
    if (!tempx || !result) {
        SYSTEM_ERROR("Can't convert argument to 'mpc'.");
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    result->rc = mpc_abs(result->f, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    MPFR_SUBNORMALIZE(result);
    MPFR_CHECK_INVALID(result, "invalid operation in 'mpc' __abs__");
    MPFR_CHECK_UNDERFLOW(result, "underflow in 'mpc' __abs__");
    MPFR_CHECK_OVERFLOW(result, "overflow in 'mpc' __abs__");
    MPFR_CHECK_INEXACT(result, "inexact result in 'mpc' __abs__");
  done:
    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}
