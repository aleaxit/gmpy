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
 *   GMPy_Context_Abs; called by context.abs()
 *
 * Public API
 * ==========
 * The following functions are availabe as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 * The first four functions check the type of the first argument and will set
 * an exception and return NULL if the check fails.
 *
 *   GMPy_Integer_Abs(Integer, context|NULL)
 *   GMPy_Rational_Abs(Rational, context|NULL)
 *   GMPy_Real_Abs(Real, context|NULL)
 *   GMPy_Complex_Abs(Complex, context|NULL)
 *   GMPy_Number_Abs(Number, context|NULL)
 *
 */

static PyObject *
GMPy_Integer_Abs(PyObject *x, GMPyContextObject *context)
{
    PympzObject *result;

    if (IS_INTEGER(x)) {
        if ((result = Pympz_From_Integer(x))) {
            mpz_abs(result->z, result->z);
        }
        return (PyObject*)result;
    }
    else {
        TYPE_ERROR("argument is not an integer number");
        return NULL;
    }
}

static PyObject *
GMPy_mpz_abs_fast(PympzObject *x)
{
    PympzObject *result;

    if ((result = (PympzObject*)Pympz_new())) {
        mpz_abs(result->z, x->z);
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_Abs(PyObject *x, GMPyContextObject *context)
{
    PympqObject *result;

    if (IS_RATIONAL(x)) {
        if ((result = Pympq_From_Rational(x))) {
            mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
        }
        return (PyObject*)result;
    }
    else {
        TYPE_ERROR("argument is not a rational number");
        return NULL;
    }
}

static PyObject *
GMPy_mpq_abs_fast(PympqObject *x)
{
    PympqObject *result;

    if ((result = (PympqObject*)Pympq_new())) {
        mpq_set(result->q, x->q);
        mpz_abs(mpq_numref(result->q), mpq_numref(result->q));
    }

    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Abs(PyObject *x, GMPyContextObject *context)
{
    PympfrObject *result;

    if (!context)
        CURRENT_CONTEXT(context);

    if (!(result = (PympfrObject*)Pympfr_new_context(context)))
        return NULL;

    if (Pympfr_CheckAndExp(x)) {
        SET_EXPONENT(context);
        mpfr_clear_flags();
        result->rc = mpfr_abs(result->f, Pympfr_AS_MPFR(x), GET_MPFR_ROUND(context));
        MERGE_FLAGS;
        CHECK_FLAGS("abs()");
        goto done;
    }
    else if (IS_REAL(x)) {
        PympfrObject *tempx;

        if (!(tempx = Pympfr_From_Real_context(x, 0, context))) {
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        SET_EXPONENT(context);
        mpfr_clear_flags();
        result->rc = mpfr_abs(result->f, tempx->f, GET_MPFR_ROUND(context));

        /* Per MPFR source, the only flag that can be set is overflow. The
         * following code should be cleaned up in the future.
         */

        Py_DECREF((PyObject*)tempx);
        MERGE_FLAGS;
        CHECK_FLAGS("abs()");
        goto done;
    }
    else {
        TYPE_ERROR("argument is not a real number");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

  done:
    MPFR_CLEANUP_RESULT("abs()");
    return (PyObject*)result;
}

static PyObject *
GMPy_mpfr_abs_fast(PympfrObject *x)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);
    return GMPy_Real_Abs((PyObject*)x, context);
}

static PyObject *
GMPy_Complex_Abs(PyObject *x, GMPyContextObject *context)
{
    PympfrObject *result;

    if (!context)
        CURRENT_CONTEXT(context);

    if (!(result = (PympfrObject*)Pympfr_new_context(context)))
        return NULL;

    if (IS_COMPLEX(x)) {
        PympcObject *tempx;

        if (!(tempx = Pympc_From_Complex_context(x, context)))  {
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        /* Per MPC, mpc_abs just calls mpfr_hypot. We will just call mpfr_hypot
         * directly.
         */
        SET_EXPONENT(context);
        mpfr_clear_flags();
        result->rc = mpfr_hypot(result->f, mpc_realref(tempx->c),
                                mpc_imagref(tempx->c), GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        MERGE_FLAGS;
        CHECK_FLAGS("abs()");
        goto done;
    }
    else {
        TYPE_ERROR("argument is not a complex number");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

  done:
    MPFR_CLEANUP_RESULT("abs()");
    return (PyObject*)result;
}

static PyObject *
GMPy_mpc_abs_fast(PympcObject *x)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);
    return GMPy_Complex_Abs((PyObject*)x, context);
}

static PyObject *
GMPy_Number_Abs(PyObject *x, GMPyContextObject *context)
{
    if (!context)
        CURRENT_CONTEXT(context);

    if (IS_INTEGER(x))
        return GMPy_Integer_Abs(x, context);

    if (IS_RATIONAL(x))
        return GMPy_Rational_Abs(x, context);

    if (IS_REAL(x))
        return GMPy_Real_Abs(x, context);

    if (IS_COMPLEX(x))
        return GMPy_Complex_Abs(x, context);

    TYPE_ERROR("argument type not supported");
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
        return GMPy_Integer_Abs(arg0, context);

    if (IS_RATIONAL(arg0))
        return GMPy_Rational_Abs(arg0, context);

    if (IS_REAL(arg0))
        return GMPy_Real_Abs(arg0, context);

    if (IS_COMPLEX(arg0))
        return GMPy_Complex_Abs(arg0, context);

    TYPE_ERROR("abs() argument type not supported");
    return NULL;
}

