/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_square.c                                                             *
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

/* Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_Square(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_Integer_Square(Integer, Integer, context|NULL)
 *   GMPy_Rational_Square(Rational, Rational, context|NULL)
 *   GMPy_Real_Square(Real, Real, context|NULL)
 *   GMPy_Complex_Square(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Square(context, args)
 *
 */

static PyObject *
GMPy_Integer_Square(PyObject *x, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *tempx = NULL;

    result = GMPy_MPZ_New(context);
    tempx = GMPy_MPZ_From_Integer(x, context);
    if (!result || !tempx) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpz_mul(result->z, tempx->z, tempx->z);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

static PyObject *
GMPy_Rational_Square(PyObject *x, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *tempx = NULL;

    result = GMPy_MPQ_New(context);
    tempx = GMPy_MPQ_From_Rational(x, context);
    if (!result || !tempx) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpq_mul(result->q, tempx->q, tempx->q);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Square(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);
    
    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    if (!result || !tempx) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();
    mpfr_mul(result->f, tempx->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPFR_CLEANUP(result, context, "square()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Square(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    result = GMPy_MPC_New(0, 0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        return NULL;
    }

    mpc_mul(result->c, tempx->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPC_CLEANUP(result, context, "square()");
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_square,
"square(x) -> number\n\n"
"Return x * x. If x is an integer, then the result is an 'mpz'.\n"
"If x is a rational, then the result is an 'mpq'. If x is a float,\n"
"then the result is an 'mpfr'. If x is a complex number, then the\n"
"result is an 'mpc'.");

PyDoc_STRVAR(GMPy_doc_context_square,
"context.square(x) -> number\n\n"
"Return x * x. If x is an integer, then the result is an 'mpz'.\n"
"If x is a rational, then the result is an 'mpq'. If x is a float,\n"
"then the result is an 'mpfr'. If x is a complex number, then the\n"
"result is an 'mpc'.");

static PyObject *
GMPy_Number_Square(PyObject *x, CTXT_Object *context)
{
    if (IS_INTEGER(x))
        return GMPy_Integer_Square(x, context);

    if (IS_RATIONAL(x))
        return GMPy_Rational_Square(x, context);

    if (IS_REAL(x))
        return GMPy_Real_Square(x, context);

    if (IS_COMPLEX(x))
        return GMPy_Complex_Square(x, context);

    TYPE_ERROR("square() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Square(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Square(other, context);
}

