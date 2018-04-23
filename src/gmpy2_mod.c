/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mod.c                                                             *
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
 *   GMPy_MPZ_Mod_Slot
 *   GMPy_MPQ_Mod_Slot
 *   GMPy_MPFR_Mod_Slot
 *   GMPy_MPC_Mod_Slot
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
    MPZ_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            int error;
            native_si temp = GMPy_Integer_AsNative_siAndError(y, &error);

            if (!error) {
                if (temp > 0) {
                    mpz_fdiv_r_ui(result->z, MPZ(x), temp);
                }
                else if (temp == 0) {
                    ZERO_ERROR("division or modulo by zero");
                    Py_DECREF((PyObject*)result);
                    return NULL;
                }
                else {
                    mpz_cdiv_r_ui(result->z, MPZ(x), -temp);
                }
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, y);
                mpz_fdiv_r(result->z, MPZ(x), global.tempz);
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
            mpz_set_PyIntOrLong(global.tempz, x);
            mpz_fdiv_r(result->z, global.tempz, MPZ(y));
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        MPZ_Object *tempx, *tempy;

        tempx = GMPy_MPZ_From_Integer(x, context);
        tempy = GMPy_MPZ_From_Integer(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("could not convert Integer to mpz");
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
GMPy_MPZ_Mod_Slot(PyObject *x, PyObject *y)
{
    if (MPZ_Check(x) && MPZ_Check(y)) {
        MPZ_Object *result = NULL;

        if ((result = GMPy_MPZ_New(NULL))) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpz_fdiv_r(result->z, MPZ(x), MPZ(y));
        }
        return (PyObject*)result;
    }

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
    MPQ_Object *tempx, *tempy, *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        tempx = GMPy_MPQ_From_Number(x, context);
        tempy = GMPy_MPQ_From_Number(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("could not convert Rational to mpq");
            goto error;
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }

        mpq_div(result->q, tempx->q, tempy->q);
        mpz_fdiv_q(global.tempz, mpq_numref(result->q), mpq_denref(result->q));
        /* Need to calculate x - tempz * y. */
        mpq_set_z(result->q, global.tempz);
        mpq_mul(result->q, result->q, tempy->q);
        mpq_sub(result->q, tempx->q, result->q);
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
GMPy_MPQ_Mod_Slot(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mod(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mod(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mod(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_Real_Mod(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx = NULL, *tempy = NULL, *temp, *result;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(0, context);
    if (!result || !temp) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)temp);
        return NULL;
    }

    if (IS_REAL(x) && IS_REAL(y)) {
        tempx = GMPy_MPFR_From_Real(x, 1, context);
        tempy = GMPy_MPFR_From_Real(y, 1, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("could not convert Real to mpfr");
            goto error;
        }
        if (mpfr_zero_p(tempy->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.traps & TRAP_DIVZERO) {
                GMPY_DIVZERO("mod() modulo by zero");
                goto error;
            }
        }

        mpfr_clear_flags();

        if (mpfr_nan_p(tempx->f) || mpfr_nan_p(tempy->f) || mpfr_inf_p(tempx->f)) {

            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("mod() invalid operation");
                goto error;
            }
            mpfr_set_nan(result->f);
        }
        else if (mpfr_inf_p(tempy->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.traps & TRAP_INVALID) {
                GMPY_INVALID("mod() invalid operation");
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
            mpfr_fmod(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context));

            if (!mpfr_zero_p(result->f)) {
                if ((mpfr_sgn(tempy->f) < 0) != (mpfr_sgn(result->f) < 0)) {
                    mpfr_add(result->f, result->f, tempy->f, GET_MPFR_ROUND(context));
                }
            }
            else {
                mpfr_copysign(result->f, result->f, tempy->f, GET_MPFR_ROUND(context));
            }

            Py_DECREF((PyObject*)temp);
        }
        _GMPy_MPFR_Cleanup(&result, context);

        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PyObject *
GMPy_MPFR_Mod_Slot(PyObject *x, PyObject *y)
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
    TYPE_ERROR("can't take mod of complex number");
    return NULL;
}

static PyObject *
GMPy_MPC_Mod_Slot(PyObject *x, PyObject *y)
{
    return GMPy_Complex_Mod(x, y, NULL);
}

PyDoc_STRVAR(GMPy_doc_mod,
"mod(x, y) -> number\n\n"
"Return mod(x, y).\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"mpfr arguments to mod().");

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
"Return mod(x, y).\n"
"Note: overflow, underflow, and inexact exceptions are not supported for\n"
"mpfr arguments to context.mod().");

static PyObject *
GMPy_Context_Mod(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mod() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Mod(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
}

