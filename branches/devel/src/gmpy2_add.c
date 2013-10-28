/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_add.c                                                              *
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

/* This file implements the + operator and context.add().
 *
 * Private API
 * ===========
 * The Python + operator calls the nb_add slot of a numeric type. This
 * file implements the following private functions:
 *
 *   GMPy_mpz_add_fast; called by + via the nb_add slot of mpz
 *   GMPy_mpq_add_fast; called by + via the nb_add slot of mpq
 *   GMPy_mpfr_add_fast; called by + via the nb_add slot of mpfr
 *   GMPy_mpc_add_fast; called by + via the nb_add slot of mpc
 *
 *   GMPy_Context_Add; called by context.add()
 *
 * Public API
 * ==========
 * The following functions are availabe as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 * The first four functions check the type of the first argument and will set
 * an exception and return NULL if the check fails.
 *
 *   GMPy_Integer_Add(Integer, Integer, context|NULL)
 *   GMPy_Rational_Add(Rational, Rational, context|NULL)
 *   GMPy_Real_Add(Real, Real, context|NULL)
 *   GMPy_Complex_Add(Complex, Complex, context|NULL)
 *   GMPy_Number_Add(Number, Number, context|NULL)
 *
 */

/* Add two Integer objects (see convert.c/isInteger). If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted
 * into an mpz, Py_NotImplemented is returned. */

static PyObject *
GMPy_Integer_Add(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    PympzObject *result;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(result = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_add(result->z, Pympz_AS_MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si >= 0) {
                mpz_add_ui(result->z, Pympz_AS_MPZ(x), temp_si);
            }
            else {
                mpz_sub_ui(result->z, Pympz_AS_MPZ(x), -temp_si);
            }
            return (PyObject*)result;
        }

        if (CHECK_MPZANY(y)) {
            mpz_add(result->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
            return (PyObject*)result;
        }
    }

    if (CHECK_MPZANY(y)) {
        if (PyIntOrLong_Check(x)) {
            temp_si = PyLong_AsSIAndOverflow(x, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, x);
                mpz_add(result->z, Pympz_AS_MPZ(y), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_add_ui(result->z, Pympz_AS_MPZ(y), temp_si);
            }
            else {
                mpz_sub_ui(result->z, Pympz_AS_MPZ(y), -temp_si);
            }
            return (PyObject*)result;
        }
    }

    if (PyIntOrLong_Check(x) && PyIntOrLong_Check(y)) {
        PympzObject *tempx, *tempy;

        tempx = Pympz_From_PyLong(x);
        tempy = Pympz_From_PyLong(y);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Integer to mpz.");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpz_add(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement __add__ for Pympz. On entry, one of the two arguments must
 * be a Pympz. If the other object is an Integer, add and return a Pympz.
 * If the other object isn't a Pympz, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_mpz_add_fast(PyObject *x, PyObject *y)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Add(x, y, NULL);
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Add(x, y, NULL);
    else if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, NULL);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Add two Rational objects (see convert.c/isRational). Returns None and
 * raises TypeError if both objects are not valid rationals. Pympq_Add_Rational
 * is intended to be called from Pympany_Add_Number. */

static PyObject *
GMPy_Rational_Add(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    PympqObject *result;

    if (!(result = (PympqObject*)Pympq_new()))
        return NULL;

    if (Pympq_Check(x) && Pympq_Check(y)) {
        mpq_add(result->q, Pympq_AS_MPQ(x), Pympq_AS_MPQ(y));
        return (PyObject*)result;
    }

    if (isRational(x) && isRational(y)) {
        PympqObject *tempx, *tempy;

        tempx = Pympq_From_Number(x);
        tempy = Pympq_From_Number(y);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Rational to mpq.");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpq_add(result->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement __add__ for Pympq. On entry, one of the two arguments must
 * be a Pympq. If the other object is a Rational, add and return a Pympq.
 * If the other object isn't a Pympq, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_mpq_add_fast(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Add(x, y, NULL);
    else if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, NULL);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Addition can be performed by the equivalent of mpfr.__add__ or by
 * gmpy2.add().
 *
 *   GMPy_Real_Add(x, y, context) returns x+y using the provided context. If
 *   provided context is NULL, then the current context is used. If an error
 *   occurs, NULL is returned and an exception is set. If either x or y can't
 *   be converted to an mpfr, then Py_NotImplemented is returned.
*    GMPy_Real_Add() will not try to promote the result to a different type
 *   (i.e. mpc).
 *
 *   GMPy_mpfr_add_fast(x, y) is the entry point for mpfr.__add__.
 */

/* Attempt to add two numbers and return an mpfr. The code path is optimized by
 * checking for mpfr objects first. Returns Py_NotImplemented if both objects
 * are not valid reals.  */

static PyObject *
GMPy_Real_Add(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    PympfrObject *result;

    if (!context)
        CURRENT_CONTEXT(context);

    SET_EXPONENT(context);

    if (!(result = (PympfrObject*)Pympfr_new_context(context)))
        return NULL;

    /* This only processes mpfr if the exponent is still in-bounds. Need
     * to handle the rare case at the end. */

    if (Pympfr_CheckAndExp(x) && Pympfr_CheckAndExp(y)) {
        mpfr_clear_flags();
        result->rc = mpfr_add(result->f, Pympfr_AS_MPFR(x), Pympfr_AS_MPFR(y),
                              GET_MPFR_ROUND(context));
        goto done;
    }

    if (Pympfr_CheckAndExp(x)) {
        if (PyIntOrLong_Check(y)) {
            mpz_t tempz;
            mpir_si temp_si;
            int overflow;

            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpfr_clear_flags();
                result->rc = mpfr_add_z(result->f, Pympfr_AS_MPFR(x),
                                        tempz, GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_add_si(result->f, Pympfr_AS_MPFR(x),
                                         temp_si, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_add_z(result->f, Pympfr_AS_MPFR(x),
                                    Pympz_AS_MPZ(y), GET_MPFR_ROUND(context));
            goto done;
        }

        if (isRational(y) || isDecimal(y)) {
            PympqObject *tempy;

            if (!(tempy = Pympq_From_Number(y))) {
                SYSTEM_ERROR("Can not convert Rational or Decimal to 'mpq'");
                Py_DECREF(result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_add_q(result->f, Pympfr_AS_MPFR(x), tempy->q,
                                    GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempy);
            goto done;
        }

        if (PyFloat_Check(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_add_d(result->f, Pympfr_AS_MPFR(x),
                                    PyFloat_AS_DOUBLE(y),
                                    GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (Pympfr_CheckAndExp(y)) {
        if (PyIntOrLong_Check(x)) {
            mpz_t tempz;
            mpir_si temp_si;
            int overflow;

            temp_si = PyLong_AsSIAndOverflow(x, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, x);
                mpfr_clear_flags();
                result->rc = mpfr_add_z(result->f, Pympfr_AS_MPFR(y),
                                        tempz, GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_add_si(result->f, Pympfr_AS_MPFR(y),
                                         temp_si, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_add_z(result->f, Pympfr_AS_MPFR(y),
                                    Pympz_AS_MPZ(x), GET_MPFR_ROUND(context));
            goto done;
        }

        if (isRational(x) || isDecimal(x)) {
            PympqObject *tempx;

            if (!(tempx = Pympq_From_Number(x))) {
                SYSTEM_ERROR("Can not convert Rational or Decimal to 'mpq'");
                Py_DECREF(result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_add_q(result->f, Pympfr_AS_MPFR(y), tempx->q,
                                    GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempx);
            goto done;
        }

        if (PyFloat_Check(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_add_d(result->f, Pympfr_AS_MPFR(y),
                                    PyFloat_AS_DOUBLE(x),
                                    GET_MPFR_ROUND(context));
            goto done;
        }
    }

    /* In addition to handling PyFloat + PyFloat, the rare case when the
     * exponent bounds have been changed is handled here. See
     * Pympfr_From_Real() for details. */

    if (IS_REAL(x) && IS_REAL(y)) {
        PympfrObject *tempx, *tempy;

        tempx = Pympfr_From_Real_context(x, 0, context);
        tempy = Pympfr_From_Real_context(y, 0, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        mpfr_clear_flags();
        result->rc = mpfr_add(result->f, Pympfr_AS_MPFR(tempx),
                              Pympfr_AS_MPFR(tempy),
                              GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    MPFR_CLEANUP_2(result, context, "addition");
    return (PyObject*)result;
}

/* Implement __add__ for Pympfr. On entry, one of the two arguments must
 * be a Pympfr. If the other object is a Real, add and return a Pympfr.
 * If the other object isn't a Pympfr, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_mpfr_add_fast(PyObject *x, PyObject *y)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, NULL);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* GMPy_Complex_Add(x, y, context) returns x+y using the provided context. If
 * context is NULL, then the current context is used. If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted to
 * an mpc, then Py_NotImplemented is returned. */

static PyObject *
GMPy_Complex_Add(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    PympcObject *result = NULL;

    if (!context)
        CURRENT_CONTEXT(context);

    SET_EXPONENT(context);

    if (!(result = (PympcObject*)Pympc_new_bits_context(0, 0, context)))
        return NULL;

    if (Pympc_CheckAndExp(x) && Pympc_CheckAndExp(y)) {
        result->rc = mpc_add(result->c, Pympc_AS_MPC(x), Pympc_AS_MPC(y),
                             GET_MPC_ROUND(context));
        goto done;
    }

    if (isComplex(x) && isComplex(y)) {
        PympcObject *tempx, *tempy;

        tempx = Pympc_From_Complex_bits_context(x, 0, 0, context);
        tempy = Pympc_From_Complex_bits_context(y, 0, 0, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Can not convert Complex to 'mpc'");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        result->rc = mpc_add(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    MPC_CLEANUP_RESULT("addition");
    return (PyObject*)result;
}

/* Pympc_add_fast() is called by mpc.__add__. It just gets a borrowed reference
 * to the current context and call Pympc_Add_Complex(). Since mpc is the last
 * step of the numeric ladder, the NotImplemented return value from
 * Pympc_Add_Complex() is correct and is just passed on. */

static PyObject *
GMPy_mpc_add_fast(PyObject *x, PyObject *y)
{
    return GMPy_Complex_Add(x, y, NULL);
}

static PyObject *
GMPy_Number_Add(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    if (!context)
        CURRENT_CONTEXT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Add(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Add(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Add(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Add(x, y, context);

    TYPE_ERROR("add(): argument type not supported");
    return NULL;
}

/* Implement context.add() and gmpy2.add(). */

PyDoc_STRVAR(GMPy_doc_add,
"add(x, y) -> number\n\n"
"Return x + y.");

PyDoc_STRVAR(GMPy_doc_context_add,
"context.add(x, y) -> number\n\n"
"Return x + y.");

static PyObject *
GMPy_Context_Add(PyObject *self, PyObject *args)
{
    PyObject *result;
    GMPyContextObject *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("add(): requires 2 arguments.");
        return NULL;
    }

    if (self && GMPyContext_Check(self)) {
        /* If we are passed a read-only context, make a copy of it before
         * proceeding. Remember to decref context when we're done. */

        if (((GMPyContextObject*)self)->ctx.readonly) {
            context = (GMPyContextObject*)GMPyContext_context_copy(self, NULL);
            if (!context)
                return NULL;
        }
        else {
            context = (GMPyContextObject*)self;
            Py_INCREF((PyObject*)context);
        }
    }
    else {
        CURRENT_CONTEXT(context);
        Py_INCREF((PyObject*)context);
    }

    result = GMPy_Number_Add(PyTuple_GET_ITEM(args, 0),
                             PyTuple_GET_ITEM(args, 1),
                             context);
    Py_DECREF((PyObject*)context);
    return result;
}

