/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mul.c                                                              *
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

/* This file implements the * operator, gmpy2.mul() and context.mul().
 *
 * Private API
 * ===========
 * The Python * operator calls the nb_add slot of a numeric type. This
 * file implements the following private functions:
 *
 *   GMPy_mpz_mul_fast; called by + via the nb_add slot of mpz
 *   GMPy_mpq_mul_fast; called by + via the nb_add slot of mpq
 *   GMPy_mpfr_mul_fast; called by + via the nb_add slot of mpfr
 *   GMPy_mpc_mul_fast; called by + via the nb_add slot of mpc
 *
 *   GMPy_Context_Mul; called by context.add()
 *
 * Public API
 * ==========
 * The following functions are availabe as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 * The first four functions check the type of the first argument and will set
 * an exception and return NULL if the check fails.
 *
 *   GMPy_Integer_Mul(Integer, Integer, context|NULL)
 *   GMPy_Rational_Mul(Rational, Rational, context|NULL)
 *   GMPy_Real_Mul(Real, Real, context|NULL)
 *   GMPy_Complex_Mul(Complex, Complex, context|NULL)
 *   GMPy_Number_Mul(Number, Number, context|NULL)
 *
 */

/* Multiply two Integer objects (see gmpy2_convert.c). If an error occurs,
 * NULL is returned and an exception is set. If either x or y can't be
 * converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
GMPy_Integer_Mul(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPZ_Object *result;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_mul(result->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(result->z, MPZ(x), temp_si);
            }
            return (PyObject*)result;
        }

        if (CHECK_MPZANY(y)) {
            mpz_mul(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
    }

    if (CHECK_MPZANY(y)) {
        if (PyIntOrLong_Check(x)) {
            temp_si = PyLong_AsSIAndOverflow(x, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, x);
                mpz_mul(result->z, MPZ(y), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(result->z, MPZ(y), temp_si);
            }
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        MPZ_Object *tempx, *tempy;

        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpz_mul(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement __mul__ for MPZ_Object. On entry, one of the two arguments must
 * be an MPZ_Object. If the other object is an Integer, add and return an
 * MPZ_Object. If the other object isn't an MPZ_Object, call the appropriate
 * function. If no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_mpz_mul_fast(PyObject *x, PyObject *y)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Mul(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mul(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Multiply two Rational objects (see convert.c/IS_RATIONAL). Returns None and
 * raises TypeError if both objects are not valid rationals. GMPy_Rational_Mul
 * is intended to be called from GMPy_Number_Mul. */

static PyObject *
GMPy_Rational_Mul(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New()))
        return NULL;

    if (MPQ_Check(x) && MPQ_Check(y)) {
        mpq_mul(result->q, MPQ(x), MPQ(y));
        return (PyObject*)result;
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        MPQ_Object *tempx, *tempy;

        tempx = GMPy_MPQ_From_Number_Temp(x);
        tempy = GMPy_MPQ_From_Number_Temp(y);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpq_mul(result->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement __mul__ for Pympq. On entry, one of the two arguments must
 * be a Pympq. If the other object is a Rational, add and return a Pympq.
 * If the other object isn't a Pympq, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_mpq_mul_fast(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mul(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Attempt to multiply two numbers and return an mpfr. The code path is
 * optimized by checking for mpfr objects first. Returns Py_NotImplemented if
 * both objects are not valid reals.  */

static PyObject *
GMPy_Real_Mul(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPFR_Object *result;

    if (!context)
        CURRENT_CONTEXT(context);

    SET_EXPONENT(context);

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    /* This only processes mpfr if the exponent is still in-bounds. Need
     * to handle the rare case at the end. */

    if (MPFR_CheckAndExp(x) && MPFR_CheckAndExp(y)) {
        mpfr_clear_flags();
        result->rc = mpfr_mul(result->f, MPFR(x), MPFR(y),
                              GET_MPFR_ROUND(context));
        goto done;
    }

    if (MPFR_CheckAndExp(x)) {
        if (PyIntOrLong_Check(y)) {
            mpz_t tempz;
            mpir_si temp_si;
            int overflow;

            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpfr_clear_flags();
                result->rc = mpfr_mul_z(result->f, MPFR(x), tempz,
                                        GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_mul_si(result->f, MPFR(x), temp_si,
                                         GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_mul_z(result->f, MPFR(x), MPZ(y),
                                    GET_MPFR_ROUND(context));
            goto done;
        }

        if (IS_RATIONAL(y) || IS_DECIMAL(y)) {
            MPQ_Object *tempy;

            if (!(tempy = GMPy_MPQ_From_Number_Temp(y))) {
                Py_DECREF(result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_mul_q(result->f, MPFR(x), tempy->q,
                                    GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempy);
            goto done;
        }

        if (PyFloat_Check(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_mul_d(result->f, MPFR(x), PyFloat_AS_DOUBLE(y),
                                    GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (MPFR_CheckAndExp(y)) {
        if (PyIntOrLong_Check(x)) {
            mpz_t tempz;
            mpir_si temp_si;
            int overflow;

            temp_si = PyLong_AsSIAndOverflow(x, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, x);
                mpfr_clear_flags();
                result->rc = mpfr_mul_z(result->f, MPFR(y), tempz,
                                        GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_mul_si(result->f, MPFR(y), temp_si,
                                         GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_mul_z(result->f, MPFR(y), MPZ(x),
                                    GET_MPFR_ROUND(context));
            goto done;
        }

        if (IS_RATIONAL(x) || IS_DECIMAL(x)) {
            MPQ_Object *tempx;

            if (!(tempx = GMPy_MPQ_From_Number_Temp(x))) {
                Py_DECREF(result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_mul_q(result->f, MPFR(y), tempx->q,
                                    GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempx);
            goto done;
        }

        if (PyFloat_Check(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_mul_d(result->f, MPFR(y), PyFloat_AS_DOUBLE(x),
                                    GET_MPFR_ROUND(context));
            goto done;
        }
    }

    /* In addition to handling PyFloat + PyFloat, the rare case when the
     * exponent bounds have been changed is handled here. See
     * Pympfr_From_Real() for details. */

    if (IS_REAL(x) && IS_REAL(y)) {
        MPFR_Object *tempx, *tempy;

        tempx = GMPy_MPFR_From_Real_Temp(x, context);
        tempy = GMPy_MPFR_From_Real_Temp(y, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        mpfr_clear_flags();
        result->rc = mpfr_mul(result->f, MPFR(tempx), MPFR(tempy),
                              GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    MPFR_CLEANUP_2(result, context, "multiplication");
    return (PyObject*)result;
}

/* Implement __mul__ for Pympfr. On entry, one of the two arguments must
 * be a Pympfr. If the other object is a Real, add and return a Pympfr.
 * If the other object isn't a Pympfr, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_mpfr_mul_fast(PyObject *x, PyObject *y)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* GMPy_Complex_Mul(x, y, context) returns x*y using the provided context. If
 * an error occurs, NULL is returned and an exception is set. If either x or
 * y can't be converted to an mpc, then Py_NotImplemented is returned. */

static PyObject *
GMPy_Complex_Mul(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPC_Object *result = NULL;

    if (!context)
        CURRENT_CONTEXT(context);

    SET_EXPONENT(context);

    if (!(result = (MPC_Object*)Pympc_new_bits_context(0, 0, context)))
        return NULL;

    if (MPC_CheckAndExp(x) && MPC_CheckAndExp(y)) {
        result->rc = mpc_mul(result->c, MPC(x), MPC(y),
                             GET_MPC_ROUND(context));
        goto done;
    }

    if (IS_COMPLEX(x) && IS_COMPLEX(y)) {
        MPC_Object *tempx, *tempy;

        tempx = GMPy_MPC_From_Complex_Temp(x, context);
        tempy = GMPy_MPC_From_Complex_Temp(y, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        result->rc = mpc_mul(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    MPC_CLEANUP_2(result, context, "multiplication");
    return (PyObject*)result;
}

/* Pympc_mul_fast() is called by mpc.__mul__. It just gets a borrowed reference
 * to the current context and call Pympc_Mul_Complex(). Since mpc is the last
 * step of the numeric ladder, the NotImplemented return value from
 * Pympc_Add_Complex() is correct and is just passed on. */

static PyObject *
GMPy_mpc_mul_fast(PyObject *x, PyObject *y)
{
    return GMPy_Complex_Mul(x, y, NULL);
}

static PyObject *
GMPy_Number_Mul(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Mul(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Mul(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Mul(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Mul(x, y, context);

    TYPE_ERROR("mul(): argument type not supported");
    return NULL;
}

/* Implement context.add() and gmpy2.add(). */

PyDoc_STRVAR(GMPy_doc_mul,
"mul(x, y) -> number\n\n"
"Return x * y.");

PyDoc_STRVAR(GMPy_doc_context_mul,
"context.mul(x, y) -> number\n\n"
"Return x * y.");

static PyObject *
GMPy_Context_Mul(PyObject *self, PyObject *args)
{
    PyObject *result;
    GMPyContextObject *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mul(): requires 2 arguments.");
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

    result = GMPy_Number_Mul(PyTuple_GET_ITEM(args, 0),
                             PyTuple_GET_ITEM(args, 1),
                             context);
    Py_DECREF((PyObject*)context);
    return result;
}

