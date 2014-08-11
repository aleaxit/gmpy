/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_sub.c                                                             *
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

/* This file implements the - operator, gmpy2.add(), and context.sub().
 *
 * Public API
 * ==========
 * The following function is available as part of GMPY2's C API. A NULL value
 * for context implies the function should use the currently active context.
 *
 *   GMPy_Number_Sub(Number, Number, context|NULL)
 *
 * Private API
 * ===========
 *   GMPy_MPZ_Sub_Slot
 *   GMPy_MPQ_Sub_Slot
 *   GMPy_MPFR_Sub_Slot
 *   GMPy_MPC_Sub_Slot
 *
 *   GMPy_Integer_Sub(Integer, Integer, context|NULL)
 *   GMPy_Rational_Sub(Rational, Rational, context|NULL)
 *   GMPy_Real_Sub(Real, Real, context|NULL)
 *   GMPy_Complex_Sub(Complex, Complex, context|NULL)
 *
 *   GMPy_Context_Sub(context, args)
 *
 */

/* Subtract two Integer objects (see gmpy2_convert.c). If an error occurs,
 * NULL is returned and an exception is set. If either x or y can't be
 * converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
GMPy_Integer_Sub(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    if (CHECK_MPZANY(x)) {
        
#ifdef PY2
        if (PyInt_CheckExact(y)) {
            long i = PyInt_AS_LONG(y);

            if (i > 0) {
                mpz_sub_ui(result->z, MPZ(x), i);
            }
            else if (i < 0) {
                mpz_add_ui(result->z, MPZ(x), -i);
            }
            else {
                mpz_set(result->z, MPZ(x));
            }
            return (PyObject*)result;
        }
#endif

        if (PyLong_CheckExact(y)) {
            long temp = 0;
            Py_ssize_t i = Py_SIZE((PyLongObject*)y);
            
            switch (i) {
            CASE_NEGATIVE(temp, y)
                mpz_add_ui(result->z, MPZ(x), temp);
                return (PyObject*)result;
                
            case 0:
                mpz_set(result->z, MPZ(x));
                return (PyObject*)result;
                
            CASE_POSITIVE(temp, y)
                mpz_sub_ui(result->z, MPZ(x), temp);
                return (PyObject*)result;
            }
        }
        
        if (CHECK_MPZANY(y)) {
            mpz_sub(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
        
        if (PyIntOrLong_Check(y)) {
            mpz_t tempz;
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, y);
            mpz_sub(result->z, MPZ(x), tempz);
            mpz_cloc(tempz);
            return (PyObject*)result;
        }
    }

    if (CHECK_MPZANY(y)) {
#ifdef PY2
        if (PyInt_CheckExact(x)) {
            long i = PyInt_AS_LONG(x);

            if (i >= 0) {
                mpz_ui_sub(result->z, i, MPZ(y));
            }
            else {
                mpz_add_ui(result->z, MPZ(y), -i);
                mpz_neg(result->z, result->z);
            }
            return (PyObject*)result;
        }
#endif

        if (PyLong_CheckExact(x)) {
            long temp = 0;
            Py_ssize_t i = Py_SIZE((PyLongObject*)x);
            
            switch (i) {
            CASE_NEGATIVE(temp, x)
                mpz_add_ui(result->z, MPZ(y), temp);
                mpz_neg(result->z, result->z);
                return (PyObject*)result;
                
            case 0:
                mpz_set(result->z, MPZ(y));
                mpz_neg(result->z, result->z);
                return (PyObject*)result;
                
            CASE_POSITIVE(temp, x)
                mpz_ui_sub(result->z, temp, MPZ(y));
                return (PyObject*)result;
            }
        }
        
        if (PyIntOrLong_Check(x)) {
            mpz_t tempz;
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, x);
            mpz_sub(result->z, tempz, MPZ(y));
            mpz_cloc(tempz);
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        MPZ_Object *tempx, *tempy;

        tempx = GMPy_MPZ_From_Integer(x, context);
        tempy = GMPy_MPZ_From_Integer(y, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpz_sub(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement __sub__ for MPZ_Object. On entry, one of the two arguments must
 * be an MPZ_Object. If the other object is an Integer, add and return an
 * MPZ_Object. If the other object isn't an MPZ_Object, call the appropriate
 * function. If no appropriate function can be found, return NotImplemented.
 */

static PyObject *
GMPy_MPZ_Sub_Slot(PyObject *x, PyObject *y)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Sub(x, y, NULL);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Sub(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Sub(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Sub(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Subtract two Rational objects (see gmpy2_convert.h). Returns None and
 * raises TypeError if both objects are not valid rationals. Pympq_Sub_Rational
 * is intended to be called from GMPy_Number_Sub(). */

static PyObject *
GMPy_Rational_Sub(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPQ_Object *result;

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (MPQ_Check(x) && MPQ_Check(y)) {
        mpq_sub(result->q, MPQ(x), MPQ(y));
        return (PyObject*)result;
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        MPQ_Object *tempx, *tempy;

        tempx = GMPy_MPQ_From_Rational(x, context);
        tempy = GMPy_MPQ_From_Rational(y, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpq_sub(result->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement __sub__ for MPQ_Object. On entry, one of the two arguments must
 * be an MPQ_Object. If the other object is a Rational, subtrace and return
 * an MPQ_Object. If the other object isn't an MPQ_Object, call the
 * appropriate function. If no appropriate function can be found, return
 * NotImplemented.
 */

static PyObject *
GMPy_MPQ_Sub_Slot(PyObject *x, PyObject *y)
{
    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Sub(x, y, NULL);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Sub(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Sub(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/*   GMPy_Real_Sub(x, y, context) returns x-y using the provided context. If
 *   provided context is NULL, then the current context is used. If an error
 *   occurs, NULL is returned and an exception is set. If either x or y can't
 *   be converted to an mpfr, then Py_NotImplemented is returned.
 *   GMPy_Real_Sub() will not try to promote the result to a different type
 *   (i.e. mpc).
 *
 *   GMPy_mpfr_sub_fast(x, y) is the entry point for mpfr.__sub__.
 */

static PyObject *
GMPy_Real_Sub(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    /* This only processes mpfr if the exponent is still in-bounds. Need
     * to handle the rare case at the end. */

    if (MPFR_Check(x) && MPFR_Check(y)) {
        mpfr_clear_flags();
        result->rc = mpfr_sub(result->f, MPFR(x), MPFR(y),
                              GET_MPFR_ROUND(context));
        goto done;
    }

    if (MPFR_Check(x)) {
        if (PyIntOrLong_Check(y)) {
            mpz_t tempz;
            long temp;
            int overflow;

            temp = PyLong_AsLongAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpfr_clear_flags();
                result->rc = mpfr_sub_z(result->f, MPFR(x), tempz, GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_sub_si(result->f, MPFR(x), temp, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_sub_z(result->f, MPFR(x), MPZ(y), GET_MPFR_ROUND(context));
            goto done;
        }

        if (IS_RATIONAL(y)) {
            MPQ_Object *tempy;

            if (!(tempy = GMPy_MPQ_From_Number(y, context))) {
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_sub_q(result->f, MPFR(x), tempy->q, GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempy);
            goto done;
        }

        if (PyFloat_Check(y)) {
            mpfr_clear_flags();
            result->rc = mpfr_sub_d(result->f, MPFR(x), PyFloat_AS_DOUBLE(y), GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (MPFR_Check(y)) {
        if (PyIntOrLong_Check(x)) {
            mpz_t tempz;
            long temp;
            int overflow;

            temp = PyLong_AsLongAndOverflow(x, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, x);
                mpfr_clear_flags();
                result->rc = mpfr_sub_z(result->f, MPFR(y), tempz, GET_MPFR_ROUND(context));
                mpfr_neg(result->f, result->f, GET_MPFR_ROUND(context));
                mpz_cloc(tempz);
                goto done;
            }
            else {
                mpfr_clear_flags();
                result->rc = mpfr_sub_si(result->f, MPFR(y), temp, GET_MPFR_ROUND(context));
                mpfr_neg(result->f, result->f, GET_MPFR_ROUND(context));
                goto done;
            }
        }

        if (CHECK_MPZANY(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_sub_z(result->f, MPFR(y), MPZ(x), GET_MPFR_ROUND(context));
            mpfr_neg(result->f, result->f, GET_MPFR_ROUND(context));
            goto done;
        }

        if (IS_RATIONAL(x)) {
            MPQ_Object *tempx;

            if (!(tempx = GMPy_MPQ_From_Number(x, context))) {
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpfr_clear_flags();
            result->rc = mpfr_sub_q(result->f, MPFR(y), tempx->q, GET_MPFR_ROUND(context));
            mpfr_neg(result->f, result->f, GET_MPFR_ROUND(context));
            Py_DECREF((PyObject*)tempx);
            goto done;
        }

        if (PyFloat_Check(x)) {
            mpfr_clear_flags();
            result->rc = mpfr_sub_d(result->f, MPFR(y), PyFloat_AS_DOUBLE(x), GET_MPFR_ROUND(context));
            mpfr_neg(result->f, result->f, GET_MPFR_ROUND(context));
            goto done;
        }
    }

    if (IS_REAL(x) && IS_REAL(y)) {
        MPFR_Object *tempx, *tempy;

        tempx = GMPy_MPFR_From_Real(x, 1, context);
        tempy = GMPy_MPFR_From_Real(y, 1, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpfr_clear_flags();
        result->rc = mpfr_sub(result->f, MPFR(tempx), MPFR(tempy), GET_MPFR_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    GMPY_MPFR_CLEANUP(result, context, "subtraction");
    return (PyObject*)result;
}

/* Implement __sub__ for Pympfr. On entry, one of the two arguments must
 * be a Pympfr. If the other object is a Real, add and return a Pympfr.
 * If the other object isn't a Pympfr, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
GMPy_MPFR_Sub_Slot(PyObject *x, PyObject *y)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Sub(x, y, NULL);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Sub(x, y, NULL);

    Py_RETURN_NOTIMPLEMENTED;
}

/* GMPy_Complex_Sub(x, y, context) returns x-y using the provided context. If
 * context is NULL, then the current context is used. If an error occurs, NULL
 * is returned and an exception is set. If either x or y can't be converted to
 * an mpc, then Py_NotImplemented is returned. */

static PyObject *
GMPy_Complex_Sub(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context)))
        return NULL;

    if (MPC_Check(x) && MPC_Check(y)) {
        result->rc = mpc_sub(result->c, MPC(x), MPC(y),
                             GET_MPC_ROUND(context));
        goto done;
    }

    if (IS_COMPLEX(x) && IS_COMPLEX(y)) {
        MPC_Object *tempx, *tempy;

        tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
        tempy = GMPy_MPC_From_Complex(y, 1, 1, context);
        if (!tempx || !tempy) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        result->rc = mpc_sub(result->c, tempx->c, tempy->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        goto done;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  done:
    GMPY_MPC_CLEANUP(result, context, "subtraction");
    return (PyObject*)result;
}

/* Pympc_sub_fast() is called by mpc.__sub__. It just gets a borrowed reference
 * to the current context and call Pympc_Add_Complex(). Since mpc is the last
 * step of the numeric ladder, the NotImplemented return value from
 * Pympc_Sub_Complex() is correct and is just passed on. */

static PyObject *
GMPy_MPC_Sub_Slot(PyObject *x, PyObject *y)
{
    return GMPy_Complex_Sub(x, y, NULL);
}

static PyObject *
GMPy_Number_Sub(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_INTEGER(x) && IS_INTEGER(y))
        return GMPy_Integer_Sub(x, y, context);

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return GMPy_Rational_Sub(x, y, context);

    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Real_Sub(x, y, context);

    if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return GMPy_Complex_Sub(x, y, context);

    TYPE_ERROR("sub() argument type not supported");
    return NULL;
}

/* Implement context.sub() and gmpy2.sub(). */

PyDoc_STRVAR(GMPy_doc_sub,
"sub(x, y) -> number\n\n"
"Return x - y.");

PyDoc_STRVAR(GMPy_doc_context_sub,
"context.sub(x, y) -> number\n\n"
"Return x - y.");

static PyObject *
GMPy_Context_Sub(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("sub() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Sub(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1),
                           context);
}

