/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpq_misc.c                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
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
GMPy_MPQ_Attrib_GetNumer(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, mpq_numref(self->q));
    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Attrib_GetReal(MPQ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPQ_Attrib_GetDenom(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, mpq_denref(self->q));
    return (PyObject*)result;
}

static PyObject *
GMPy_MPQ_Attrib_GetImag(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set_ui(result->z, 0);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_function_numer,
"numer(x) -> mpz\n\n"
"Return the numerator of x.");

static PyObject *
GMPy_MPQ_Function_Numer(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    MPQ_Object *tempq;
    CTXT_Object *context = NULL;

    if (!(result = (MPZ_Object*)GMPy_MPZ_New(context)))
        return NULL;

    if (!(tempq = GMPy_MPQ_From_Rational(other, context))) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_set(result->z, mpq_numref(tempq->q));
    Py_DECREF((PyObject*)tempq);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_function_denom,
"denom(x) -> mpz\n\n"
"Return the denominator of x.");

static PyObject *
GMPy_MPQ_Function_Denom(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    MPQ_Object *tempq;
    CTXT_Object *context = NULL;

    if (!(result = (MPZ_Object*)GMPy_MPZ_New(context)))
        return NULL;

    if (!(tempq = GMPy_MPQ_From_Rational(other, context))) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_set(result->z, mpq_denref(tempq->q));
    Py_DECREF((PyObject*)tempq);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_qdiv,
"qdiv(x[, y=1]) -> number\n\n"
"Return x/y as 'mpz' if possible, or as 'mpq' if x is not exactly\n"
"divisible by y.");

static PyObject *
GMPy_MPQ_Function_Qdiv(PyObject *self, PyObject *args)
{
    Py_ssize_t argc;
    PyObject *result = NULL, *x, *y;
    MPQ_Object *tempx = NULL, *tempy = NULL, *tempr = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    /* Validate the argument(s). */

    /* If there is only one argument, it should be either an integer or
     * rational type. If it is an integer, then immediately return an mpz().
     * If it is a rational type, convert it to an mpq and check the denominator.
     */

    argc = PyTuple_GET_SIZE(args);
    if (argc == 1) {
        x = PyTuple_GET_ITEM(args, 0);

        if (!IS_RATIONAL(x)) {
            goto arg_error;
        }

        if (IS_INTEGER(x)) {
            return (PyObject*)GMPy_MPZ_From_Integer(x, context);
        }

        if (!(tempx = GMPy_MPQ_From_Rational(x, context))) {
            return NULL;
        }

        if (mpz_cmp_ui(mpq_denref(tempx->q), 1) == 0) {
            if ((result = (PyObject*)GMPy_MPZ_New(context))) {
                mpz_set(MPZ(result), mpq_numref(tempx->q));
            }
            Py_DECREF((PyObject*)tempx);
            return result;
        }
        else {
            return (PyObject*)tempx;
        }
    }

    /* If there are two rational arguments, just convert them both to mpq,
     * divide, and then check the denominator.
     */

    if (argc == 2) {
        x = PyTuple_GET_ITEM(args, 0);
        y = PyTuple_GET_ITEM(args, 1);

        if (!IS_RATIONAL(x) || !IS_RATIONAL(y)) {
            goto arg_error;
        }

        if (!(tempx = GMPy_MPQ_From_Rational(x, context)) ||
            !(tempy = GMPy_MPQ_From_Rational(y, context))) {
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }

        if (mpq_sgn(tempy->q) == 0) {
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            ZERO_ERROR("qdiv() division by zero");
            return NULL;
        }

        /* tempr contains the result of the division. */

        if (!(tempr = GMPy_MPQ_New(context))) {
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            return NULL;
        }

        mpq_div(tempr->q, tempx->q, tempy->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);

        if (mpz_cmp_ui(mpq_denref(tempr->q), 1) == 0) {
            if ((result = (PyObject*)GMPy_MPZ_New(context))) {
                mpz_set(MPZ(result), mpq_numref(tempr->q));
            }
            Py_DECREF((PyObject*)tempr);
            return result;
        }
        else {
            return (PyObject*)tempr;
        };
    }

  arg_error:
    TYPE_ERROR("qdiv() requires 1 or 2 integer or rational arguments");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_floor,
"Return greatest integer less than or equal to an mpq.");

static PyObject *
GMPy_MPQ_Method_Floor(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        mpz_fdiv_q(result->z, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_ceil,
"Return least integer greater than or equal to an mpq.");

static PyObject *
GMPy_MPQ_Method_Ceil(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        mpz_cdiv_q(result->z, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_trunc,
"Return integer portion of an mpq.");

static PyObject *
GMPy_MPQ_Method_Trunc(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPZ_New(context))) {
        mpz_tdiv_q(result->z, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_round, "Round an mpq to power of 10.");

static PyObject *
GMPy_MPQ_Method_Round(PyObject *self, PyObject *args)
{
    Py_ssize_t round_digits = 0;
    MPQ_Object *resultq;
    MPZ_Object *resultz;
    mpz_t temp, rem;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    /* If args is NULL or the size of args is 0, we just return an mpz. */

    if (!args || PyTuple_GET_SIZE(args) == 0) {
        if (!(resultz = GMPy_MPZ_New(context))) {
            return NULL;
        }

        mpz_init(rem);
        mpz_fdiv_qr(resultz->z, rem, mpq_numref(MPQ(self)), mpq_denref(MPQ(self)));
        mpz_mul_2exp(rem, rem, 1);
        if (mpz_cmp(rem, mpq_denref(MPQ(self))) > 0) {
            mpz_add_ui(resultz->z, resultz->z, 1);
        }
        else if (mpz_cmp(rem, mpq_denref(MPQ(self))) == 0) {
            if (mpz_odd_p(resultz->z)) {
                mpz_add_ui(resultz->z, resultz->z, 1);
            }
        }
        mpz_clear(rem);
        return (PyObject*)resultz;
    }

    if (PyTuple_GET_SIZE(args) > 1) {
        TYPE_ERROR("Too many arguments for __round__()");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        round_digits = PyIntOrLong_AsSsize_t(PyTuple_GET_ITEM(args, 0));
        if (round_digits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("__round__() requires 'int' argument");
            return NULL;
        }
    }

    if (!(resultq = GMPy_MPQ_New(context))) {
        return NULL;
    }

    mpz_init(temp);
    mpz_ui_pow_ui(temp, 10, round_digits > 0 ? round_digits : -round_digits);

    mpq_set(resultq->q, MPQ(self));
    if (round_digits > 0) {
        mpz_mul(mpq_numref(resultq->q), mpq_numref(resultq->q), temp);
        mpq_canonicalize(resultq->q);
        if (!(resultz = (MPZ_Object*)GMPy_MPQ_Method_Round((PyObject*)resultq, NULL))) {
            mpz_clear(temp);
            return NULL;
        }
        mpz_set(mpq_numref(resultq->q), resultz->z);
        Py_DECREF((PyObject*)resultz);
        mpz_set(mpq_denref(resultq->q), temp);
        mpz_clear(temp);
        mpq_canonicalize(resultq->q);
    }
    else {
        mpz_mul(mpq_denref(resultq->q), mpq_denref(resultq->q), temp);
        mpq_canonicalize(resultq->q);
        if (!(resultz = (MPZ_Object*)GMPy_MPQ_Method_Round((PyObject*)resultq, NULL))) {
            mpz_clear(temp);
            return NULL;
        }
        mpq_set_ui(resultq->q, 0, 1);
        mpz_mul(mpq_numref(resultq->q), resultz->z, temp);
        Py_DECREF((PyObject*)resultz);
        mpz_clear(temp);
        mpq_canonicalize(resultq->q);
    }
    return (PyObject*)resultq;
}

static int
GMPy_MPQ_NonZero_Slot(MPQ_Object *self)
{
    return mpq_sgn(self->q) != 0;
}

PyDoc_STRVAR(GMPy_doc_mpq_method_sizeof,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x. Note: deleted mpq objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
GMPy_MPQ_Method_Sizeof(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPQ_Object) + \
        (mpq_numref(MPQ(self))->_mp_alloc * sizeof(mp_limb_t)) + \
        (mpq_denref(MPQ(self))->_mp_alloc * sizeof(mp_limb_t)));
}

