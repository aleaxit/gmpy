/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr_misc.c                                                       *
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

PyDoc_STRVAR(GMPy_doc_function_f2q,
"f2q(x,[err]) -> mpq\n\n"
"Return the 'best' mpq approximating x to within relative error 'err'.\n"
"Default is the precision of x. Uses Stern-Brocot tree to find the\n"
"'best' approximation. An 'mpz' is returned if the the denominator\n"
"is 1. If 'err'<0, relative error is 2.0 ** err.");

static PyObject *
GMPy_Real_F2Q(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx, *tempy = NULL;
    PyObject *result;

    CHECK_CONTEXT(context);

    if (y) {
        if (!(tempy = GMPy_MPFR_From_Real(y, 1, context))) {
            return NULL;
        }
    }

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context))) {
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }
        
    result = stern_brocot(tempx, tempy, 0, 1, context);
    Py_DECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    return result;
}

static PyObject *
GMPy_Number_F2Q(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && (!y || IS_REAL(y)))
        return GMPy_Real_F2Q(x, y, context);

    TYPE_ERROR("f2q() argument types not supported");
    return NULL;
}

static PyObject *
GMPy_Context_F2Q(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) < 1 || PyTuple_GET_SIZE(args) > 2) {
        TYPE_ERROR("f2q() requires 1 or 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        return GMPy_Number_F2Q(PyTuple_GET_ITEM(args, 0), NULL, context);
    }
    else {
        return GMPy_Number_F2Q(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
    }
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_emin_min,
"get_emin_min() -> integer\n\n"
"Return the minimum possible exponent that can be set for 'mpfr'.");

static PyObject *
GMPy_MPFR_get_emin_min(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_emin_min());
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_emax_max,
"get_emax_max() -> integer\n\n"
"Return the maximum possible exponent that can be set for 'mpfr'.");

static PyObject *
GMPy_MPFR_get_emax_max(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_emax_max());
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_max_precision,
"get_max_precision() -> integer\n\n"
"Return the maximum bits of precision that can be used for calculations.\n"
"Note: to allow extra precision for intermediate calculations, avoid\n"
"setting precision close the maximum precision.");

static PyObject *
GMPy_MPFR_get_max_precision(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)MPFR_PREC_MAX);
}

PyDoc_STRVAR(GMPy_doc_mpfr_get_exp,
"get_exp(mpfr) -> integer\n\n"
"Return the exponent of an mpfr. Returns 0 for NaN or Infinity and\n"
"sets the erange flag and will raise an exception if trap_erange\n"
"is set.");

static PyObject *
GMPy_MPFR_get_exp(PyObject *self, PyObject *other)
{
    PyObject *result = NULL;
    Py_ssize_t exp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(MPFR_Check(other))) {
        TYPE_ERROR("get_exp() requires 'mpfr' argument");
        return NULL;
    }

    if (mpfr_regular_p(MPFR(other))) {
        exp = (Py_ssize_t)mpfr_get_exp(MPFR(other));
        result = PyIntOrLong_FromSsize_t(exp);
    }
    else if (mpfr_zero_p(MPFR(other))) {
        result = PyIntOrLong_FromSsize_t(0);
    }
    else {
        context->ctx.erange = 1;
        if (context->ctx.traps & TRAP_ERANGE) {
            GMPY_ERANGE("Can not get exponent from NaN or Infinity.");
        }
        else {
            result = PyIntOrLong_FromSsize_t(0);
        }
    }
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_exp,
"set_exp(mpfr, n) -> mpfr\n\n"
"Set the exponent of an mpfr to n. If n is outside the range of\n"
"valid exponents, set_exp() will set the erange flag and either\n"
"return the original value or raise an exception if trap_erange\n"
"is set.");

static PyObject *
GMPy_MPFR_set_exp(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *temp;
    mpfr_exp_t _oldemin, _oldemax, exp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2 ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 0)) ||
        !PyIntOrLong_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("set_exp() requires 'mpfr', 'integer' arguments");
        return NULL;
    }

    temp = PyTuple_GET_ITEM(args, 0);
    exp = (mpfr_exp_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(args, 1));
    if (exp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("exponent too large");
        return NULL;
    }
    
    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(temp)), context))) {
        return NULL;
    }
    
    _oldemin = mpfr_get_emin();
    _oldemax = mpfr_get_emax();
    mpfr_set_emin(context->ctx.emin);
    mpfr_set_emax(context->ctx.emax);

    mpfr_set(MPFR(result), MPFR(temp), GET_MPFR_ROUND(context));
    result->rc = mpfr_set_exp(MPFR(result), exp);
    
    mpfr_set_emin(_oldemin);
    mpfr_set_emax(_oldemax);

    if (result->rc) {
        context->ctx.erange = 1;
        if (context->ctx.traps & TRAP_ERANGE) {
            GMPY_ERANGE("new exponent is out-of-bounds");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_sign,
"set_sign(mpfr, bool) -> mpfr\n\n"
"If 'bool' is True, then return an 'mpfr' with the sign bit set.");

static PyObject *
GMPy_MPFR_set_sign(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2 ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 0)) ||
        !PyIntOrLong_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("set_sign() requires 'mpfr', 'boolean' arguments");
        return NULL;
    }
    
    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    result->rc = mpfr_setsign(MPFR(result), MPFR(PyTuple_GET_ITEM(args, 0)),
                              PyObject_IsTrue(PyTuple_GET_ITEM(args, 1)),
                              GET_MPFR_ROUND(context));

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_copy_sign,
"copy_sign(mpfr, mpfr) -> mpfr\n\n"
"Return an 'mpfr' composed of the first argument with the sign of the\n"
"second argument.");

static PyObject *
GMPy_MPFR_copy_sign(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 2 ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 0)) ||
        !MPFR_Check(PyTuple_GET_ITEM(args, 1))) {
        TYPE_ERROR("copy_sign() requires 'mpfr', 'boolean' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    result->rc = mpfr_copysign(MPFR(result), MPFR(PyTuple_GET_ITEM(args, 0)),
                               MPFR(PyTuple_GET_ITEM(args, 1)),
                               GET_MPFR_ROUND(context));

    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_nan,
"nan() -> mpfr\n\n"
"Return an 'mpfr' initialized to NaN (Not-A-Number).");

static PyObject *
GMPy_MPFR_set_nan(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_set_nan(result->f);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_inf,
"inf(n) -> mpfr\n\n"
"Return an 'mpfr' initialized to Infinity with the same sign as n.\n"
"If n is not given, +Infinity is returned.");

static PyObject *
GMPy_MPFR_set_inf(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long s = 1;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_Size(args) == 1) {
        s = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (s == -1 && PyErr_Occurred()) {
            TYPE_ERROR("inf() requires 'int' argument");
            return NULL;
        }
    }

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_set_inf(result->f, s < 0 ? -1 : 1);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_set_zero,
"zero(n) -> mpfr\n\n"
"Return an 'mpfr' inialized to 0.0 with the same sign as n.\n"
"If n is not given, +0.0 is returned.");

static PyObject *
GMPy_MPFR_set_zero(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long s = 1;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_Size(args) == 1) {
        s = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (s == -1 && PyErr_Occurred()) {
            TYPE_ERROR("zero() requires 'int' argument");
            return NULL;
        }
    }

    if ((result = GMPy_MPFR_New(0, context))) {
        mpfr_set_zero(result->f, s < 0 ? -1 : 1);
    }
    return (PyObject*)result;
}

#define MPFR_TEST_OTHER(NAME, msg) \
static PyObject * \
Pympfr_is_##NAME(PyObject *self, PyObject *other) \
{ \
    int res; \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT_SET_EXPONENT(context); \
    if(self && MPFR_Check(self)) { \
        Py_INCREF(self); \
    } \
    else if(MPFR_Check(other)) { \
        self = other; \
        Py_INCREF((PyObject*)self); \
    } \
    else if (!(self = (PyObject*)GMPy_MPFR_From_Real(other, 1, context))) { \
        PyErr_SetString(PyExc_TypeError, msg); \
        return NULL; \
    } \
    res = mpfr_##NAME##_p(MPFR(self)); \
    Py_DECREF(self); \
    if (res) \
        Py_RETURN_TRUE; \
    else\
        Py_RETURN_FALSE; \
}

PyDoc_STRVAR(doc_mpfr_integer_ratio,
"x.as_integer_ratio() -> (num, den)\n\n"
"Return the exact rational equivalent of an mpfr. Value is a tuple\n"
"for compatibility with Python's float.as_integer_ratio().");

static PyObject *
Pympfr_integer_ratio(PyObject *self, PyObject *args)
{
    MPZ_Object *num = 0, *den = 0;
    mpfr_exp_t temp, twocount;
    PyObject *result;
    CTXT_Object *context = NULL;

    if (mpfr_nan_p(MPFR(self))) {
        VALUE_ERROR("Cannot pass NaN to mpfr.as_integer_ratio.");
        return NULL;
    }

    if (mpfr_inf_p(MPFR(self))) {
        OVERFLOW_ERROR("Cannot pass Infinity to mpfr.as_integer_ratio.");
        return NULL;
    }

    num = GMPy_MPZ_New(context);
    den = GMPy_MPZ_New(context);
    if (!num || !den) {
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        return NULL;
    }

    if (mpfr_zero_p(MPFR(self))) {
        mpz_set_ui(num->z, 0);
        mpz_set_ui(den->z, 1);
    }
    else {
        temp = mpfr_get_z_2exp(num->z, MPFR(self));
        twocount = (mpfr_exp_t)mpz_scan1(num->z, 0);
        if (twocount) {
            temp += twocount;
            mpz_div_2exp(num->z, num->z, twocount);
        }
        mpz_set_ui(den->z, 1);
        if (temp > 0)
            mpz_mul_2exp(num->z, num->z, temp);
        else if (temp < 0)
            mpz_mul_2exp(den->z, den->z, -temp);
    }
    result = Py_BuildValue("(NN)", (PyObject*)num, (PyObject*)den);
    if (!result) {
        Py_DECREF((PyObject*)num);
        Py_DECREF((PyObject*)den);
    }
    return result;
}

PyDoc_STRVAR(doc_mpfr_mantissa_exp,
"x.as_mantissa_exp() -> (mantissa,exponent)\n\n"
"Return the mantissa and exponent of an mpfr.");

static PyObject *
Pympfr_mantissa_exp(PyObject *self, PyObject *args)
{
    MPZ_Object *mantissa = 0, *exponent = 0;
    mpfr_exp_t temp;
    PyObject *result;
    CTXT_Object *context = NULL;

    if (mpfr_nan_p(MPFR(self))) {
        VALUE_ERROR("Cannot pass NaN to mpfr.as_mantissa_exp.");
        return NULL;
    }

    if (mpfr_inf_p(MPFR(self))) {
        OVERFLOW_ERROR("Cannot pass Infinity to mpfr.as_mantissa_exp.");
        return NULL;
    }

    mantissa = GMPy_MPZ_New(context);
    exponent = GMPy_MPZ_New(context);
    if (!mantissa || !exponent) {
        Py_XDECREF((PyObject*)mantissa);
        Py_XDECREF((PyObject*)exponent);
        return NULL;
    }

    if (mpfr_zero_p(MPFR(self))) {
        mpz_set_ui(mantissa->z, 0);
        mpz_set_ui(exponent->z, 1);
    }
    else {
        temp = mpfr_get_z_2exp(mantissa->z, MPFR(self));
        mpz_set_si(exponent->z, temp);
    }
    result = Py_BuildValue("(NN)", (PyObject*)mantissa, (PyObject*)exponent);
    if (!result) {
        Py_DECREF((PyObject*)mantissa);
        Py_DECREF((PyObject*)exponent);
    }
    return result;
}

PyDoc_STRVAR(doc_mpfr_simple_fraction,
"x.as_simple_fraction([precision=0]) -> mpq\n\n"
"Return a simple rational approximation to x. The result will be\n"
"accurate to 'precision' bits. If 'precision' is 0, the precision\n"
"of 'x' will be used.");

static PyObject *
Pympfr_simple_fraction(PyObject *self, PyObject *args, PyObject *keywds)
{
    mpfr_prec_t prec = 0;
    static char *kwlist[] = {"precision", NULL};
    CTXT_Object *context = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &prec))
        return NULL;

    return (PyObject*)stern_brocot((MPFR_Object*)self, 0, prec, 0, context);
}

