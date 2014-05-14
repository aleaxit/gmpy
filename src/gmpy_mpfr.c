/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpfr.c                                                             *
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

PyDoc_STRVAR(GMPy_doc_mpfr_factory,
"mpfr() -> mpfr(0.0)\n\n"
"      If no argument is given, return mpfr(0.0).\n\n"
"mpfr(n [, precision=0]) -> mpfr\n\n"
"      Return an 'mpfr' object after converting a numeric value. See\n"
"      below for the interpretation of precision.\n\n"
"mpfr(s [, precision=0 [, base=0]]) -> mpfr\n\n"
"      Return a new 'mpfr' object by converting a string s made of\n"
"      digits in the given base, possibly with fraction-part (with a\n"
"      period as a separator) and/or exponent-part (with an exponent\n"
"      marker 'e' for base<=10, else '@'). The base of the string\n"
"      representation must be 0 or in the interval [2,62]. If the base\n"
"      is 0, the leading digits of the string are used to identify the\n"
"      base: 0b implies base=2, 0x implies base=16, otherwise base=10\n"
"      is assumed.\n\n"
"Note: If a precision greater than or equal to 2 is specified, then it\n"
"      is used.\n\n"
"      A precision of 0 (the default) implies the precision of the\n"
"      current context is used.\n\n"
"      A precision of 1 minimizes the loss of precision by following\n"
"      these rules:\n"
"        1) If n is a radix-2 floating point number, then the full\n"
"           precision of n is retained.\n"
"        2) For all other n, the precision of the result is the context\n"
"           precision + guard_bits.\n" );

static PyObject *
GMPy_MPFR_Factory(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    PyObject *arg0 = NULL;
    int base = 10;
    Py_ssize_t argc, keywdc = 0;
    CTXT_Object *context = NULL;

    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t prec = 0;

    static char *kwlist_s[] = {"s", "precision", "base", NULL};
    static char *kwlist_n[] = {"n", "precision", NULL};

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    argc = PyTuple_Size(args);
    if (keywds) {
        keywdc = PyDict_Size(keywds);
    }

    if (argc + keywdc > 3) {
        TYPE_ERROR("mpfr() takes at most 3 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPFR_New(0, context))) {
            mpfr_set_ui(result->f, 0, MPFR_RNDN);
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpfr() requires at least one non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    /* A string can have both precision and base additional arguments. */
    if (PyStrOrUnicode_Check(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|li", kwlist_s,
                                              &arg0, &prec, &base)))
                return NULL;
        }

        if (prec < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
            return NULL;
        }

        if (base != 0 && (base < 2 || base > 62)) {
            VALUE_ERROR("base for mpfr() must be 0 or in the interval [2, 62]");
            return NULL;
        }

        return (PyObject*)GMPy_MPFR_From_PyStr(arg0, base, prec, context);
    }

    /* A number can only have precision additional argument. */
    if (IS_REAL(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist_n,
                                              &arg0, &prec)))
                return NULL;
        }

        if (prec < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
            return NULL;
        }

        return (PyObject*)GMPy_MPFR_From_Real(arg0, prec, context);
    }

    TYPE_ERROR("mpfr() requires numeric or string argument");
    return NULL;
}

PyDoc_STRVAR(doc_g_mpfr_f2q,
"f2q(x,[err]) -> mpq\n\n"
"Return the 'best' mpq approximating x to within relative error 'err'.\n"
"Default is the precision of x. Uses Stern-Brocot tree to find the\n"
"'best' approximation. An 'mpz' is returned if the the denominator\n"
"is 1. If 'err'<0, error sought is 2.0 ** err.");

static PyObject *
Pympfr_f2q(PyObject *self, PyObject *args)
{
    MPFR_Object *err = 0;
    PyObject *result;
    CTXT_Object *context = NULL;

    if (!PyArg_ParseTuple(args, "O&|O&", GMPy_MPFR_convert_arg, &self,
                          GMPy_MPFR_convert_arg, &err)) {
        TYPE_ERROR("f2q() requires 'mpfr', ['mpfr'] arguments");
        return NULL;
    }

    result = (PyObject*)stern_brocot((MPFR_Object*)self, err, 0, 1, context);
    Py_DECREF(self);
    Py_XDECREF((PyObject*)err);
    return result;
}

/* Implement the .precision attribute of an mpfr. */

static PyObject *
Pympfr_getprec_attrib(MPFR_Object *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_prec(self->f));
}

/* Implement the .rc attribute of an mpfr. */

static PyObject *
Pympfr_getrc_attrib(MPFR_Object *self, void *closure)
{
    return PyIntOrLong_FromLong((long)self->rc);
}

/* Implement the .imag attribute of an mpfr. */

static PyObject *
Pympfr_getimag_attrib(MPFR_Object *self, void *closure)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPFR_New(0, context)))
        mpfr_set_zero(result->f, 1);
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpfr. */

static PyObject *
Pympfr_getreal_attrib(MPFR_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

/* Implement the nb_bool slot. */

static int
Pympfr_nonzero(MPFR_Object *self)
{
    return !mpfr_zero_p(self->f);
}

/* Implement the conjugate() method. */

PyDoc_STRVAR(doc_mpfr_conjugate,
"x.conjugate() -> mpfr\n\n"
"Return the conjugate of x (which is just a copy of x since x is\n"
"not a complex number).");

static PyObject *
Pympfr_conjugate(PyObject *self, PyObject *args)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

PyDoc_STRVAR(doc_g_mpfr_get_emin_min,
"get_emin_min() -> integer\n\n"
"Return the minimum possible exponent that can be set for 'mpfr'.");

static PyObject *
Pympfr_get_emin_min(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_emin_min());
}

PyDoc_STRVAR(doc_g_mpfr_get_emax_max,
"get_emax_max() -> integer\n\n"
"Return the maximum possible exponent that can be set for 'mpfr'.");

static PyObject *
Pympfr_get_emax_max(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_emax_max());
}

PyDoc_STRVAR(doc_g_mpfr_get_max_precision,
"get_max_precision() -> integer\n\n"
"Return the maximum bits of precision that can be used for calculations.\n"
"Note: to allow extra precision for intermediate calculations, avoid\n"
"setting precision close the maximum precision.");

static PyObject *
Pympfr_get_max_precision(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)MPFR_PREC_MAX);
}

PyDoc_STRVAR(doc_g_mpfr_get_exp,
"get_exp(mpfr) -> integer\n\n"
"Return the exponent of an mpfr. Returns 0 for NaN or Infinity and\n"
"sets the erange flag and will raise an exception if trap_erange\n"
"is set.");

static PyObject *
Pympfr_get_exp(PyObject *self, PyObject *other)
{
    PyObject *result = 0;
    Py_ssize_t exp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("get_exp() requires 'mpfr' argument");

    if (mpfr_regular_p(MPFR(self))) {
        exp = (Py_ssize_t)mpfr_get_exp(MPFR(self));
        result = PyIntOrLong_FromSsize_t((Py_ssize_t)exp);
    }
    else if (mpfr_zero_p(MPFR(self))) {
        Py_DECREF(self);
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
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_set_exp,
"set_exp(mpfr, n) -> mpfr\n\n"
"Set the exponent of an mpfr to n. If n is outside the range of\n"
"valid exponents, set_exp() will set the erange flag and either\n"
"return the original value or raise an exception if trap_erange\n"
"is set.");

static PyObject *
Pympfr_set_exp(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *temp;
    long exp = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!PyArg_ParseTuple(args, "O!l", &MPFR_Type, &temp, &exp)) {
        TYPE_ERROR("set_exp() requires 'mpfr', 'integer' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(temp)), context)))
        return NULL;

    mpfr_set(MPFR(result), MPFR(temp), GET_MPFR_ROUND(context));
    result->rc = mpfr_set_exp(MPFR(result), exp);

    if (result->rc) {
        context->ctx.erange = 1;
        if (context->ctx.traps & TRAP_ERANGE) {
            GMPY_ERANGE("New exponent is out-of-bounds.");
            Py_DECREF(result);
            return NULL;
        }
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_sign,
"set_sign(mpfr, bool) -> mpfr\n\n"
"If 'bool' is True, then return an 'mpfr' with the sign bit set.");

static PyObject *
Pympfr_set_sign(PyObject *self, PyObject *args)
{
    MPFR_Object *result = 0;
    PyObject *boolean = 0;
    int s;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!PyArg_ParseTuple(args, "O&O", GMPy_MPFR_convert_arg, &self, &boolean)) {
        TYPE_ERROR("set_sign() requires 'mpfr', 'boolean' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    s = PyObject_IsTrue(boolean);
    if (s == -1) {
        TYPE_ERROR("set_sign() requires 'mpfr', 'boolean' arguments");
        Py_DECREF(self);
        Py_DECREF(boolean);
        Py_DECREF(result);
        return NULL;
    }

    result->rc = mpfr_setsign(MPFR(result), MPFR(self),
                              s, context->ctx.mpfr_round);

    Py_DECREF(self);
    Py_DECREF(boolean);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_copy_sign,
"copy_sign(mpfr, mpfr) -> mpfr\n\n"
"Return an 'mpfr' composed of the first argument with the sign of the\n"
"second argument.");

static PyObject *
Pympfr_copy_sign(PyObject *self, PyObject *args)
{
    MPFR_Object *result = 0;
    PyObject *other = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!PyArg_ParseTuple(args, "O&O&", GMPy_MPFR_convert_arg, &self,
                          GMPy_MPFR_convert_arg, &other)) {
        TYPE_ERROR("copy_sign() requires 'mpfr', 'mpfr' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    result->rc = mpfr_copysign(MPFR(result), MPFR(self),
                              MPFR(other), context->ctx.mpfr_round);

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_nan,
"nan() -> mpfr\n\n"
"Return an 'mpfr' initialized to NaN (Not-A-Number).");

static PyObject *
Pympfr_set_nan(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPFR_New(0, context)))
        mpfr_set_nan(result->f);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_inf,
"inf(n) -> mpfr\n\n"
"Return an 'mpfr' initialized to Infinity with the same sign as n.\n"
"If n is not given, +Infinity is returned.");

static PyObject *
Pympfr_set_inf(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long s = 1;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (PyTuple_Size(args) == 1) {
        s = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (s == -1 && PyErr_Occurred()) {
            TYPE_ERROR("inf() requires 'int' argument");
            return NULL;
        }
    }

    if ((result = GMPy_MPFR_New(0, context)))
        mpfr_set_inf(result->f, s<0?-1:1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_zero,
"zero(n) -> mpfr\n\n"
"Return an 'mpfr' inialized to 0.0 with the same sign as n.\n"
"If n is not given, +0.0 is returned.");

static PyObject *
Pympfr_set_zero(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long s = 1;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (PyTuple_Size(args) == 1) {
        s = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (s == -1 && PyErr_Occurred()) {
            TYPE_ERROR("zero() requires 'int' argument");
            return NULL;
        }
    }

    if ((result = GMPy_MPFR_New(0, context)))
        mpfr_set_zero(result->f, s<0?-1:1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_is_signed,
"is_signed(x) -> boolean\n\n"
"Return True if the sign bit of x is set.");

static PyObject *
Pympfr_is_signed(PyObject *self, PyObject *other)
{
    int res;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if(self && MPFR_Check(self)) {
        Py_INCREF(self);
    }
    else if(MPFR_Check(other)) {
        self = other;
        Py_INCREF((PyObject*)self);
    }
    else if (!(self = (PyObject*)GMPy_MPFR_From_Real(other, 1, context))) {
        TYPE_ERROR("is_signed() requires 'mpfr' argument");
        return NULL;
    }
    res = mpfr_signbit(MPFR(self));
    Py_DECREF(self);
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
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

PyDoc_STRVAR(doc_g_mpfr_is_regular,
"is_regular(x) -> boolean\n\n"
"Return True if x is not zero, NaN, or Infinity; False otherwise.");

MPFR_TEST_OTHER(regular, "is_regular() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_integer,
"x.is_integer() -> boolean\n\n"
"Return True if x is an integer; False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_integer,
"is_integer(x) -> boolean\n\n"
"Return True if x is an integer; False otherwise.");

MPFR_TEST_OTHER(integer, "is_integer() requires 'mpfr' argument");

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

PyDoc_STRVAR(doc_mpfr_root,
"root(x, n) -> mpfr\n\n"
"Return n-th root of x. The result always an 'mpfr'.");

static PyObject *
Pympfr_root(PyObject *self, PyObject *args)
{
    long n;
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_REQ_CLONG(&n, "root() requires 'mpfr','int' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    if (n <= 0) {
        VALUE_ERROR("n must be > 0");
        goto done;
    }

    mpfr_clear_flags();
    result->rc = mpfr_root(result->f, MPFR(self), n,
                           context->ctx.mpfr_round);

    MPFR_CLEANUP_SELF("root()");
}

PyDoc_STRVAR(doc_g_mpfr_round2,
"round2(x[, n]) -> mpfr\n\n"
"Return x rounded to n bits. Uses default precision if n is not\n"
"specified. See round_away() to access the mpfr_round() function.");

static PyObject *
Pympfr_round2(PyObject *self, PyObject *args)
{
    mpfr_prec_t prec;
    MPFR_Object *result = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);
    prec = context->ctx.mpfr_prec;

    PARSE_ONE_MPFR_OPT_CLONG(&prec,
            "round2() requires 'mpfr',['int'] arguments");

    if (prec < MPFR_PREC_MIN || prec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid precision");
        goto done;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(self)), context))) {
        goto done;
    }

    mpfr_clear_flags();
    /* Duplicate the code from Pympfr_pos. */
    mpfr_set(result->f, MPFR(self), context->ctx.mpfr_round);
    result->round_mode = ((MPFR_Object*)self)->round_mode;
    result->rc = ((MPFR_Object*)self)->rc;
    result->rc = mpfr_check_range(result->f, result->rc, result->round_mode);
    result->rc = mpfr_prec_round(result->f, prec, context->ctx.mpfr_round);

    MPFR_CLEANUP_SELF("round2()");
}

PyDoc_STRVAR(doc_g_mpfr_round10,
"__round__(x[, n = 0]) -> mpfr\n\n"
"Return x rounded to n decimal digits before (n < 0) or after (n > 0)\n"
"the decimal point. Rounds to an integer if n is not specified.");

static PyObject *
Pympfr_round10(PyObject *self, PyObject *args)
{
    Py_ssize_t digits = 0;
    mpz_t temp;
    MPFR_Object *resultf = 0;
    MPZ_Object *resultz;
    CTXT_Object *context = NULL;

    /* If the size of args is 0, we just return an mpz. */

    if (PyTuple_GET_SIZE(args) == 0) {
        if ((resultz = GMPy_MPZ_New(context))) {
            if (mpfr_nan_p(MPFR(self))) {
                Py_DECREF((PyObject*)resultz);
                VALUE_ERROR("'mpz' does not support NaN");
                return NULL;
            }
            if (mpfr_inf_p(MPFR(self))) {
                Py_DECREF((PyObject*)resultz);
                OVERFLOW_ERROR("'mpz' does not support Infinity");
                return NULL;
            }
            /* return code is ignored */
            mpfr_get_z(resultz->z, MPFR(self), MPFR_RNDN);
        }
        return (PyObject*)resultz;
    }

    /* Now we need to return an mpfr, so handle the simple cases. */

    if (!mpfr_regular_p(MPFR(self))) {
        Py_INCREF(self);
        return self;
    }

    if (PyTuple_GET_SIZE(args) > 1) {
        TYPE_ERROR("Too many arguments for __round__().");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        digits = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (digits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("__round__() requires 'int' argument");
            return NULL;
        }
    }

    /* TODO: better error analysis, or else convert the mpfr to an exact
     * fraction, round the fraction, and then convert back to an mpfr.
     */

    resultf = GMPy_MPFR_New(mpfr_get_prec(MPFR(self))+100, context);
    if (!resultf)
        return NULL;

    mpz_inoc(temp);
    mpz_ui_pow_ui(temp, 10, digits > 0 ? digits : -digits);
    if (digits >= 0) {
        mpfr_mul_z(resultf->f, MPFR(self), temp, MPFR_RNDN);
    }
    else {
        mpfr_div_z(resultf->f, MPFR(self), temp, MPFR_RNDN);
    }

    mpfr_rint(resultf->f, resultf->f, MPFR_RNDN);

    if (digits >= 0) {
        mpfr_div_z(resultf->f, resultf->f, temp, MPFR_RNDN);
    }
    else {
        mpfr_mul_z(resultf->f, resultf->f, temp, MPFR_RNDN);
    }
    mpfr_prec_round(resultf->f, mpfr_get_prec(MPFR(self)), MPFR_RNDN);

    mpz_cloc(temp);
    return((PyObject*)resultf);
}

PyDoc_STRVAR(doc_g_mpfr_reldiff,
"reldiff(x, y) -> mpfr\n\n"
"Return the relative difference between x and y. Result is equal to\n"
"abs(x-y)/x.");

static PyObject *
Pympfr_reldiff(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "reldiff() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    /* mpfr_reldiff doesn't guarantee correct rounding and doesn't appear
     * to set any exceptions.
     */
    mpfr_reldiff(result->f, MPFR(self), MPFR(other),
                 context->ctx.mpfr_round);
    result->rc = 0;
    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympfr_sign(PyObject *self, PyObject *other)
{
    long sign;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("sign() requires 'mpfr' argument");

    mpfr_clear_flags();
    sign = mpfr_sgn(MPFR(self));

    MERGE_FLAGS;
    CHECK_ERANGE("range error in 'mpfr' sign(), NaN argument");

  done:
    Py_DECREF((PyObject*)self);
    if (PyErr_Occurred())
        return NULL;
    else
        return PyIntOrLong_FromLong(sign);
}

#define MPFR_MONOP(NAME) \
static PyObject * \
Py##NAME(MPFR_Object *x) \
{ \
    MPFR_Object *r; \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT_SET_EXPONENT(context); \
    if (!(r = GMPy_MPFR_New(0, context))) \
        return NULL; \
    if (MPFR_Check(x)) { \
        r->rc = NAME(r->f, x->f, context->ctx.mpfr_round); \
    } \
    else { \
        mpfr_set(r->f, x->f, context->ctx.mpfr_round); \
        r->round_mode = x->round_mode; \
        r->rc = x->rc; \
        mpfr_clear_flags(); \
        mpfr_check_range(r->f, r->rc, r->round_mode); \
        r->rc = NAME(r->f, r->f, context->ctx.mpfr_round); \
        MERGE_FLAGS; \
        CHECK_FLAGS(#NAME "()"); \
    } \
  done: \
    return (PyObject *) r; \
}


#define MPFR_UNIOP_NOROUND(NAME) \
static PyObject * \
Pympfr_##NAME(PyObject* self, PyObject *other) \
{ \
    MPFR_Object *result; \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT_SET_EXPONENT(context); \
    PARSE_ONE_MPFR_OTHER(#NAME "() requires 'mpfr' argument"); \
    if (!(result = GMPy_MPFR_New(0, context))) goto done; \
    mpfr_clear_flags(); \
    result->rc = mpfr_##NAME(result->f, MPFR(self)); \
    MPFR_CLEANUP_SELF(#NAME "()"); \
}

PyDoc_STRVAR(doc_mpfr_ceil,
"x.__ceil__() -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

PyDoc_STRVAR(doc_g_mpfr_ceil,
"ceil(x) ->mpfr\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

MPFR_UNIOP_NOROUND(ceil)

PyDoc_STRVAR(doc_mpfr_floor,
"x.__floor__() -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer <= x.");

PyDoc_STRVAR(doc_g_mpfr_floor,
"floor(x) -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer <= x.");

MPFR_UNIOP_NOROUND(floor);

PyDoc_STRVAR(doc_mpfr_trunc,
"x.__trunc__() -> mpfr\n\n"
"Return an 'mpfr' that is truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

PyDoc_STRVAR(doc_g_mpfr_trunc,
"trunc(x) -> mpfr\n\n"
"Return an 'mpfr' that is x truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

MPFR_UNIOP_NOROUND(trunc)

PyDoc_STRVAR(doc_g_mpfr_round_away,
"round_away(x) -> mpfr\n\n"
"Return an 'mpfr' that is x rounded to the nearest integer,\n"
"with ties rounded away from 0.");

static PyObject *
Pympfr_round_away(PyObject* self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("round_away() requires 'mpfr' argument");
    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;
    mpfr_clear_flags();
    result->rc = mpfr_round(result->f, MPFR(self));
    MPFR_CLEANUP_SELF("round_away()");
}

#define MPFR_UNIOP(NAME) \
static PyObject * \
Pympfr_##NAME(PyObject* self, PyObject *other) \
{ \
    MPFR_Object *result; \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT_SET_EXPONENT(context); \
    PARSE_ONE_MPFR_OTHER(#NAME "() requires 'mpfr' argument"); \
    if (!(result = GMPy_MPFR_New(0, context))) goto done; \
    mpfr_clear_flags(); \
    result->rc = mpfr_##NAME(result->f, MPFR(self), context->ctx.mpfr_round); \
    MPFR_CLEANUP_SELF(#NAME "()"); \
}

PyDoc_STRVAR(doc_g_mpfr_modf,
"modf(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the integer and fractional portions\n"
"of x.");

static PyObject *
Pympfr_modf(PyObject *self, PyObject *other)
{
    MPFR_Object *s, *c;
    PyObject *result;
    int code;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("modf() requires 'mpfr' argument");

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result)
        goto done;

    mpfr_clear_flags();
    code = mpfr_modf(s->f, c->f, MPFR(self),
                     context->ctx.mpfr_round);
    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;
    SUBNORMALIZE(s);
    SUBNORMALIZE(c);
    MERGE_FLAGS;
    CHECK_FLAGS("modf()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)s);
        PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    }
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_lgamma,
"lgamma(x) -> (mpfr, int)\n\n"
"Return a tuple containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

static PyObject *
Pympfr_lgamma(PyObject* self, PyObject *other)
{
    PyObject *result;
    MPFR_Object *value;
    int signp = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("lgamma() requires 'mpfr' argument");

    value = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!value || !result)
        goto done;

    mpfr_clear_flags();
    value->rc = mpfr_lgamma(value->f, &signp, MPFR(self),
                            context->ctx.mpfr_round);
    SUBNORMALIZE(value);
    MERGE_FLAGS;
    CHECK_FLAGS("lgamma()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)value);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)value);
        PyTuple_SET_ITEM(result, 1, PyIntOrLong_FromLong((long)signp));
    }
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_jn,
"jn(x,n) -> mpfr\n\n"
"Return the first kind Bessel function of order n of x.");

static PyObject *
Pympfr_jn(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long n = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_REQ_CLONG(&n, "jn() requires 'mpfr','int' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_jn(result->f, n, MPFR(self),
                         context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF("jn()");
}

PyDoc_STRVAR(doc_g_mpfr_yn,
"yn(x,n) -> mpfr\n\n"
"Return the second kind Bessel function of order n of x.");

static PyObject *
Pympfr_yn(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    long n = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_REQ_CLONG(&n, "yn() requires 'mpfr','int' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_yn(result->f, n, MPFR(self),
                         context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF("yn()");
}

PyDoc_STRVAR(doc_g_mpfr_fmod,
"fmod(x, y) -> mpfr\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to 0.");

static PyObject *
Pympfr_fmod(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "fmod() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_fmod(result->f, MPFR(self),
                           MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("fmod()");
}

PyDoc_STRVAR(doc_g_mpfr_remainder,
"remainder(x, y) -> mpfr\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

static PyObject *
Pympfr_remainder(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "remainder() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_remainder(result->f, MPFR(self),
                                MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("remainder()");
}

PyDoc_STRVAR(doc_g_mpfr_remquo,
"remquo(x, y) -> (mpfr, int)\n\n"
"Return a tuple containing the remainder(x,y) and the low bits of the\n"
"quotient.");

static PyObject *
Pympfr_remquo(PyObject* self, PyObject *args)
{
    PyObject *result, *other;
    MPFR_Object *value;
    long quobits = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "remquo() requires 'mpfr', 'mpfr' argument");

    value = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!value || !result)
        goto done;

    mpfr_clear_flags();
    value->rc = mpfr_remquo(value->f, &quobits, MPFR(self),
                            MPFR(other), context->ctx.mpfr_round);
    SUBNORMALIZE(value);
    MERGE_FLAGS;
    CHECK_FLAGS("remquo()");

  done:
    Py_DECREF(self);
    Py_DECREF(other);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)value);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)value);
        PyTuple_SET_ITEM(result, 1, PyIntOrLong_FromLong(quobits));
    }
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_frexp,
"frexp(x) -> (int, mpfr)\n\n"
"Return a tuple containing the exponent and mantissa of x.");

static PyObject *
Pympfr_frexp(PyObject *self, PyObject *other)
{
    PyObject *result;
    MPFR_Object *value;
    mpfr_exp_t exp = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("frexp() requires 'mpfr' argument");

    value = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!value || !result)
        goto done;

    mpfr_clear_flags();
    value->rc = mpfr_frexp(&exp, value->f, MPFR(self),
                           context->ctx.mpfr_round);
    MERGE_FLAGS;
    CHECK_FLAGS("frexp()");

  done:
    Py_DECREF(self);
    Py_DECREF(other);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)value);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, PyIntOrLong_FromSsize_t((Py_ssize_t)exp));
        PyTuple_SET_ITEM(result, 1, (PyObject*)value);
    }
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_agm,
"agm(x, y) -> mpfr\n\n"
"Return arithmetic-geometric mean of x and y.");

static PyObject *
Pympfr_agm(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "agm() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_agm(result->f, MPFR(self),
                          MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("agm()");
}

PyDoc_STRVAR(doc_g_mpfr_max2,
"max2(x, y) -> mpfr\n\n"
"Return the maximum number of x and y. This function is deprecated.\n"
"Please use maxnum() instead.");

PyDoc_STRVAR(doc_g_mpfr_maxnum,
"maxnum(x, y) -> mpfr\n\n"
"Return the maximum number of x and y. If x and y are not 'mpfr', they are\n"
"converted to 'mpfr'. The result is rounded to match the current\n"
"context. If only one of x or y is a number, then that number is returned.");

static PyObject *
Pympfr_max2(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "max2() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_max(result->f, MPFR(self),
                          MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("max2()");
}

PyDoc_STRVAR(doc_g_mpfr_min2,
"min2(x, y) -> mpfr\n\n"
"Return the minimum of x and y. This function is deprecated.\n"
"Please use minnum() instead.");

PyDoc_STRVAR(doc_g_mpfr_minnum,
"minnum(x, y) -> mpfr\n\n"
"Return the minimum of x and y. If x and y are not 'mpfr', they are\n"
"converted to 'mpfr'. The result is rounded to match the current\n"
"context. If only one of x or y is a number, then that number is returned.");

static PyObject *
Pympfr_min2(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "min2() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_min(result->f, MPFR(self),
                          MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("min2()");
}

PyDoc_STRVAR(doc_g_mpfr_nexttoward,
"next_toward(y, x) -> mpfr\n\n"
"Return the next 'mpfr' from x in the direction of y.");

static PyObject *
Pympfr_nexttoward(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "next_toward() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(self)), context)))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, MPFR(self), context->ctx.mpfr_round);
    mpfr_nexttoward(result->f, MPFR(other));
    result->rc = 0;
    MPFR_CLEANUP_SELF_OTHER("next_toward()");
}

PyDoc_STRVAR(doc_g_mpfr_nextabove,
"next_above(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward +Infinity.");

static PyObject *
Pympfr_nextabove(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("next_above() requires 'mpfr' argument");

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(self)), context)))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, MPFR(self), context->ctx.mpfr_round);
    mpfr_nextabove(result->f);
    result->rc = 0;
    MPFR_CLEANUP_SELF("next_above()");
}

PyDoc_STRVAR(doc_g_mpfr_nextbelow,
"next_below(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward -Infinity.");

static PyObject *
Pympfr_nextbelow(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("next_below() requires 'mpfr' argument");

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(MPFR(self)), context)))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, MPFR(self), context->ctx.mpfr_round);
    mpfr_nextbelow(result->f);
    result->rc = 0;
    MPFR_CLEANUP_SELF("next_below()");
}

PyDoc_STRVAR(doc_g_mpfr_factorial,
"factorial(n) -> mpfr\n\n"
"Return the floating-point approximation to the factorial of n.\n\n"
"See fac(n) to get the exact integer result.");

static PyObject *
Pympfr_factorial(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    long n;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    n = clong_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("factorial() requires 'int' argument");
        return NULL;
    }

    if (n < 0) {
        VALUE_ERROR("factorial() of negative number");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    mpfr_clear_flags();
    mpfr_fac_ui(result->f, n, context->ctx.mpfr_round);

    MERGE_FLAGS;
    CHECK_FLAGS("factorial()");
  done:
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_is_lessgreater,
"is_lessgreater(x,y) -> boolean\n\n"
"Return True if x > y or x < y. Return False if x == y or either x\n"
"and/or y is NaN.");

static PyObject *
Pympfr_is_lessgreater(PyObject *self, PyObject *args)
{
    PyObject *other;
    int temp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "is_lessgreater() requires 'mpfr','mpfr' arguments");

    temp = mpfr_lessgreater_p(MPFR(self), MPFR(other));
    Py_DECREF(self);
    Py_DECREF(other);
    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_g_mpfr_is_unordered,
"is_unordered(x,y) -> boolean\n\n"
"Return True if either x and/or y is NaN.");

static PyObject *
Pympfr_is_unordered(PyObject *self, PyObject *args)
{
    PyObject *other;
    int temp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "unordered() requires 'mpfr','mpfr' arguments");

    temp = mpfr_unordered_p(MPFR(self), MPFR(other));
    Py_DECREF(self);
    Py_DECREF(other);
    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_g_mpfr_check_range,
"check_range(x) -> mpfr\n\n"
"Return a new 'mpfr' with exponent that lies within the current range\n"
"of emin and emax.");

static PyObject *
Pympfr_check_range(PyObject *self, PyObject *other)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (self && MPFR_Check(self)) {
        if ((result = GMPy_MPFR_New(mpfr_get_prec(MPFR(self)), context))) {
            mpfr_set(result->f, MPFR(self), context->ctx.mpfr_round);
            result->round_mode = ((MPFR_Object*)self)->round_mode;
            result->rc = ((MPFR_Object*)self)->rc;
            mpfr_clear_flags();
            result->rc = mpfr_check_range(result->f, result->rc,
                                          result->round_mode);
        }
    }
    else if (MPFR_Check(other)) {
        if ((result = GMPy_MPFR_New(mpfr_get_prec(MPFR(other)), context))) {
            mpfr_set(result->f, MPFR(other), context->ctx.mpfr_round);
            result->round_mode = ((MPFR_Object*)other)->round_mode;
            result->rc = ((MPFR_Object*)other)->rc;
            mpfr_clear_flags();
            result->rc = mpfr_check_range(result->f, result->rc,
                                          result->round_mode);
        }
    }
    else {
        TYPE_ERROR("check_range() requires 'mpfr' argument");
    }
    MERGE_FLAGS;
    CHECK_FLAGS("check_range()");
  done:
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_fsum,
"fsum(iterable) -> mpfr\n\n"
"Return an accurate sum of the values in the iterable.");

static PyObject *
Pympfr_fsum(PyObject *self, PyObject *other)
{
    MPFR_Object *temp, *result;
    mpfr_ptr *tab;
    int errcode;
    Py_ssize_t i, seq_length = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(result = GMPy_MPFR_New(0, context)))
        return NULL;

    if (!(other = PySequence_List(other))) {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }

    /* other contains a new list containing all the values from the
     * iterable. Now make sure each item in the list is an mpfr.
     */

    seq_length = PyList_GET_SIZE(other);
    for (i=0; i < seq_length; i++) {
        if (!(temp = GMPy_MPFR_From_Real(PyList_GET_ITEM(other, i), 1, context))) {
            Py_DECREF(other);
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("all items in iterable must be real numbers");
            return NULL;
        }

        errcode = PyList_SetItem(other, i,(PyObject*)temp);
        if (errcode < 0) {
            Py_DECREF(other);
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("all items in iterable must be real numbers");
            return NULL;
        }
    }

    /* create an array of pointers to the mpfr_t field of a Pympfr object */

    if (!(tab = (mpfr_ptr *)GMPY_MALLOC((sizeof(mpfr_srcptr) * seq_length)))) {
        Py_DECREF(other);
        Py_DECREF((PyObject*)result);
        return PyErr_NoMemory();
    }
    for (i=0; i < seq_length; i++) {
        temp = (MPFR_Object*)PyList_GET_ITEM(other, i);
        tab[i] = temp->f;
    }
    result->rc = mpfr_sum(result->f, tab, seq_length, context->ctx.mpfr_round);
    Py_DECREF(other);
    GMPY_FREE(tab);

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpfr_sizeof,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x.");

static PyObject *
Pympfr_sizeof(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPFR_Object) + \
        (((MPFR(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t));
}

#ifdef PY3
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) GMPy_MPFR_Add_Slot,         /* nb_add                  */
    (binaryfunc) GMPy_MPFR_Sub_Slot,         /* nb_subtract             */
    (binaryfunc) GMPy_MPFR_Mul_Slot,         /* nb_multiply             */
    (binaryfunc) GMPy_MPFR_Mod_Slot,         /* nb_remainder            */
    (binaryfunc) GMPy_MPFR_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,       /* nb_power                */
    (unaryfunc) GMPy_MPFR_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPFR_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPFR_Abs_Slot,          /* nb_absolute             */
    (inquiry) Pympfr_nonzero,                /* nb_bool                 */
        0,                                   /* nb_invert               */
        0,                                   /* nb_lshift               */
        0,                                   /* nb_rshift               */
        0,                                   /* nb_and                  */
        0,                                   /* nb_xor                  */
        0,                                   /* nb_or                   */
    (unaryfunc) GMPy_MPFR_Int_Slot,          /* nb_int                  */
        0,                                   /* nb_reserved             */
    (unaryfunc) GMPy_MPFR_Float_Slot,        /* nb_float                */
        0,                                   /* nb_inplace_add          */
        0,                                   /* nb_inplace_subtract     */
        0,                                   /* nb_inplace_multiply     */
        0,                                   /* nb_inplace_remainder    */
        0,                                   /* nb_inplace_power        */
        0,                                   /* nb_inplace_lshift       */
        0,                                   /* nb_inplace_rshift       */
        0,                                   /* nb_inplace_and          */
        0,                                   /* nb_inplace_xor          */
        0,                                   /* nb_inplace_or           */
    (binaryfunc) GMPy_MPFR_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_MPFR_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                   /* nb_inplace_floor_divide */
        0,                                   /* nb_inplace_true_divide  */
        0,                                   /* nb_index                */
};
#else
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) GMPy_MPFR_Add_Slot,         /* nb_add                  */
    (binaryfunc) GMPy_MPFR_Sub_Slot,         /* nb_subtract             */
    (binaryfunc) GMPy_MPFR_Mul_Slot,         /* nb_multiply             */
    (binaryfunc) GMPy_MPFR_TrueDiv_Slot,     /* nb_divide               */
    (binaryfunc) GMPy_MPFR_Mod_Slot,         /* nb_remainder            */
    (binaryfunc) GMPy_MPFR_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,       /* nb_power                */
    (unaryfunc) GMPy_MPFR_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPFR_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPFR_Abs_Slot,          /* nb_absolute             */
    (inquiry) Pympfr_nonzero,                /* nb_bool                 */
        0,                                   /* nb_invert               */
        0,                                   /* nb_lshift               */
        0,                                   /* nb_rshift               */
        0,                                   /* nb_and                  */
        0,                                   /* nb_xor                  */
        0,                                   /* nb_or                   */
        0,                                   /* nb_coerce               */
    (unaryfunc) GMPy_MPFR_Int_Slot,          /* nb_int                  */
    (unaryfunc) GMPy_MPFR_Long_Slot,         /* nb_long                 */
    (unaryfunc) GMPy_MPFR_Float_Slot,        /* nb_float                */
        0,                                   /* nb_oct                  */
        0,                                   /* nb_hex                  */
        0,                                   /* nb_inplace_add          */
        0,                                   /* nb_inplace_subtract     */
        0,                                   /* nb_inplace_multiply     */
        0,                                   /* nb_inplace_divide       */
        0,                                   /* nb_inplace_remainder    */
        0,                                   /* nb_inplace_power        */
        0,                                   /* nb_inplace_lshift       */
        0,                                   /* nb_inplace_rshift       */
        0,                                   /* nb_inplace_and          */
        0,                                   /* nb_inplace_xor          */
        0,                                   /* nb_inplace_or           */
    (binaryfunc) GMPy_MPFR_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_MPFR_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                   /* nb_inplace_floor_divide */
        0,                                   /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympfr_getseters[] =
{
    {"precision", (getter)Pympfr_getprec_attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)Pympfr_getrc_attrib, NULL, "return code", NULL},
    {"imag", (getter)Pympfr_getimag_attrib, NULL, "imaginary component", NULL},
    {"real", (getter)Pympfr_getreal_attrib, NULL, "real component", NULL},
    {NULL}
};

static PyMethodDef Pympfr_methods [] =
{
    { "__ceil__", Pympfr_ceil, METH_NOARGS, doc_mpfr_ceil },
    { "__floor__", Pympfr_floor, METH_NOARGS, doc_mpfr_floor },
    { "__format__", GMPy_MPFR_Format, METH_VARARGS, GMPy_doc_mpfr_format },
    { "__round__", Pympfr_round10, METH_VARARGS, doc_g_mpfr_round10 },
    { "__sizeof__", Pympfr_sizeof, METH_NOARGS, doc_mpfr_sizeof },
    { "__trunc__", Pympfr_trunc, METH_NOARGS, doc_mpfr_trunc },
    { "as_integer_ratio", Pympfr_integer_ratio, METH_NOARGS, doc_mpfr_integer_ratio },
    { "as_mantissa_exp", Pympfr_mantissa_exp, METH_NOARGS, doc_mpfr_mantissa_exp },
    { "as_simple_fraction", (PyCFunction)Pympfr_simple_fraction, METH_VARARGS | METH_KEYWORDS, doc_mpfr_simple_fraction },
    { "conjugate", Pympfr_conjugate, METH_NOARGS, doc_mpfr_conjugate },
    { "digits", GMPy_MPFR_Digits_Method, METH_VARARGS, GMPy_doc_mpfr_digits_method },
    { "is_finite", GMPy_MPFR_Is_Finite_Method, METH_NOARGS, GMPy_doc_method_is_finite },
    { "is_infinite", GMPy_MPFR_Is_Infinite_Method, METH_NOARGS, GMPy_doc_method_is_infinite },
    { "is_integer", Pympfr_is_integer, METH_NOARGS, doc_mpfr_is_integer },
    { "is_nan", GMPy_MPFR_Is_NAN_Method, METH_NOARGS, GMPy_doc_method_is_nan },
    { "is_zero", GMPy_MPFR_Is_Zero_Method, METH_NOARGS, GMPy_doc_method_is_zero },
    { NULL, NULL, 1 }
};

static PyTypeObject MPFR_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpfr",                                 /* tp_name          */
    sizeof(MPFR_Object),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPFR_Dealloc,         /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPFR_Repr_Slot,         /* tp_repr          */
    &mpfr_number_methods,                   /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) GMPy_MPFR_Hash_Slot,         /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPFR_Str_Slot,          /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "Multiple precision real",              /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympfr_methods,                         /* tp_methods       */
        0,                                  /* tp_members       */
    Pympfr_getseters,                       /* tp_getset        */
};

