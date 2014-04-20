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

/* NAME is used as part of the GMPy function name. It usually uses an upper-
 *      case first character.
 * FUNC is the component of the actual name used by MPFR and MPC.
 */

#define GMPY_MPFR_MPC_UNIOP(NAME, FUNC) \
static PyObject * \
GMPy_Real_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPFR_Object *result = NULL, *tempx = NULL; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPFR_New(0, context); \
    tempx = GMPy_MPFR_From_Real(x, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    mpfr_clear_flags(); \
    result->rc = mpfr_##FUNC(result->f, tempx->f, GET_MPFR_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    GMPY_MPFR_CLEANUP(result, context, #FUNC "()"); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Complex_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    MPC_Object *result = NULL, *tempx = NULL; \
    CHECK_CONTEXT(context); \
    result = GMPy_MPC_New(0, 0, context); \
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context); \
    if (!result || !tempx) { \
        Py_XDECREF((PyObject*)result); \
        Py_XDECREF((PyObject*)tempx); \
        return NULL; \
    } \
    result->rc = mpc_##FUNC(result->c, tempx->c, GET_MPC_ROUND(context)); \
    Py_DECREF((PyObject*)tempx); \
    GMPY_MPC_CLEANUP(result, context, #FUNC "()"); \
    return (PyObject*)result; \
} \
static PyObject * \
GMPy_Number_##NAME(PyObject *x, CTXT_Object *context) \
{ \
    if (IS_REAL(x)) \
        return GMPy_Real_##NAME(x, context); \
    if (IS_COMPLEX(x)) \
        return GMPy_Complex_##NAME(x, context); \
    TYPE_ERROR(#FUNC"() argument type not supported"); \
    return NULL; \
} \
static PyObject * \
GMPy_Function_##NAME(PyObject *self, PyObject *other) \
{ \
    return GMPy_Number_##NAME(other, NULL); \
} \
static PyObject * \
GMPy_Context_##NAME(PyObject *self, PyObject *args) \
{ \
    CTXT_Object *context = NULL; \
    if (PyTuple_GET_SIZE(args) != 1) { \
        TYPE_ERROR(#FUNC "() requires 1 argument"); \
        return NULL; \
    } \
    if (self && CTXT_Check(self)) { \
        context = (CTXT_Object*)self; \
    } \
    else { \
        CHECK_CONTEXT(context); \
    } \
    return GMPy_Number_##NAME(PyTuple_GET_ITEM(args, 0), context); \
}

PyDoc_STRVAR(GMPy_doc_context_sin,
"context.sin(x) -> number\n\n"
"Return sine of x; x in radians.");
PyDoc_STRVAR(GMPy_doc_function_sin,
"sin(x) -> number\n\n"
"Return sine of x; x in radians.");
GMPY_MPFR_MPC_UNIOP(Sin, sin)

PyDoc_STRVAR(GMPy_doc_context_cos,
"context.cos(x) -> number\n\n"
"Return cosine of x; x in radians.");
PyDoc_STRVAR(GMPy_doc_function_cos,
"cos(x) -> number\n\n"
"Return cosine of x; x in radians.");
GMPY_MPFR_MPC_UNIOP(Cos, cos)

PyDoc_STRVAR(GMPy_doc_context_tan,
"context.tan(x) -> number\n\n"
"Return tangent of x; x in radians.");
PyDoc_STRVAR(GMPy_doc_function_tan,
"tan(x) -> number\n\n"
"Return tangent of x; x in radians.");
GMPY_MPFR_MPC_UNIOP(Tan, tan)

PyDoc_STRVAR(GMPy_doc_context_atan,
"context.atan(x) -> number\n\n"
"Return inverse tangent of x; x in radians.");
PyDoc_STRVAR(GMPy_doc_function_atan,
"atan(x) -> number\n\n"
"Return inverse tangent of x; x in radians.");
GMPY_MPFR_MPC_UNIOP(Atan, atan)

PyDoc_STRVAR(GMPy_doc_context_sinh,
"context.sinh(x) -> number\n\n"
"Return hyperbolic sine of x.");
PyDoc_STRVAR(GMPy_doc_function_sinh,
"sinh(x) -> number\n\n"
"Return hyperbolic sine of x.");
GMPY_MPFR_MPC_UNIOP(Sinh, sinh)

PyDoc_STRVAR(GMPy_doc_context_cosh,
"context.cosh(x) -> number\n\n"
"Return hyperbolic cosine of x.");
PyDoc_STRVAR(GMPy_doc_function_cosh,
"cosh(x) -> number\n\n"
"Return hyperbolic cosine of x.");
GMPY_MPFR_MPC_UNIOP(Cosh, cosh)

PyDoc_STRVAR(GMPy_doc_context_tanh,
"context.tanh(x) -> number\n\n"
"Return hyperbolic tangent of x.");
PyDoc_STRVAR(GMPy_doc_function_tanh,
"tanh(x) -> number\n\n"
"Return hyperbolic tangent of x.");
GMPY_MPFR_MPC_UNIOP(Tanh, tanh)

PyDoc_STRVAR(GMPy_doc_context_asinh,
"context.asinh(x) -> number\n\n"
"Return inverse hyperbolic sine of x.");
PyDoc_STRVAR(GMPy_doc_function_asinh,
"asinh(x) -> number\n\n"
"Return inverse hyperbolic sine of x.");
GMPY_MPFR_MPC_UNIOP(Asinh, asinh)

PyDoc_STRVAR(GMPy_doc_context_acosh,
"context.acosh(x) -> number\n\n"
"Return inverse hyperbolic cosine of x.");
PyDoc_STRVAR(GMPy_doc_function_acosh,
"acosh(x) -> number\n\n"
"Return inverse hyperbolic cosine of x.");
GMPY_MPFR_MPC_UNIOP(Acosh, acosh)



PyDoc_STRVAR(doc_g_mpfr_sec,
"sec(x) -> mpfr\n\n"
"Return secant of x; x in radians.");

MPFR_UNIOP(sec)

PyDoc_STRVAR(doc_g_mpfr_csc,
"csc(x) -> mpfr\n\n"
"Return cosecant of x; x in radians.");

MPFR_UNIOP(csc)

PyDoc_STRVAR(doc_g_mpfr_cot,
"cot(x) -> mpfr\n\n"
"Return cotangent of x; x in radians.");

MPFR_UNIOP(cot)

static PyObject *
Pympfr_acos(PyObject* self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("acos() requires 'mpfr' argument");

    if (!mpfr_nan_p(MPFR(self)) &&
            (mpfr_cmp_si(MPFR(self), 1) > 0 ||
            mpfr_cmp_si(MPFR(self), -1) < 0) &&
            context->ctx.allow_complex) {
        Py_DECREF(self);
        return Pympc_acos(self, other);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_clear_flags();
    result->rc = mpfr_acos(result->f, MPFR(self),
                           context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF("acos()");
}

static PyObject *
Pympfr_asin(PyObject* self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("asin() requires 'mpfr' argument");

    if (!mpfr_nan_p(MPFR(self)) &&
            (mpfr_cmp_si(MPFR(self), 1) > 0 ||
            mpfr_cmp_si(MPFR(self), -1) < 0) &&
            context->ctx.allow_complex) {
        Py_DECREF(self);
        return Pympc_asin(self, other);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_clear_flags();
    result->rc = mpfr_asin(result->f, MPFR(self),
                           context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF("asin()");
}

PyDoc_STRVAR(doc_g_mpfr_sech,
"sech(x) -> mpfr\n\n"
"Returns hyperbolic secant of x.");

MPFR_UNIOP(sech)

PyDoc_STRVAR(doc_g_mpfr_csch,
"csch(x) -> mpfr\n\n"
"Return hyperbolic cosecant of x.");

MPFR_UNIOP(csch)

PyDoc_STRVAR(doc_g_mpfr_coth,
"coth(x) -> mpfr\n\n"
"Return hyperbolic cotangent of x.");

MPFR_UNIOP(coth)

static PyObject *
Pympfr_atanh(PyObject* self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("atanh() requires 'mpfr' argument");

    if (!mpfr_nan_p(MPFR(self)) &&
            (mpfr_cmp_si(MPFR(self), 1) > 0 ||
            mpfr_cmp_si(MPFR(self), -1) < 0) &&
            context->ctx.allow_complex) {
        Py_DECREF(self);
        return Pympc_atanh(self, other);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_clear_flags();
    result->rc = mpfr_atanh(result->f, MPFR(self),
                           context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF("atanh()");
}


PyDoc_STRVAR(doc_g_mpfr_atan2,
"atan2(y, x) -> mpfr\n\n"
"Return arc-tangent of (y/x).");

static PyObject *
Pympfr_atan2(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "atan2() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_atan2(result->f, MPFR(self),
                            MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("atan2()");
}

PyDoc_STRVAR(doc_g_mpfr_hypot,
"hypot(y, x) -> mpfr\n\n"
"Return square root of (x**2 + y**2).");

static PyObject *
Pympfr_hypot(PyObject *self, PyObject *args)
{
    MPFR_Object *result;
    PyObject *other;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "hypot() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPFR_New(0, context)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_hypot(result->f, MPFR(self),
                            MPFR(other), context->ctx.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("hypot()");
}

static PyObject *
Pympfr_sin_cos(PyObject *self, PyObject *other)
{
    MPFR_Object *s, *c;
    PyObject *result;
    int code;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("sin_cos() requires 'mpfr' argument");

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result)
        goto done;

    mpfr_clear_flags();
    code = mpfr_sin_cos(s->f, c->f, MPFR(self),
                        context->ctx.mpfr_round);
    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;
    SUBNORMALIZE(s);
    SUBNORMALIZE(c);
    MERGE_FLAGS;
    CHECK_FLAGS("sin_cos()");

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

PyDoc_STRVAR(doc_g_mpfr_sinh_cosh,
"sinh_cosh(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

static PyObject *
Pympfr_sinh_cosh(PyObject *self, PyObject *other)
{
    MPFR_Object *s, *c;
    PyObject *result;
    int code;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("sinh_cosh() requires 'mpfr' argument");

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result)
        goto done;

    mpfr_clear_flags();
    code = mpfr_sinh_cosh(s->f, c->f, MPFR(self),
                          context->ctx.mpfr_round);
    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;
    SUBNORMALIZE(s);
    SUBNORMALIZE(c);
    MERGE_FLAGS;
    CHECK_FLAGS("sin_cos()");

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

PyDoc_STRVAR(doc_g_mpfr_degrees,
"degrees(x) -> mpfr\n\n"
"Convert angle x from radians to degrees.");

static PyObject *
Pympfr_degrees(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *temp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("degrees() requires 'mpfr' argument");

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(context->ctx.mpfr_prec + 20, context);
    if (!result || !temp) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_clear_flags();
    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_ui_div(temp->f, 180, temp->f, MPFR_RNDN);
    mpfr_mul(result->f, temp->f, MPFR(self), MPFR_RNDN);
    Py_DECREF((PyObject*)temp);
    MPFR_CLEANUP_SELF("degrees()");
}

PyDoc_STRVAR(doc_g_mpfr_radians,
"radians(x) -> mpfr\n\n"
"Convert angle x from degrees to radians.");

static PyObject *
Pympfr_radians(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *temp;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPFR_OTHER("radians() requires 'mpfr' argument");

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(context->ctx.mpfr_prec + 20, context);
    if (!result || !temp) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_clear_flags();
    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_div_ui(temp->f, temp->f, 180, MPFR_RNDN);
    mpfr_mul(result->f, MPFR(self), temp->f, MPFR_RNDN);
    Py_DECREF((PyObject*)temp);
    MPFR_CLEANUP_SELF("radians()");
}

