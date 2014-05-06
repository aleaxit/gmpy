/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_math.c                                                            *
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
"Return inverse tangent of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_atan,
"atan(x) -> number\n\n"
"Return inverse tangent of x; result in radians.");

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

PyDoc_STRVAR(GMPy_doc_context_sec,
"context.sec(x) -> number\n\n"
"Return secant of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sec,
"sec(x) -> number\n\n"
"Return secant of x; x in radians.");

GMPY_MPFR_UNIOP(Sec, sec)

PyDoc_STRVAR(GMPy_doc_context_csc,
"context.csc(x) -> number\n\n"
"Return cosecant of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_csc,
"csc(x) -> number\n\n"
"Return cosecant of x; x in radians.");

GMPY_MPFR_UNIOP(Csc, csc)

PyDoc_STRVAR(GMPy_doc_context_cot,
"context.cot(x) -> number\n\n"
"Return cotangent of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_cot,
"cot(x) -> number\n\n"
"Return cotangent of x; x in radians.");

GMPY_MPFR_UNIOP(Cot, cot)

PyDoc_STRVAR(GMPy_doc_context_sech,
"context.sech(x) -> number\n\n"
"Return hyperbolic secant of x.");

PyDoc_STRVAR(GMPy_doc_function_sech,
"sech(x) -> number\n\n"
"Return hyperbolic secant of x.");

GMPY_MPFR_UNIOP(Sech, sech)

PyDoc_STRVAR(GMPy_doc_context_csch,
"context.csch(x) -> number\n\n"
"Return hyperbolic cosecant of x.");

PyDoc_STRVAR(GMPy_doc_function_csch,
"csch(x) -> number\n\n"
"Return hyperbolic cosecant of x.");

GMPY_MPFR_UNIOP(Csch, csch)

PyDoc_STRVAR(GMPy_doc_context_coth,
"context.coth(x) -> number\n\n"
"Return hyperbolic cotangent of x.");

PyDoc_STRVAR(GMPy_doc_function_coth,
"coth(x) -> number\n\n"
"Return hyperbolic cotangent of x.");

GMPY_MPFR_UNIOP(Coth, coth)

PyDoc_STRVAR(GMPy_doc_context_acos,
"context.acos(x) -> number\n\n"
"Return inverse cosine of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_acos,
"acos(x) -> number\n\n"
"Return inverse cosine of x; result in radians.");

static PyObject *
GMPy_Real_Acos(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)))
        return NULL;

    if (!mpfr_nan_p(tempx->f) &&
            (mpfr_cmp_si(tempx->f, 1) > 0 || mpfr_cmp_si(tempx->f, -1) < 0) &&
            context->ctx.allow_complex
       ) {
        Py_DECREF((PyObject*)tempx);
        return GMPy_Complex_Acos(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_acos(result->f, tempx->f, GET_MPFR_ROUND(context));
    GMPY_MPFR_CLEANUP(result, context, "acos()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Acos(PyObject *x, CTXT_Object *context)
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

    result->rc = mpc_acos(result->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPC_CLEANUP(result, context, "acos()");
    return (PyObject*)result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE(Acos, acos)

PyDoc_STRVAR(GMPy_doc_context_asin,
"context.asin(x) -> number\n\n"
"Return inverse sine of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_asin,
"asin(x) -> number\n\n"
"Return inverse sine of x; result in radians.");

static PyObject *
GMPy_Real_Asin(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)))
        return NULL;

    if (!mpfr_nan_p(tempx->f) &&
            (mpfr_cmp_si(tempx->f, 1) > 0 || mpfr_cmp_si(tempx->f, -1) < 0) &&
            context->ctx.allow_complex
       ) {
        Py_DECREF((PyObject*)tempx);
        return GMPy_Complex_Asin(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_asin(result->f, tempx->f, GET_MPFR_ROUND(context));
    GMPY_MPFR_CLEANUP(result, context, "asin()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Asin(PyObject *x, CTXT_Object *context)
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

    result->rc = mpc_asin(result->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPC_CLEANUP(result, context, "asin()");
    return (PyObject*)result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE(Asin, asin)

PyDoc_STRVAR(GMPy_doc_context_atanh,
"context.atanh(x) -> number\n\n"
"Return inverse hyperbolic tanget of x.");

PyDoc_STRVAR(GMPy_doc_function_atanh,
"atanh(x) -> number\n\n"
"Return inverse hyperbolic tangent of x.");

static PyObject *
GMPy_Real_Atanh(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)))
        return NULL;

    if (!mpfr_nan_p(tempx->f) &&
            (mpfr_cmp_si(tempx->f, 1) > 0 || mpfr_cmp_si(tempx->f, -1) < 0) &&
            context->ctx.allow_complex
       ) {
        Py_DECREF((PyObject*)tempx);
        return GMPy_Complex_Atanh(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_atanh(result->f, tempx->f, GET_MPFR_ROUND(context));
    GMPY_MPFR_CLEANUP(result, context, "atanh()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Atanh(PyObject *x, CTXT_Object *context)
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

    result->rc = mpc_atanh(result->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPC_CLEANUP(result, context, "atanh()");
    return (PyObject*)result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE(Atanh, atanh)

PyDoc_STRVAR(GMPy_doc_function_atan2,
"atan2(y, x) -> number\n\n"
"Return arc-tangent of (y/x); result in radians.");

PyDoc_STRVAR(GMPy_doc_context_atan2,
"context.atan2(y, x) -> number\n\n"
"Return arc-tangent of (y/x); result in radians.");

GMPY_MPFR_BINOP(Atan2, atan2)

PyDoc_STRVAR(GMPy_doc_function_hypot,
"hypot(x, y) -> number\n\n"
"Return square root of (x**2 + y**2).");

PyDoc_STRVAR(GMPy_doc_context_hypot,
"context.hypot(x, y) -> number\n\n"
"Return square root of (x**2 + y**2).");

GMPY_MPFR_BINOP(Hypot, hypot)

static PyObject *
GMPy_Real_Sin_Cos(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *s, *c, *tempx;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPFR_From_Real(x, 1, context);
    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!tempx || !s || !c || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();
    code = mpfr_sin_cos(s->f, c->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    GMPY_MPFR_CLEANUP(s, context, "sin_cos()");
    GMPY_MPFR_CLEANUP(c, context, "sin_cos()");

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

static PyObject *
GMPy_Complex_Sin_Cos(PyObject *x, CTXT_Object *context)
{
    MPC_Object *s, *c, *tempx;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    s = GMPy_MPC_New(0, 0, context);
    c = GMPy_MPC_New(0, 0, context);
    result = PyTuple_New(2);
    if (!tempx || !s || !c || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    code = mpc_sin_cos(s->c, c->c, tempx->c, GET_MPC_ROUND(context), GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    s->rc = MPC_INEX1(code);
    c->rc = MPC_INEX2(code);

    GMPY_MPC_CLEANUP(s, context, "sin_cos()");
    GMPY_MPC_CLEANUP(c, context, "sin_cos()");

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_sin_cos,
"context.sin_cos(x) -> (number, number)\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sin_cos,
"sin_cos(x) -> (number, number)\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_TEMPLATE(Sin_Cos, sin_cos)

static PyObject *
GMPy_Real_Sinh_Cosh(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *s, *c, *tempx;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPFR_From_Real(x, 1, context);
    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!tempx || !s || !c || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();
    code = mpfr_sinh_cosh(s->f, c->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    GMPY_MPFR_CLEANUP(s, context, "sinh_cosh()");
    GMPY_MPFR_CLEANUP(c, context, "sinh_cosh()");

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_sinh_cosh,
"context.sinh_cosh(x) -> (number, number)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_sinh_cosh,
"sinh_cosh(x) -> (number, number)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

GMPY_MPFR_UNIOP_TEMPLATE(Sinh_Cosh, sinh_cosh)

PyDoc_STRVAR(GMPy_doc_function_degrees,
"degrees(x) -> mpfr\n\n"
"Convert angle x from radians to degrees.\n"
"Note: In rare cases the result may not be correctly rounded.");

PyDoc_STRVAR(GMPy_doc_context_degrees,
"context.degrees(x) -> mpfr\n\n"
"Convert angle x from radians to degrees.\n"
"Note: In rare cases the result may not be correctly rounded.");

static PyObject *
GMPy_Context_Degrees(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx, *temp;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(context->ctx.mpfr_prec + 100, context);
    tempx = GMPy_MPFR_From_Real(other, 1, context);
    if (!result || !temp || !tempx) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_ui_div(temp->f, 180, temp->f, MPFR_RNDN);

    mpfr_clear_flags();
    mpfr_mul(result->f, temp->f, tempx->f, MPFR_RNDN);

    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)tempx);
    GMPY_MPFR_CLEANUP(result, context, "degrees()");
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_radians,
"radians(x) -> mpfr\n\n"
"Convert angle x from degrees to radians.\n"
"Note: In rare cases the result may not be correctly rounded.");

PyDoc_STRVAR(GMPy_doc_context_radians,
"context.radians(x) -> mpfr\n\n"
"Convert angle x from degrees to radians.\n"
"Note: In rare cases the result may not be correctly rounded.");

static PyObject *
GMPy_Context_Radians(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx, *temp;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    result = GMPy_MPFR_New(0, context);
    temp = GMPy_MPFR_New(context->ctx.mpfr_prec + 100, context);
    tempx = GMPy_MPFR_From_Real(other, 1, context);
    if (!result || !temp || !tempx) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_div_ui(temp->f, temp->f, 180, MPFR_RNDN);
    mpfr_mul(result->f, MPFR(self), temp->f, MPFR_RNDN);

    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)tempx);
    GMPY_MPFR_CLEANUP(result, context, "radians()");
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_log10,
"context.log10(x) -> number\n\n"
"Return the base-10 logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log10,
"log10(x) -> number\n\n"
"Return the base-10 logarithm of x.");

GMPY_MPFR_MPC_UNIOP(Log10, log10)

PyDoc_STRVAR(GMPy_doc_context_log,
"context.log(x) -> number\n\n"
"Return the natural logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log,
"log(x) -> number\n\n"
"Return the natural logarithm of x.");

GMPY_MPFR_MPC_UNIOP(Log, log)

PyDoc_STRVAR(GMPy_doc_context_exp,
"context.exp(x) -> number\n\n"
"Return the exponential of x.");

PyDoc_STRVAR(GMPy_doc_function_exp,
"exp(x) -> number\n\n"
"Return the exponential of x.");

GMPY_MPFR_MPC_UNIOP(Exp, exp)

PyDoc_STRVAR(GMPy_doc_context_sqrt,
"context.sqrt(x) -> number\n\n"
"Return the square root of x.");

PyDoc_STRVAR(GMPy_doc_function_sqrt,
"sqrt(x) -> number\n\n"
"Return the square root of x.");

static PyObject *
GMPy_Real_Sqrt(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context)))
        return NULL;

    if (mpfr_sgn(tempx->f) < 0 && context->ctx.allow_complex) {
        Py_DECREF((PyObject*)tempx);
        return GMPy_Complex_Sqrt(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_sqrt(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPFR_CLEANUP(result, context, "sqrt()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Sqrt(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL, *tempx = NULL;

    CHECK_CONTEXT(context);

    if (!(tempx = GMPy_MPC_From_Complex(x, 1, 1, context)))
        return NULL;

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    result->rc = mpc_sqrt(result->c, tempx->c, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    GMPY_MPC_CLEANUP(result, context, "sqrt()");
    return (PyObject*)result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE(Sqrt, sqrt)


