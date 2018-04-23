/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_math.c                                                            *
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

PyDoc_STRVAR(GMPy_doc_context_sin,
"context.sin(x) -> number\n\n"
"Return sine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sin,
"sin(x) -> number\n\n"
"Return sine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_EX(Sin, sin)

PyDoc_STRVAR(GMPy_doc_context_cos,
"context.cos(x) -> number\n\n"
"Return cosine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_cos,
"cos(x) -> number\n\n"
"Return cosine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_EX(Cos, cos)

PyDoc_STRVAR(GMPy_doc_context_tan,
"context.tan(x) -> number\n\n"
"Return tangent of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_tan,
"tan(x) -> number\n\n"
"Return tangent of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_EX(Tan, tan)

PyDoc_STRVAR(GMPy_doc_context_atan,
"context.atan(x) -> number\n\n"
"Return inverse tangent of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_atan,
"atan(x) -> number\n\n"
"Return inverse tangent of x; result in radians.");

GMPY_MPFR_MPC_UNIOP_EX(Atan, atan)

PyDoc_STRVAR(GMPy_doc_context_sinh,
"context.sinh(x) -> number\n\n"
"Return hyperbolic sine of x.");

PyDoc_STRVAR(GMPy_doc_function_sinh,
"sinh(x) -> number\n\n"
"Return hyperbolic sine of x.");

GMPY_MPFR_MPC_UNIOP_EX(Sinh, sinh)

PyDoc_STRVAR(GMPy_doc_context_cosh,
"context.cosh(x) -> number\n\n"
"Return hyperbolic cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_cosh,
"cosh(x) -> number\n\n"
"Return hyperbolic cosine of x.");

GMPY_MPFR_MPC_UNIOP_EX(Cosh, cosh)

PyDoc_STRVAR(GMPy_doc_context_tanh,
"context.tanh(x) -> number\n\n"
"Return hyperbolic tangent of x.");

PyDoc_STRVAR(GMPy_doc_function_tanh,
"tanh(x) -> number\n\n"
"Return hyperbolic tangent of x.");

GMPY_MPFR_MPC_UNIOP_EX(Tanh, tanh)

PyDoc_STRVAR(GMPy_doc_context_asinh,
"context.asinh(x) -> number\n\n"
"Return inverse hyperbolic sine of x.");

PyDoc_STRVAR(GMPy_doc_function_asinh,
"asinh(x) -> number\n\n"
"Return inverse hyperbolic sine of x.");

GMPY_MPFR_MPC_UNIOP_EX(Asinh, asinh)

PyDoc_STRVAR(GMPy_doc_context_acosh,
"context.acosh(x) -> number\n\n"
"Return inverse hyperbolic cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_acosh,
"acosh(x) -> number\n\n"
"Return inverse hyperbolic cosine of x.");

GMPY_MPFR_MPC_UNIOP_EX(Acosh, acosh)

/* Section 2:
 * These functions accept a single argument and return an mpfr result.
 *
 * GMPY_MPFR_UNIOP(NAME, FUNC) creates the following functions:
 *     GMPy_Real_NAME(x, context)
 *     GMPy_Number_NAME(x, context)
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 */

PyDoc_STRVAR(GMPy_doc_context_sec,
"context.sec(x) -> number\n\n"
"Return secant of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sec,
"sec(x) -> number\n\n"
"Return secant of x; x in radians.");

GMPY_MPFR_UNIOP_EX(Sec, sec)

PyDoc_STRVAR(GMPy_doc_context_csc,
"context.csc(x) -> number\n\n"
"Return cosecant of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_csc,
"csc(x) -> number\n\n"
"Return cosecant of x; x in radians.");

GMPY_MPFR_UNIOP_EX(Csc, csc)

PyDoc_STRVAR(GMPy_doc_context_cot,
"context.cot(x) -> number\n\n"
"Return cotangent of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_cot,
"cot(x) -> number\n\n"
"Return cotangent of x; x in radians.");

GMPY_MPFR_UNIOP_EX(Cot, cot)

PyDoc_STRVAR(GMPy_doc_context_sech,
"context.sech(x) -> number\n\n"
"Return hyperbolic secant of x.");

PyDoc_STRVAR(GMPy_doc_function_sech,
"sech(x) -> number\n\n"
"Return hyperbolic secant of x.");

GMPY_MPFR_UNIOP_EX(Sech, sech)

PyDoc_STRVAR(GMPy_doc_context_csch,
"context.csch(x) -> number\n\n"
"Return hyperbolic cosecant of x.");

PyDoc_STRVAR(GMPy_doc_function_csch,
"csch(x) -> number\n\n"
"Return hyperbolic cosecant of x.");

GMPY_MPFR_UNIOP_EX(Csch, csch)

PyDoc_STRVAR(GMPy_doc_context_coth,
"context.coth(x) -> number\n\n"
"Return hyperbolic cotangent of x.");

PyDoc_STRVAR(GMPy_doc_function_coth,
"coth(x) -> number\n\n"
"Return hyperbolic cotangent of x.");

GMPY_MPFR_UNIOP_EX(Coth, coth)

PyDoc_STRVAR(GMPy_doc_context_rec_sqrt,
"context.rec_sqrt(x) -> number\n\n"
"Return the reciprocal of the square root of x.");

PyDoc_STRVAR(GMPy_doc_function_rec_sqrt,
"rec_sqrt(x) -> number\n\n"
"Return the reciprocal of the square root of x.");

GMPY_MPFR_UNIOP_EX(RecSqrt, rec_sqrt)

PyDoc_STRVAR(GMPy_doc_context_rint,
"context.rint(x) -> number\n\n"
"Return x rounded to the nearest integer using the context rounding\n"
"mode.");

PyDoc_STRVAR(GMPy_doc_function_rint,
"rint(x) -> number\n\n"
"Return x rounded to the nearest integer using the current rounding\n"
"mode.");

GMPY_MPFR_UNIOP_EX(Rint, rint)

PyDoc_STRVAR(GMPy_doc_context_rint_ceil,
"context.rint_ceil(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next higher or equal integer and then, if needed, using the context\n"
"rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_ceil,
"rint_ceil(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next higher or equal integer and then, if needed, using the current\n"
"rounding mode.");

GMPY_MPFR_UNIOP_EX(RintCeil, rint_ceil)

PyDoc_STRVAR(GMPy_doc_context_rint_floor,
"context.rint_floor(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next lower or equal integer and then, if needed, using the context\n"
"rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_floor,
"rint_floor(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next lower or equal integer and then, if needed, using the current\n"
"rounding mode.");

GMPY_MPFR_UNIOP_EX(RintFloor, rint_floor)

PyDoc_STRVAR(GMPy_doc_context_rint_round,
"context.rint_round(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the context rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_round,
"rint_round(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the current rounding mode.");

GMPY_MPFR_UNIOP_EX(RintRound, rint_round)

PyDoc_STRVAR(GMPy_doc_context_rint_trunc,
"context.rint_trunc(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the context rounding mode.");

PyDoc_STRVAR(GMPy_doc_function_rint_trunc,
"rint_trunc(x) -> number\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the current rounding mode.");

GMPY_MPFR_UNIOP_EX(RintTrunc, rint_trunc)

PyDoc_STRVAR(GMPy_doc_context_frac,
"context.frac(x) -> number\n\n"
"Return fractional part of x.");

PyDoc_STRVAR(GMPy_doc_function_frac,
"frac(x) -> number\n\n"
"Return fractional part of x.");

GMPY_MPFR_UNIOP_EX(Frac, frac)

PyDoc_STRVAR(GMPy_doc_context_cbrt,
"context.cbrt(x) -> number\n\n"
"Return the cube root of x.");

PyDoc_STRVAR(GMPy_doc_function_cbrt,
"cbrt(x) -> number\n\n"
"Return the cube root of x.");

GMPY_MPFR_UNIOP_EX(Cbrt, cbrt)

PyDoc_STRVAR(GMPy_doc_context_log2,
"context.log2(x) -> number\n\n"
"Return base-2 logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log2,
"log2(x) -> number\n\n"
"Return base-2 logarithm of x.");

GMPY_MPFR_UNIOP_EX(Log2, log2)

PyDoc_STRVAR(GMPy_doc_context_exp2,
"context.exp2(x) -> number\n\n"
"Return 2**x.");

PyDoc_STRVAR(GMPy_doc_function_exp2,
"exp2(x) -> number\n\n"
"Return 2**x.");

GMPY_MPFR_UNIOP_EX(Exp2, exp2)

PyDoc_STRVAR(GMPy_doc_context_exp10,
"context.exp10(x) -> number\n\n"
"Return 10**x.");

PyDoc_STRVAR(GMPy_doc_function_exp10,
"exp10(x) -> number\n\n"
"Return 10**x.");

GMPY_MPFR_UNIOP_EX(Exp10, exp10)

PyDoc_STRVAR(GMPy_doc_context_log1p,
"context.log1p(x) -> number\n\n"
"Return natural logarithm of (1+x).");

PyDoc_STRVAR(GMPy_doc_function_log1p,
"log1p(x) -> number\n\n"
"Return natural logarithm of (1+x).");

GMPY_MPFR_UNIOP_EX(Log1p, log1p)

PyDoc_STRVAR(GMPy_doc_context_expm1,
"context.expm1(x) -> number\n\n"
"Return exp(x) - 1.");

PyDoc_STRVAR(GMPy_doc_function_expm1,
"expm1(x) -> number\n\n"
"Return exp(x) - 1.");

GMPY_MPFR_UNIOP_EX(Expm1, expm1)

PyDoc_STRVAR(GMPy_doc_context_eint,
"context.eint(x) -> number\n\n"
"Return exponential integral of x.");

PyDoc_STRVAR(GMPy_doc_function_eint,
"eint(x) -> number\n\n"
"Return exponential integral of x.");

GMPY_MPFR_UNIOP_EX(Eint, eint)

PyDoc_STRVAR(GMPy_doc_context_li2,
"context.li2(x) -> number\n\n"
"Return real part of dilogarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_li2,
"li2(x) -> number\n\n"
"Return real part of dilogarithm of x.");

GMPY_MPFR_UNIOP_EX(Li2, li2)

PyDoc_STRVAR(GMPy_doc_context_gamma,
"context.gamma(x) -> number\n\n"
"Return gamma of x.");

PyDoc_STRVAR(GMPy_doc_function_gamma,
"gamma(x) -> number\n\n"
"Return gamma of x.");

GMPY_MPFR_UNIOP_EX(Gamma, gamma)

PyDoc_STRVAR(GMPy_doc_context_lngamma,
"context.lngamma(x) -> number\n\n"
"Return natural logarithm of gamma(x).");

PyDoc_STRVAR(GMPy_doc_function_lngamma,
"lngamma(x) -> number\n\n"
"Return natural logarithm of gamma(x).");

GMPY_MPFR_UNIOP_EX(Lngamma, lngamma)

PyDoc_STRVAR(GMPy_doc_context_digamma,
"context.digamma(x) -> number\n\n"
"Return digamma of x.");

PyDoc_STRVAR(GMPy_doc_function_digamma,
"digamma(x) -> number\n\n"
"Return digamma of x.");

GMPY_MPFR_UNIOP_EX(Digamma, digamma)

PyDoc_STRVAR(GMPy_doc_context_zeta,
"context.zeta(x) -> number\n\n"
"Return Riemann zeta of x.");

PyDoc_STRVAR(GMPy_doc_function_zeta,
"zeta(x) -> number\n\n"
"Return Riemann zeta of x.");

GMPY_MPFR_UNIOP_EX(Zeta, zeta)

PyDoc_STRVAR(GMPy_doc_context_erf,
"context.erf(x) -> number\n\n"
"Return error function of x.");

PyDoc_STRVAR(GMPy_doc_function_erf,
"erf(x) -> number\n\n"
"Return error function of x.");

GMPY_MPFR_UNIOP_EX(Erf, erf)

PyDoc_STRVAR(GMPy_doc_context_erfc,
"context.erfc(x) -> number\n\n"
"Return complementary error function of x.");

PyDoc_STRVAR(GMPy_doc_function_erfc,
"erfc(x) -> number\n\n"
"Return complementary error function of x.");

GMPY_MPFR_UNIOP_EX(Erfc, erfc)

PyDoc_STRVAR(GMPy_doc_context_j0,
"context.j0(x) -> number\n\n"
"Return first kind Bessel function of order 0 of x.");

PyDoc_STRVAR(GMPy_doc_function_j0,
"j0(x) -> number\n\n"
"Return first kind Bessel function of order 0 of x.");

GMPY_MPFR_UNIOP_EX(J0, j0)

PyDoc_STRVAR(GMPy_doc_context_j1,
"context.j1(x) -> number\n\n"
"Return first kind Bessel function of order 1 of x.");

PyDoc_STRVAR(GMPy_doc_function_j1,
"j1(x) -> number\n\n"
"Return first kind Bessel function of order 1 of x.");

GMPY_MPFR_UNIOP_EX(J1, j1)

PyDoc_STRVAR(GMPy_doc_context_y0,
"context.y0(x) -> number\n\n"
"Return second kind Bessel function of order 0 of x.");

PyDoc_STRVAR(GMPy_doc_function_y0,
"y0(x) -> number\n\n"
"Return second kind Bessel function of order 0 of x.");

GMPY_MPFR_UNIOP_EX(Y0, y0)

PyDoc_STRVAR(GMPy_doc_context_y1,
"context.y1(x) -> number\n\n"
"Return second kind Bessel function of order 1 of x.");

PyDoc_STRVAR(GMPy_doc_function_y1,
"y1(x) -> number\n\n"
"Return second kind Bessel function of order 1 of x.");

GMPY_MPFR_UNIOP_EX(Y1, y1)

PyDoc_STRVAR(GMPy_doc_context_ai,
"context.ai(x) -> number\n\n"
"Return Airy function of x.");

PyDoc_STRVAR(GMPy_doc_function_ai,
"ai(x) -> number\n\n"
"Return Airy function of x.");

GMPY_MPFR_UNIOP_EX(Ai, ai)

/* Section 3:
 * The following functions may return an mpc result for certain mpfr arguments.
 * Since the expectional values vary between functions, the 'Real' and 'Complex'
 * functions do not use macros. However, they do use a macro to create the
 * higher-level functions.
 *
 * GMPY_MPFR_MPC_UNIOP_TEMPLATE(NAME, FUNC) creates the following functions:
 *     GMPy_Number_NAME(x, context)
 *     - assumes GMPy_Real_NAME & GMPy_Complex_NAME exist
 *     GMPy_Context_NAME(self, other)
 *     - called with METH_O
 */

PyDoc_STRVAR(GMPy_doc_context_acos,
"context.acos(x) -> number\n\n"
"Return inverse cosine of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_acos,
"acos(x) -> number\n\n"
"Return inverse cosine of x; result in radians.");

static PyObject *
_GMPy_MPFR_Acos(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    if (!mpfr_nan_p(MPFR(x)) &&
            (mpfr_cmp_si(MPFR(x), 1) > 0 || mpfr_cmp_si(MPFR(x), -1) < 0) &&
            context->ctx.allow_complex
       ) {
        return GMPy_Complex_Acos(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_acos(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Acos(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Acos(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPC_Acos(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_acos(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Acos(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Acos(tempx, context);
    Py_DECREF(tempx);
    return result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EX(Acos, acos)

PyDoc_STRVAR(GMPy_doc_context_asin,
"context.asin(x) -> number\n\n"
"Return inverse sine of x; result in radians.");

PyDoc_STRVAR(GMPy_doc_function_asin,
"asin(x) -> number\n\n"
"Return inverse sine of x; result in radians.");

static PyObject *
_GMPy_MPFR_Asin(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!mpfr_nan_p(MPFR(x)) &&
            (mpfr_cmp_si(MPFR(x), 1) > 0 || mpfr_cmp_si(MPFR(x), -1) < 0) &&
            context->ctx.allow_complex
       ) {
        return GMPy_Complex_Asin(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_asin(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Asin(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Asin(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPC_Asin(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_asin(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Asin(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Asin(tempx, context);
    Py_DECREF(tempx);
    return result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EX(Asin, asin)

PyDoc_STRVAR(GMPy_doc_context_atanh,
"context.atanh(x) -> number\n\n"
"Return inverse hyperbolic tanget of x.");

PyDoc_STRVAR(GMPy_doc_function_atanh,
"atanh(x) -> number\n\n"
"Return inverse hyperbolic tangent of x.");

static PyObject *
_GMPy_MPFR_Atanh(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (!mpfr_nan_p(MPFR(x)) &&
            (mpfr_cmp_si(MPFR(x), 1) > 0 || mpfr_cmp_si(MPFR(x), -1) < 0) &&
            context->ctx.allow_complex
       ) {
        return GMPy_Complex_Atanh(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_atanh(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Atanh(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Atanh(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPC_Atanh(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_atanh(result->c, MPC(x), GET_MPC_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Atanh(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Atanh(tempx, context);
    Py_DECREF(tempx);
    return result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EX(Atanh, atanh)

PyDoc_STRVAR(GMPy_doc_function_atan2,
"atan2(y, x) -> number\n\n"
"Return arc-tangent of (y/x); result in radians.");

PyDoc_STRVAR(GMPy_doc_context_atan2,
"context.atan2(y, x) -> number\n\n"
"Return arc-tangent of (y/x); result in radians.");

GMPY_MPFR_BINOP_EX(Atan2, atan2)

PyDoc_STRVAR(GMPy_doc_function_hypot,
"hypot(x, y) -> number\n\n"
"Return square root of (x**2 + y**2).");

PyDoc_STRVAR(GMPy_doc_context_hypot,
"context.hypot(x, y) -> number\n\n"
"Return square root of (x**2 + y**2).");

GMPY_MPFR_BINOP_EX(Hypot, hypot)

static PyObject *
_GMPy_MPFR_Sin_Cos(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *s, *c;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    code = mpfr_sin_cos(s->f, c->f, MPFR(x), GET_MPFR_ROUND(context));

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    _GMPy_MPFR_Cleanup(&s, context);
    _GMPy_MPFR_Cleanup(&c, context);

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

static PyObject *
GMPy_Real_Sin_Cos(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Sin_Cos(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPC_Sin_Cos(PyObject *x, CTXT_Object *context)
{
    MPC_Object *s, *c;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    s = GMPy_MPC_New(0, 0, context);
    c = GMPy_MPC_New(0, 0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    code = mpc_sin_cos(s->c, c->c, MPC(x), GET_MPC_ROUND(context), GET_MPC_ROUND(context));

    s->rc = MPC_INEX1(code);
    c->rc = MPC_INEX2(code);

    _GMPy_MPC_Cleanup(&s, context);
    _GMPy_MPC_Cleanup(&c, context);

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
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Sin_Cos(tempx, context);
    Py_DECREF(tempx);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_sin_cos,
"context.sin_cos(x) -> (number, number)\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.");

PyDoc_STRVAR(GMPy_doc_function_sin_cos,
"sin_cos(x) -> (number, number)\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.");

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EX(Sin_Cos, sin_cos)

static PyObject *
_GMPy_MPFR_Sinh_Cosh(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *s, *c;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    code = mpfr_sinh_cosh(s->f, c->f, MPFR(x), GET_MPFR_ROUND(context));

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    _GMPy_MPFR_Cleanup(&s, context);
    _GMPy_MPFR_Cleanup(&c, context);

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
GMPy_Real_Sinh_Cosh(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Sinh_Cosh(tempx, context);
    Py_DECREF(tempx);
    return result;
}

PyDoc_STRVAR(GMPy_doc_context_sinh_cosh,
"context.sinh_cosh(x) -> (number, number)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

PyDoc_STRVAR(GMPy_doc_function_sinh_cosh,
"sinh_cosh(x) -> (number, number)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

GMPY_MPFR_UNIOP_TEMPLATE_EX(Sinh_Cosh, sinh_cosh)

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
    _GMPy_MPFR_Cleanup(&result, context);
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

    mpfr_clear_flags();

    mpfr_mul(result->f, MPFR(self), temp->f, MPFR_RNDN);

    Py_DECREF((PyObject*)temp);
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_context_log10,
"context.log10(x) -> number\n\n"
"Return the base-10 logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log10,
"log10(x) -> number\n\n"
"Return the base-10 logarithm of x.");

GMPY_MPFR_MPC_UNIOP_EX(Log10, log10)

PyDoc_STRVAR(GMPy_doc_context_log,
"context.log(x) -> number\n\n"
"Return the natural logarithm of x.");

PyDoc_STRVAR(GMPy_doc_function_log,
"log(x) -> number\n\n"
"Return the natural logarithm of x.");

GMPY_MPFR_MPC_UNIOP_EX(Log, log)

PyDoc_STRVAR(GMPy_doc_context_exp,
"context.exp(x) -> number\n\n"
"Return the exponential of x.");

PyDoc_STRVAR(GMPy_doc_function_exp,
"exp(x) -> number\n\n"
"Return the exponential of x.");

GMPY_MPFR_MPC_UNIOP_EX(Exp, exp)

PyDoc_STRVAR(GMPy_doc_context_sqrt,
"context.sqrt(x) -> number\n\n"
"Return the square root of x.");

PyDoc_STRVAR(GMPy_doc_function_sqrt,
"sqrt(x) -> number\n\n"
"Return the square root of x.");

static PyObject *
_GMPy_MPFR_Sqrt(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;

    CHECK_CONTEXT(context);

    if (mpfr_sgn(MPFR(x)) < 0 && context->ctx.allow_complex) {
        return GMPy_Complex_Sqrt(x, context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    result->rc = mpfr_sqrt(result->f, MPFR(x), GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Real_Sqrt(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPFR_Sqrt(tempx, context);
    Py_DECREF(tempx);
    return result;
}

static PyObject *
_GMPy_MPC_Sqrt(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_sqrt(result->c, MPC(x), GET_MPFR_ROUND(context));
    _GMPy_MPC_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Complex_Sqrt(PyObject *x, CTXT_Object *context)
{
    PyObject *result, *tempx;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    result = _GMPy_MPC_Sqrt(tempx, context);
    Py_DECREF(tempx);
    return result;
}

GMPY_MPFR_MPC_UNIOP_TEMPLATE_EX(Sqrt, sqrt)

PyDoc_STRVAR(GMPy_doc_function_root,
"root(x, n) -> mpfr\n\n"
"Return n-th root of x. The result always an 'mpfr'.");

PyDoc_STRVAR(GMPy_doc_context_root,
"context.root(x, n) -> mpfr\n\n"
"Return n-th root of x. The result always an 'mpfr'.");

GMPY_MPFR_BINOP_REAL_ULONG(Root, root)

PyDoc_STRVAR(GMPy_doc_function_jn,
"jn(x,n) -> mpfr\n\n"
"Return the first kind Bessel function of order n of x.");

PyDoc_STRVAR(GMPy_doc_context_jn,
"context.jn(x,n) -> mpfr\n\n"
"Return the first kind Bessel function of order n of x.");

GMPY_MPFR_BINOP_REAL_LONG(Jn, jn)

PyDoc_STRVAR(GMPy_doc_function_yn,
"yn(x,n) -> mpfr\n\n"
"Return the second kind Bessel function of order n of x.");

PyDoc_STRVAR(GMPy_doc_context_yn,
"context.yn(x,n) -> mpfr\n\n"
"Return the second kind Bessel function of order n of x.");

GMPY_MPFR_BINOP_REAL_LONG(Yn, yn)

PyDoc_STRVAR(GMPy_doc_function_agm,
"agm(x, y) -> mpfr\n\n"
"Return arithmetic-geometric mean of x and y.");

PyDoc_STRVAR(GMPy_doc_context_agm,
"context.agm(x, y) -> mpfr\n\n"
"Return arithmetic-geometric mean of x and y.");

GMPY_MPFR_BINOP(AGM, agm)

PyDoc_STRVAR(GMPy_doc_function_maxnum,
"maxnum(x, y) -> mpfr\n\n"
"Return the maximum number of x and y. If x and y are not 'mpfr', they are\n"
"converted to 'mpfr'. The result is rounded to match the current context.\n"
"If only one of x or y is a number, then that number is returned.");

PyDoc_STRVAR(GMPy_doc_context_maxnum,
"context.maxnum(x, y) -> mpfr\n\n"
"Return the maximum number of x and y. If x and y are not 'mpfr', they are\n"
"converted to 'mpfr'. The result is rounded to match the specified context.\n"
"If only one of x or y is a number, then that number is returned.");

GMPY_MPFR_BINOP(Maxnum, max)

PyDoc_STRVAR(GMPy_doc_function_minnum,
"minnum(x, y) -> mpfr\n\n"
"Return the minimum number of x and y. If x and y are not 'mpfr', they are\n"
"converted to 'mpfr'. The result is rounded to match the current context.\n"
"If only one of x or y is a number, then that number is returned.");

PyDoc_STRVAR(GMPy_doc_context_minnum,
"context.minnum(x, y) -> mpfr\n\n"
"Return the minimum number of x and y. If x and y are not 'mpfr', they are\n"
"converted to 'mpfr'. The result is rounded to match the specified context.\n"
"If only one of x or y is a number, then that number is returned.");

GMPY_MPFR_BINOP(Minnum, min)

PyDoc_STRVAR(GMPy_doc_function_remainder,
"remainder(x, y) -> mpfr\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

PyDoc_STRVAR(GMPy_doc_context_remainder,
"context.remainder(x, y) -> mpfr\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

GMPY_MPFR_BINOP(Remainder, remainder)

PyDoc_STRVAR(GMPy_doc_function_fmod,
"fmod(x, y) -> mpfr\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to 0.");

PyDoc_STRVAR(GMPy_doc_context_fmod,
"context.fmod(x, y) -> mpfr\n\n"
"Return x - n*y where n is the integer quotient of x/y, rounded to 0.");

GMPY_MPFR_BINOP(Fmod, fmod)

PyDoc_STRVAR(GMPy_doc_function_round2,
"round2(x[, n]) -> mpfr\n\n"
"Return x rounded to n bits. Uses default precision if n is not\n"
"specified. See round_away() to access the mpfr_round() function.");

PyDoc_STRVAR(GMPy_doc_context_round2,
"context.round2(x[, n]) -> mpfr\n\n"
"Return x rounded to n bits. Uses default precision if n is not\n"
"specified. See round_away() to access the mpfr_round() function.");

static PyObject *
GMPy_Real_Round2(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *result, *tempx;
    long n = 0;

    CHECK_CONTEXT(context);

    if (y) {
        n = PyIntOrLong_AsLong(y);
        if ( (n == -1 && PyErr_Occurred()) || n < MPFR_PREC_MIN || n > MPFR_PREC_MAX) {
            VALUE_ERROR("invalid precision");
            return NULL;
        }
    }

    if (!(tempx = GMPy_MPFR_From_Real(x, 1, context))) {
        return NULL;
    }
    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    mpfr_clear_flags();

    result->rc = mpfr_prec_round(result->f, n, GET_MPFR_ROUND(context));
    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_Round2(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && (!y || PyIntOrLong_Check(y)))
        return GMPy_Real_Round2(x, y, context);

    TYPE_ERROR("round2() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Round2(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) < 1 || PyTuple_GET_SIZE(args) > 2) {
        TYPE_ERROR("round2() requires 1 or 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        return GMPy_Number_Round2(PyTuple_GET_ITEM(args, 0), NULL, context);
    }
    else {
        return GMPy_Number_Round2(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
    }
}

PyDoc_STRVAR(GMPy_doc_function_reldiff,
"reldiff(x, y) -> mpfr\n\n"
"Return the relative difference between x and y. Result is equal to\n"
"abs(x-y)/x.");

PyDoc_STRVAR(GMPy_doc_context_reldiff,
"context.reldiff(x, y) -> mpfr\n\n"
"Return the relative difference between x and y. Result is equal to\n"
"abs(x-y)/x.");

static PyObject *
GMPy_Real_RelDiff(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx, *tempy, *result;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    tempy = GMPy_MPFR_From_Real(y, 1, context);
    if (!result || !tempx || !tempy) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_reldiff(result->f, tempx->f, tempy->f, GET_MPFR_ROUND(context));
    result->rc = 0;
    _GMPy_MPFR_Cleanup(&result, context);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)result;
}

GMPY_MPFR_BINOP_TEMPLATE(RelDiff, reldiff)

PyDoc_STRVAR(GMPy_doc_mpfr_ceil_method,
"x.__ceil__() -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

PyDoc_STRVAR(GMPy_doc_function_ceil,
"ceil(x) ->mpfr\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

PyDoc_STRVAR(GMPy_doc_context_ceil,
"context.ceil(x) ->mpfr\n\n"
"Return an 'mpfr' that is the smallest integer >= x.");

GMPY_MPFR_UNIOP_NOROUND(Ceil, ceil)

PyDoc_STRVAR(GMPy_doc_mpfr_floor_method,
"x.__floor__() -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer <= x.");

PyDoc_STRVAR(GMPy_doc_function_floor,
"floor(x) -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer <= x.");

PyDoc_STRVAR(GMPy_doc_context_floor,
"context.floor(x) -> mpfr\n\n"
"Return an 'mpfr' that is the smallest integer <= x.");

GMPY_MPFR_UNIOP_NOROUND(Floor, floor);

PyDoc_STRVAR(GMPy_doc_mpfr_trunc_method,
"x.__trunc__() -> mpfr\n\n"
"Return an 'mpfr' that is truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

PyDoc_STRVAR(GMPy_doc_function_trunc,
"trunc(x) -> mpfr\n\n"
"Return an 'mpfr' that is x truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

PyDoc_STRVAR(GMPy_doc_context_trunc,
"context.trunc(x) -> mpfr\n\n"
"Return an 'mpfr' that is x truncated towards 0. Same as\n"
"x.floor() if x>=0 or x.ceil() if x<0.");

GMPY_MPFR_UNIOP_NOROUND(Trunc, trunc)

PyDoc_STRVAR(GMPy_doc_function_round_away,
"round_away(x) -> mpfr\n\n"
"Return an 'mpfr' that is x rounded to the nearest integer,\n"
"with ties rounded away from 0.");

PyDoc_STRVAR(GMPy_doc_context_round_away,
"context.round_away(x) -> mpfr\n\n"
"Return an 'mpfr' that is x rounded to the nearest integer,\n"
"with ties rounded away from 0.");

GMPY_MPFR_UNIOP_NOROUND_NOMETHOD(RoundAway, round)

PyDoc_STRVAR(GMPy_doc_function_modf,
"modf(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the integer and fractional portions\n"
"of x.");

PyDoc_STRVAR(GMPy_doc_context_modf,
"context.modf(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the integer and fractional portions\n"
"of x.");

static PyObject *
GMPy_Real_Modf(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *s, *c, *tempx;
    PyObject *result;
    int code;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPFR_From_Real(x, 1, context);
    s = GMPy_MPFR_New(0, context);
    c = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (! tempx || !s || !c || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    code = mpfr_modf(s->f, c->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;

    _GMPy_MPFR_Cleanup(&s, context);
    _GMPy_MPFR_Cleanup(&c, context);

    if (!s || !c) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    return result;
}

GMPY_MPFR_UNIOP_TEMPLATE(Modf, modf)

PyDoc_STRVAR(GMPy_doc_function_lgamma,
"lgamma(x) -> (mpfr, int)\n\n"
"Return a tuple containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

PyDoc_STRVAR(GMPy_doc_context_lgamma,
"context.lgamma(x) -> (mpfr, int)\n\n"
"Return a tuple containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

static PyObject *
GMPy_Real_Lgamma(PyObject *x, CTXT_Object *context)
{
    PyObject *result;
    MPFR_Object *value, *tempx;
    int signp = 0;

    CHECK_CONTEXT(context)

    tempx = GMPy_MPFR_From_Real(x, 1, context);
    value = GMPy_MPFR_New(0, context);
    result = PyTuple_New(2);
    if (!tempx || !value || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)value);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    value->rc = mpfr_lgamma(value->f, &signp, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    _GMPy_MPFR_Cleanup(&value, context);

    if (!value) {
        Py_DECREF(result);
        return NULL;
    }

    PyTuple_SET_ITEM(result, 0, (PyObject*)value);
    PyTuple_SET_ITEM(result, 1, PyIntOrLong_FromLong((long)signp));
    return result;
}

GMPY_MPFR_UNIOP_TEMPLATE(Lgamma, lgamma)

PyDoc_STRVAR(GMPy_doc_function_remquo,
"remquo(x, y) -> (mpfr, int)\n\n"
"Return a tuple containing the remainder(x,y) and the low bits of the\n"
"quotient.");

PyDoc_STRVAR(GMPy_doc_context_remquo,
"context.remquo(x, y) -> (mpfr, int)\n\n"
"Return a tuple containing the remainder(x,y) and the low bits of the\n"
"quotient.");

static PyObject *
GMPy_Real_RemQuo(PyObject *x, PyObject *y, CTXT_Object *context)
{
    PyObject *result;
    MPFR_Object *value, *tempx, *tempy;
    long quobits = 0;

    CHECK_CONTEXT(context);

    value = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    tempy = GMPy_MPFR_From_Real(y, 1, context);
    result = PyTuple_New(2);
    if (!value || !tempx || !tempx || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)value);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    value->rc = mpfr_remquo(value->f, &quobits, tempx->f, tempy->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    _GMPy_MPFR_Cleanup(&value, context);

    PyTuple_SET_ITEM(result, 0, (PyObject*)value);
    PyTuple_SET_ITEM(result, 1, PyIntOrLong_FromLong(quobits));
    return result;
}

GMPY_MPFR_BINOP_TEMPLATE(RemQuo, remquo);

PyDoc_STRVAR(GMPy_doc_function_frexp,
"frexp(x) -> (int, mpfr)\n\n"
"Return a tuple containing the exponent and mantissa of x.");

PyDoc_STRVAR(GMPy_doc_context_frexp,
"context.frexp(x) -> (int, mpfr)\n\n"
"Return a tuple containing the exponent and mantissa of x.");

static PyObject *
GMPy_Real_Frexp(PyObject *x, CTXT_Object *context)
{
    PyObject *result;
    MPFR_Object *value, *tempx;
    mpfr_exp_t exp = 0;

    CHECK_CONTEXT(context);

    value = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPFR_From_Real(x, 1, context);
    result = PyTuple_New(2);
    if (!value || !result || !tempx) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)value);
        Py_XDECREF(result);
        return NULL;
    }

    mpfr_clear_flags();

    value->rc = mpfr_frexp(&exp, value->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    _GMPy_MPFR_Cleanup(&value, context);

    PyTuple_SET_ITEM(result, 0, PyIntOrLong_FromSsize_t((Py_ssize_t)exp));
    PyTuple_SET_ITEM(result, 1, (PyObject*)value);
    return result;
}

GMPY_MPFR_UNIOP_TEMPLATE(Frexp, frexp)

PyDoc_STRVAR(GMPy_doc_function_next_toward,
"next_toward(x, y) -> mpfr\n\n"
"Return the next 'mpfr' from x in the direction of y. The result has\n"
"the same precision as x.");

PyDoc_STRVAR(GMPy_doc_context_next_toward,
"context.next_toward(x, y) -> mpfr\n\n"
"Return the next 'mpfr' from x in the direction of y. The result has\n"
"the same precision as x.");

static PyObject *
GMPy_Context_NextToward(PyObject *self, PyObject *args)
{
    MPFR_Object *result, *tempx, *tempy;
    CTXT_Object *context = NULL;
    int direction;
    mpfr_rnd_t temp_round;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("next_toward() requires 2 arguments");
        return NULL;
    }

    tempx = GMPy_MPFR_From_Real(PyTuple_GET_ITEM(args, 0), 1, context);
    tempy = GMPy_MPFR_From_Real(PyTuple_GET_ITEM(args, 1), 1, context);
    if (!tempx || !tempy) {
        TYPE_ERROR("next_toward() argument type not supported");
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    mpfr_nexttoward(result->f, tempy->f);
    result->rc = 0;
    direction = mpfr_signbit(tempy->f);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    temp_round = GET_MPFR_ROUND(context);
    if (direction)
        context->ctx.mpfr_round = MPFR_RNDD;
    else
         context->ctx.mpfr_round = MPFR_RNDU;
    _GMPy_MPFR_Cleanup(&result, context);
    context->ctx.mpfr_round = temp_round;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_next_above,
"next_above(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward +Infinity.");

PyDoc_STRVAR(GMPy_doc_context_next_above,
"context.next_above(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward +Infinity.");

static PyObject *
GMPy_Context_NextAbove(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx;
    CTXT_Object *context = NULL;
    mpfr_rnd_t temp_round;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(tempx = GMPy_MPFR_From_Real(other, 1, context))) {
        TYPE_ERROR("next_above() argument type not supported");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    mpfr_nextabove(result->f);
    result->rc = 0;
    temp_round = GET_MPFR_ROUND(context);
    context->ctx.mpfr_round = MPFR_RNDU;
    _GMPy_MPFR_Cleanup(&result, context);
    context->ctx.mpfr_round = temp_round;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_next_below,
"next_below(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward -Infinity.");

PyDoc_STRVAR(GMPy_doc_context_next_below,
"context.next_below(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward -Infinity.");

static PyObject *
GMPy_Context_NextBelow(PyObject *self, PyObject *other)
{
    MPFR_Object *result, *tempx;
    CTXT_Object *context = NULL;
    mpfr_rnd_t temp_round;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(tempx = GMPy_MPFR_From_Real(other, 1, context))) {
        TYPE_ERROR("next_below() argument type not supported");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(mpfr_get_prec(tempx->f), context))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_set(result->f, tempx->f, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    mpfr_nextbelow(result->f);
    result->rc = 0;
    temp_round = GET_MPFR_ROUND(context);
    context->ctx.mpfr_round = MPFR_RNDD;
    _GMPy_MPFR_Cleanup(&result, context);
    context->ctx.mpfr_round = temp_round;
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_factorial,
"factorial(n) -> mpfr\n\n"
"Return the floating-point approximation to the factorial of n.\n\n"
"See fac(n) to get the exact integer result.");

PyDoc_STRVAR(GMPy_doc_context_factorial,
"context.factorial(n) -> mpfr\n\n"
"Return the floating-point approximation to the factorial of n.\n\n"
"See fac(n) to get the exact integer result.");

static PyObject *
GMPy_Context_Factorial(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    long n;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    n = PyLong_AsLong(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("factorial() requires 'int' argument");
        return NULL;
    }

    if (n < 0) {
        VALUE_ERROR("factorial() of negative number");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    mpfr_clear_flags();

    mpfr_fac_ui(result->f, n, GET_MPFR_ROUND(context));

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_function_fsum,
"fsum(iterable) -> mpfr\n\n"
"Return an accurate sum of the values in the iterable.");

PyDoc_STRVAR(GMPy_doc_context_fsum,
"fsum(iterable) -> mpfr\n\n"
"Return an accurate sum of the values in the iterable.");

static PyObject *
GMPy_Context_Fsum(PyObject *self, PyObject *other)
{
    MPFR_Object *temp, *result;
    mpfr_ptr *tab;
    int errcode;
    Py_ssize_t i, seq_length = 0;
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    if (!(result = GMPy_MPFR_New(0, context))) {
        return NULL;
    }

    if (!(other = PySequence_List(other))) {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }

    /* other contains a new list containing all the values from the
     * iterable. Now make sure each item in the list is an mpfr.
     */

    seq_length = PyList_GET_SIZE(other);
    if (seq_length > LONG_MAX) {
        OVERFLOW_ERROR("temporary array is too large");
	Py_DECREF(other);
	Py_DECREF((PyObject*)result);
	return NULL;
    }
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

    if (!(tab = (mpfr_ptr *)malloc((sizeof(mpfr_srcptr) * seq_length)))) {
        Py_DECREF(other);
        Py_DECREF((PyObject*)result);
        return PyErr_NoMemory();
    }
    for (i=0; i < seq_length; i++) {
        temp = (MPFR_Object*)PyList_GET_ITEM(other, i);
        tab[i] = temp->f;
    }

    mpfr_clear_flags();

    /* The cast is safe since we have compared seq_length to LONG_MAX. */
    result->rc = mpfr_sum(result->f, tab, (unsigned long)seq_length, GET_MPFR_ROUND(context));
    Py_DECREF(other);
    free(tab);

    _GMPy_MPFR_Cleanup(&result, context);
    return (PyObject*)result;
}
