/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_math.h                                                            *
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

#ifndef GMPY_MATH_H
#define GMPY_MATH_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_Real_Sin(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Sin(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sin(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sin(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Cos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Cos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Cos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Cos(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Tan(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Tan(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Tan(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Tan(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Atan(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Atan(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Atan(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Atan(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Sinh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Sinh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sinh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sinh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Cosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Cosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Cosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Cosh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Tanh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Tanh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Tanh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Tanh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Asinh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Asinh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Asinh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Asinh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Acosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Acosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Acosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Acosh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Sec(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sec(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sec(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Csc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Csc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Csc(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Cot(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Cot(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Cot(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Sech(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sech(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sech(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Csch(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Csch(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Csch(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Coth(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Coth(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Coth(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Acos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Acos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Acos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Acos(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Asin(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Asin(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Asin(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Asin(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Atanh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Atanh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Atanh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Atanh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Sin_Cos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Sin_Cos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sin_Cos(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sin_Cos(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Sinh_Cosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sinh_Cosh(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sinh_Cosh(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Atan2(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Atan2(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Atan2(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Hypot(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Hypot(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Hypot(PyObject *self, PyObject *args);

static PyObject * GMPy_Context_Degrees(PyObject *self, PyObject *other);
static PyObject * GMPy_Context_Radians(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Log(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Log(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Log(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Log(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Log10(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Log10(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Log10(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Log10(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Exp(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Exp(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Exp(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Exp(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Sqrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Complex_Sqrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Sqrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Sqrt(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_RecSqrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_RecSqrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_RecSqrt(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Rint(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Rint(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Rint(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_RintCeil(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_RintCeil(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_RintCeil(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_RintFloor(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_RintFloor(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_RintFloor(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_RintRound(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_RintRound(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_RintRound(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_RintTrunc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_RintTrunc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_RintTrunc(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Frac(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Frac(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Frac(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Cbrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Cbrt(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Cbrt(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Log2(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Log2(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Log2(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Exp2(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Exp2(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Exp2(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Exp10(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Exp10(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Exp10(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Log1p(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Log1p(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Log1p(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Expm1(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Expm1(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Expm1(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Eint(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Eint(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Eint(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Li2(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Li2(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Li2(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Lngamma(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Lngamma(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Lngamma(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Digamma(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Digamma(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Digamma(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Zeta(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Zeta(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Zeta(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Erf(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Erf(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Erf(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Erfc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Erfc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Erfc(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_J0(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_J0(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_J0(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_J1(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_J1(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_J1(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Y0(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Y0(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Y0(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Y1(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Y1(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Y1(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Ai(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Ai(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Ai(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);
static PyObject * GMPy_Number_FMA(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);
static PyObject * GMPy_Context_FMA(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);
static PyObject * GMPy_Number_FMS(PyObject *x, PyObject *y, PyObject *z, CTXT_Object *context);
static PyObject * GMPy_Context_FMS(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Root(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Root(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Root(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Jn(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Jn(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Jn(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Yn(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Yn(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Yn(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_AGM(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_AGM(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_AGM(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Maxnum(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Maxnum(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Maxnum(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Minnum(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Minnum(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Minnum(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Remainder(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Remainder(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Remainder(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Fmod(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Fmod(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Fmod(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_RelDiff(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_RelDiff(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_RelDiff(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Ceil(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Ceil(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_MPFR_Method_Ceil(PyObject *self, PyObject *args);
static PyObject * GMPy_Context_Ceil(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Floor(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Floor(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_MPFR_Method_Floor(PyObject *self, PyObject *args);
static PyObject * GMPy_Context_Floor(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Trunc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Trunc(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_MPFR_Method_Trunc(PyObject *self, PyObject *args);
static PyObject * GMPy_Context_Trunc(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_Round2(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Round2(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Round2(PyObject *self, PyObject *args);

static PyObject * GMPy_Real_RoundAway(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_RoundAway(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_RoundAway(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Modf(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Modf(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Modf(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Lgamma(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Lgamma(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Lgamma(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_RemQuo(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_RemQuo(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_RemQuo(PyObject *self, PyObject *other);

static PyObject * GMPy_Real_Frexp(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Frexp(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Frexp(PyObject *self, PyObject *other);

static PyObject * GMPy_Context_NextToward(PyObject *self, PyObject *args);

static PyObject * GMPy_Context_NextAbove(PyObject *self, PyObject *other);

static PyObject * GMPy_Context_NextBelow(PyObject *self, PyObject *other);

static PyObject * GMPy_Context_Factorial(PyObject *self, PyObject *other);

static PyObject * GMPy_Context_Fsum(PyObject *self, PyObject *other);

#if 0
static PyObject * GMPy_Context_Vector(PyObject *self, PyObject *other);
#endif

#ifdef __cplusplus
}
#endif
#endif
