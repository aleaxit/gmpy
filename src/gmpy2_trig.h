/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr_trig.h                                                       *
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

#ifndef GMPY_MPFR_TRIG_H
#define GMPY_MPFR_TRIG_H

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

static PyObject * Pympfr_atan2(PyObject *self, PyObject *args);
static PyObject * Pympfr_hypot(PyObject *self, PyObject *args);
static PyObject * Pympfr_sin_cos(PyObject *self, PyObject *other);
static PyObject * Pympfr_sinh_cosh(PyObject *self, PyObject *other);
static PyObject * GMPy_Context_Degrees(PyObject *self, PyObject *other);
static PyObject * GMPy_Context_Radians(PyObject *self, PyObject *other);

#ifdef __cplusplus
}
#endif
#endif
