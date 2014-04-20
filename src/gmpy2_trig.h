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
static PyObject * GMPy_Context_Sin(PyObject *self, PyObject *args);
static PyObject * GMPy_Function_Sin(PyObject *self, PyObject *other);

static PyObject * Pympfr_cos(PyObject* self, PyObject *other);
static PyObject * Pympfr_tan(PyObject* self, PyObject *other);
static PyObject * Pympfr_sec(PyObject* self, PyObject *other);
static PyObject * Pympfr_csc(PyObject* self, PyObject *other);
static PyObject * Pympfr_cot(PyObject* self, PyObject *other);
static PyObject * Pympfr_acos(PyObject* self, PyObject *other);
static PyObject * Pympfr_asin(PyObject* self, PyObject *other);
static PyObject * Pympfr_atan(PyObject* self, PyObject *other);
static PyObject * Pympfr_cosh(PyObject* self, PyObject *other);
static PyObject * Pympfr_sinh(PyObject* self, PyObject *other);
static PyObject * Pympfr_tanh(PyObject* self, PyObject *other);
static PyObject * Pympfr_sech(PyObject* self, PyObject *other);
static PyObject * Pympfr_csch(PyObject* self, PyObject *other);
static PyObject * Pympfr_coth(PyObject* self, PyObject *other);
static PyObject * Pympfr_acosh(PyObject* self, PyObject *other);
static PyObject * Pympfr_asinh(PyObject* self, PyObject *other);
static PyObject * Pympfr_atanh(PyObject* self, PyObject *other);
static PyObject * Pympfr_atan2(PyObject *self, PyObject *args);
static PyObject * Pympfr_hypot(PyObject *self, PyObject *args);
static PyObject * Pympfr_sin_cos(PyObject *self, PyObject *other);
static PyObject * Pympfr_sinh_cosh(PyObject *self, PyObject *other);
static PyObject * Pympfr_degrees(PyObject *self, PyObject *other);
static PyObject * Pympfr_radians(PyObject *self, PyObject *other);

#ifdef __cplusplus
}
#endif
#endif
