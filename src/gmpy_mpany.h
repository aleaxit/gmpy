/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_context.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013 Case Van Horsen            *
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

#ifndef GMPY_MPANY_H
#define GMPY_MPANY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations begin here. */

static PyObject * Pympany_square(PyObject *self, PyObject *other);
static PyObject * Pympany_digits(PyObject *self, PyObject *args);
static PyObject * Pympany_sign(PyObject *self, PyObject *other);
static PyObject * Pympany_add(PyObject *self, PyObject *args);
static PyObject * Pympany_sub(PyObject *self, PyObject *args);
static PyObject * Pympany_mul(PyObject *self, PyObject *args);
static PyObject * Pympany_div(PyObject *self, PyObject *args);
static PyObject * Pympany_to_binary(PyObject *self, PyObject *other);
static PyObject * Pympany_pow(PyObject *base, PyObject *exp, PyObject *mod);
static PyObject * Pympany_printf(PyObject *self, PyObject *args);
#ifdef WITHMPFR
static PyObject * Pympany_is_nan(PyObject *self, PyObject *other);
static PyObject * Pympany_is_inf(PyObject *self, PyObject *other);
static PyObject * Pympany_is_infinite(PyObject *self, PyObject *other);
static PyObject * Pympany_is_zero(PyObject *self, PyObject *other);
static PyObject * Pympany_log(PyObject *self, PyObject *other);
static PyObject * Pympany_exp(PyObject *self, PyObject *other);
static PyObject * Pympany_cos(PyObject *self, PyObject *other);
static PyObject * Pympany_sin(PyObject *self, PyObject *other);
static PyObject * Pympany_tan(PyObject *self, PyObject *other);
static PyObject * Pympany_acos(PyObject *self, PyObject *other);
static PyObject * Pympany_asin(PyObject *self, PyObject *other);
static PyObject * Pympany_atan(PyObject *self, PyObject *other);
static PyObject * Pympany_cosh(PyObject *self, PyObject *other);
static PyObject * Pympany_sinh(PyObject *self, PyObject *other);
static PyObject * Pympany_tanh(PyObject *self, PyObject *other);
static PyObject * Pympany_acosh(PyObject *self, PyObject *other);
static PyObject * Pympany_asinh(PyObject *self, PyObject *other);
static PyObject * Pympany_atanh(PyObject *self, PyObject *other);
static PyObject * Pympany_sqrt(PyObject *self, PyObject *other);
static PyObject * Pympany_sin_cos(PyObject *self, PyObject *other);
static PyObject * Pympany_fma(PyObject *self, PyObject *other);
static PyObject * Pympany_fms(PyObject *self, PyObject *other);
static PyObject * Pympany_div_2exp(PyObject *self, PyObject *other);
static PyObject * Pympany_mul_2exp(PyObject *self, PyObject *other);
#endif
static PyObject * mpany_richcompare(PyObject *a, PyObject *b, int op);

#ifdef __cplusplus
}
#endif
#endif
