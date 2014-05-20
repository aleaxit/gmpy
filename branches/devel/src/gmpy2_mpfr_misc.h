/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpfr_misc.h                                                       *
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

#ifndef GMPY_MPFR_MISC_H
#define GMPY_MPFR_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_Real_F2Q(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_F2Q(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_F2Q(PyObject *self, PyObject *args);

static PyObject * GMPy_MPFR_get_emin_min(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_get_emax_max(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_get_max_precision(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_get_exp(PyObject *self, PyObject *other);
static PyObject * GMPy_MPFR_set_exp(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_set_sign(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_copy_sign(PyObject *self, PyObject *args);
static PyObject * Pympfr_integer_ratio(PyObject *self, PyObject *args);
static PyObject * Pympfr_mantissa_exp(PyObject *self, PyObject *args);
static PyObject * Pympfr_simple_fraction(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * GMPy_MPFR_set_nan(PyObject *self, PyObject *other);
static PyObject * GMPy_MPFR_set_inf(PyObject *self, PyObject *args);
static PyObject * GMPy_MPFR_set_zero(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
