/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mul.h                                                              *
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

#ifndef GMPY2_MUL_H
#define GMPY2_MUL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Public API */

static PyObject * GMPy_Integer_Mul(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Rational_Mul(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Real_Mul(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Complex_Mul(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Mul(PyObject *x, PyObject *y, CTXT_Object *context);

/* Private API */

static PyObject * GMPy_mpz_mul_fast(PyObject *x, PyObject *y);
static PyObject * GMPy_mpq_mul_fast(PyObject *x, PyObject *y);
static PyObject * GMPy_mpfr_mul_fast(PyObject *x, PyObject *y);
static PyObject * GMPy_mpc_mul_fast(PyObject *x, PyObject *y);

static PyObject * GMPy_Context_Mul(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
