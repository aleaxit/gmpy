/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_xmpz_inplace.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015 Case Van Horsen                                          *
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

#ifndef GMPY_XMPZ_INPLACE_H
#define GMPY_XMPZ_INPLACE_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * Pyxmpz_inplace_add(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_sub(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_mul(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_floordiv(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_rem(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_rshift(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_lshift(PyObject *a, PyObject *b);
static PyObject * Pyxmpz_inplace_pow(PyObject *in_b, PyObject *in_e, PyObject *in_m);
static PyObject * Pyxmpz_inplace_and(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_inplace_xor(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_inplace_ior(PyObject *self, PyObject *other);

#ifdef __cplusplus
}
#endif
#endif
