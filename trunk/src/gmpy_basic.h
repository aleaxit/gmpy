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
#ifndef GMPY_BASIC_H
#define GMPY_BASIC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Forward declarations begin here. */

static PyObject * Pybasic_add(PyObject *a, PyObject *b);
static PyObject * Pybasic_sub(PyObject *a, PyObject *b);
static PyObject * Pybasic_mul(PyObject *a, PyObject *b);
static PyObject * Pybasic_floordiv(PyObject *a, PyObject *b);
static PyObject * Pybasic_truediv(PyObject *a, PyObject *b);
#ifdef PY2
static PyObject * Pybasic_div2(PyObject *a, PyObject *b);
#endif
static PyObject * Pybasic_rem(PyObject *a, PyObject *b);
static PyObject * Pybasic_divmod(PyObject *a, PyObject *b);

#ifdef __cplusplus
}
#endif
#endif








