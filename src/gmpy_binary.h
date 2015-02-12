/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_binary.c                                                           *
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

#ifndef GMPY_BINARY_H
#define GMPY_BINARY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Conversion routines between GMPY2 objects and a compact, portable
 * binary representation. The binary format of GMPY2 is not compatible
 * with GMPY 1.x. Methods to read the old format are provided.
 */

static PyObject * Pympz_From_Old_Binary(PyObject *self, PyObject *other);
static PyObject * Pympq_From_Old_Binary(PyObject *self, PyObject *other);

static PyObject * Pympany_From_Binary(PyObject *self, PyObject *other);
static PyObject * Pympz_To_Binary(PympzObject *self);
static PyObject * Pyxmpz_To_Binary(PyxmpzObject *self);
static PyObject * Pympq_To_Binary(PympqObject *self);

#ifdef WITHMPFR
static PyObject * Pympfr_From_Old_Binary(PyObject *self, PyObject *other);
static PyObject * Pympfr_To_Binary(PympfrObject *self);
#endif

#ifdef WITHMPC
static PyObject * Pympc_To_Binary(PympcObject *self);
#endif


#ifdef __cplusplus
}
#endif
#endif
