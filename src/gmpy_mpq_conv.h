/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq_conv.c                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen                  *
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

#ifndef GMPY_MPQ_CONV_H
#define GMPY_MPQ_CONV_H

#ifdef __cplusplus
extern "C" {
#endif

static PympqObject * Pympz2Pympq(PyObject *self);
static PympqObject * Pyxmpz2Pympq(PyObject * obj);
static PympzObject * Pympq2Pympz(PyObject *self);
static PyxmpzObject * Pympq2Pyxmpz(PyObject *self);
static PympqObject * PyLong2Pympq(PyObject *self);
static PympqObject * PyFloat2Pympq(PyObject *self);
static PympqObject * PyStr2Pympq(PyObject *stringarg, long base);
static PyObject * Pympq2PyLong(PympqObject *self);
#ifdef PY2
static PympqObject * PyInt2Pympq(PyObject *self);
static PyObject * Pympq2PyInt(PympqObject *self);
#endif
static PyObject * Pympq2PyFloat(PympqObject *self);
static PyObject * Pympq_ascii(PympqObject *self, int base, int option);

static int isRational(PyObject* obj);

/* NOTE: Pympq_From_Decimal returns an invalid mpq object when attempting to
 *       convert a NaN or inifinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 */

static PympqObject * Pympq_From_Decimal(PyObject* obj);
static PympqObject * Pympq_From_Real(PyObject* obj);
static PympqObject * Pympq_From_Rational(PyObject* obj);
int Pympq_convert_arg(PyObject *arg, PyObject **ptr);
static PyObject * Pympq2str(PympqObject *self);
static PyObject * Pympq2repr(PympqObject *self);

#ifdef __cplusplus
}
#endif
#endif



