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

static int isRational(PyObject* obj);

int Pympq_convert_arg(PyObject *arg, PyObject **ptr);

#ifdef PY2
static PympqObject *   Pympq_From_PyInt(PyObject *self);
static PyObject *      Pympq_To_PyInt(PympqObject *self);
#endif

static PympqObject *   Pympq_From_PyLong(PyObject *self);
static PympqObject *   Pympq_From_Pympz(PyObject *self);
static PympqObject *   Pympq_From_Pyxmpz(PyObject * obj);
static PympqObject *   Pympq_From_PyFloat(PyObject *self);
static PympqObject *   Pympq_From_PyStr(PyObject *stringarg, long base);
static PympqObject *   Pympq_From_Decimal(PyObject* obj);
static PyObject *      Pympq_To_PyLong(PympqObject *self);
static PympzObject *   Pympq_To_Pympz(PyObject *self);
static PyxmpzObject *  Pympq_To_Pyxmpz(PyObject *self);
static PyObject *      Pympq_To_PyFloat(PympqObject *self);
static PyObject *      Pympq_To_PyStr(PympqObject *self, int base, int option);

/* NOTE: Pympq_From_Decimal returns an invalid mpq object when attempting to
 *       convert a NaN or inifinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 */

static PympqObject *   Pympq_From_Number(PyObject* obj);

/* support str() and repr() */
static PyObject *      Pympq_To_Str(PympqObject *self);
static PyObject *      Pympq_To_Repr(PympqObject *self);

#ifdef __cplusplus
}
#endif
#endif



