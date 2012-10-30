/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_conv.h                                                         *
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

#ifndef GMPY_MPZ_CONV_H
#define GMPY_MPZ_CONV_H

#ifdef __cplusplus
extern "C" {
#endif

static int isFraction(PyObject* obj);
static int isDecimal(PyObject* obj);
static int isInteger(PyObject* obj);

static int mpz_set_PyStr(mpz_ptr z, PyObject *s, long base);
static PyObject * mpz_ascii(mpz_t z, int base, int option);
static PyObject * xmpz_ascii(mpz_t z, int base, int option);
#if 0
static int Pympz_convert_arg(PyObject *arg, PyObject **ptr);
#endif
static Py_ssize_t ssize_t_From_Integer(PyObject *obj);

/* conversions with Pympz */
#ifdef PY2
static PympzObject *   Pympz_From_PyInt(PyObject *self);
#endif
static PympzObject *   Pympz_From_PyLong(PyObject * obj);
static PympzObject *   Pympz_From_Pyxmpz(PyObject *self);
static PympzObject *   Pympz_From_PyFloat(PyObject *self);
static PympzObject *   Pympz_From_PyStr(PyObject *s, long base);
static PyObject *      Pympz_To_PyLong(PympzObject *self);
static PyObject *      Pympz_To_PyIntOrLong(PympzObject *self);
static PyObject *      Pympz_To_PyFloat(PympzObject *self);
static PyObject *      Pympz_To_PyStr(PympzObject *self, int base, int option);

static PympzObject *   Pympz_From_Number(PyObject* obj);
static PympzObject *   Pympz_From_Integer(PyObject* obj);

/* conversions with Pyxmpz */
#ifdef PY2
static PyxmpzObject *  Pyxmpz_From_PyInt(PyObject *self);
#endif
static PyxmpzObject *  Pyxmpz_From_PyLong(PyObject * obj);
static PyxmpzObject *  Pyxmpz_From_Pyxmpz(PyObject *self);
static PyxmpzObject *  Pyxmpz_From_Pympz(PyObject *self);
static PyxmpzObject *  Pyxmpz_From_PyFloat(PyObject *self);
static PyxmpzObject *  Pyxmpz_From_PyStr(PyObject *s, long base);
static PyObject *      Pyxmpz_To_PyLong(PyxmpzObject *self);
static PyObject *      Pyxmpz_To_PyIntOrLong(PyxmpzObject *self);
static PyObject *      Pyxmpz_To_PyStr(PyxmpzObject *self, int base, int option);

static PyxmpzObject *  Pyxmpz_From_Number(PyObject* obj);

/* for str() and repr() */
static PyObject *      Pympz_To_Str(PympzObject *self);
static PyObject *      Pympz_To_Repr(PympzObject *self);
static PyObject *      Pyxmpz_To_Str(PyxmpzObject *self);
static PyObject *      Pyxmpz_To_Repr(PyxmpzObject *self);

/* Should only be used by MPFR/MPC related code. */
static long clong_From_Integer(PyObject *obj);

/* Should be used by all MPIR/GMP related code. */
static mpir_si SI_From_Integer(PyObject *obj);
static mpir_ui UI_From_Integer(PyObject *obj);

#ifdef __cplusplus
}
#endif
#endif



