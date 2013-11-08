/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_mpfr.h                                                    *
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

#ifndef GMPY2_CONVERT_MPFR_H
#define GMPY2_CONVERT_MPFR_H

#ifdef __cplusplus
extern "C" {
#endif

/********* Pympfr Conversions *********/

/* Conversions with Pympfr */
#ifdef PY2
static MPFR_Object *    Pympfr_From_PyInt(PyObject *self, mpfr_prec_t bits);
static PyObject *       Pympfr_To_PyInt(MPFR_Object *self);
#endif
static MPFR_Object *    Pympfr_From_Pympfr(PyObject *self, mpfr_prec_t bits);
#if 0
static MPFR_Object *    Pympfr_From_Pympfr_bits_context(PyObject *self, mpfr_prec_t bits, GMPyContextObject *context);
#endif
static MPFR_Object *    Pympfr_From_Pympfr_context(PyObject *self, GMPyContextObject *context);
static MPFR_Object *    Pympfr_From_PyFloat(PyObject *self, mpfr_prec_t bits);
static MPFR_Object *    Pympfr_From_PyFloat_bits_context(PyObject *self, mpfr_prec_t bits, GMPyContextObject *context);
static MPFR_Object *    Pympfr_From_PyLong(PyObject *self, mpfr_prec_t bits);
static MPFR_Object *    Pympfr_From_PyLong_context(PyObject *self, mpfr_prec_t bits, GMPyContextObject *context);
static MPFR_Object *    Pympfr_From_Pympz(PyObject *self, mpfr_prec_t bits);
static MPFR_Object *    Pympfr_From_Pympz_context(PyObject *self, mpfr_prec_t bits, GMPyContextObject *context);
static MPFR_Object *    Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits);
static MPFR_Object *    GMPy_MPFR_From_Real_Temp(PyObject* obj, GMPyContextObject *context);
static MPFR_Object *    Pympfr_From_Pympq(PyObject *self, mpfr_prec_t bits);
static MPFR_Object *    Pympfr_From_Pympq_bits_context(PyObject *self, mpfr_prec_t bits, GMPyContextObject *context);
static MPFR_Object *    Pympfr_From_PyStr(PyObject *s, int base, mpfr_prec_t bits);
static MPFR_Object *    Pympfr_From_Decimal(PyObject *obj, mpfr_prec_t bits);
static MPFR_Object *    Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits);

static MPZ_Object   *   Pympfr_To_Pympz(PyObject *self);
static XMPZ_Object *    Pympfr_To_Pyxmpz(PyObject *self);
static MPQ_Object  *    Pympfr_To_Pympq(PyObject *self);
static PyObject *       Pympfr_To_PyLong(MPFR_Object *self);
static PyObject *       Pympfr_To_PyFloat(MPFR_Object *self);
static PyObject *       Pympfr_To_PyStr(MPFR_Object *self, int base, int digits);

/* support str() and repr() */
static PyObject *       Pympfr_To_Str(MPFR_Object *self);
static PyObject *       Pympfr_To_Repr(MPFR_Object *self);

/* Miscellaneous */
static int Pympfr_convert_arg(PyObject *arg, PyObject **ptr);
static MPQ_Object * stern_brocot(MPFR_Object* self, MPFR_Object *err, mpfr_prec_t prec, int mayz);

static PyObject * raw_mpfr_ascii(mpfr_t self, int base, int digits, int round);

#ifdef __cplusplus
}
#endif
#endif
