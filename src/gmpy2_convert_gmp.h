/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert.h                                                          *
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

#ifndef GMPY2_CONVERT_GMP_H
#define GMPY2_CONVERT_GMP_H

#ifdef __cplusplus
extern "C" {
#endif

/* ======================================================================== *
 * Conversion between native Python objects and MPZ.                        *
 * ======================================================================== */

static MPZ_Object *    GMPy_MPZ_From_PyIntOrLong(PyObject *obj);
static MPZ_Object *    GMPy_MPZ_From_PyStr(PyObject *s, int base);
static MPZ_Object *    GMPy_MPZ_From_PyFloat(PyObject *obj);

static PyObject *      GMPy_PyLong_From_MPZ(MPZ_Object *obj);
static PyObject *      GMPy_PyIntOrLong_From_MPZ(MPZ_Object *obj);
static PyObject *      GMPy_PyFloat_From_MPZ(MPZ_Object *obj);
static PyObject *      GMPy_PyStr_From_MPZ(MPZ_Object *obj, int base, int option);

static MPZ_Object *    GMPy_MPZ_From_Number_New(PyObject *obj);
static MPZ_Object *    GMPy_MPZ_From_Number_Temp(PyObject *obj);
static MPZ_Object *    GMPy_MPZ_From_Integer_New(PyObject *obj);
static MPZ_Object *    GMPy_MPZ_From_Integer_Temp(PyObject *obj);

static PyObject *      GMPy_MPZ_Str_Slot(MPZ_Object *self);
static PyObject *      GMPy_MPZ_Repr_Slot(MPZ_Object *self);

/* Helper function not currently used. */
#if 0
static int             GMPy_MPZ_convert_arg(PyObject *arg, PyObject **ptr);
#endif

/* ======================================================================== *
 * Conversion between native Python objects/MPZ and XMPZ.                   *
 * ======================================================================== */

static XMPZ_Object *   GMPy_XMPZ_From_PyIntOrLong(PyObject *obj);
static XMPZ_Object *   GMPy_XMPZ_From_PyStr(PyObject *s, int base);
static XMPZ_Object *   GMPy_XMPZ_From_PyFloat(PyObject *self);
static XMPZ_Object *   GMPy_XMPZ_From_MPZ(MPZ_Object *obj);
static XMPZ_Object *   GMPy_XMPZ_From_XMPZ(XMPZ_Object *obj);

static PyObject *      GMPy_PyLong_From_XMPZ(XMPZ_Object *self);
static PyObject *      GMPy_PyIntOrLong_From_XMPZ(XMPZ_Object *self);
static PyObject *      GMPy_PyStr_From_XMPZ(XMPZ_Object *self, int base, int option);
static MPZ_Object *    GMPy_MPZ_From_XMPZ(XMPZ_Object *obj);

static XMPZ_Object *   GMPy_XMPZ_From_Number_New(PyObject *obj);
static XMPZ_Object *   GMPy_XMPZ_From_Number_Temp(PyObject *obj);

static PyObject *      GMPy_XMPZ_Str_Slot(XMPZ_Object *self);
static PyObject *      GMPy_XMPZ_Repr_Slot(XMPZ_Object *self);

/* ======================================================================== *
 * Conversion between Integer objects and C types.                          *
 * ======================================================================== */

/* Should only be used by MPFR/MPC related code. */
static long            clong_From_Integer(PyObject *obj);

/* MPIR 2.6.0 introduces two new data types: mpir_si and mpir_ui. On all
 * platforms except 64-bit Windows, those data types correspond to "long" and
 * "unsigned long". On 64-bit Windows, those data type correspond to
 * "long long" and "unsigned long long".
 */

/* Should be used by all MPIR/GMP related code. */
static mpir_si         SI_From_Integer(PyObject *obj);
static mpir_ui         UI_From_Integer(PyObject *obj);

static Py_ssize_t      ssize_t_From_Integer(PyObject *obj);

/* ======================================================================== *
 * Conversion between native Python objects/MPZ/XMPZ and MPQ.               *
 * ======================================================================== */

static MPQ_Object *    GMPy_MPQ_From_PyIntOrLong(PyObject *obj);
static MPQ_Object *    GMPy_MPQ_From_PyStr(PyObject *s, int base);
static MPQ_Object *    GMPy_MPQ_From_MPZ(MPZ_Object *obj);
static MPQ_Object *    GMPy_MPQ_From_XMPZ(XMPZ_Object *obj);

static MPQ_Object *    Pympq_From_PyFloat(PyObject *obj);
static MPQ_Object *    Pympq_From_Fraction(PyObject *obj);

/* NOTE: Pympq_From_Decimal returns an invalid mpq object when attempting to
 *       convert a NaN or Infinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 *
 *       If the numerator is 0 and the denominator is not 0, then the sign of
 *       the denominator is the sign of the 0.
 */
static MPQ_Object *    Pympq_From_Decimal(PyObject* obj);

static MPQ_Object *    Pympq_From_Number(PyObject* obj);

static PyObject *      Pympq_To_PyLong(MPQ_Object *self);
static MPZ_Object *    GMPy_MPZ_From_MPQ(PyObject *obj);
static XMPZ_Object *   GMPy_XMPZ_From_MPQ(PyObject *obj);
static PyObject *      Pympq_To_PyFloat(MPQ_Object *self);
static PyObject *      Pympq_To_PyStr(MPQ_Object *self, int base, int option);

/* support str() and repr() */
static PyObject *      Pympq_To_Str(MPQ_Object *self);
static PyObject *      Pympq_To_Repr(MPQ_Object *self);

/* Miscellaneous rational conversion functions. */
int Pympq_convert_arg(PyObject *arg, PyObject **ptr);


#ifdef __cplusplus
}
#endif
#endif
