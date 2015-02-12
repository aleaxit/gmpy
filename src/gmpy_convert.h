/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert.h                                                          *
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

#ifndef GMPY_CONVERT_H
#define GMPY_CONVERT_H

#ifdef __cplusplus
extern "C" {
#endif

/* The following functions identify and classify the numeric types that are
 * supported by gmpy2.
 *
 * These checks are currently implemented as functions but may be
 * implemented as macros in the future.
 */

/* Checks for mpz, xmpz, and the integer types included with Python. */
static int isInteger(PyObject* obj);

/* Checks for the Fraction type included with Python. */
static int isFraction(PyObject* obj);

/* Combined mpq, isInteger() and isFraction() check. */
static int isRational(PyObject* obj);

#ifdef WITHMPFR
/* Checks for the Decimal type included with Python. */
static int isDecimal(PyObject* obj);

/* Combined mpfr, PyFloat, isDecimal() and isRational() check. */
static int isReal(PyObject* obj);
#endif

#ifdef WITHMPC
/* Combined mpc, PyComplex, and isReal() check. */
static int isComplex(PyObject* obj);
#endif

/********* Integer Conversions *********/

/* conversions with Pympz */
#ifdef PY2
static PympzObject *   Pympz_From_PyInt(PyObject *self);
#endif
static PympzObject *   Pympz_From_PyStr(PyObject *s, int base);
static PympzObject *   Pympz_From_PyLong(PyObject *obj);
static PympzObject *   Pympz_From_Pyxmpz(PyObject *self);
static PympzObject *   Pympz_From_PyFloat(PyObject *self);
static PympzObject *   Pympz_From_Number(PyObject *obj);

static PyObject *      Pympz_To_PyLong(PympzObject *self);
static PyObject *      Pympz_To_PyIntOrLong(PympzObject *self);
static PyObject *      Pympz_To_PyFloat(PympzObject *self);
static PyObject *      Pympz_To_PyStr(PympzObject *self, int base, int option);

/* conversions with Pyxmpz */
#ifdef PY2
static PyxmpzObject *  Pyxmpz_From_PyInt(PyObject *self);
#endif
static PyxmpzObject *  Pyxmpz_From_PyStr(PyObject *s, int base);
static PyxmpzObject *  Pyxmpz_From_PyLong(PyObject *obj);
static PyxmpzObject *  Pyxmpz_From_Pyxmpz(PyObject *self);
static PyxmpzObject *  Pyxmpz_From_Pympz(PyObject *self);
static PyxmpzObject *  Pyxmpz_From_PyFloat(PyObject *self);
static PyxmpzObject *  Pyxmpz_From_Number(PyObject *obj);

static PyObject *      Pyxmpz_To_PyLong(PyxmpzObject *self);
static PyObject *      Pyxmpz_To_PyIntOrLong(PyxmpzObject *self);
static PyObject *      Pyxmpz_To_PyStr(PyxmpzObject *self, int base, int option);

/* for str() and repr() */
static PyObject *      Pympz_To_Str(PympzObject *self);
static PyObject *      Pympz_To_Repr(PympzObject *self);
static PyObject *      Pyxmpz_To_Str(PyxmpzObject *self);
static PyObject *      Pyxmpz_To_Repr(PyxmpzObject *self);

/* Miscellaneous integer conversion functions. */
static int mpz_set_PyStr(mpz_ptr z, PyObject *s, int base);
static PyObject * mpz_ascii(mpz_t z, int base, int option);
static PyObject * xmpz_ascii(mpz_t z, int base, int option);
#if 0
static int Pympz_convert_arg(PyObject *arg, PyObject **ptr);
#endif

/* The following functions convert isInteger() objects into C types. */

#ifdef WITHMPFR
/* Should only be used by MPFR/MPC related code. */
static long clong_From_Integer(PyObject *obj);
#endif

/* MPIR 2.6.0 introduces two new data types: mpir_si and mpir_ui. On all
 * platforms except 64-bit Windows, those data types correspond to "long" and
 * "unsigned long". On 64-bit Windows, those data type correspond to
 * "long long" and "unsigned long long".
 */

/* Should be used by all MPIR/GMP related code. */
static mpir_si SI_From_Integer(PyObject *obj);
static mpir_ui UI_From_Integer(PyObject *obj);

static Py_ssize_t ssize_t_From_Integer(PyObject *obj);

/********* Pympq Conversions *********/

/* Conversions with Pympq */
#ifdef PY2
static PympqObject *   Pympq_From_PyInt(PyObject *self);
static PyObject *      Pympq_To_PyInt(PympqObject *self);
#endif
static PympqObject *   Pympq_From_PyLong(PyObject *self);
static PympqObject *   Pympq_From_Pympz(PyObject *self);
static PympqObject *   Pympq_From_Pyxmpz(PyObject *obj);
static PympqObject *   Pympq_From_PyFloat(PyObject *self);
static PympqObject *   Pympq_From_Fraction(PyObject *obj);
static PympqObject *   Pympq_From_PyStr(PyObject *stringarg, int base);

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
static PympqObject *   Pympq_From_Decimal(PyObject* obj);

static PympqObject *   Pympq_From_Number(PyObject* obj);

static PyObject *      Pympq_To_PyLong(PympqObject *self);
static PympzObject *   Pympq_To_Pympz(PyObject *self);
static PyxmpzObject *  Pympq_To_Pyxmpz(PyObject *self);
static PyObject *      Pympq_To_PyFloat(PympqObject *self);
static PyObject *      Pympq_To_PyStr(PympqObject *self, int base, int option);

/* support str() and repr() */
static PyObject *      Pympq_To_Str(PympqObject *self);
static PyObject *      Pympq_To_Repr(PympqObject *self);

/* Miscellaneous rational conversion functions. */
int Pympq_convert_arg(PyObject *arg, PyObject **ptr);

#ifdef WITHMPFR
/********* Pympfr Conversions *********/

/* Conversions with Pympfr */
#ifdef PY2
static PympfrObject *   Pympfr_From_PyInt(PyObject *self, mpfr_prec_t bits);
static PyObject *       Pympfr_To_PyInt(PympfrObject *self);
#endif
static PympfrObject *   Pympfr_From_Pympfr(PyObject *self, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_PyFloat(PyObject *self, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_PyLong(PyObject *self, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_Pympz(PyObject *self, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_Pympq(PyObject *self, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_PyStr(PyObject *s, int base, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_Decimal(PyObject *obj, mpfr_prec_t bits);
static PympfrObject *   Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits);

static PympzObject *    Pympfr_To_Pympz(PyObject *self);
static PyxmpzObject *   Pympfr_To_Pyxmpz(PyObject *self);
static PympqObject *    Pympfr_To_Pympq(PyObject *self);
static PyObject *       Pympfr_To_PyLong(PympfrObject *self);
static PyObject *       Pympfr_To_PyFloat(PympfrObject *self);
static PyObject *       Pympfr_To_PyStr(PympfrObject *self, int base, int digits);

/* support str() and repr() */
static PyObject *       Pympfr_To_Str(PympfrObject *self);
static PyObject *       Pympfr_To_Repr(PympfrObject *self);

/* Miscellaneous */
static int Pympfr_convert_arg(PyObject *arg, PyObject **ptr);
static PympqObject * stern_brocot(PympfrObject* self, PympfrObject *err, mpfr_prec_t prec, int mayz);
#endif

#ifdef WITHMPC
/********* Pympc Conversions *********/

/* Conversions with Pympc */
#ifdef PY2
static PympcObject *   Pympc_From_PyInt(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PyObject *      Pympc_To_PyIntOrLong(PyObject *self);
#endif
static PympcObject *   Pympc_From_Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_PyComplex(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_Pympfr(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_PyFloat(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_Pympz(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_Pympq(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_PyLong(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static PympcObject *   Pympc_From_PyStr(PyObject *s, int base, mpfr_prec_t rbits, mpfr_prec_t ibits);
static PympcObject *   Pympc_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec);

static PyObject *      Pympc_To_PyFloat(PyObject *self);
static PyObject *      Pympc_To_PyLong(PyObject *self);
static PyObject *      Pympc_To_PyStr(PympcObject *self, int base, int digits);
static PyObject *      Pympc_To_PyComplex(PyObject *self, PyObject *other);
/* support str() and repr() */
static PyObject *      Pympc_To_Str(PympcObject *self);
static PyObject *      Pympc_To_Repr(PympcObject *self);

/* Miscellaneous */
static PyObject * raw_mpfr_ascii(mpfr_t self, int base, int digits, int round);
int Pympc_convert_arg(PyObject *arg, PyObject **ptr);
#endif

#ifdef __cplusplus
}
#endif
#endif
