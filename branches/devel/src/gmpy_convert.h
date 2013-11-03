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
#ifdef PY2
#define IS_INTEGER_ONLY(x) (MPZ_Check(x) || PyInt_Check(x) || PyLong_Check(x) || XMPZ_Check(x))
#else
#define IS_INTEGER_ONLY(x) (MPZ_Check(x) || PyLong_Check(x) || XMPZ_Check(x))
#endif
#define IS_INTEGER(x) IS_INTEGER_ONLY(x)

#define IS_RATIONAL_ONLY(x) (MPQ_Check(x) || (!strcmp(Py_TYPE(x)->tp_name, "Fraction")))
#define IS_RATIONAL(x) (IS_INTEGER(x) || IS_RATIONAL_ONLY(x))

#if PY_VERSION_HEX < 0x03030000
#define IS_DECIMAL_ONLY(x) (!strcmp(Py_TYPE(x)->tp_name, "Decimal"))
#else
#define IS_DECIMAL_ONLY(x) (!strcmp(Py_TYPE(x)->tp_name, "decimal.Decimal"))
#endif

#define IS_REAL_ONLY(x) (MPFR_Check(x) || PyFloat_Check(x) || IS_DECIMAL_ONLY(x))
#define IS_REAL(x) (IS_RATIONAL(x) || IS_REAL_ONLY(x))

#define IS_COMPLEX_ONLY(x) (MPC_Check(x) || PyComplex_Check(x))
#define IS_COMPLEX(x) (IS_REAL(x) || IS_COMPLEX_ONLY(x))

/* Checks for mpz, xmpz, and the integer types included with Python. */
static int isInteger(PyObject* obj);

/* Checks for the Fraction type included with Python. */
static int isFraction(PyObject* obj);

/* Combined mpq, isInteger() and isFraction() check. */
static int isRational(PyObject* obj);

/* Checks for the Decimal type included with Python. */
static int isDecimal(PyObject* obj);

/* Combined mpfr, PyFloat, isDecimal() and isRational() check. */
static int isReal(PyObject* obj);

/* Combined mpc, PyComplex, and isReal() check. */
static int isComplex(PyObject* obj);

/********* Integer Conversions *********/

/* conversions with MPZ_Object */
#ifdef PY2
static MPZ_Object *    Pympz_From_PyInt(PyObject *self);
#endif
static MPZ_Object *    Pympz_From_PyStr(PyObject *s, int base);
static MPZ_Object *    GMPy_MPZ_From_PyLong(PyObject *obj);
static MPZ_Object *    Pympz_From_Pyxmpz(PyObject *self);
static MPZ_Object *    Pympz_From_PyFloat(PyObject *self);
static MPZ_Object *    GMPy_MPZ_From_Integer(PyObject *obj);
static MPZ_Object *    Pympz_From_Number(PyObject *obj);

static PyObject *      Pympz_To_PyLong(MPZ_Object *self);
static PyObject *      Pympz_To_PyIntOrLong(MPZ_Object *self);
static PyObject *      Pympz_To_PyFloat(MPZ_Object *self);
static PyObject *      Pympz_To_PyStr(MPZ_Object *self, int base, int option);

/* conversions with XMPZ_Object */
#ifdef PY2
static XMPZ_Object *   Pyxmpz_From_PyInt(PyObject *self);
#endif
static XMPZ_Object *   Pyxmpz_From_PyStr(PyObject *s, int base);
static XMPZ_Object *   GMPy_XMPZ_From_PyLong(PyObject *obj);
static XMPZ_Object *   Pyxmpz_From_Pyxmpz(PyObject *self);
static XMPZ_Object *   Pyxmpz_From_Pympz(PyObject *self);
static XMPZ_Object *   Pyxmpz_From_PyFloat(PyObject *self);
static XMPZ_Object *   Pyxmpz_From_Number(PyObject *obj);

static PyObject *      Pyxmpz_To_PyLong(XMPZ_Object *self);
static PyObject *      Pyxmpz_To_PyIntOrLong(XMPZ_Object *self);
static PyObject *      Pyxmpz_To_PyStr(XMPZ_Object *self, int base, int option);

/* for str() and repr() */
static PyObject *      Pympz_To_Str(MPZ_Object *self);
static PyObject *      Pympz_To_Repr(MPZ_Object *self);
static PyObject *      Pyxmpz_To_Str(XMPZ_Object *self);
static PyObject *      Pyxmpz_To_Repr(XMPZ_Object *self);

/* Miscellaneous integer conversion functions. */
static int mpz_set_PyStr(mpz_ptr z, PyObject *s, int base);
static PyObject * mpz_ascii(mpz_t z, int base, int option);
static PyObject * xmpz_ascii(mpz_t z, int base, int option);
#if 0
static int Pympz_convert_arg(PyObject *arg, PyObject **ptr);
#endif

/* The following functions convert isInteger() objects into C types. */

/* Should only be used by MPFR/MPC related code. */
static long clong_From_Integer(PyObject *obj);

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
static MPQ_Object *    Pympq_From_PyInt(PyObject *self);
static PyObject *      Pympq_To_PyInt(MPQ_Object *self);
#endif
static MPQ_Object *    Pympq_From_PyLong(PyObject *self);
static MPQ_Object *    Pympq_From_Pympz(PyObject *self);
static MPQ_Object *    Pympq_From_Pyxmpz(PyObject *obj);
static MPQ_Object *    Pympq_From_PyFloat(PyObject *self);
static MPQ_Object *    Pympq_From_Fraction(PyObject *obj);
static MPQ_Object *    Pympq_From_PyStr(PyObject *stringarg, int base);

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
static MPZ_Object *    Pympq_To_Pympz(PyObject *self);
static XMPZ_Object *   Pympq_To_Pyxmpz(PyObject *self);
static PyObject *      Pympq_To_PyFloat(MPQ_Object *self);
static PyObject *      Pympq_To_PyStr(MPQ_Object *self, int base, int option);

/* support str() and repr() */
static PyObject *      Pympq_To_Str(MPQ_Object *self);
static PyObject *      Pympq_To_Repr(MPQ_Object *self);

/* Miscellaneous rational conversion functions. */
int Pympq_convert_arg(PyObject *arg, PyObject **ptr);

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
static XMPZ_Object *   Pympfr_To_Pyxmpz(PyObject *self);
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

/********* Pympc Conversions *********/

/* Conversions with Pympc */
#ifdef PY2
static MPC_Object *   Pympc_From_PyInt_bits_context(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec, GMPyContextObject *context);
static PyObject *      Pympc_To_PyIntOrLong(PyObject *self);
#endif
static MPC_Object *   Pympc_From_Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_PyComplex(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_Pympfr(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_PyFloat(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_Pympz(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_Pympq(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_PyLong(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   Pympc_From_PyStr(PyObject *s, int base, mpfr_prec_t rbits, mpfr_prec_t ibits);
static MPC_Object *   Pympc_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec);
static MPC_Object *   GMPy_MPC_From_Complex_Temp(PyObject* obj, GMPyContextObject *context);

static PyObject *      Pympc_To_PyFloat(PyObject *self);
static PyObject *      Pympc_To_PyLong(PyObject *self);
static PyObject *      Pympc_To_PyStr(MPC_Object *self, int base, int digits);
static PyObject *      Pympc_To_PyComplex(PyObject *self, PyObject *other);
/* support str() and repr() */
static PyObject *      Pympc_To_Str(MPC_Object *self);
static PyObject *      Pympc_To_Repr(MPC_Object *self);

/* Miscellaneous */
static PyObject * raw_mpfr_ascii(mpfr_t self, int base, int digits, int round);
int Pympc_convert_arg(PyObject *arg, PyObject **ptr);

#ifdef __cplusplus
}
#endif
#endif
