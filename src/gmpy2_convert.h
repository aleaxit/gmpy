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
 *           2015, 2016, 2017, 2018, 2019, 2020 Case Van Horsen            *
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

#ifndef GMPY2_CONVERT_H
#define GMPY2_CONVERT_H

#ifdef __cplusplus
extern "C" {
#endif

/* The following macros classify the numeric types that are supported by
 * gmpy2.
 */
#define HAS_MPZ_CONVERSION(x) (PyObject_HasAttrString(x, "__mpz__"))
#define HAS_MPQ_CONVERSION(x) (PyObject_HasAttrString(x, "__mpq__"))
#define HAS_MPFR_CONVERSION(x) (PyObject_HasAttrString(x, "__mpfr__"))
#define HAS_MPC_CONVERSION(x) (PyObject_HasAttrString(x, "__mpc__"))

#define HAS_STRICT_MPZ_CONVERSION(x) (HAS_MPZ_CONVERSION(x) && \
                                     !HAS_MPQ_CONVERSION(x))
#define HAS_STRICT_MPFR_CONVERSION(x) (HAS_MPFR_CONVERSION(x) && \
                                      !HAS_MPC_CONVERSION(x))

#define IS_FRACTION(x) (!strcmp(Py_TYPE(x)->tp_name, "Fraction"))

#define IS_RATIONAL_ONLY(x) (MPQ_Check(x) || IS_FRACTION(x) || HAS_MPQ_CONVERSION(x))

#ifdef PY2
#define IS_INTEGER(x) (MPZ_Check(x) || PyInt_Check(x) || \
                        PyLong_Check(x) || XMPZ_Check(x) || \
                        HAS_STRICT_MPZ_CONVERSION(x))
#define IS_RATIONAL(x) (MPQ_Check(x) || IS_FRACTION(x) || \
                        MPZ_Check(x) || PyInt_Check(x) || PyLong_Check(x) || \
                        XMPZ_Check(x) || HAS_MPQ_CONVERSION(x) || \
                        HAS_MPZ_CONVERSION(x))
#else
#define IS_INTEGER(x) (MPZ_Check(x) || PyLong_Check(x) || \
                        XMPZ_Check(x) || HAS_STRICT_MPZ_CONVERSION(x))
#define IS_RATIONAL(x) (MPQ_Check(x) || IS_FRACTION(x) || \
                        MPZ_Check(x) || PyLong_Check(x) || \
                        XMPZ_Check(x) || HAS_MPQ_CONVERSION(x) || \
                        HAS_MPZ_CONVERSION(x))
#endif

#define IS_REAL_ONLY(x) (MPFR_Check(x) || PyFloat_Check(x) || HAS_STRICT_MPFR_CONVERSION(x))
#define IS_REAL(x) (IS_RATIONAL(x) || IS_REAL_ONLY(x))

#define IS_COMPLEX_ONLY(x) (MPC_Check(x) || PyComplex_Check(x) || HAS_MPC_CONVERSION(x))
#define IS_COMPLEX(x) (IS_REAL(x) || IS_COMPLEX_ONLY(x))

/* Define constants used in gmpy2_convert.c->GMPy_ObjectType. */

#define OBJ_TYPE_MPZ        1
#define OBJ_TYPE_PyInteger  2
/* 3 TO 14 reserved for additional integer types. */
#define OBJ_TYPE_HAS_MPZ    15

#define OBJ_TYPE_MPQ        16
#define OBJ_TYPE_PyFraction 17
/* 18 to 30 reserved for additional rational types. */
#define OBJ_TYPE_HAS_MPQ    31

#define OBJ_TYPE_MPFR       32
#define OBJ_TYPE_PyFloat    33
/* 34 to 46 reserved for additional real types. */
#define OBJ_TYPE_HAS_MPFR   47

#define OBJ_TYPE_MPC        48
#define OBJ_TYPE_PyComplex  49
/* 50 to 62 reserved for additional complex types. */
#define OBJ_TYPE_HAS_MPC    63

#define OBJ_TYPE_MAX        64

/* Since the macros are used in gmpy2's codebase, these functions are skipped
 * until they are needed for the C API in the future.
 */

#if 0
/* Checks for mpz, xmpz, and the integer types included with Python. */
static int GMPy_isInteger(PyObject *obj);

/* Checks for the Fraction type included with Python. */
static int GMPy_isFraction(PyObject *obj);

/* Combined mpq, isInteger() and isFraction() check. */
static int GMPy_isRational(PyObject *obj);

/* Combined mpfr, PyFloat, and isRational() check. */
static int GMPy_isReal(PyObject *obj);

/* Combined mpc, PyComplex, and isReal() check. */
static int GMPy_isComplex(PyObject *obj);
#endif

/* ======== C helper routines ======== */
static int             mpz_set_PyStr(mpz_ptr z, PyObject *s, int base);
static PyObject *      mpz_ascii(mpz_t z, int base, int option, int which);

#ifdef __cplusplus
}
#endif
#endif
