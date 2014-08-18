/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_misc.h                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Case Van Horsen      *
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

#ifndef GMPY_MPZ_MISC_H
#define GMPY_MPZ_MISC_H

#ifdef __cplusplus
extern "C" {
#endif


static int        GMPy_MPZ_NonZero_Slot(MPZ_Object *self);
static PyObject * GMPy_MPZ_Attrib_GetNumer(MPZ_Object *self, void *closure);
static PyObject * GMPy_MPZ_Attrib_GetDenom(MPQ_Object *self, void *closure);
static PyObject * GMPy_MPZ_Method_Ceil(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_Floor(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_Trunc(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_Round(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_NumDigits(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_NumDigits(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Iroot(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_IrootRem(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Bincoef(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_GCD(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_LCM(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_GCDext(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Divm(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Fac(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Primorial(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_DoubleFac(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_MultiFac(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Fib(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Fib2(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Lucas(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Lucas2(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Isqrt(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_IsqrtRem(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Remove(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Invert(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Divexact(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_IsSquare(PyObject *self, PyObject *other);

static PyObject * GMPy_MPZ_Function_IsDivisible(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Method_IsDivisible(PyObject *self, PyObject *other);

static PyObject * GMPy_MPZ_Function_IsCongruent(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Method_IsCongruent(PyObject *self, PyObject *args);

static PyObject * GMPy_MPZ_Function_IsPower(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_IsPrime(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_NextPrime(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_Jacobi(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Legendre(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_Kronecker(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Function_IsEven(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Function_IsOdd(PyObject *self, PyObject *other);
static Py_ssize_t GMPy_MPZ_Method_Length(MPZ_Object *self);
static PyObject * GMPy_MPZ_Method_SubScript(MPZ_Object *self, PyObject *item);

#if PY_MAJOR_VERSION < 3
static PyObject * GMPy_MPZ_Oct_Slot(MPZ_Object *self);
static PyObject * GMPy_MPZ_Hex_Slot(MPZ_Object *self);
#endif

#ifdef __cplusplus
}
#endif
#endif
