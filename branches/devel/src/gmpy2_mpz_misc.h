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


static PyObject * GMPy_MPZ_NumDigits(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Iroot(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_IrootRem(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Bincoef(PyObject *self, PyObject *args);
static PyObject * GMPy_MPZ_Method_Ceil(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_Floor(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_Trunc(PyObject *self, PyObject *other);
static PyObject * GMPy_MPZ_Method_Round(PyObject *self, PyObject *other);
static int        GMPy_MPZ_NonZero_Slot(MPZ_Object *self);
static PyObject * GMPy_MPZ_Attrib_GetNumer(MPZ_Object *self, void *closure);
static PyObject * GMPy_MPZ_Attrib_GetDenom(MPQ_Object *self, void *closure);
static PyObject * GMPy_MPZ_Function_GCD(PyObject *self, PyObject *args);

static PyObject * Pygmpy_lcm(PyObject *self, PyObject *args);
static PyObject * Pygmpy_gcdext(PyObject *self, PyObject *args);
static PyObject * Pygmpy_divm(PyObject *self, PyObject *args);
static PyObject * Pygmpy_fac(PyObject *self, PyObject *other);
static PyObject * Pygmpy_fib(PyObject *self, PyObject *other);
static PyObject * Pygmpy_fib2(PyObject *self, PyObject *other);
static PyObject * Pygmpy_lucas(PyObject *self, PyObject *other);
static PyObject * Pygmpy_lucas2(PyObject *self, PyObject *other);
static PyObject * Pympz_isqrt(PyObject *self, PyObject *other);
static PyObject * Pympz_isqrt_rem(PyObject *self, PyObject *args);
static PyObject * Pympz_remove(PyObject *self, PyObject *args);
static PyObject * Pygmpy_invert(PyObject *self, PyObject *args);
static PyObject * Pygmpy_divexact(PyObject *self, PyObject *args);
static PyObject * Pympz_is_square(PyObject *self, PyObject *other);
static PyObject * Pympz_is_power(PyObject *self, PyObject *other);
static PyObject * Pympz_is_prime(PyObject *self, PyObject *args);
static PyObject * Pympz_next_prime(PyObject *self, PyObject *other);
static PyObject * Pympz_jacobi(PyObject *self, PyObject *args);
static PyObject * Pympz_legendre(PyObject *self, PyObject *args);
static PyObject * Pympz_kronecker(PyObject *self, PyObject *args);
static PyObject * Pympz_is_even(PyObject *self, PyObject *other);
static PyObject * Pympz_is_odd(PyObject *self, PyObject *other);
static Py_ssize_t Pympz_nbits(MPZ_Object *self);
static PyObject * Pympz_subscript(MPZ_Object *self, PyObject *item);

#if PY_MAJOR_VERSION < 3
static PyObject * GMPy_MPZ_Oct_Slot(MPZ_Object *self);
static PyObject * GMPy_MPZ_Hex_Slot(MPZ_Object *self);
#endif

#ifdef __cplusplus
}
#endif
#endif
