/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz.h                                                              *
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

#ifndef GMPY_MPZ_H
#define GMPY_MPZ_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    PyObject_HEAD
    mpz_t z;
    Py_hash_t hash_cache;
} PympzObject;

#define Pympz_AS_MPZ(obj) (((PympzObject *)(obj))->z)

static PyTypeObject Pympz_Type;

#define Pympz_Check(v) (((PyObject*)v)->ob_type == &Pympz_Type)

static PyObject * Pygmpy_mpz(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * Pympz_digits(PyObject *self, PyObject *args);
static PyObject * Pympz_num_digits(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_length(PyObject *self, PyObject *other);
static PyObject * Pympz_bit_mask(PyObject *self, PyObject *other);
static PyObject * Pympz_bit_scan0(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_scan1(PyObject *self, PyObject *args);
static PyObject * Pympz_popcount(PyObject *self, PyObject *other);
static PyObject * Pygmpy_bit_test(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_test(PyObject *self, PyObject *other);
static PyObject * Pygmpy_bit_clear(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_clear(PyObject *self, PyObject *other);
static PyObject * Pygmpy_bit_set(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_set(PyObject *self, PyObject *other);
static PyObject * Pygmpy_bit_flip(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_flip(PyObject *self, PyObject *other);
static PyObject * Pympz_iroot(PyObject *self, PyObject *args);
static PyObject * Pympz_iroot_rem(PyObject *self, PyObject *args);
static PyObject * Pympz_sign(PyObject *self, PyObject *other);
static PyObject * Pympz_abs(PympzObject *self);
static PyObject * Pympz_neg(PympzObject *self);
static PyObject * Pympz_pos(PympzObject *self);
static PyObject * Pympz_ceil(PyObject *self, PyObject *other);
static PyObject * Pympz_floor(PyObject *self, PyObject *other);
static PyObject * Pympz_round(PyObject *self, PyObject *other);
static PyObject * Pympz_trunc(PyObject *self, PyObject *other);
static PyObject * Pympz_square(PyObject *self, PyObject *other);
static PyObject * Pympz_pow(PyObject *b, PyObject *e, PyObject *m);
static int Pympz_nonzero(PympzObject *self);
static PyObject * Pympz_com(PympzObject *self);
static PyObject * Pympz_and(PyObject *self, PyObject *other);
static PyObject * Pympz_ior(PyObject *self, PyObject *other);
static PyObject * Pympz_xor(PyObject *self, PyObject *other);
static PyObject * Pympz_rshift(PyObject *self, PyObject *other);
static PyObject * Pympz_lshift(PyObject *self, PyObject *other);
#if PY_MAJOR_VERSION < 3
/* hex/oct formatting (mpz-only) */
static PyObject * Pympz_oct(PympzObject *self);
static PyObject * Pympz_hex(PympzObject *self);
#endif
static Py_hash_t Pympz_hash(PympzObject *self);
static PyObject * Pygmpy_gcd(PyObject *self, PyObject *args);
static PyObject * Pygmpy_lcm(PyObject *self, PyObject *args);
static PyObject * Pygmpy_gcdext(PyObject *self, PyObject *args);
static PyObject * Pygmpy_divm(PyObject *self, PyObject *args);
static PyObject * Pygmpy_fac(PyObject *self, PyObject *other);
static PyObject * Pygmpy_fib(PyObject *self, PyObject *other);
static PyObject * Pygmpy_fib2(PyObject *self, PyObject *other);
static PyObject * Pygmpy_lucas(PyObject *self, PyObject *other);
static PyObject * Pygmpy_lucas2(PyObject *self, PyObject *other);
static PyObject * Pympz_bincoef(PyObject *self, PyObject *args);
static PyObject * Pympz_isqrt(PyObject *self, PyObject *other);
static PyObject * Pympz_isqrt_rem(PyObject *self, PyObject *args);
static PyObject * Pympz_remove(PyObject *self, PyObject *args);
static PyObject * Pygmpy_invert(PyObject *self, PyObject *args);
static PyObject * Pympz_hamdist(PyObject *self, PyObject *args);
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
static Py_ssize_t Pympz_nbits(PympzObject *self);
static PyObject * Pympz_subscript(PympzObject *self, PyObject *item);
static PyObject * Pympz_format(PyObject *self, PyObject *args);
static PyObject * Pympz_add(PyObject *self, PyObject *args);
static PyObject * Pympz_sub(PyObject *self, PyObject *args);
static PyObject * Pympz_mul(PyObject *self, PyObject *args);
static PyObject * Pympz_div(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
