/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz.h                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen             *
 *                                                                         *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2.1 of the License, or  *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This library is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * Lesser General Public License for more details.                         *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with this library; if not, write to the Free Software     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA           *
 * 02110-1301  USA                                                         *
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

typedef struct {
    PyObject_HEAD
    mpz_t z;
} PyxmpzObject;

#define Pympz_AS_MPZ(obj) (((PympzObject *)(obj))->z)

#define Pyxmpz_AS_MPZ(obj) (((PyxmpzObject *)(obj))->z)

static PyTypeObject Pympz_Type;

#define Pympz_Check(v) (((PyObject*)v)->ob_type == &Pympz_Type)

static PyTypeObject Pyxmpz_Type;

#define Pyxmpz_Check(v) (((PyObject*)v)->ob_type == &Pyxmpz_Type)

#define CHECK_MPZANY(v) (Pympz_Check(v) || Pyxmpz_Check(v))

static PyObject * Pygmpy_mpz(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * Pygmpy_xmpz(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * Pympz_digits(PyObject *self, PyObject *args);
static PyObject * Pyxmpz_digits(PyObject *self, PyObject *args);
static PyObject * Pympz_numdigits(PyObject *self, PyObject *args);
static PyObject * Pympz_bit_length(PyObject *self, PyObject *other);
static PyObject * Pympz_bit_mask(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_xbit_mask(PyObject *self, PyObject *other);
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
static PyObject * Pympz_abs(PympzObject *x);
static PyObject * Pyxmpz_abs(PyxmpzObject *x);
static PyObject * Pympz_neg(PympzObject *x);
static PyObject * Pyxmpz_neg(PyxmpzObject *x);
static PyObject * Pympz_pos(PympzObject *x);
static PyObject * Pyxmpz_pos(PyxmpzObject *x);
static PyObject * Pympz_square(PyObject *self, PyObject *other);
static PyObject * Pympz_pow(PyObject *b, PyObject *e, PyObject *m);
static int Pympz_nonzero(PympzObject *x);
static int Pyxmpz_nonzero(PyxmpzObject *x);
static PyObject * Pympz_com(PympzObject *x);
static PyObject * Pyxmpz_com(PyxmpzObject *x);
static PyObject * Pympz_and(PyObject *a, PyObject *b);
static PyObject * Pympz_ior(PyObject *a, PyObject *b);
static PyObject * Pympz_xor(PyObject *a, PyObject *b);
static PyObject * Pympz_rshift(PyObject *a, PyObject *b);
static PyObject * Pympz_lshift(PyObject *a, PyObject *b);
#if PY_MAJOR_VERSION < 3
/* hex/oct formatting (mpz-only) */
static PyObject * Pympz_oct(PympzObject *self);
static PyObject * Pyxmpz_oct(PyxmpzObject *self);
static PyObject * Pympz_hex(PympzObject *self);
static PyObject * Pyxmpz_hex(PyxmpzObject *self);
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
static PyObject * Pyxmpz_make_mpz(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_copy(PyObject *self, PyObject *other);
static Py_ssize_t Pyxmpz_nbits(PyxmpzObject *obj);
static PyObject * Pyxmpz_subscript(PyxmpzObject* self, PyObject* item);
static int Pyxmpz_assign_subscript(PyxmpzObject* self, PyObject* item, PyObject* value);
static Py_ssize_t Pympz_nbits(PyxmpzObject *obj);
static PyObject * Pympz_subscript(PyxmpzObject* self, PyObject* item);
static PyObject * Pympz_format(PyObject *self, PyObject *args);
static PyObject * Pympz_add(PyObject *self, PyObject *args);
static PyObject * Pympz_sub(PyObject *self, PyObject *args);
static PyObject * Pympz_mul(PyObject *self, PyObject *args);
static PyObject * Pympz_div(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
