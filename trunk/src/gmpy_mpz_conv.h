/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_conv.h                                                         *
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

#ifndef GMPY_MPZ_CONV_H
#define GMPY_MPZ_CONV_H

#ifdef __cplusplus
extern "C" {
#endif

static int isFraction(PyObject* obj);
static int isDecimal(PyObject* obj);
static int isInteger(PyObject* obj);

static PympzObject * Pyxmpz2Pympz(PyObject *self);
#ifdef PY2
static PympzObject * PyInt2Pympz(PyObject *self);
#endif
static PympzObject * PyFloat2Pympz(PyObject *self);
static PympzObject * PyLong2Pympz(PyObject * obj);
static PyxmpzObject * PyLong2Pyxmpz(PyObject * obj);
static int mpz_set_PyStr(mpz_ptr z, PyObject *s, long base);
static PympzObject * PyStr2Pympz(PyObject *s, long base);
static PyObject * Pympz2PyLong(PympzObject *self);
static PyObject * Pympz_To_Integer(PympzObject *self);
static PyObject * Pympz2PyFloat(PympzObject *self);
static PyObject * mpz2binary(mpz_ptr z);
static PyObject * Pympz2binary(PympzObject *self);
static PyObject * mpz_ascii(mpz_t z, int base, int option);
static PyObject * Pympz_ascii(PympzObject *self, int base, int option);
static PympzObject * anynum2Pympz(PyObject* obj);
static PympzObject * Pympz_From_Integer(PyObject* obj);
static PyObject * Pympz2str(PympzObject *self);
static PyObject * Pympz2repr(PympzObject *self);
int Pympz_convert_arg(PyObject *arg, PyObject **ptr);

/* Should only be used by MPFR/MPC related code. */
static long clong_From_Integer(PyObject *obj);

/* Should be used by all MPIR/GMP related code. */
static gmp_si gmp_si_From_Integer(PyObject *obj);
static gmp_ui gmp_ui_From_Integer(PyObject *obj);

static Py_ssize_t ssize_t_From_Integer(PyObject *obj);

static PyxmpzObject * Pyxmpz2Pyxmpz(PyObject *self);
static PyxmpzObject * Pympz2Pyxmpz(PyObject *self);
#ifdef PY2
static PyxmpzObject * PyInt2Pyxmpz(PyObject *self);
#endif
static PyxmpzObject * PyFloat2Pyxmpz(PyObject *self);
static PyxmpzObject * PyStr2Pyxmpz(PyObject *s, long base);
static PyObject * Pyxmpz2PyLong(PyxmpzObject *self);
static PyObject * Pyxmpz_To_Integer(PyxmpzObject *self);
static PyObject * Pyxmpz2binary(PyxmpzObject *self);
static PyObject * xmpz_ascii(mpz_t z, int base, int option);
static PyObject * Pyxmpz_ascii(PyxmpzObject *self, int base, int option);
static PyxmpzObject * anynum2Pyxmpz(PyObject* obj);
static PyObject * Pyxmpz2str(PyxmpzObject *self);
static PyObject * Pyxmpz2repr(PyxmpzObject *self);

#ifdef __cplusplus
}
#endif
#endif



