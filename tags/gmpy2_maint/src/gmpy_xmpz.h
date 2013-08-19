/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_xmpz.h                                                             *
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

#ifndef GMPY_XMPZ_H
#define GMPY_XMPZ_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    PyObject_HEAD
    mpz_t z;
} PyxmpzObject;

#define Pyxmpz_AS_MPZ(obj) (((PyxmpzObject *)(obj))->z)

static PyTypeObject Pyxmpz_Type;

#define Pyxmpz_Check(v) (((PyObject*)v)->ob_type == &Pyxmpz_Type)

#define CHECK_MPZANY(v) (Pympz_Check(v) || Pyxmpz_Check(v))

typedef struct {
    PyObject_HEAD
    PyxmpzObject *bitmap;
    Py_ssize_t start, stop;
    int iter_type;
} GMPYIterObject;

static PyTypeObject GMPYIter_Type;

#define GMPyIter_Check(v) (((PyObject*)v)->ob_type == &GMPYIter_Type)

static PyObject * Pygmpy_xmpz(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * Pyxmpz_digits(PyObject *self, PyObject *args);
static PyObject * Pyxmpz_xbit_mask(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_abs(PyxmpzObject *x);
static PyObject * Pyxmpz_neg(PyxmpzObject *x);
static PyObject * Pyxmpz_pos(PyxmpzObject *x);
static int Pyxmpz_nonzero(PyxmpzObject *x);
static PyObject * Pyxmpz_com(PyxmpzObject *x);
#if PY_MAJOR_VERSION < 3
/* hex/oct formatting (mpz-only) */
static PyObject * Pyxmpz_oct(PyxmpzObject *self);
static PyObject * Pyxmpz_hex(PyxmpzObject *self);
#endif
static PyObject * Pyxmpz_make_mpz(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_copy(PyObject *self, PyObject *other);
static Py_ssize_t Pyxmpz_nbits(PyxmpzObject *obj);
static PyObject * Pyxmpz_subscript(PyxmpzObject* self, PyObject* item);
static int Pyxmpz_assign_subscript(PyxmpzObject* self, PyObject* item, PyObject* value);

#ifdef __cplusplus
}
#endif
#endif
