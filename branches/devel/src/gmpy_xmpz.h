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
} XMPZ_Object;

static PyTypeObject XMPZ_Type;

#define XMPZ_Check(v) (((PyObject*)v)->ob_type == &XMPZ_Type)

#define CHECK_MPZANY(v) (MPZ_Check(v) || XMPZ_Check(v))

typedef struct {
    PyObject_HEAD
    XMPZ_Object *bitmap;
    Py_ssize_t start, stop;
    int iter_type;
} GMPyIterObject;

static PyTypeObject GMPyIter_Type;

#define GMPyIter_Check(v) (((PyObject*)v)->ob_type == &GMPyIter_Type)

static PyObject * Pygmpy_xmpz(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * Pyxmpz_digits(PyObject *self, PyObject *args);
static PyObject * Pyxmpz_xbit_mask(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_abs(XMPZ_Object *x);
static PyObject * Pyxmpz_neg(XMPZ_Object *x);
static PyObject * Pyxmpz_pos(XMPZ_Object *x);
static int Pyxmpz_nonzero(XMPZ_Object *x);
static PyObject * Pyxmpz_com(XMPZ_Object *x);
#if PY_MAJOR_VERSION < 3
/* hex/oct formatting (mpz-only) */
static PyObject * Pyxmpz_oct(XMPZ_Object *self);
static PyObject * Pyxmpz_hex(XMPZ_Object *self);
#endif
static PyObject * Pyxmpz_make_mpz(PyObject *self, PyObject *other);
static PyObject * Pyxmpz_copy(PyObject *self, PyObject *other);
static Py_ssize_t Pyxmpz_nbits(XMPZ_Object *obj);
static PyObject * Pyxmpz_subscript(XMPZ_Object* self, PyObject* item);
static int Pyxmpz_assign_subscript(XMPZ_Object* self, PyObject* item, PyObject* value);

#ifdef __cplusplus
}
#endif
#endif
