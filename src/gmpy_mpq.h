/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq.h                                                              *
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

#ifndef GMPY_MPQ_H
#define GMPY_MPQ_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    PyObject_HEAD
    mpq_t q;
    Py_hash_t  hash_cache;
} MPQ_Object;

#define MPQ(obj) (((MPQ_Object *)(obj))->q)

static PyTypeObject MPQ_Type;
#define MPQ_Check(v) (((PyObject*)v)->ob_type == &MPQ_Type)

static PyObject * Pygmpy_mpq(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * Pympq_digits(PyObject *self, PyObject *args);
static PyObject * Pympq_sign(PyObject *self, PyObject *other);
static PyObject * Pympq_numer(PyObject *self, PyObject *args);
static PyObject * Pympq_getnumer(MPQ_Object *self, void *closure);
static PyObject * Pympq_denom(PyObject *self, PyObject *args);
static PyObject * Pympq_getdenom(MPQ_Object *self, void *closure);
static PyObject * Pympq_qdiv(PyObject *self, PyObject *args);
static PyObject * Pympq_ceil(PyObject *self, PyObject *other);
static PyObject * Pympq_floor(PyObject *self, PyObject *other);
static PyObject * Pympq_trunc(PyObject *self, PyObject *other);
static int Pympq_nonzero(MPQ_Object *x);
static Py_hash_t Pympq_hash(MPQ_Object *self);

#ifdef __cplusplus
}
#endif
#endif
