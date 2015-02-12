/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_random.h                                                           *
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

#ifndef GMPY_RANDOM_H
#define GMPY_RANDOM_H

#ifdef __cplusplus
extern "C" {
#endif

/* gmpy_random C API extension header file.
 *
 * Provide support random number state.
 *
 * Version 2.00, December 2011 (created) casevh
 *
 * This file is expected to be included from gmpy.h
 */

typedef struct {
    PyObject_HEAD
    gmp_randstate_t state;
} GMPYRandomStateObject;

static PyTypeObject GMPYRandomState_Type;
#define PyObj_AS_STATE(obj) (((GMPYRandomStateObject *)(obj))->state)
#define GMPYRandomState_Check(v) (((PyObject*)v)->ob_type == &GMPYRandomState_Type)

static GMPYRandomStateObject * GMPYRandomState_New(void);
static void GMPYRandomState_Dealloc(GMPYRandomStateObject *self);
static PyObject * GMPYRandomState_Repr(GMPYRandomStateObject *self);
static PyObject * GMPY_random_state(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_urandomb(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_rrandomb(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_random(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
