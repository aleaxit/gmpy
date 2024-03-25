/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_add.h                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2024 Case Van Horsen                                   *
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

#ifndef GMPY2_ADD_H
#define GMPY2_ADD_H

#ifdef __cplusplus
extern "C" {
#endif

/* Public API */

static PyObject * GMPy_Number_Add(PyObject *x, PyObject *y, CTXT_Object *context);

/* Private API */

static PyObject * GMPy_Integer_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context);
static PyObject * GMPy_Rational_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context);
static PyObject * GMPy_Real_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context);
static PyObject * GMPy_Complex_AddWithType(PyObject *x, int xtype, PyObject *y, int ytype, CTXT_Object *context);
static PyObject * GMPy_Number_Add_Slot(PyObject *x, PyObject *y);
static PyObject * GMPy_Context_Add(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
