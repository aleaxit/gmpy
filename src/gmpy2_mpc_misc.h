/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpc_misc.h                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018 Case Van Horsen                        *
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

#ifndef GMPY_MPC_MISC_H
#define GMPY_MPC_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPy_Complex_Phase(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Phase(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Phase(PyObject *self, PyObject *other);

static PyObject * GMPy_Complex_Norm(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Norm(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Norm(PyObject *self, PyObject *other);

static PyObject * GMPy_Complex_Polar(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Polar(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Polar(PyObject *self, PyObject *other);

static PyObject * GMPy_Complex_Rect(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Number_Rect(PyObject *x, PyObject *y, CTXT_Object *context);
static PyObject * GMPy_Context_Rect(PyObject *self, PyObject *args);

static PyObject * GMPy_Complex_Proj(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Number_Proj(PyObject *x, CTXT_Object *context);
static PyObject * GMPy_Context_Proj(PyObject *self, PyObject *other);

static PyObject * GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure);
static int        GMPy_MPC_NonZero_Slot(MPC_Object *self);
static PyObject * GMPy_MPC_SizeOf_Method(PyObject *self, PyObject *other);

static PyObject * GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure);
static int        GMPy_MPC_NonZero_Slot(MPC_Object *self);

#ifdef __cplusplus
}
#endif
#endif
