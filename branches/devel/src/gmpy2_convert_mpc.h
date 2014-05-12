/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_mpc.h                                                     *
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

#ifndef GMPY2_CONVERT_MPC_H
#define GMPY2_CONVERT_MPC_H

#ifdef __cplusplus
extern "C" {
#endif

/* The following functions identify and classify the numeric types that are
 * supported by gmpy2.
 *
 * These checks are currently implemented as functions but may be
 * implemented as macros in the future.
 */

    /* Conversions with Pympc */
static MPC_Object *   GMPy_MPC_From_MPC(MPC_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyComplex(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_MPFR(MPFR_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyFloat(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_MPZ(MPZ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_MPQ(MPQ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_Fraction(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_Decimal(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyIntOrLong(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_PyStr(PyObject *s, int base, mpfr_prec_t rbits, mpfr_prec_t ibits, CTXT_Object *context);
static MPC_Object *   GMPy_MPC_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);

static PyObject *     GMPy_MPC_Float_Slot(PyObject *self);
#ifdef PY2
static PyObject *     GMPy_MPC_Long_Slot(PyObject *self);
#endif
static PyObject *     GMPy_MPC_Int_Slot(PyObject *self);
static PyObject *     GMPy_PyStr_From_MPC(MPC_Object *self, int base, int digits, CTXT_Object *context);
static PyObject *     GMPy_PyComplex_From_MPC(PyObject *self, PyObject *other);
/* support str() and repr() */
static PyObject *     GMPy_MPC_Str_Slot(MPC_Object *self);
static PyObject *     GMPy_MPC_Repr_Slot(MPC_Object *self);

/* Miscellaneous */
#if 0
static int            GMPy_MPC_convert_arg(PyObject *arg, PyObject **ptr);
#endif

#ifdef __cplusplus
}
#endif
#endif
