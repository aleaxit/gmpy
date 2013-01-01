/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_context.h                                                          *
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

#ifndef GMPY_CONTEXT_H
#define GMPY_CONTEXT_H

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    mpfr_prec_t mpfr_prec;   /* current precision in bits, for MPFR */
    mpfr_rnd_t mpfr_round;   /* current rounding mode for float (MPFR) */
    mpfr_exp_t emax;         /* maximum exponent */
    mpfr_exp_t emin;         /* minimum exponent */
    int subnormalize;        /* if 1, subnormalization is performed */
    int underflow;           /* did an underflow occur? */
    int overflow;            /* did an overflow occur? */
    int inexact;             /* was the result inexact? */
    int invalid;             /* invalid operation (i.e. NaN)? */
    int erange;              /* did a range error occur? */
    int divzero;             /* divided by zero? */
    int trap_underflow;      /* if 1, raise exception for underflow */
    int trap_overflow;       /* if 1, raise exception for overflow */
    int trap_inexact;        /* if 1, raise exception for inexact */
    int trap_invalid;        /* if 1, raise exception for invalid (NaN) */
    int trap_erange;         /* if 1, raise exception for range error */
    int trap_divzero;        /* if 1, raise exception for divide by zero */
    int trap_expbound;       /* if 1, raise exception if mpfr/mpc exponents */
                             /*       are out of bounds */
#ifdef WITHMPC
    mpfr_prec_t real_prec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t imag_prec;   /* current precision in bits, for Im(MPC) */
    mpfr_rnd_t real_round;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t imag_round;   /* current rounding mode for Im(MPC) */
    int allow_complex;       /* if 1, allow mpfr functions to return an mpc */
#endif
} gmpy_context;

typedef struct {
    PyObject_HEAD
    gmpy_context ctx;
} GMPyContextObject;

typedef struct {
    PyObject_HEAD
    gmpy_context new_ctx;    /* Context that will be returned when __enter__
                              * is called. */
    gmpy_context old_ctx;    /* Context that will restored when __exit__ is
                              * is called. */
} GMPyContextManagerObject;


static PyTypeObject GMPyContext_Type;
static PyTypeObject GMPyContextManager_Type;

#define GMPyContext_Check(v) (((PyObject*)v)->ob_type == &GMPyContext_Type)
#define GMPyContextManager_Check(v) (((PyObject*)v)->ob_type == &GMPyContextManager_Type)

#define GET_MPFR_PREC(c) (c->ctx.mpfr_prec)
#define GET_REAL_PREC(c) ((c->ctx.real_prec==GMPY_DEFAULT)?GET_MPFR_PREC(c):c->ctx.real_prec)
#define GET_IMAG_PREC(c) ((c->ctx.imag_prec==GMPY_DEFAULT)?GET_REAL_PREC(c):c->ctx.imag_prec)
#define GET_MPFR_ROUND(c) (c->ctx.mpfr_round)
#define GET_REAL_ROUND(c) ((c->ctx.real_round==GMPY_DEFAULT)?GET_MPFR_ROUND(c):c->ctx.real_round)
#define GET_IMAG_ROUND(c) ((c->ctx.imag_round==GMPY_DEFAULT)?GET_REAL_ROUND(c):c->ctx.imag_round)
#define GET_MPC_ROUND(c) (MPC_RND(GET_REAL_ROUND(c), GET_IMAG_ROUND(c)))

static PyObject * GMPyContextManager_new(void);
static void GMPyContextManager_dealloc(GMPyContextManagerObject *self);
static PyObject * GMPyContextManager_repr(GMPyContextManagerObject *self);
static PyObject * GMPyContextManager_enter(PyObject *self, PyObject *args);
static PyObject * GMPyContextManager_exit(PyObject *self, PyObject *args);

static PyObject * GMPyContext_new(void);
static void GMPyContext_dealloc(GMPyContextObject *self);
static PyObject * GMPyContext_repr(GMPyContextObject *self);
static PyObject * GMPyContext_get_context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * GMPyContext_local_context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * GMPyContext_context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * GMPyContext_set_context(PyObject *self, PyObject *other);
static PyObject * GMPyContext_clear_flags(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif

