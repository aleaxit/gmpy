/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_context.h                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen                  *
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
#ifdef WITHMPC
    mpfr_prec_t mpc_rprec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t mpc_iprec;   /* current precision in bits, for Im(MPC) */
#endif
    mpfr_rnd_t mpfr_round;   /* current rounding mode for float (MPFR) */
#ifdef WITHMPC
    mpfr_rnd_t mpc_rround;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t mpc_iround;   /* current rounding mode for Im(MPC) */
#endif
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
    int allow_complex;       /* if 1, allow mpfr functions to return an mpc */
#endif
} gmpy_context;

typedef struct {
    PyObject_HEAD
    gmpy_context now;        /* The "new" values, used by __enter__ */
    PyObject *orig;          /* Original context, restored by __exit__*/
} GMPyContextObject;

static PyTypeObject GMPyContext_Type;

#define GMPyContext_Check(v) (((PyObject*)v)->ob_type == &GMPyContext_Type)

#define GET_MPFR_PREC(c) (c->now.mpfr_prec)
#define GET_MPC_RPREC(c) ((c->now.mpc_rprec==GMPY_DEFAULT)?GET_MPFR_PREC(c):c->now.mpc_rprec)
#define GET_MPC_IPREC(c) ((c->now.mpc_iprec==GMPY_DEFAULT)?GET_MPC_RPREC(c):c->now.mpc_iprec)
#define GET_MPFR_ROUND(c) (c->now.mpfr_round)
#define GET_MPC_RROUND(c) ((c->now.mpc_rround==GMPY_DEFAULT)?GET_MPFR_ROUND(c):c->now.mpc_rround)
#define GET_MPC_IROUND(c) ((c->now.mpc_iround==GMPY_DEFAULT)?GET_MPC_RROUND(c):c->now.mpc_iround)
#define GET_MPC_ROUND(c) (RNDC(GET_MPC_RROUND(c), GET_MPC_IROUND(c)))

static GMPyContextObject * GMPyContext_new(void);
static void GMPyContext_dealloc(GMPyContextObject *self);
static PyObject * GMPyContext_repr(GMPyContextObject *self);
static PyObject * Pygmpy_get_context(PyObject *self, PyObject *args);
static PyObject * Pygmpy_local_context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * Pygmpy_context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * Pygmpy_set_context(PyObject *self, PyObject *other);
static PyObject * GMPyContext_enter(PyObject *self, PyObject *args);
static PyObject * GMPyContext_exit(PyObject *self, PyObject *args);
static PyObject * GMPyContext_clear_flags(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif

