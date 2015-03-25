/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_context.h                                                         *
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

#define TRAP_NONE      0
#define TRAP_UNDERFLOW 1
#define TRAP_OVERFLOW  2
#define TRAP_INEXACT   4
#define TRAP_INVALID   8
#define TRAP_ERANGE    16
#define TRAP_DIVZERO   32

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
    int traps;               /* if 0, do not trap any exceptions */
                             /* if not 0, then raise traps per bits above  */
    mpfr_prec_t real_prec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t imag_prec;   /* current precision in bits, for Im(MPC) */
    mpfr_rnd_t real_round;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t imag_round;   /* current rounding mode for Im(MPC) */
    int allow_complex;       /* if 1, allow mpfr functions to return an mpc */
    int rational_division;   /* if 1, mpz/mpz returns an mpq result */
    int guard_bits;          /* number of additional guard bits to use */
                             /*   must be less than MAX_GUARD_BITS     */
    int convert_exact;       /* if 1, str -> mpfr via mpq */
    int mpfr_divmod_exact;   /* if 1, divmod(mpfr, mpfr) uses mpq */
} gmpy_context;

typedef struct {
    PyObject_HEAD
    gmpy_context ctx;
#ifndef WITHOUT_THREADS
    PyThreadState *tstate;
#endif
} CTXT_Object;

typedef struct {
    PyObject_HEAD
    CTXT_Object *new_context; /* Context that will be returned when
                               * __enter__ is called. */
    CTXT_Object *old_context; /* Context that will restored when
                               * __exit__ is called. */
} CTXT_Manager_Object;


static PyTypeObject CTXT_Type;
static PyTypeObject CTXT_Manager_Type;

#ifdef WITHOUT_THREADS
#define CURRENT_CONTEXT(obj) obj = module_context;
#else
#define CURRENT_CONTEXT(obj) obj = GMPy_current_context();
#endif

#define CHECK_CONTEXT(context) \
    if (!context) CURRENT_CONTEXT(context);

#define CTXT_Check(v) (((PyObject*)v)->ob_type == &CTXT_Type)
#define CTXT_Manager_Check(v) (((PyObject*)v)->ob_type == &CTXT_Manager_Type)

#define GET_MPFR_PREC(c) (c->ctx.mpfr_prec)
#define GET_REAL_PREC(c) ((c->ctx.real_prec==GMPY_DEFAULT)?GET_MPFR_PREC(c):c->ctx.real_prec)
#define GET_IMAG_PREC(c) ((c->ctx.imag_prec==GMPY_DEFAULT)?GET_REAL_PREC(c):c->ctx.imag_prec)
#define GET_MPFR_ROUND(c) (c->ctx.mpfr_round)
#define GET_REAL_ROUND(c) ((c->ctx.real_round==GMPY_DEFAULT)?GET_MPFR_ROUND(c):c->ctx.real_round)
#define GET_IMAG_ROUND(c) ((c->ctx.imag_round==GMPY_DEFAULT)?GET_REAL_ROUND(c):c->ctx.imag_round)
#define GET_MPC_ROUND(c) (MPC_RND(GET_REAL_ROUND(c), GET_IMAG_ROUND(c)))

#define GET_DIV_MODE(c) (c->ctx.rational_division)
#define GET_GUARD_BITS(c) (c->ctx.guard_bits)
#define MAX_GUARD_BITS 1000

#define GET_CONVERT_EXACT(c) (c->ctx.convert_exact)
#define GET_DIVMOD_EXACT(c) (c->ctx.mpfr_divmod_exact)

static PyObject *    GMPy_CTXT_Manager_New(void);
static void          GMPy_CTXT_Manager_Dealloc(CTXT_Manager_Object *self);
static PyObject *    GMPy_CTXT_Manager_Repr_Slot(CTXT_Manager_Object *self);
static PyObject *    GMPy_CTXT_Manager_Enter(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Manager_Exit(PyObject *self, PyObject *args);

static PyObject *    GMPy_CTXT_New(void);
static void          GMPy_CTXT_Dealloc(CTXT_Object *self);
static PyObject *    GMPy_CTXT_Repr_Slot(CTXT_Object *self);
static PyObject *    GMPy_CTXT_Get(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Local(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *    GMPy_CTXT_Context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *    GMPy_CTXT_Set(PyObject *self, PyObject *other);
static PyObject *    GMPy_CTXT_Clear_Flags(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Copy(PyObject *self, PyObject *other);
static PyObject *    GMPy_CTXT_ieee(PyObject *self, PyObject *other);
static PyObject *    GMPy_CTXT_Enter(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Exit(PyObject *self, PyObject *args);

#ifndef WITHOUT_THREADS
static CTXT_Object * GMPy_current_context(void);
#endif

#ifdef __cplusplus
}
#endif
#endif

