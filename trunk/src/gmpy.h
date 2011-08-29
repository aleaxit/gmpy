/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy.h                                                                  *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008 - 2011 Case Van Horsen                              *
 *                                                                         *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2.1 of the License, or  *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This library is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * Lesser General Public License for more details.                         *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with this library; if not, write to the Free Software     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA           *
 * 02110-1301  USA                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/*
  gmpy C API extension header file.
  Part of Python's gmpy module since version 0.4

  Created by Pearu Peterson <pearu@cens.ioc.ee>, November 2000.
  Edited by A. Martelli <aleaxit@yahoo.com>, December 2000.
  Edited by Case Van Horsen <casevh@gmail.com>, 2009, 2010, 2011.

  Version 1.02, February 2007.
  Version 1.03, June 2008
  Version 1.04, June 2008 (no changes)
  Version 1.05, February 2009 (support MPIR)
  Version 1.20, January 2010 (remove obsolete MS hacks) casevh
  Version 2.00, April 2010 (change to gmpy2) casevh
                October 2010 (added Py_hash_t) casevh
                December 2010 (added mpfr, mpc) casevh
                January 2011 (add Pygmpy_context) casevh
                April 2011 (split into multiple files) casevh
 */

#ifndef Py_GMPYMODULE_H
#define Py_GMPYMODULE_H

#ifdef __cplusplus
extern "C" {
#endif

#if PY_VERSION_HEX < 0x02060000
#  error "GMPY2 requires Python 2.6 or later."
#endif

#if PY_VERSION_HEX < 0x030200A4
typedef long Py_hash_t;
typedef unsigned long Py_uhash_t;
#endif

#if defined(MS_WIN32) && defined(_MSC_VER)
   /* so one won't need to link explicitly to gmp.lib...: */
#  if defined(MPIR)
#    pragma comment(lib,"mpir.lib")
#  else
#    pragma comment(lib,"gmp.lib")
#  endif
#  define isnan _isnan
#  define isinf !_finite
#  define USE_ALLOCA 1
#  define inline __inline
#endif

#ifdef MPIR
#  include "mpir.h"
#else
#  include "gmp.h"
#endif

/* Header file for gmpy2 */
typedef struct {
    PyObject_HEAD
    /* PyObject* callable; */
    /* long flags; */
} mpob;

typedef struct {
    mpob ob;
    mpz_t z;
    Py_hash_t hash_cache;
} PympzObject;

typedef struct {
    mpob ob;
    mpq_t q;
    Py_hash_t  hash_cache;
} PympqObject;

typedef struct {
    mpob ob;
    mpz_t z;
} PyxmpzObject;

#ifdef WITHMPFR
#  include "mpfr.h"
#  include "gmpy_mpfr.h"
#endif

#ifdef WITHMPC
#  include "mpc.h"
#  include "gmpy_mpc.h"
#endif

#ifdef __GNUC__
#define USE_ALLOCA 1
#endif

#ifndef alloca
# ifdef __GNUC__
#  define alloca __builtin_alloca
# else
#   ifdef _MSC_VER
#    include <malloc.h>
#    define alloca _alloca
#   else
#    if HAVE_ALLOCA_H
#     include <alloca.h>
#    else
       char *alloca ();
#    endif
#   endif
# endif
#endif

#define ALLOC_THRESHOLD 8192

#ifdef DEBUG
#  define TRACE(msg) if(options.debug) fprintf(stderr, msg)
#else
#  define TRACE(msg)
#endif

#define TYPE_ERROR(msg) PyErr_SetString(PyExc_TypeError, msg)
#define VALUE_ERROR(msg) PyErr_SetString(PyExc_ValueError, msg)
#define ZERO_ERROR(msg) PyErr_SetString(PyExc_ZeroDivisionError, msg)
#define SYSTEM_ERROR(msg) PyErr_SetString(PyExc_SystemError, msg)
#define OVERFLOW_ERROR(msg) PyErr_SetString(PyExc_OverflowError, msg)

#define GMPY_DEFAULT -1

#ifdef USE_ALLOCA
#  define TEMP_ALLOC(B, S) \
    if(S < ALLOC_THRESHOLD) { \
        B = alloca(S); \
    } else { \
        if(!(B = PyMem_Malloc(S))) { \
            PyErr_NoMemory(); \
            return NULL; \
        } \
    }
#  define TEMP_FREE(B, S) if(S >= ALLOC_THRESHOLD) PyMem_Free(B)
#else
#  define TEMP_ALLOC(B, S) \
    if(!(B = PyMem_Malloc(S)))  { \
        PyErr_NoMemory(); \
        return NULL; \
    }
#  define TEMP_FREE(B, S) PyMem_Free(B)
#endif

/* Various defs to mask differences between Python versions. */

#define Py_RETURN_NOTIMPLEMENTED\
    return Py_INCREF(Py_NotImplemented), Py_NotImplemented

#ifndef Py_SIZE
#define Py_SIZE(ob)     (((PyVarObject*)(ob))->ob_size)
#endif

#ifndef Py_TYPE
#define Py_TYPE(ob)     (((PyObject*)(ob))->ob_type)
#endif

#define Pympz_AS_MPZ(obj) (((PympzObject *)(obj))->z)
#define Pympq_AS_MPQ(obj) (((PympqObject *)(obj))->q)
#define Pyxmpz_AS_MPZ(obj) (((PyxmpzObject *)(obj))->z)

static PyTypeObject Pympz_Type;
#define Pympz_Check(v) (((PyObject*)v)->ob_type == &Pympz_Type)

static PyTypeObject Pympq_Type;
#define Pympq_Check(v) (((PyObject*)v)->ob_type == &Pympq_Type)

static PyTypeObject Pyxmpz_Type;
#define Pyxmpz_Check(v) (((PyObject*)v)->ob_type == &Pyxmpz_Type)

#define CHECK_MPZANY(v) (Pympz_Check(v) || Pyxmpz_Check(v))

#ifdef WITHMPFR
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
    mpc_rnd_t mpc_round;     /* current rounding mode for complex (MPC)*/
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
#endif

/* Choose which memory manager is used: Python or C */
#define USE_PYMEM 1
#ifdef USE_PYMEM
#define GMPY_FREE(NAME) PyMem_Free(NAME)
#else
#define GMPY_FREE(NAME) free(NAME)
#endif

#ifdef __cplusplus
}
#endif
#endif /* !defined(Py_GMPYMODULE_H */
