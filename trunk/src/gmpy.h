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

#if !defined(FLT_RADIX) || (FLT_RADIX!=2)
#   error "FLT_RADIX undefined or != 2, GMPY2 is confused. :("
#endif

#if defined(MS_WIN32) && defined(_MSC_VER)
   /* so one won't need to link explicitly to gmp.lib...: */
#  if defined(MPIR)
#    pragma comment(lib,"mpir.lib")
#  else
#    pragma comment(lib,"gmp.lib")
#  endif
#  pragma comment(lib,"mpfr.lib")
#  pragma comment(lib,"mpc.lib")
#  define isnan _isnan
#  define isinf !_finite
#  define USE_ALLOCA 1
#  define inline __inline
#endif

#if defined MPIR
#  include "mpir.h"
#else
#  include "gmp.h"
#endif

#include "mpfr.h"
#include "mpc.h"

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

#define TYPE_ERROR(msg) PyErr_SetString(PyExc_TypeError, msg)
#define VALUE_ERROR(msg) PyErr_SetString(PyExc_ValueError, msg)
#define ZERO_ERROR(msg) PyErr_SetString(PyExc_ZeroDivisionError, msg)
#define SYSTEM_ERROR(msg) PyErr_SetString(PyExc_SystemError, msg)
#define OVERFLOW_ERROR(msg) PyErr_SetString(PyExc_OverflowError, msg)

#define GMPY_DIVZERO(msg) PyErr_SetString(GMPyExc_DivZero, msg)
#define GMPY_INEXACT(msg) PyErr_SetString(GMPyExc_Inexact, msg)
#define GMPY_INVALID(msg) PyErr_SetString(GMPyExc_Invalid, msg)
#define GMPY_OVERFLOW(msg) PyErr_SetString(GMPyExc_Overflow, msg)
#define GMPY_UNDERFLOW(msg) PyErr_SetString(GMPyExc_Underflow, msg)
#define GMPY_ERANGE(msg) PyErr_SetString(GMPyExc_Erange, msg)

#define CHECK_UNDERFLOW(msg) \
    if (mpfr_underflow_p() && context->now.trap_underflow) { \
        GMPY_UNDERFLOW(msg); \
        goto done; \
    }
#define CHECK_OVERFLOW(msg) \
    if (mpfr_overflow_p() && context->now.trap_overflow) { \
        GMPY_OVERFLOW(msg); \
        goto done; \
    }
#define CHECK_INVALID(msg) \
    if (mpfr_nanflag_p() && context->now.trap_invalid) { \
        GMPY_INVALID(msg); \
        goto done; \
    }
#define CHECK_INEXACT(msg) \
    if (mpfr_inexflag_p() && context->now.trap_inexact) { \
        GMPY_INEXACT(msg); \
        goto done; \
    }
#define CHECK_ERANGE(msg) \
    if (mpfr_erangeflag_p() && context->now.trap_erange) { \
        GMPY_ERANGE(msg); \
        goto done; \
    }

#define MERGE_FLAGS \
    context->now.underflow |= mpfr_underflow_p(); \
    context->now.overflow |= mpfr_overflow_p(); \
    context->now.invalid |= mpfr_nanflag_p(); \
    context->now.inexact |= mpfr_inexflag_p(); \
    context->now.erange |= mpfr_erangeflag_p();

#define CHECK_FLAGS(NAME) \
    CHECK_INVALID("invalid operation in 'mpfr' "NAME); \
    CHECK_UNDERFLOW("underflow in 'mpfr' "NAME); \
    CHECK_OVERFLOW("overflow in 'mpfr' "NAME); \
    CHECK_INEXACT("inexact result in 'mpfr' "NAME); \

#define SUBNORMALIZE(NAME) \
    if (context->now.subnormalize) \
        NAME->rc = mpfr_subnormalize(NAME->f, NAME->rc, context->now.mpfr_round);


#define MPFR_CLEANUP_SELF(NAME) \
    SUBNORMALIZE(result); \
    MERGE_FLAGS; \
    CHECK_INVALID("invalid operation in 'mpfr' "NAME); \
    CHECK_UNDERFLOW("underflow in 'mpfr' "NAME); \
    CHECK_OVERFLOW("overflow in 'mpfr' "NAME); \
    CHECK_INEXACT("inexact result in 'mpfr' "NAME); \
  done: \
    Py_DECREF(self); \
    if (PyErr_Occurred()) { \
        Py_XDECREF((PyObject*)result); \
        result = NULL; \
    } \
    return (PyObject*)result;

#define MPFR_CLEANUP_SELF_OTHER(NAME) \
    SUBNORMALIZE(result); \
    MERGE_FLAGS; \
    CHECK_INVALID("invalid operation in 'mpfr' "NAME); \
    CHECK_UNDERFLOW("underflow in 'mpfr' "NAME); \
    CHECK_OVERFLOW("overflow in 'mpfr' "NAME); \
    CHECK_INEXACT("inexact result in 'mpfr' "NAME); \
  done: \
    Py_DECREF(self); \
    Py_DECREF(other); \
    if (PyErr_Occurred()) { \
        Py_XDECREF((PyObject*)result); \
        result = NULL; \
    } \
    return (PyObject*)result;

#define MPFR_CLEANUP_RF(NAME) \
    SUBNORMALIZE(rf); \
    MERGE_FLAGS; \
    if (mpfr_underflow_p() && context->now.trap_underflow) { \
        GMPY_UNDERFLOW("underflow in 'mpfr' " #NAME); \
        Py_DECREF((PyObject*)rf); \
        return NULL; \
    } \
    if (mpfr_overflow_p() && context->now.trap_overflow) { \
        GMPY_OVERFLOW("overflow in 'mpfr' " #NAME); \
        Py_DECREF((PyObject*)rf); \
        return NULL; \
    } \
    if (mpfr_inexflag_p() && context->now.trap_inexact) { \
        GMPY_INEXACT("inexact result in 'mpfr' " #NAME); \
        Py_DECREF((PyObject*)rf); \
        return NULL; \
    } \
    return (PyObject*)rf;


#ifdef DEBUG
#define TRACE(msg) if(options.debug) fprintf(stderr, msg)
#else
#define TRACE(msg)
#endif

#define GMPY_DEFAULT -1

#ifdef USE_ALLOCA
#define TEMP_ALLOC(B, S) \
    if(S < ALLOC_THRESHOLD) { \
        B = alloca(S); \
    } else { \
        if(!(B = PyMem_Malloc(S))) { \
            PyErr_NoMemory(); \
            return NULL; \
        } \
    }
#define TEMP_FREE(B, S) if(S >= ALLOC_THRESHOLD) PyMem_Free(B)
#else
#define TEMP_ALLOC(B, S) \
    if(!(B = PyMem_Malloc(S)))  { \
        PyErr_NoMemory(); \
        return NULL; \
    }
#define TEMP_FREE(B, S) PyMem_Free(B)
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
    mpfr_t f;
    Py_hash_t hash_cache;
    int rc;
    int round_mode;
} PympfrObject;

typedef struct {
    mpob ob;
    mpc_t c;
    Py_hash_t hash_cache;
    int rc;
    int round_mode;
} PympcObject;

typedef struct {
    mpob ob;
    mpz_t z;
} PyxmpzObject;

#define Pympz_AS_MPZ(obj) (((PympzObject *)(obj))->z)
#define Pympq_AS_MPQ(obj) (((PympqObject *)(obj))->q)
#define Pympfr_AS_MPFR(obj) (((PympfrObject *)(obj))->f)
#define Pympc_AS_MPC(obj) (((PympcObject *)(obj))->c)
#define Pyxmpz_AS_MPZ(obj) (((PyxmpzObject *)(obj))->z)

static PyTypeObject Pympz_Type;
#define Pympz_Check(v) (((PyObject*)v)->ob_type == &Pympz_Type)

static PyTypeObject Pympq_Type;
#define Pympq_Check(v) (((PyObject*)v)->ob_type == &Pympq_Type)

static PyTypeObject Pympfr_Type;
#define Pympfr_Check(v) (((PyObject*)v)->ob_type == &Pympfr_Type)
#define Pympfr_CheckAndExp(v) (Pympfr_Check(v) && \
(mpfr_zero_p(Pympfr_AS_MPFR(v)) || \
(mpfr_regular_p(Pympfr_AS_MPFR(v)) && \
(Pympfr_AS_MPFR(v)->_mpfr_exp >= context->now.emin) && \
(Pympfr_AS_MPFR(v)->_mpfr_exp <= context->now.emax))))

static PyTypeObject Pympc_Type;
#define Pympc_Check(v) (((PyObject*)v)->ob_type == &Pympc_Type)

static PyTypeObject Pyxmpz_Type;
#define Pyxmpz_Check(v) (((PyObject*)v)->ob_type == &Pyxmpz_Type)

#define CHECK_MPZANY(v) (Pympz_Check(v) || Pyxmpz_Check(v))

typedef struct {
    int subnormalize;        /* if 1, subnormalization is performed */
    mpfr_prec_t mpfr_prec;   /* current precision in bits, for MPFR */
    mpfr_prec_t mpc_rprec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t mpc_iprec;   /* current precision in bits, for Im(MPC) */
    mpfr_rnd_t mpfr_round;   /* current rounding mode for float (MPFR) */
    mpfr_rnd_t mpc_rround;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t mpc_iround;   /* current rounding mode for Im(MPC) */
    mpc_rnd_t mpc_round;     /* current rounding mode for complex (MPC)*/
    mpfr_exp_t emax;         /* maximum exponent */
    mpfr_exp_t emin;         /* minimum exponent */
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
} gmpy_context;

typedef struct {
    PyObject_HEAD;
    gmpy_context now;        /* The "new" values, used by __enter__ */
    PyObject *orig;          /* Original context, restored by __exit__*/
} GMPyContextObject;

static PyTypeObject GMPyContext_Type;
#define GMPyContext_Check(v) (((PyObject*)v)->ob_type == &GMPyContext_Type)

#ifdef __cplusplus
}
#endif
#endif /* !defined(Py_GMPYMODULE_H */
