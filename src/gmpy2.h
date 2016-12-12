/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy.h                                                                  *
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
  Version 2.10  August 2014 (reflect major rewrite during 2013/2014) casevh
 */

#ifndef Py_GMPYMODULE_H
#define Py_GMPYMODULE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Check for Python version requirements. */

#if PY_VERSION_HEX < 0x02060000
#  error "GMPY2 requires Python 2.6 or later."
#endif

#if PY_VERSION_HEX < 0x030200A4
typedef long Py_hash_t;
typedef unsigned long Py_uhash_t;
#  define _PyHASH_IMAG 1000003
#endif

/* Define various macros to deal with differences between Python 2 and 3. */

#if (PY_MAJOR_VERSION == 3)
#define PY3
#define Py2or3String_FromString     PyUnicode_FromString
#define Py2or3String_FromFormat     PyUnicode_FromFormat
#define Py2or3String_Check          PyUnicode_Check
#define Py2or3String_Format         PyUnicode_Format
#define Py2or3String_AsString       PyUnicode_AS_DATA
#define PyStrOrUnicode_Check(op)    (PyBytes_Check(op) || PyUnicode_Check(op))
#define PyIntOrLong_FromLong        PyLong_FromLong
#define PyIntOrLong_Check(op)       PyLong_Check(op)
#define PyIntOrLong_FromSize_t      PyLong_FromSize_t
#define PyIntOrLong_FromSsize_t     PyLong_FromSsize_t
#define PyIntOrLong_AsSsize_t       PyLong_AsSsize_t
#define PyIntOrLong_AsLong          PyLong_AsLong
#else
#define PY2
#define Py2or3String_FromString     PyString_FromString
#define Py2or3String_FromFormat     PyString_FromFormat
#define Py2or3String_Check          PyString_Check
#define Py2or3String_Format         PyString_Format
#define Py2or3String_AsString       PyString_AsString
#define PyStrOrUnicode_Check(op)    (PyString_Check(op) || PyUnicode_Check(op))
#define PyIntOrLong_FromLong        PyInt_FromLong
#define PyIntOrLong_Check(op)       (PyInt_Check(op) || PyLong_Check(op))
#define PyIntOrLong_FromSize_t      PyInt_FromSize_t
#define PyIntOrLong_FromSsize_t     PyInt_FromSsize_t
#define PyIntOrLong_AsSsize_t       PyInt_AsSsize_t
#define PyIntOrLong_AsLong          PyInt_AsLong
#endif

/* Support MPIR, if requested. */

#ifdef MPIR
#  include <mpir.h>
#else
#  include <gmp.h>
#endif

#include <mpfr.h>
#include <mpc.h>

#ifndef ABS
#define ABS(a)  (((a) < 0) ? -(a) : (a))
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

#ifdef __GNUC__
#  define USE_ALLOCA 1
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

#define INDEX_ERROR(msg) PyErr_SetString(PyExc_IndexError, msg)
#define TYPE_ERROR(msg) PyErr_SetString(PyExc_TypeError, msg)
#define VALUE_ERROR(msg) PyErr_SetString(PyExc_ValueError, msg)
#define ZERO_ERROR(msg) PyErr_SetString(PyExc_ZeroDivisionError, msg)
#define SYSTEM_ERROR(msg) PyErr_SetString(PyExc_SystemError, msg)
#define OVERFLOW_ERROR(msg) PyErr_SetString(PyExc_OverflowError, msg)
#define RUNTIME_ERROR(msg) PyErr_SetString(PyExc_RuntimeError, msg)

#define GMPY_DEFAULT -1

/* To prevent excessive memory usage, we don't want to save very large
 * numbers in the cache. The default value specified in the options
 * structure is 128 words (512 bytes on 32-bit platforms, 1024 bytes on
 * 64-bit platforms).
 */
#define MAX_CACHE_LIMBS 16384

/* The maximum number of objects that can be saved in a cache is specified
 * here. The default value is 100.*/
#define MAX_CACHE 1000

#ifdef USE_ALLOCA
#  define TEMP_ALLOC(B, S) \
    if(S < ALLOC_THRESHOLD) { \
        B = alloca(S); \
    } else { \
        if(!(B = malloc(S))) { \
            PyErr_NoMemory(); \
            return NULL; \
        } \
    }
#  define TEMP_FREE(B, S) if(S >= ALLOC_THRESHOLD) free(B)
#else
#  define TEMP_ALLOC(B, S) \
    if(!(B = malloc(S)))  { \
        PyErr_NoMemory(); \
        return NULL; \
    }
#  define TEMP_FREE(B, S) free(B)
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

#ifdef FAST

/* Very bad code ahead. I've copied portions of mpfr-impl.h and
 * hacked them so they work. This code will only be enabled if you
 * specify the --fast option.
 */

#define MPFR_THREAD_ATTR __thread
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR unsigned int __gmpfr_flags;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_exp_t   __gmpfr_emin;
__MPFR_DECLSPEC extern MPFR_THREAD_ATTR mpfr_exp_t   __gmpfr_emax;

/* Flags of __gmpfr_flags */
#define MPFR_FLAGS_UNDERFLOW 1
#define MPFR_FLAGS_OVERFLOW 2
#define MPFR_FLAGS_NAN 4
#define MPFR_FLAGS_INEXACT 8
#define MPFR_FLAGS_ERANGE 16
#define MPFR_FLAGS_DIVBY0 32
#define MPFR_FLAGS_ALL 63

/* Replace some common functions for direct access to the global vars */
#define mpfr_get_emin() (__gmpfr_emin + 0)
#define mpfr_get_emax() (__gmpfr_emax + 0)

#define mpfr_clear_flags()      ((void) (__gmpfr_flags = 0))
#define mpfr_clear_underflow()  ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_UNDERFLOW))
#define mpfr_clear_overflow()   ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_OVERFLOW))
#define mpfr_clear_nanflag()    ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_NAN))
#define mpfr_clear_inexflag()   ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_INEXACT))
#define mpfr_clear_erangeflag() ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_ERANGE))
#define mpfr_clear_divby0()     ((void) (__gmpfr_flags &= MPFR_FLAGS_ALL ^ MPFR_FLAGS_DIVBY0))
#define mpfr_underflow_p()      ((int) (__gmpfr_flags & MPFR_FLAGS_UNDERFLOW))
#define mpfr_overflow_p()       ((int) (__gmpfr_flags & MPFR_FLAGS_OVERFLOW))
#define mpfr_nanflag_p()        ((int) (__gmpfr_flags & MPFR_FLAGS_NAN))
#define mpfr_inexflag_p()       ((int) (__gmpfr_flags & MPFR_FLAGS_INEXACT))
#define mpfr_erangeflag_p()     ((int) (__gmpfr_flags & MPFR_FLAGS_ERANGE))
#define mpfr_divby0_p()         ((int) (__gmpfr_flags & MPFR_FLAGS_DIVBY0))

#define mpfr_check_range(x,t,r) \
 ((((x)->_mpfr_exp) >= __gmpfr_emin && ((x)->_mpfr_exp) <= __gmpfr_emax) \
  ? ((t) ? (__gmpfr_flags |= MPFR_FLAGS_INEXACT, (t)) : 0)                   \
  : mpfr_check_range(x,t,r))

/* End of the really bad code. Hopefully.
 */

#endif

#include "gmpy2_macros.h"

#include "gmpy2_context.h"

#include "gmpy2_mpz.h"
#include "gmpy2_xmpz.h"
#include "gmpy2_mpz_divmod.h"
#include "gmpy2_mpz_divmod2exp.h"
#include "gmpy2_mpz_pack.h"
#include "gmpy2_mpz_bitops.h"
#include "gmpy2_mpz_inplace.h"
#include "gmpy2_xmpz_inplace.h"

#include "gmpy2_mpq.h"

#include "gmpy2_mpfr.h"
#if MPFR_VERSION < 0x030100
#  error gmpy2 requires MPFR 3.1.0 or later
#endif

#include "gmpy2_mpc.h"
#if MPC_VERSION < 0x010000
#  error gmpy2 requires MPC 1.0.0 or later
#endif

#include "gmpy2_convert.h"
#include "gmpy2_convert_utils.h"
#include "gmpy2_convert_gmp.h"
#include "gmpy2_convert_mpfr.h"
#include "gmpy2_convert_mpc.h"

/* Support object caching, creation, and deletion. */

#include "gmpy2_cache.h"

/* Suport for miscellaneous functions (ie. version, license, etc.). */

#include "gmpy2_misc.h"

/* Support conversion to/from binary format. */

#include "gmpy2_binary.h"

/* Support random number generators. */

#include "gmpy2_random.h"

/* Support Lucas sequences. */

#include "gmpy_mpz_lucas.h"

/* Support probable-prime tests. */

#include "gmpy_mpz_prp.h"

/* Begin includes for refactored code. */

#include "gmpy2_abs.h"
#include "gmpy2_add.h"
#include "gmpy2_divmod.h"
#include "gmpy2_floordiv.h"
#include "gmpy2_minus.h"
#include "gmpy2_mod.h"
#include "gmpy2_mul.h"
#include "gmpy2_plus.h"
#include "gmpy2_pow.h"
#include "gmpy2_sub.h"
#include "gmpy2_truediv.h"
#include "gmpy2_math.h"
#include "gmpy2_const.h"
#include "gmpy2_square.h"
#include "gmpy2_format.h"
#include "gmpy2_hash.h"
#include "gmpy2_fused.h"
#include "gmpy2_muldiv_2exp.h"
#include "gmpy2_predicate.h"
#include "gmpy2_sign.h"
#include "gmpy2_richcompare.h"
#include "gmpy2_mpc_misc.h"
#include "gmpy2_mpfr_misc.h"
#include "gmpy2_mpq_misc.h"
#include "gmpy2_mpz_misc.h"
#include "gmpy2_xmpz_misc.h"

#ifdef VECTOR
#include "gmpy2_vector.h"
#endif

#ifdef __cplusplus
}
#endif
#endif /* !defined(Py_GMPYMODULE_H */
