/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2.c                                                                 *
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

 /* Todo list
  * ---------
  * Add all MPFR and MPC functions as context methods.
  * All MPFR and MPC functions need to set exponent range on entry. The
  *    current approach where only set_context() and context.__enter__ set
  *    the exponent range fails for context methods.
  * Should a read-only (or template) context prevent the setting of
  *    exception flags?
  * Add context option to control the result of integer division:
  *    integer (mpz), exact (mpq), or true (mpfr).
  * Add modular arithmetic functions.
  * Implement Chinese Remainder Theorem.
  * Update PRP code.
  */

/*
 * originally written for GMP-2.0 (by AMK...?)
 * Rewritten by Niels M�ller, May 1996
 *
 * Version for GMP-4, Python 2.X, with support for MSVC++6,
 * addition of mpf's, &c: Alex Martelli (now aleaxit@gmail.com, Nov 2000).
 * cleanups & reorgs leading to 1.0: Alex Martelli (until Aug 2003)
 * further cleanups and bugfixes leading to 1.01, Alex Martelli (Nov 2005)
 * minor bugfixes+new decimal (&c) support to 1.02, Alex Martelli (Feb 2006)
 * various bugfixes for 64-bit platforms, 1.03, aleaxit and casevh (Jun 2008)
 * rich comparisons, 1.04, aleaxit and casevh (Jan 2009)
 * support for Python 3.x, 1.10, casevh (Oct 2009)
 *
 * Some hacks by Gustavo Niemeyer <niemeyer@conectiva.com>.
 *
 *   0.1, pre-alpha; date: 2000-11-06  first placed on sourceforge
 *
 *   0.2, still pre-alpha: 2000-11-15: bugfixes re formatting (tx, Peanu!)
 *   no tags on oct() and hex() of mpz's
 *   insert 'tagoff' in options (gmpy.mpz() vs mpz() in repr) (for Peanu!)
 *   speedups for _nonzero & _cmp (tx, Peanu!)
 *   slight speedup (7/8%?) for excess reallocs 4<->8 bytes (Peanu's help!)
 *   added copy/fcopy; bin; fib; remove
 *
 *   0.3, still pre-alpha, but...:
 *   performance tweaks via mpz-caching & fixed-constants
 *   added get/set functions for zcache, zco min/max
 *   added get-only function for versions (of gmp, and of gmpy)
 *   removed all 'traces' of mutability (to be re-done... much later!)
 *   cleaned up all of the mpz_cmp_ui(X,0) to mpz_sgn(X)
 *   cleaned up Py_BuildValue usage (N vs O, explicit-() for tuples)
 *   added numdigits, lowbits, root, next_prime, invert, popcount,
 *      hamdist, scan0, scan1
 *   renamed bin to bincoef
 *
 *   0.4:
 *   split gmpy.c/gmpy.h introducing C-API interface (Pearu's suggestion)
 *   cleanup some casts using Pearu's new macros
 *   further cache-tweaks at Pearu's suggestion (macros introduced)
 *   added sign (Pearu's request), getbit, setbit
 *   added docstrings
 *   renamed copy functions to start with _ ('internal, private')
 *   added .comb as a synonym of .bincoef
 *
 *   0.5:
 *   added jacobi, legendre, kronecker
 *   added random-number generation, seed set/save, shuffling
 *   added mpq (at last!-)
 *
 *   0.6: (lots of good ideas from Pearu once more!-):
 *   fixed silly bugs in kronecker and mpq_abs
 *   gmpy-level workaround for scan0/scan1 bugs (?) in GMP 3.1.1
 *   added qdiv; anynum->mpq substituted for all such conversions
 *       (also anynum->mpz and anynum->mpf by analogy, with care!)
 *   added global.fcoform for optional use of intermediate string in
 *       float2mpf (used for any float->mpf conversion)
 *   added set_fcoform function for global.fcoform access
 *   general cleanup of sources; added alloca for MSVC++;
 *       many sundry minor bugfixes & uniformization;
 *       a little useful refactoring (more would be good...)
 *   added caching of mpq objects
 *   power for mpq
 *   Stern-Brocot algorithm for mpf->mpq (also exposed as f2q)
 *       also used for float->mpq
 *       with stricter tracking of mpf's requested-precision
 *       added getrprec method to mpf, getrprec module-function
 *   exposed ceil, floor and trunc methods/functions for mpf's
 *   changed a couple exceptions from Value to ZeroDivision
 *   added 'qual' and 'floa' options to gmpy.rand
 *
 *   0.7: (good feedback from Keith Briggs, some advice from Tim Peters
 *      and Fred Lundh -- thanks all!):
 *   fixed bug of '"%d" where "%ld" was meant' in many places
 *      and other sundry minor warnings given by gcc
 *   fixed hash (delegating to Python) so mp[nqz](x) will
 *      produce the same value as hash(x) for any Python number x
 *   workaround for GMP 3.1.1 bug, mpz_root wrongly returning
 *      'exact' for non-exact root if dest==source, which stopped
 *      needed value-error for inexact mpq**mpq operations
 *   determined correct 'actual precision' of floats
 *   explicitly stored precision with binary-form mpf's
 *   extended explicit-bits request to all ->mpf operations
 *     (good in itself, plus, preparing for future MPFR)
 *   removed the limitation of no binary-form for <0 mpz
 *   introduced macros to parse args, for conciseness
 *
 *   0.8: (again, requests & suggestions by great Pearu!)
 *   raise test coverage 72.5% -> 90.0%
 *   introduced callbacks (not documented/tested for now;
 *       Pearu will test/support/document in PySymbolic)
 *   some errors went undiagnosed, caused crash: now fixed
 *   workaround for GMP bug(?s?) in mpz_fits_... (?)
 *   added exposure of mpf_ sqrt and pow_ui
 *
 *   0.9: (ditto)
 *   change ValueError to OverflowError for 'too-large' errors
 *   fix bug in mpq_pow (negative base, exp. with odd denominator)
 *       (fix now corrected -- _even_ denominator is the error!)
 *   fixed gcc warnings reported by K. Briggs
 *
 *   0.9b:
 *   support GMP 4 (but added no GMP4-only functionality yet)
 *
 *   0.9c:
 *   updated tests to 0.9, better coverage
 *
 *   1.0:
 *   minor cleanups, ensure support for Python 2.3
 *   fixed misdiagnosis of some argument counts in macro
 *     SELF_ONE_ARG_CONVERTED (tx to Paul Rubin!)
 *
 *   1.01:
 *   cleanups, ensure support for Python 2.4.1 on MacOSX 10.4/XCode 2.1
 *     as well as Python 2.2 and 2.3 (on MacOSX and Linux)
 *   fixed memory leak on divm (thanks to mensanator@aol.com)
 *   fixed bug on mpq('123') [[str2mpq on string w/o a slash]]
 *   added floordiv and truediv operators, and tests for them
 *   NOT tested on GMP 3 (have none left around...), ONLY on GMP 4.*
 *
 *   1.02:
 *   fix warning in comparison of mpq's
 *   added support of mpq('12.34') [[string w/o a slash, but with a dot]]
 *   fixes for 64-bit build (thanks to a patch by dmcooke)
 *   added experimental support for decimal.Decimal (and user-coded types)
 *      via wider use of special conversion methods (if present) and their
 *      sly insertion on-the-fly into the decimal.Decimal class (!)
 *   two bugfixes, thanks to Simon Burton
 *   Brought back into C89 compliance (thanks to Chip Turner), had
 *      drifted to C99 (declarations in the middle of the code).
 *   Python 2.5 support (Py_ssize_t, __index__) thanks to Chip Turner
 *   Pushed coverage to 93.3% (missing only "sanity check" level error
 *      tests [mostly for out-of-memory conditions], output to stderr
 *      conditioned by global.debug, & a couple of very obscure cases)
 *
 *   1.03:
 *   Fixed the bug that caused crashes on gmpy.mpf(float('inf')) and
 *      other such conversions, implicit and explicit
 *   Fixed a bug in get_zconst's prototype affecting 64-bit machines,
 *      thanks to Gary Bunting
 *   Fixed a bug in hashing on 64-bit systems. hash(long) now equals
 *      hash(mpz) for large values. (casevh)
 *   Changed int() to return a long value instead of OverFlowError.
 *      Complies with PEP 237. (casevh)
 *   Added support in setup.py for darwinports/macports build of GMP
 *      on MacOSX. (aleaxit)
 *
 *   1.04:
 *   Avoid GMP/mingw32 bug when converting very small floats to mpz.
 *      (casevh)
 *   Significant performance improvement for long->mpz and mpz->long.
 *      (casevh)
 *   Added "rich comparisons" to mpz, mpq and mpf types (aleaxit)
 *   Added additional tests (casevh, aleaxit)
 *   Fixed bug when converting very large mpz to str (casevh)
 *   Faster conversion from mpz->binary and binary->mpz (casevh)
 *   Added support for pickling (casevh)
 *   Added divexact (casevh)
 *   Fixed mpf comparisons by rounding mpf results when GMP returns
 *      a longer result. Added fround() (casevh)
 *   Added bit_length (Thanks Mario Pernici)
 *   Added helper functions for mpmath (casevh)
 *   Faster conversion from mpq->binary and binary->mpq (casevh)
 *   Recognize MPIR, mpir_version() (casevh)
 *
 *   1.10:
 *   Remove dependancy on pymemcompat.h (casevh)
 *   Remove callback (casevh)
 *   Added support for -DMPIR to include MPIR instead of GMP (casevh)
 *   Major code revisions to add support for Python 3.x (casevh)
 *   Fixed bug in binary() and qbinary() (casevh)
 *   Fixed bug in rich comparisons (casevh)
 *   Added % and divmod support to mpq and mpf (casevh)
 *   Changed memory allocation functions to use PyMem (casevh)
 *   Removed small number interning (casevh)
 *   Added tdivmod, cdivmod, and fdivmod (casevh)
 *   Added more helper functions for mpmath (casevh)
 *   Faster mpz<>PyLong conversion (casevh)
 *   Faster hash(mpz) (casevh)
 *
 *   1.11:
 *   Recognize True/False (bug in 1.10) (casevh)
 *   Optimize argument handling (casevh)
 *   Added caching for mpz (casevh)
 *
 ************************************************************************
 *
 *   2.0.0 alpha and b1:
 *   Added caching for mpq (casevh)
 *   Added rootrem, fib2, lucas, lucas2 (casevh)
 *   Removed mpf.setprec(), use mpf.round() (casevh)
 *   Fix test compatibility with Python 3.1.2 and 3.2 (casevh)
 *   Support changed hash function in Python 3.2 (casevh)
 *   Added is_even, is_odd (casevh)
 *   Rename to gmpy2 to allow backwards incompatible changes (casevh)
 *   Remove old random number functions, to be replaced later (casevh)
 *   Add caching of the calculated hash value (casevh)
 *   Add xmpz (mutable mpz) type (casevh)
 *   Fix mpq formatting issue (casevh)
 *   Add read/write bit access using slices to xmpz (casevh)
 *   Add read-only bit access using slices to mpz (casevh)
 *   Add pack()/unpack() methods (casevh)
 *   Remove tagoff option (casevh)
 *   Add support for MPFR (casevh)
 *   Debug messages only available if compiled with -DDEBUG (casevh)
 *   Removed fcoform float conversion modifier (casevh)
 *   Add support for MPC (casevh)
 *   Renamed 'mpf' to 'mpfr' to reflect use of MPFR (casevh)
 *   Added context manager (casevh)
 *   Allow building with just GMP/MPIR if MPFR not available (casevh)
 *   Allow building with GMP/MPIR and MPFR if MPC not available (casevh)
 *   Removed most instance methods in favor of gmpy2.method (casevh)
 *   Added __ceil__, __floor__, and __trunc__ methods (casevh)
 *   Removed gmpy2.pow to avoid conflicts (casevh)
 *   Removed gmpy2._copy and added xmpz.copy (casevh)
 *   Added support for __format__ (casevh)
 *   Completed support for MPC (casevh)
 *   Added as_integer_ratio, as_mantissa_exp, as_simple_fraction (casevh)
 *   Update rich_compare (casevh)
 *   Require MPFR 3.1.0+ to get divby0 support (casevh)
 *   Added fsum(), degrees(), radians() (casevh)
 *   Renamed context() -> local_context(), new_context() -> context() (casevh)
 *   Added get_context() (casevh)
 *   Added random number generation support (casevh)
 *   Changed license to LGPL 3+ (casevh)
 *   Added lucasu, lucasu_mod, lucasv, and lucasv_mod (casevh)
 *      (Based on code contributed by David Cleaver.)
 *   Added probable-prime tests (casevh)
 *      (Based on code contributed by David Cleaver.)
 *   Added to_binary()/from_binary (casevh)
 *   Renamed numdigits() to num_digits() (casevh)
 *   Added keyword precision to constants (casevh)
 *   Added addmul() and submul() (casevh)
 *   Added __round__(), round2(), round_away() for mpfr (casevh)
 *   round() is no longer a module level function (casevh)
 *   pow() is no longer a module level function (casevh)
 *   Renamed module functions min()/max() to min2()/max2() (casevh)
 *      No longer conflicts with builtin min() and max()
 *   Removed set_debug() and related functionality (casevh)
 *   Released as 2.0.0b1
 *
 *   2.0.0b2
 *   Allow xmpz slice assignment to increase length of xmpz instance by
 *      specifying a value for stop (casevh)
 *   Fixed ref-count bug in several is_xxx_prp tests (casevh)
 *   Added iter_bits, iter_clear, iter_set methods to xmpz (casevh)
 *   Added powmod() for easy access to three argument pow() (casevh)
 *   Removed addmul() and submul() since they are slower than (casevh)
 *      just using Python code
 *   Bug fix in gcd_ext when both arguments are not mpz (casevh)
 *   Added ieee() to create contexts for 32, 64, or 128 bit floats (casevh)
 *   Bug fix in context() not setting emax/emin correctly if they (casevh)
 *      had been changed earlier
 *   Contexts can be directly used in with statement without (casevh)
 *      requiring set_context()/local_context() sequence
 *   local_context() now accepts an optional context (casevh)
 *
 *   2.0.0b3
 *   mp_version(), mpc_version(), and mpfr_version() shouldn't (casevh)
 *      return Unicode on Python 2.x
 *   Fix warnings when shifting 32-bit integer by 32 bits (casevh)
 *   Faster conversion of Fraction to gmpy2 types (casevh)
 *   Fix conversion with Decimal, especially on Python 3.3 (casevh)
 *   Consistently return OverflowError when converting "inf" (casevh)
 *   Fix mpz.__format__() with # code (casevh)
 *   Add is_infinite(), deprecate is_inf() (casevh)
 *   Add is_finite(), deprecate is_number() (casevh)
 *   Fixed issues with mpc() and various is_XXX() functions (casevh)
 *   Fixed error handling with mpc(); mpc(1,"nan") is properly (casevh)
 *      handled
 *   Added caching for mpc objects; faster when real and  (casevh)
 *      imaginary precisions are equal
 *   Add optimal path for mpfr/mpc + - * / when both operands (casevh)
 *      have the same type
 *   Fix mpfr + float segmentation fault (casevh)
 *
 *   2.0.0b4
 *   Add __ceil__, __floor__, __trunc__, __round__ to mpz & mpq (casevh)
 *   Add __complex__ to mpc (casevh)
 *   round(mpfr) now correctly returns an mpz (casevh)
 *   Add mpz.denominator and mpz.numerator (casevh)
 *   mpz() returns mpz(0); also xmpz, mpq, mpfr, and mpc (casevh)
 *   Fix bug when comparing mpz to mpq (with mpz on left) (casevh)
 *   Add __sizeof__ (casevh)
 *
 *   2.0.0
 *   Fix segfault in _mpmath_normalize if rnd not specified (casevh)
 *   Improved setup.py (casevh)
 *   Fix issues encountered when compiled without MPFR support (casevh)
 *   Conversion of too large an mpz to float now raises OverflowError (casevh)
 *   Renamed module functions min2()/max2() to minnum()/maxnum() (casevh)
 *   Added copy() method to contexts (casevh)
 *   get_context() no longer supports keyword arguments (casevh)
 *
 *   2.1.0
 *   Improvements to setup.py (casevh)
 *   Add thread-safe contexts (casevh)
 *   MPFR and MPC are now required (casevh)
 *   Invalid Operation exception now raised for addition, etc. (casevh)
 *   inverse() now raises exception if inverse does not exist (casevh)
 *   Add context methods (casevh) (in progress)
 *
 ************************************************************************
 *
 * Discussion on sizes of C integer types.
 *
 * GMP, MPIR, MPFR, and MPC use typedef to create integer objects with
 * different sizes. It can become confusing to map the different types
 * onto the standard C types used by Python's C API. Below are external
 * types and how the map to C types. The assumptions are verified when
 * the module is initialized.
 *
 * mp_limb_t: This is usually an 'unsigned long' but is an 'unsigned
 *     long long' on MPIR/64-bit Windows.
 *
 * mp_bitcnt_t: This is usually an 'unsigned long' but is an 'unsigned
 *     long long' on MPIR/64-bit Windows. 'size_t' is probably the best
 *     match since that is the return type of mpz_sizeinbase().
 *
 * mp_size_t: This is a 'long'.
 *
 * mp_exp_t: This is defined by GMP to be a 'long'.
 *
 * mpfr_rnd_t: This is an 'int' ???
 *
 * mpfr_prec_t: This can be 'short', 'int', or 'long'. GMPY assumes it
 *     will be a 'long'.
 *
 * mpfr_sign_t: This is an 'int'.
 *
 * mpfr_exp_t: This is currently the same as mp_exp_t but will change
 *     to a signed 64-bit integer in the future.
 *
 * mpc_rnd_t: This is an 'int'.
 *
 * mpc_prec_t: See mpfr_exp_t.
 *
 ************************************************************************
 */
#define PY_SSIZE_T_CLEAN
#include "Python.h"

/*
 * we do have a dependence on Python's internals, specifically:
 * how Python "long int"s are internally represented.
 */
#include "longintrepr.h"

#include <assert.h>
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <ctype.h>

/* Include functions missing from the Python 2.6 C-API that are include
 * in later versions of Python.
 */

#if ((PY_VERSION_HEX < 0x02070000) || ((PY_MAJOR_VERSION == 3) && (PY_MINOR_VERSION < 2)))
#include "py3intcompat.c"
#endif

#include "gmpy.h"

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Global data declarations begin here.                                    *
 * NOTE: Because of these global declarations, GMPY2 is not thread-safe!   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* The following global strings are used by gmpy_misc.c. */

char gmpy_version[] = "2.1.0a0";

char _gmpy_cvs[] = "$Id$";

char gmpy_license[] = "\
The GMPY2 source code is licensed under LGPL 3 or later. The supported \
versions of the GMP/MPIR, MPFR, and MPC libraries are also licensed under \
LGPL 3 or later.";

/* The following global structures are used by gmpy_cache.c.
 */

static struct gmpy_global {
    int cache_size;          /* size of cache, for all caches */
    int cache_obsize;        /* maximum size of the objects that are cached */
} global = {
    100,                     /* cache_size */
    128,                     /* cache_obsize */
};

static mpz_t* zcache;
static int in_zcache;

static MPZ_Object **gmpympzcache;
static int in_gmpympzcache;

static XMPZ_Object **gmpyxmpzcache;
static int in_gmpyxmpzcache;

static MPQ_Object **gmpympqcache;
static int in_gmpympqcache;

static MPFR_Object **gmpympfrcache;
static int in_gmpympfrcache;

static MPC_Object **gmpympccache;
static int in_gmpympccache;

/* Support for context manager. */

#ifdef WITHOUT_THREADS
/* Use a module-level context. */
static CTXT_Object *module_context = NULL;
#else
/* Key for thread state dictionary */
static PyObject *tls_context_key = NULL;
/* Invariant: NULL or the most recently accessed thread local context */
static CTXT_Object *cached_context = NULL;
#endif


/* Define gmpy2 specific errors for mpfr and mpc data types. No change will
 * be made the exceptions raised by mpz, xmpz, and mpq.
 */

static PyObject *GMPyExc_GmpyError = NULL;
static PyObject *GMPyExc_DivZero = NULL;
static PyObject *GMPyExc_Inexact = NULL;
static PyObject *GMPyExc_Invalid = NULL;
static PyObject *GMPyExc_Overflow = NULL;
static PyObject *GMPyExc_Underflow = NULL;
static PyObject *GMPyExc_Erange = NULL;


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * End of global data declarations.                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* The code for object creation, deletion, and caching is in gmpy_cache.c. */

#include "gmpy2_cache.c"

/* Miscellaneous helper functions and simple methods are in gmpy_misc.c. */

#include "gmpy2_misc.c"

/* Include fast mpz to/from PyLong conversion from sage. */

#include "mpz_pylong.c"

/* Support for conversion to/from binary representation. */

#include "gmpy2_binary.c"

/* Support for conversions to/from numeric types. */

#include "gmpy2_convert.c"
#include "gmpy2_convert_gmp.c"
#include "gmpy2_convert_mpfr.c"
#include "gmpy2_convert_mpc.c"

/* Support for random numbers. */

#include "gmpy2_random.c"

/* Support for Lucas sequences. */

#include "gmpy_mpz_lucas.c"

/* Support for probable-prime tests. */

#include "gmpy_mpz_prp.c"

/* Include helper functions for mpmath. */

#include "gmpy_mpmath.c"

#include "gmpy2_mpz_divmod.c"
#include "gmpy2_mpz_divmod2exp.c"
#include "gmpy2_mpz_pack.c"
#include "gmpy2_mpz_bitops.c"
#include "gmpy2_mpz_inplace.c"
#include "gmpy2_xmpz_inplace.c"

/* Begin includes of refactored code. */

#include "gmpy2_abs.c"
#include "gmpy2_add.c"
#include "gmpy2_divmod.c"
#include "gmpy2_floordiv.c"
#include "gmpy2_minus.c"
#include "gmpy2_mod.c"
#include "gmpy2_mul.c"
#include "gmpy2_plus.c"
#include "gmpy2_pow.c"
#include "gmpy2_sub.c"
#include "gmpy2_truediv.c"
#include "gmpy2_math.c"
#include "gmpy2_const.c"
#include "gmpy2_square.c"
#include "gmpy2_format.c"
#include "gmpy2_hash.c"
#include "gmpy2_fused.c"
#include "gmpy2_muldiv_2exp.c"
#include "gmpy2_predicate.c"
#include "gmpy2_sign.c"
#include "gmpy2_richcompare.c"
#include "gmpy2_mpc_misc.c"
#include "gmpy2_mpfr_misc.c"

/* Include gmpy_context last to avoid adding doc names to .h files. */

#include "gmpy_mpz.c"
#include "gmpy_xmpz.c"
#include "gmpy_mpq.c"
#include "gmpy_mpfr.c"
#include "gmpy2_mpc.c"

#include "gmpy2_context.c"

static PyMethodDef Pygmpy_methods [] =
{
    { "_cvsid", GMPy_get_cvsid, METH_NOARGS, GMPy_doc_cvsid },
    { "_printf", GMPy_printf, METH_VARARGS, GMPy_doc_function_printf },
    { "add", GMPy_Context_Add, METH_VARARGS, GMPy_doc_function_add },
    { "bit_clear", GMPy_MPZ_bit_clear_function, METH_VARARGS, doc_bit_clear_function },
    { "bit_flip", GMPy_MPZ_bit_flip_function, METH_VARARGS, doc_bit_flip_function },
    { "bit_length", GMPy_MPZ_bit_length_function, METH_O, doc_bit_length_function },
    { "bit_mask", GMPy_MPZ_bit_mask, METH_O, doc_bit_mask },
    { "bit_scan0", GMPy_MPZ_bit_scan0_function, METH_VARARGS, doc_bit_scan0_function },
    { "bit_scan1", GMPy_MPZ_bit_scan1_function, METH_VARARGS, doc_bit_scan1_function },
    { "bit_set", GMPy_MPZ_bit_set_function, METH_VARARGS, doc_bit_set_function },
    { "bit_test", GMPy_MPZ_bit_test_function, METH_VARARGS, doc_bit_test_function },
    { "bincoef", Pympz_bincoef, METH_VARARGS, doc_bincoefg },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combg },
    { "c_div", GMPy_MPZ_c_div, METH_VARARGS, doc_c_div },
    { "c_div_2exp", GMPy_MPZ_c_div_2exp, METH_VARARGS, doc_c_div_2exp },
    { "c_divmod", GMPy_MPZ_c_divmod, METH_VARARGS, doc_c_divmod },
    { "c_divmod_2exp", GMPy_MPZ_c_divmod_2exp, METH_VARARGS, doc_c_divmod_2exp },
    { "c_mod", GMPy_MPZ_c_mod, METH_VARARGS, doc_c_mod },
    { "c_mod_2exp", GMPy_MPZ_c_mod_2exp, METH_VARARGS, doc_c_mod_2exp },
    { "denom", Pympq_denom, METH_VARARGS, doc_denomg },
    { "digits", GMPy_Context_Digits, METH_VARARGS, GMPy_doc_context_digits },
    { "div", GMPy_Context_TrueDiv, METH_VARARGS, GMPy_doc_truediv },
    { "divexact", Pygmpy_divexact, METH_VARARGS, doc_divexactg },
    { "divm", Pygmpy_divm, METH_VARARGS, doc_divm },
    { "div_mod", GMPy_Context_DivMod, METH_VARARGS, GMPy_doc_divmod },
    { "fac", Pygmpy_fac, METH_O, doc_fac },
    { "fib", Pygmpy_fib, METH_O, doc_fib },
    { "fib2", Pygmpy_fib2, METH_O, doc_fib2 },
    { "floor_div", GMPy_Context_FloorDiv, METH_VARARGS, GMPy_doc_floordiv },
    { "from_binary", GMPy_MPANY_From_Binary, METH_O, doc_from_binary },
    { "f_div", GMPy_MPZ_f_div, METH_VARARGS, doc_f_div },
    { "f_div_2exp", GMPy_MPZ_f_div_2exp, METH_VARARGS, doc_f_div_2exp },
    { "f_divmod", GMPy_MPZ_f_divmod, METH_VARARGS, doc_f_divmod },
    { "f_divmod_2exp", GMPy_MPZ_f_divmod_2exp, METH_VARARGS, doc_f_divmod_2exp },
    { "f_mod", GMPy_MPZ_f_mod, METH_VARARGS, doc_f_mod },
    { "f_mod_2exp", GMPy_MPZ_f_mod_2exp, METH_VARARGS, doc_f_mod_2exp },
    { "gcd", Pygmpy_gcd, METH_VARARGS, doc_gcd },
    { "gcdext", Pygmpy_gcdext, METH_VARARGS, doc_gcdext },
    { "get_cache", GMPy_get_cache, METH_NOARGS, GMPy_doc_get_cache },
    { "hamdist", GMPy_MPZ_hamdist, METH_VARARGS, doc_hamdist },
    { "invert", Pygmpy_invert, METH_VARARGS, doc_invertg },
    { "isqrt", Pympz_isqrt, METH_O, doc_mpz_isqrt },
    { "isqrt_rem", Pympz_isqrt_rem, METH_VARARGS, doc_mpz_isqrt_rem },
    { "is_bpsw_prp", GMPY_mpz_is_bpsw_prp, METH_VARARGS, doc_mpz_is_bpsw_prp },
    { "is_even", Pympz_is_even, METH_O, doc_is_eveng },
    { "is_euler_prp", GMPY_mpz_is_euler_prp, METH_VARARGS, doc_mpz_is_euler_prp },
    { "is_extra_strong_lucas_prp", GMPY_mpz_is_extrastronglucas_prp, METH_VARARGS, doc_mpz_is_extrastronglucas_prp },
    { "is_fermat_prp", GMPY_mpz_is_fermat_prp, METH_VARARGS, doc_mpz_is_fermat_prp },
    { "is_fibonacci_prp", GMPY_mpz_is_fibonacci_prp, METH_VARARGS, doc_mpz_is_fibonacci_prp },
    { "is_lucas_prp", GMPY_mpz_is_lucas_prp, METH_VARARGS, doc_mpz_is_lucas_prp },
    { "is_odd", Pympz_is_odd, METH_O, doc_is_oddg },
    { "is_power", Pympz_is_power, METH_O, doc_is_powerg },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primeg },
    { "is_selfridge_prp", GMPY_mpz_is_selfridge_prp, METH_VARARGS, doc_mpz_is_selfridge_prp },
    { "is_square", Pympz_is_square, METH_O, doc_is_squareg },
    { "is_strong_prp", GMPY_mpz_is_strong_prp, METH_VARARGS, doc_mpz_is_strong_prp },
    { "is_strong_bpsw_prp", GMPY_mpz_is_strongbpsw_prp, METH_VARARGS, doc_mpz_is_strongbpsw_prp },
    { "is_strong_lucas_prp", GMPY_mpz_is_stronglucas_prp, METH_VARARGS, doc_mpz_is_stronglucas_prp },
    { "is_strong_selfridge_prp", GMPY_mpz_is_strongselfridge_prp, METH_VARARGS, doc_mpz_is_strongselfridge_prp },
    { "jacobi", Pympz_jacobi, METH_VARARGS, doc_jacobig },
    { "kronecker", Pympz_kronecker, METH_VARARGS, doc_kroneckerg },
    { "lcm", Pygmpy_lcm, METH_VARARGS, doc_lcm },
    { "legendre", Pympz_legendre, METH_VARARGS, doc_legendreg },
    { "license", GMPy_get_license, METH_NOARGS, GMPy_doc_license },
    { "lucas", Pygmpy_lucas, METH_O, doc_lucas },
    { "lucasu", GMPY_mpz_lucasu, METH_VARARGS, doc_mpz_lucasu },
    { "lucasu_mod", GMPY_mpz_lucasu_mod, METH_VARARGS, doc_mpz_lucasu_mod },
    { "lucasv", GMPY_mpz_lucasv, METH_VARARGS, doc_mpz_lucasv },
    { "lucasv_mod", GMPY_mpz_lucasv_mod, METH_VARARGS, doc_mpz_lucasv_mod },
    { "lucas2", Pygmpy_lucas2, METH_O, doc_lucas2 },
    { "mod", GMPy_Context_Mod, METH_VARARGS, GMPy_doc_mod },
    { "mp_version", GMPy_get_mp_version, METH_NOARGS, GMPy_doc_mp_version },
    { "mp_limbsize", GMPy_get_mp_limbsize, METH_NOARGS, GMPy_doc_mp_limbsize },
    { "mpc_version", GMPy_get_mpc_version, METH_NOARGS, GMPy_doc_mpc_version },
    { "mpfr_version", GMPy_get_mpfr_version, METH_NOARGS, GMPy_doc_mpfr_version },
    { "mpq", (PyCFunction)Pygmpy_mpq, METH_VARARGS | METH_KEYWORDS, doc_mpq },
    { "mpq_from_old_binary", GMPy_MPQ_From_Old_Binary, METH_O, doc_mpq_from_old_binary },
    { "mpz", (PyCFunction)Pygmpy_mpz, METH_VARARGS | METH_KEYWORDS, doc_mpz },
    { "mpz_from_old_binary", GMPy_MPZ_From_Old_Binary, METH_O, doc_mpz_from_old_binary },
    { "mpz_random", GMPy_MPZ_random_Function, METH_VARARGS, GMPy_doc_mpz_random_function },
    { "mpz_rrandomb", GMPy_MPZ_rrandomb_Function, METH_VARARGS, GMPy_doc_mpz_rrandomb_function },
    { "mpz_urandomb", GMPy_MPZ_urandomb_Function, METH_VARARGS, GMPy_doc_mpz_urandomb_function },
    { "mul", GMPy_Context_Mul, METH_VARARGS, GMPy_doc_function_mul },
    { "next_prime", Pympz_next_prime, METH_O, doc_next_primeg },
    { "numer", Pympq_numer, METH_VARARGS, doc_numerg },
    { "num_digits", Pympz_num_digits, METH_VARARGS, doc_num_digitsg },
    { "pack", GMPy_MPZ_pack, METH_VARARGS, doc_pack },
    { "popcount", GMPy_MPZ_popcount, METH_O, doc_popcount },
    { "powmod", GMPy_Integer_PowMod, METH_VARARGS, GMPy_doc_integer_powmod },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivg },
    { "remove", Pympz_remove, METH_VARARGS, doc_removeg },
    { "iroot", Pympz_iroot, METH_VARARGS, doc_mpz_iroot },
    { "iroot_rem", Pympz_iroot_rem, METH_VARARGS, doc_mpz_iroot_rem },
    { "random_state", GMPy_RandomState_Factory, METH_VARARGS, GMPy_doc_random_state_factory },
    { "set_cache", GMPy_set_cache, METH_VARARGS, GMPy_doc_set_cache },
    { "sign", GMPy_Context_Sign, METH_O, GMPy_doc_function_sign },
    { "square", GMPy_Context_Square, METH_O, GMPy_doc_function_square },
    { "sub", GMPy_Context_Sub, METH_VARARGS, GMPy_doc_sub },
    { "to_binary", GMPy_MPANY_To_Binary, METH_O, doc_to_binary },
    { "t_div", GMPy_MPZ_t_div, METH_VARARGS, doc_t_div },
    { "t_div_2exp", GMPy_MPZ_t_div_2exp, METH_VARARGS, doc_t_div_2exp },
    { "t_divmod", GMPy_MPZ_t_divmod, METH_VARARGS, doc_t_divmod },
    { "t_divmod_2exp", GMPy_MPZ_t_divmod_2exp, METH_VARARGS, doc_t_divmod_2exp },
    { "t_mod", GMPy_MPZ_t_mod, METH_VARARGS, doc_t_mod },
    { "t_mod_2exp", GMPy_MPZ_t_mod_2exp, METH_VARARGS, doc_t_mod_2exp },
    { "unpack", GMPy_MPZ_unpack, METH_VARARGS, doc_unpack },
    { "version", GMPy_get_version, METH_NOARGS, GMPy_doc_version },
    { "xbit_mask", Pyxmpz_xbit_mask, METH_O, doc_xbit_maskg },
    { "xmpz", (PyCFunction)Pygmpy_xmpz, METH_VARARGS | METH_KEYWORDS, doc_xmpz },
    { "_mpmath_normalize", Pympz_mpmath_normalize, METH_VARARGS, doc_mpmath_normalizeg },
    { "_mpmath_create", Pympz_mpmath_create, METH_VARARGS, doc_mpmath_createg },

    { "acos", GMPy_Context_Acos, METH_O, GMPy_doc_function_acos },
    { "acosh", GMPy_Context_Acosh, METH_O, GMPy_doc_function_acosh },
    { "ai", GMPy_Context_Ai, METH_O, GMPy_doc_function_ai },
    { "agm", GMPy_Context_AGM, METH_VARARGS, GMPy_doc_function_agm },
    { "asin", GMPy_Context_Asin, METH_O, GMPy_doc_function_asin },
    { "asinh", GMPy_Context_Asinh, METH_O, GMPy_doc_function_asinh },
    { "atan", GMPy_Context_Atan, METH_O, GMPy_doc_function_atan },
    { "atanh", GMPy_Context_Atanh, METH_O, GMPy_doc_function_atanh },
    { "atan2", GMPy_Context_Atan2, METH_VARARGS, GMPy_doc_function_atan2 },
    { "cbrt", GMPy_Context_Cbrt, METH_O, GMPy_doc_function_cbrt },
    { "ceil", Pympfr_ceil, METH_O, doc_g_mpfr_ceil },
    { "check_range", Pympfr_check_range, METH_O, doc_g_mpfr_check_range },
    { "const_catalan", (PyCFunction)GMPy_Function_Const_Catalan, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_catalan },
    { "const_euler", (PyCFunction)GMPy_Function_Const_Euler, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_euler },
    { "const_log2", (PyCFunction)GMPy_Function_Const_Log2, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_log2 },
    { "const_pi", (PyCFunction)GMPy_Function_Const_Pi, METH_VARARGS | METH_KEYWORDS, GMPy_doc_function_const_pi },
    { "context", (PyCFunction)GMPy_CTXT_Context, METH_VARARGS | METH_KEYWORDS, GMPy_doc_context },
    { "copy_sign", GMPy_MPFR_copy_sign, METH_VARARGS, GMPy_doc_mpfr_copy_sign },
    { "cos", GMPy_Context_Cos, METH_O, GMPy_doc_function_cos },
    { "cosh", GMPy_Context_Cosh, METH_O, GMPy_doc_function_cosh },
    { "cot", GMPy_Context_Cot, METH_O, GMPy_doc_function_cot },
    { "coth", GMPy_Context_Coth, METH_O, GMPy_doc_function_coth },
    { "csc", GMPy_Context_Csc, METH_O, GMPy_doc_function_csc },
    { "csch", GMPy_Context_Csch, METH_O, GMPy_doc_function_csch },
    { "degrees", GMPy_Context_Degrees, METH_O, GMPy_doc_function_degrees },
    { "digamma", GMPy_Context_Digamma, METH_O, GMPy_doc_function_digamma },
    { "div_2exp", GMPy_Context_Div_2exp, METH_VARARGS, GMPy_doc_function_div_2exp },
    { "eint", GMPy_Context_Eint, METH_O, GMPy_doc_function_eint },
    { "erf", GMPy_Context_Erf, METH_O, GMPy_doc_function_erf },
    { "erfc", GMPy_Context_Erfc, METH_O, GMPy_doc_function_erfc },
    { "exp", GMPy_Context_Exp, METH_O, GMPy_doc_function_exp },
    { "expm1", GMPy_Context_Expm1, METH_O, GMPy_doc_function_expm1 },
    { "exp10", GMPy_Context_Exp10, METH_O, GMPy_doc_function_exp10 },
    { "exp2", GMPy_Context_Exp2, METH_O, GMPy_doc_function_exp2 },
    { "f2q", GMPy_Context_F2Q, METH_VARARGS, GMPy_doc_function_f2q },
    { "factorial", Pympfr_factorial, METH_O, doc_g_mpfr_factorial },
    { "floor", Pympfr_floor, METH_O, doc_g_mpfr_floor },
    { "fma", GMPy_Context_FMA, METH_VARARGS, GMPy_doc_function_fma },
    { "fms", GMPy_Context_FMS, METH_VARARGS, GMPy_doc_function_fms },
    { "fmod", GMPy_Context_Fmod, METH_VARARGS, GMPy_doc_function_fmod },
    { "frac", GMPy_Context_Frac, METH_O, GMPy_doc_function_frac },
    { "frexp", Pympfr_frexp, METH_O, doc_g_mpfr_frexp },
    { "fsum", Pympfr_fsum, METH_O, doc_g_mpfr_fsum },
    { "gamma", GMPy_Context_Gamma, METH_O, GMPy_doc_function_gamma },
    { "get_context", GMPy_CTXT_Get, METH_NOARGS, GMPy_doc_get_context },
    { "get_emax_max", GMPy_MPFR_get_emax_max, METH_NOARGS, GMPy_doc_mpfr_get_emax_max },
    { "get_emin_min", GMPy_MPFR_get_emin_min, METH_NOARGS, GMPy_doc_mpfr_get_emin_min },
    { "get_exp", GMPy_MPFR_get_exp, METH_O, GMPy_doc_mpfr_get_exp },
    { "get_max_precision", GMPy_MPFR_get_max_precision, METH_NOARGS, GMPy_doc_mpfr_get_max_precision },
    { "hypot", GMPy_Context_Hypot, METH_VARARGS, GMPy_doc_function_hypot },
    { "ieee", GMPy_CTXT_ieee, METH_O, GMPy_doc_context_ieee },
    { "inf", GMPy_MPFR_set_inf, METH_VARARGS, GMPy_doc_mpfr_set_inf },
    { "is_finite", GMPy_Context_Is_Finite, METH_O, GMPy_doc_function_is_finite },
    { "is_infinite", GMPy_Context_Is_Infinite, METH_O, GMPy_doc_function_is_infinite },
    { "is_integer", GMPy_Context_Is_Integer, METH_O, GMPy_doc_function_is_integer },
    { "is_lessgreater", GMPy_Context_Is_LessGreater, METH_VARARGS, GMPy_doc_function_is_lessgreater },
    { "is_nan", GMPy_Context_Is_NAN, METH_O, GMPy_doc_function_is_nan },
    { "is_regular", GMPy_Context_Is_Regular, METH_O, GMPy_doc_function_is_regular },
    { "is_signed", GMPy_Context_Is_Signed, METH_O, GMPy_doc_function_is_signed },
    { "is_unordered", GMPy_Context_Is_Unordered, METH_VARARGS, GMPy_doc_function_is_unordered },
    { "is_zero", GMPy_Context_Is_Zero, METH_O, GMPy_doc_function_is_zero },
    { "jn", GMPy_Context_Jn, METH_VARARGS, GMPy_doc_function_jn },
    { "j0", GMPy_Context_J0, METH_O, GMPy_doc_function_j0 },
    { "j1", GMPy_Context_J1, METH_O, GMPy_doc_function_j1 },
    { "lgamma", Pympfr_lgamma, METH_O, doc_g_mpfr_lgamma },
    { "li2", GMPy_Context_Li2, METH_O, GMPy_doc_function_li2 },
    { "lngamma", GMPy_Context_Lngamma, METH_O, GMPy_doc_function_lngamma },
    { "local_context", (PyCFunction)GMPy_CTXT_Local, METH_VARARGS | METH_KEYWORDS, GMPy_doc_local_context },
    { "log", GMPy_Context_Log, METH_O, GMPy_doc_function_log },
    { "log1p", GMPy_Context_Log1p, METH_O, GMPy_doc_function_log1p },
    { "log10", GMPy_Context_Log10, METH_O, GMPy_doc_function_log10 },
    { "log2", GMPy_Context_Log2, METH_O, GMPy_doc_function_log2 },
    { "maxnum", GMPy_Context_Maxnum, METH_VARARGS, GMPy_doc_function_maxnum },
    { "minnum", GMPy_Context_Minnum, METH_VARARGS, GMPy_doc_function_minnum },
    { "modf", Pympfr_modf, METH_O, doc_g_mpfr_modf },
    { "mpfr", (PyCFunction)GMPy_MPFR_Factory, METH_VARARGS | METH_KEYWORDS, GMPy_doc_mpfr_factory },
    { "mpfr_from_old_binary", GMPy_MPFR_From_Old_Binary, METH_O, doc_mpfr_from_old_binary },
    { "mpfr_random", GMPy_MPFR_random_Function, METH_VARARGS, GMPy_doc_mpfr_random_function },
    { "mpfr_grandom", GMPy_MPFR_grandom_Function, METH_VARARGS, GMPy_doc_mpfr_grandom_function },
    { "mul_2exp", GMPy_Context_Mul_2exp, METH_VARARGS, GMPy_doc_function_mul_2exp },
    { "nan", GMPy_MPFR_set_nan, METH_NOARGS, GMPy_doc_mpfr_set_nan },
    { "next_above", Pympfr_nextabove, METH_O, doc_g_mpfr_nextabove },
    { "next_below", Pympfr_nextbelow, METH_O, doc_g_mpfr_nextbelow },
    { "next_toward", Pympfr_nexttoward, METH_VARARGS, doc_g_mpfr_nexttoward },
    { "radians", GMPy_Context_Radians, METH_O, GMPy_doc_function_radians },
    { "rec_sqrt", GMPy_Context_RecSqrt, METH_O, GMPy_doc_function_rec_sqrt },
    { "reldiff", Pympfr_reldiff, METH_VARARGS, doc_g_mpfr_reldiff },
    { "remainder", GMPy_Context_Remainder, METH_VARARGS, GMPy_doc_function_remainder },
    { "remquo", Pympfr_remquo, METH_VARARGS, doc_g_mpfr_remquo },
    { "rint", GMPy_Context_Rint, METH_O, GMPy_doc_function_rint },
    { "rint_ceil", GMPy_Context_RintCeil, METH_O, GMPy_doc_function_rint_ceil },
    { "rint_floor", GMPy_Context_RintFloor, METH_O, GMPy_doc_function_rint_floor },
    { "rint_round", GMPy_Context_RintRound, METH_O, GMPy_doc_function_rint_round },
    { "rint_trunc", GMPy_Context_RintTrunc, METH_O, GMPy_doc_function_rint_trunc },
    { "root", GMPy_Context_Root, METH_VARARGS, GMPy_doc_function_root },
    { "round_away", Pympfr_round_away, METH_O, doc_g_mpfr_round_away },
    { "round2", Pympfr_round2, METH_VARARGS, doc_g_mpfr_round2 },
    { "sec", GMPy_Context_Sec, METH_O, GMPy_doc_function_sec },
    { "sech", GMPy_Context_Sech, METH_O, GMPy_doc_function_sech },
    { "set_context", GMPy_CTXT_Set, METH_O, GMPy_doc_set_context },
    { "set_exp", GMPy_MPFR_set_exp, METH_VARARGS, GMPy_doc_mpfr_set_exp },
    { "set_sign", GMPy_MPFR_set_sign, METH_VARARGS, GMPy_doc_mpfr_set_sign },
    { "sin", GMPy_Context_Sin, METH_O, GMPy_doc_function_sin },
    { "sin_cos", GMPy_Context_Sin_Cos, METH_O, GMPy_doc_function_sin_cos },
    { "sinh", GMPy_Context_Sinh, METH_O, GMPy_doc_function_sinh },
    { "sinh_cosh", GMPy_Context_Sinh_Cosh, METH_O, GMPy_doc_function_sinh_cosh },
    { "sqrt", GMPy_Context_Sqrt, METH_O, GMPy_doc_function_sqrt },
    { "tan", GMPy_Context_Tan, METH_O, GMPy_doc_function_tan },
    { "tanh", GMPy_Context_Tanh, METH_O, GMPy_doc_function_tanh },
    { "trunc", Pympfr_trunc, METH_O, doc_g_mpfr_trunc},
    { "yn", GMPy_Context_Yn, METH_VARARGS, GMPy_doc_function_yn },
    { "y0", GMPy_Context_Y0, METH_O, GMPy_doc_function_y0 },
    { "y1", GMPy_Context_Y1, METH_O, GMPy_doc_function_y1 },
    { "zero", GMPy_MPFR_set_zero, METH_VARARGS, GMPy_doc_mpfr_set_zero },
    { "zeta", GMPy_Context_Zeta, METH_O, GMPy_doc_function_zeta },

    { "mpc", (PyCFunction)GMPy_MPC_Factory, METH_VARARGS | METH_KEYWORDS, GMPy_doc_mpc_factory },
    { "mpc_random", GMPy_MPC_random_Function, METH_VARARGS, GMPy_doc_mpc_random_function },
    { "norm", GMPy_Context_Norm, METH_O, GMPy_doc_function_norm },
    { "polar", GMPy_Context_Polar, METH_O, GMPy_doc_function_polar },
    { "phase", GMPy_Context_Phase, METH_O, GMPy_doc_function_phase },
    { "proj", GMPy_Context_Proj, METH_O, GMPy_doc_function_proj },
    { "rect", GMPy_Context_Rect, METH_VARARGS, GMPy_doc_function_rect },
    { NULL, NULL, 1}
};

/* The custom memory allocation routines either use PyMem_* or the standard
 * libraries. See gmpy.h for defines.
 */

static void *
gmpy_allocate(size_t size)
{
    void *res;

    if (!(res = GMPY_MALLOC(size)))
        Py_FatalError("Insufficient memory");

    return res;
}

static void *
gmpy_reallocate(void *ptr, size_t old_size, size_t new_size)
{
    void *res;

    if (!(res = GMPY_REALLOC(ptr, new_size)))
        Py_FatalError("Insufficient memory");

    return res;
}

static void
gmpy_free( void *ptr, size_t size)
{
    GMPY_FREE(ptr);
}

static char _gmpy_docs[] =
"gmpy2 2.1.0a0 - General Multiple-precision arithmetic for Python\n"
"\n"
"gmpy2 supports several multiple-precision libraries. Integer and\n"
"rational arithmetic is provided by either the GMP or MPIR libraries.\n"
"Real floating-point arithmetic is provided by the MPFR library.\n"
"Complex floating-point arithmetic is provided by the MPC library.\n"
"\n"
"The integer type 'mpz' has comparable functionality to Python's\n"
"builtin integers, but is faster for operations on large numbers.\n"
"A wide variety of additional functions are provided:\n"
"      - bit manipulations\n"
"      - GCD, Extended GCD, LCM\n"
"      - Fibonacci and Lucas sequences\n"
"      - primality testing\n"
"      - powers and integer Nth roots\n"
"\n"
"The rational type 'mpq' is equivalent to Python's fractions\n"
"module, but is faster.\n"
"\n"
"The real type 'mpfr' and complex type 'mpc' provide multiple-\n"
"precision real and complex numbers with user-definable precision,\n"
"rounding, and exponent range. All the advanced functions from the\n"
"MPFR and MPC libraries are available.\n\
";

/* Notes on Python 3.x support: Full support for PEP-3121 has not been
 * implemented. No per-module state has been defined.
 */

#ifdef PY3
#define INITERROR return NULL
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "gmpy2",
        _gmpy_docs,
        -1, /*sizeof(struct module_state) */
        Pygmpy_methods,
        NULL,
        NULL, /* gmpy_traverse */
        NULL, /* gmpy_clear */
        NULL
};
#else
#define INITERROR return
#endif

#ifdef PY3
PyMODINIT_FUNC PyInit_gmpy2(void)
#else
PyMODINIT_FUNC initgmpy2(void)
#endif
{
    PyObject* gmpy_module = NULL;
    PyObject* copy_reg_module = NULL;
    PyObject *temp = NULL;

    /* Validate the sizes of the various typedef'ed integer types. */
    if (sizeof(mp_limb_t) != sizeof(mpir_si)) {
        SYSTEM_ERROR("Size of mp_limb_t and mpir_si not compatible");
        INITERROR;
    }
    if (sizeof(mp_bitcnt_t) != sizeof(size_t)) {
        SYSTEM_ERROR("Size of mp_bitcnt_t and size_t not compatible");
        INITERROR;
    }
    if (sizeof(mp_size_t) != sizeof(size_t)) {
        SYSTEM_ERROR("Size of mp_size_t and size_t not compatible");
        INITERROR;
    }
    if (sizeof(mpfr_prec_t) != sizeof(long)) {
        SYSTEM_ERROR("Size of mpfr_prec_t and long not compatible");
        INITERROR;
    }
    if (sizeof(mpfr_exp_t) != sizeof(long)) {
        SYSTEM_ERROR("Size of mpfr_exp_t and long not compatible");
        INITERROR;
    }

    /* Initialize the types. */
    if (PyType_Ready(&MPZ_Type) < 0)
        INITERROR;
    if (PyType_Ready(&MPQ_Type) < 0)
        INITERROR;
    if (PyType_Ready(&XMPZ_Type) < 0)
        INITERROR;
    if (PyType_Ready(&GMPyIter_Type) < 0)
        INITERROR;
    if (PyType_Ready(&MPFR_Type) < 0)
        INITERROR;
    if (PyType_Ready(&CTXT_Type) < 0)
        INITERROR;
    if (PyType_Ready(&CTXT_Manager_Type) < 0)
        INITERROR;
    if (PyType_Ready(&MPC_Type) < 0)
        INITERROR;

    /* Initialize the custom memory handlers. */
    mp_set_memory_functions(gmpy_allocate, gmpy_reallocate, gmpy_free);

    /* Initialize object caching. */
    set_zcache();
    set_gmpympzcache();
    set_gmpympqcache();
    set_gmpyxmpzcache();
    set_gmpympfrcache();
    set_gmpympccache();

    /* Initialize exceptions. */
    GMPyExc_GmpyError = PyErr_NewException("gmpy2.gmpyError",
                                           PyExc_ArithmeticError, NULL);
    if (!GMPyExc_GmpyError)
        INITERROR;

    GMPyExc_Erange = PyErr_NewException("gmpy2.RangeError",
                                        GMPyExc_GmpyError, NULL);
    if (!GMPyExc_Erange)
        INITERROR;

    GMPyExc_Inexact = PyErr_NewException("gmpy2.InexactResultError",
                                         GMPyExc_GmpyError, NULL);
    if (!GMPyExc_Inexact)
        INITERROR;

    GMPyExc_Overflow = PyErr_NewException("gmpy2.OverflowResultError",
                                          GMPyExc_Inexact, NULL);
    if (!GMPyExc_Overflow)
        INITERROR;

    GMPyExc_Underflow = PyErr_NewException("gmpy2.UnderflowResultError",
                                           GMPyExc_Inexact, NULL);
    if (!GMPyExc_Underflow)
        INITERROR;

    temp = PyTuple_Pack(2, GMPyExc_GmpyError, PyExc_ValueError);
    if (!temp)
        INITERROR;
    GMPyExc_Invalid = PyErr_NewException("gmpy2.InvalidOperationError",
                                         temp, NULL);
    Py_DECREF(temp);
    if (!GMPyExc_Invalid)
        INITERROR;

    temp = PyTuple_Pack(2, GMPyExc_GmpyError, PyExc_ZeroDivisionError);
    if (!temp)
        INITERROR;
    GMPyExc_DivZero = PyErr_NewException("gmpy2.DivisionByZeroError",
                                         temp, NULL);
    Py_DECREF(temp);
    if (!GMPyExc_DivZero)
        INITERROR;


#ifdef PY3
    gmpy_module = PyModule_Create(&moduledef);
#else
    gmpy_module = Py_InitModule3("gmpy2", Pygmpy_methods, _gmpy_docs);
#endif

    if (gmpy_module == NULL)
        INITERROR;

    /* Initialize thread local contexts. */
#ifdef WITHOUT_THREADS
    module_context = (CTXT_Object*)GMPy_CTXT_New();
    if (!module_context)
        INITERROR;
    Py_INCREF(Py_False);
    if (PyModule_AddObject(gmpy_module, "HAVE_THREADS", Py_False) < 0) {
        Py_DECREF(Py_False);
        INITERROR;
    }
#else
    tls_context_key = PyUnicode_FromString("__GMPY2_CTX__");
    Py_INCREF(Py_True);
    if (PyModule_AddObject(gmpy_module, "HAVE_THREADS", Py_True) < 0) {
        Py_DECREF(Py_True);
        INITERROR;
    }
#endif

    /* Add the constants for defining rounding modes. */
    if (PyModule_AddIntConstant(gmpy_module, "RoundToNearest", MPFR_RNDN) < 0)
        INITERROR;
    if (PyModule_AddIntConstant(gmpy_module, "RoundToZero", MPFR_RNDZ) < 0)
        INITERROR;
    if (PyModule_AddIntConstant(gmpy_module, "RoundUp", MPFR_RNDU) < 0)
        INITERROR;
    if (PyModule_AddIntConstant(gmpy_module, "RoundDown", MPFR_RNDD) < 0)
        INITERROR;
    if (PyModule_AddIntConstant(gmpy_module, "RoundAwayZero", MPFR_RNDA) < 0)
        INITERROR;
    if (PyModule_AddIntConstant(gmpy_module, "Default", GMPY_DEFAULT) < 0)
        INITERROR;

    /* Add the exceptions. */
    Py_INCREF(GMPyExc_DivZero);
    if (PyModule_AddObject(gmpy_module, "DivisionByZeroError", GMPyExc_DivZero) < 0) {
        Py_DECREF(GMPyExc_DivZero);
        INITERROR;
    }
    Py_INCREF(GMPyExc_Inexact);
    if (PyModule_AddObject(gmpy_module, "InexactResultError", GMPyExc_Inexact) < 0) {
        Py_DECREF(GMPyExc_Inexact);
        INITERROR;
    }
    Py_INCREF(GMPyExc_Invalid);
    if (PyModule_AddObject(gmpy_module, "InvalidOperationError", GMPyExc_Invalid) < 0 ) {
        Py_DECREF(GMPyExc_Invalid);
        INITERROR;
    }
    Py_INCREF(GMPyExc_Overflow);
    if (PyModule_AddObject(gmpy_module, "OverflowResultError", GMPyExc_Overflow) < 0) {
        Py_DECREF(GMPyExc_Overflow);
        INITERROR;
    }
    Py_INCREF(GMPyExc_Underflow);
    if (PyModule_AddObject(gmpy_module, "UnderflowResultError", GMPyExc_Underflow) < 0) {
        Py_DECREF(GMPyExc_Underflow);
        INITERROR;
    }
    Py_INCREF(GMPyExc_Erange);
    if (PyModule_AddObject(gmpy_module, "RangeError", GMPyExc_Erange) < 0) {
        Py_DECREF(GMPyExc_Erange);
        INITERROR;
    }

    /* Add support for pickling. */
#ifdef PY3
    copy_reg_module = PyImport_ImportModule("copyreg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def gmpy2_reducer(x): return (gmpy2.from_binary, (gmpy2.to_binary(x),))\n"
            "copyreg.pickle(type(gmpy2.mpz(0)), gmpy2_reducer)\n"
            "copyreg.pickle(type(gmpy2.xmpz(0)), gmpy2_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpq(0)), gmpy2_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpfr(0)), gmpy2_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpc(0,0)), gmpy2_reducer)\n";
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;
        PyDict_SetItemString(namespace, "copyreg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (!result)
            PyErr_Clear();
        Py_DECREF(namespace);
        Py_DECREF(copy_reg_module);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
    }
#else
    copy_reg_module = PyImport_ImportModule("copy_reg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def gmpy2_reducer(x): return (gmpy2.from_binary, (gmpy2.to_binary(x),))\n"
            "copy_reg.pickle(type(gmpy2.mpz(0)), gmpy2_reducer)\n"
            "copy_reg.pickle(type(gmpy2.xmpz(0)), gmpy2_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpq(0)), gmpy2_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpfr(0)), gmpy2_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpc(0,0)), gmpy2_reducer)\n";
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;

        PyDict_SetItemString(namespace, "copy_reg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (!result)
            PyErr_Clear();
        Py_DECREF(namespace);
        Py_DECREF(copy_reg_module);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
    }
#endif

#ifdef PY3
    return gmpy_module;
#endif
}
