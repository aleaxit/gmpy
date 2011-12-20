/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2.c                                                                 *
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
 * originally written for GMP-2.0 (by AMK...?)
 * Rewritten by Niels Möller, May 1996
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
 *   2.00:
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

#define GMPY_MODULE
#include "gmpy.h"

/* The gmpy_args.h file includes macros that are used for argument
 * processing.
 */
#include "gmpy_args.h"

/* Define the minimum memory amount allocated. 8 has historically been
 * used, but 16 might be better for some applications or 64-bit systems.
 */
#define GMPY_ALLOC_MIN (2 * (GMP_NUMB_BITS >> 3))

/* To prevent excessive memory usage, we don't want to save very large
 * numbers in the cache. The default value specified in the options
 * structure is 128 words (512 bytes on 32-bit platforms, 1024 bytes on
 * 64-bit platforms).
 */
#define MAX_CACHE_LIMBS 16384

/* The maximum number of objects that can be saved in a cache is specified
 * here. The default value is 100.*/
#define MAX_CACHE 1000

/* Add support for PyLong_AsLongAndOverflow to older versions of Python */

#include "py3intcompat.c"

/* Define various macros to deal with differences between Python 2 and 3. */

#if (PY_MAJOR_VERSION == 3)
#define PY3
#define Py2or3String_FromString     PyUnicode_FromString
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

char gmpy_version[] = "2.0.0a3";

char _gmpy_cvs[] = "$Id$";

/*
 * global data declarations
 */

static PyObject *gmpy_module = NULL;

static struct gmpy_global {
    int debug;               /* != 0 if debug messages desired on stderr */
    int cache_size;          /* size of cache, for all caches */
    int cache_obsize;        /* maximum size of the objects that are cached */
} global = {
    0,                       /* debug */
    100,                     /* cache_size */
    128,                     /* cache_obsize */
};

#ifdef WITHMPFR
/* The context manager should really be thread-specific. Until then, gmpy2
 * is NOT thread-safe for mpfr and mpc calculations.
 */

static GMPyContextObject *context = NULL;

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
static PyObject *GMPyExc_ExpBound = NULL;
#endif

/* Support for context manager */

#ifdef WITHMPFR
#  include "gmpy_context.c"
#endif

/* Include fast mpz to/from PyLong conversion from sage. */
#include "mpz_pylong.c"

/* The code for object creation, deletion, and caching is in gmpy_cache.c. */

#include "gmpy_cache.c"

/* Miscellaneous helper functions and simple methods are in gmpy_misc.c. */

#include "gmpy_misc.c"

/* Support for conversion to/from binary representation. */

#include "gmpy_binary.c"

/* CONVERSIONS AND COPIES */

static PympzObject *
Pympz2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pympz_new()))
        mpz_set(newob->z, Pympz_AS_MPZ(self));
    return newob;
}

static PyxmpzObject *
Pyxmpz2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pyxmpz_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PympzObject *
Pyxmpz2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pyxmpz_Check(self));
    if ((newob = Pympz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PyxmpzObject *
Pympz2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PympqObject *
Pympq2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(Pympq_Check(self));
    if ((newob = Pympq_new()))
        mpq_set(newob->q, Pympq_AS_MPQ(self));
    return newob;
}

#ifdef PY2
static PympzObject *
PyInt2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympz_new()))
        mpz_set_si(newob->z, PyInt_AS_LONG(self));
    return newob;
}

static PyxmpzObject *
PyInt2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set_si(newob->z, PyInt_AsLong(self));
    return newob;
}

static PympqObject *
PyInt2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympq_new()))
        mpq_set_si(newob->q, PyInt_AsLong(self), 1);
    return newob;
}
#endif

static PympzObject *
PyFloat2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pympz_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(newob->z, d);
    }
    return newob;
}

static PyxmpzObject *
PyFloat2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pyxmpz_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(newob->z, d);
    }
    return newob;
}

static PympqObject *
Pympz2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pympq_new()))
        mpq_set_z(newob->q, Pympz_AS_MPZ(self));
    return newob;
}

static PympqObject *
Pyxmpz2Pympq(PyObject * obj)
{
    PympqObject *newob;

    assert(Pyxmpz_Check(obj));
    if ((newob = Pympq_new()))
        mpq_set_z(newob->q, Pyxmpz_AS_MPZ(obj));
    return newob;
}

static PympzObject *
Pympq2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pympq_Check(self));
    if ((newob = Pympz_new()))
        mpz_set_q(newob->z, Pympq_AS_MPQ(self));
    return newob;
}

static PyxmpzObject *
Pympq2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pympq_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set_q(newob->z, Pympq_AS_MPQ(self));
    return newob;
}

/* For fast conversion between PyLong and mpz, we use code located in
 * mpz_pylong.c.
 */
static PympzObject *
PyLong2Pympz(PyObject * obj)
{
    PympzObject *newob;
    if (!(newob = Pympz_new()))
        return NULL;
    mpz_set_PyLong(Pympz_AS_MPZ(newob), obj);
    return newob;
}

static PyxmpzObject *
PyLong2Pyxmpz(PyObject * obj)
{
    PyxmpzObject *newob;
    if (!(newob = Pyxmpz_new()))
        return NULL;
    mpz_set_PyLong(Pyxmpz_AS_MPZ(newob), obj);
    return newob;
}

static PympqObject *
PyLong2Pympq(PyObject *self)
{
    PympqObject *newob;
    PyObject *temp = (PyObject*)PyLong2Pympz(self);

    if (!temp)
        return NULL;
    newob = Pympz2Pympq(temp);
    Py_DECREF(temp);
    return newob;
}

static PympqObject *
PyFloat2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pympq_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpq' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpq' does not support Infinity");
            return NULL;
        }
        mpq_set_d(newob->q, d);
    }
    return newob;
}

/* mpz conversion from string includes from-binary (base-256 LSB string
 * of bytes) and 'true' from-string (bases 2 to 62; bases 8 and 16 are
 * special -- decorations of leading 0/0x are allowed (not required).
 *
 * Binary form was previously (0.6) limited to >=0 values; now (0.7)
 * extended to <0 values, by adding one last sign-byte, 0xFF for <0,
 * 0x00 for >0 (the latter only if the #bits is exact multiple of 8).
 */

/* mpz_set_PyStr returns -1 on error, 1 if successful. */

static int
mpz_set_PyStr(mpz_ptr z, PyObject *s, long base)
{
    unsigned char *cp;
    Py_ssize_t len;
    int i;
    PyObject *ascii_str = NULL;

    assert(PyStrOrUnicode_Check(s));

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = (unsigned char*)PyBytes_AsString(s);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return -1;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }

    /* Don't allow NULL characters */
    for (i=0; i<len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            Py_XDECREF(ascii_str);
            return -1;
        }
    }
    /* delegate rest to GMP's _set_str function */
    if (base==0) {
        if (cp[0]=='0') {
            if (cp[1]=='b') {
                base = 2;
                cp+=2;
            }
            else if (cp[1]=='o') {
                base = 8;
                cp+=2;
            }
            else if (cp[1]=='x') {
                base = 16;
                cp+=2;
            }
            else {
                base = 10;
            }
        }
        else {
            base = 10;
        }
    }
    if (-1 == mpz_set_str(z, (char*)cp, base)) {
        VALUE_ERROR("invalid digits");
        Py_XDECREF(ascii_str);
        return -1;
    }
    Py_XDECREF(ascii_str);
    return 1;
}

static PympzObject *
PyStr2Pympz(PyObject *s, long base)
{
    PympzObject *newob;

    assert(PyStrOrUnicode_Check(s));
    if (!(newob = Pympz_new()))
        return NULL;
    if (mpz_set_PyStr(newob->z, s, base) == -1) {
        Py_DECREF((PyObject*)newob);
        return NULL;
    }
    return newob;
}

static PyxmpzObject *
PyStr2Pyxmpz(PyObject *s, long base)
{
    PyxmpzObject *newob;

    assert(PyStrOrUnicode_Check(s));
    if (!(newob = Pyxmpz_new()))
        return NULL;
    if (mpz_set_PyStr(newob->z, s, base) == -1) {
        Py_DECREF((PyObject*)newob);
        return NULL;
    }
    return newob;
}

/*
 * mpq conversion from string includes from-binary (base-256 LSB string
 * of bytes) and 'true' from-string (bases 2 to 62; bases 8 and 16 are
 * special -- decorations of leading 0/0x are allowed (but not required).
 * For 'true-bases' 2..62 there is a '/' separator between numerator and
 * denominator (if none, just numerator!); decimal point NOT allowed.
 *
 * Added in gmpy 1.02: also support a string of the form '12.34', i.e.,
 * WITH a decimal point and WITHOUT a slash
 *
 * Binary-form: 4-byte numerator length (upper bit set if <0), then
 * numerator (as above for mpz), then denominator (ditto).
 */
static PympqObject *
PyStr2Pympq(PyObject *stringarg, long base)
{
    PympqObject *newob;
    unsigned char *cp;
    Py_ssize_t len;
    int i;
    PyObject *ascii_str = NULL;

    assert(PyStrOrUnicode_Check(stringarg));
    if (!(newob = Pympq_new()))
        return NULL;
    if (PyBytes_Check(stringarg)) {
        len = PyBytes_Size(stringarg);
        cp = (unsigned char*)PyBytes_AsString(stringarg);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(stringarg);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            Py_DECREF((PyObject*)newob);
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }

    /* Don't allow NULL characters */
    for (i=0; i<len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }
    }
    /* trickily delegate the rest to GMP avoiding allocations/copies */
    {
        char *whereslash = strchr((char*)cp, '/');
        char *wheredot = strchr((char*)cp, '.');
        if (whereslash && wheredot) {
            VALUE_ERROR("illegal string: both . and / found");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }

        if (wheredot) {
            char *counter;
            unsigned long digits = 0;
            if (base != 10) {
                VALUE_ERROR("illegal string: embedded . requires base=10");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }

            counter = wheredot;
            digits = 0;
            *wheredot = ' ';
            while (*++counter != '\0') {
                if (isdigit(*counter))
                    digits++;
            }
            if (-1 == mpz_set_str(mpq_numref(newob->q), (char*)cp, base)) {
                if (wheredot)
                    *wheredot = '.';
                VALUE_ERROR("invalid digits");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }
            mpz_ui_pow_ui(mpq_denref(newob->q), 10, digits);
            mpq_canonicalize(newob->q);
            *wheredot = '.';
            return (PympqObject*)newob;
        }

        if (whereslash)
            *whereslash = 0;
        if (-1 == mpz_set_str(mpq_numref(newob->q), (char*)cp, base)) {
            if (whereslash)
                *whereslash = '/';
            VALUE_ERROR("invalid digits");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }
        if (whereslash) {
            *whereslash = '/';
            if (-1==mpz_set_str(mpq_denref(newob->q), whereslash+1, base)) {
                VALUE_ERROR("invalid digits");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }
            if (0==mpz_sgn(mpq_denref(newob->q))) {
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                ZERO_ERROR("zero denominator in 'mpq'");
                return NULL;
            }
            mpq_canonicalize(newob->q);
        }
        else {
            mpz_set_ui(mpq_denref (newob->q), 1);
        }
    }
    Py_XDECREF(ascii_str);
    return newob;
}

/* For fast mpz to PyLong conversion, we use code located in mpz_pylong.
 */
static PyObject *
Pympz2PyLong(PympzObject *self)
{
    return mpz_get_PyLong(Pympz_AS_MPZ(self));
}

static PyObject *
Pyxmpz2PyLong(PyxmpzObject *self)
{
    return mpz_get_PyLong(Pyxmpz_AS_MPZ(self));
}

static PyObject *
Pympq2PyLong(PympqObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympq2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz2PyLong(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}

static PyObject *
Pympz_To_Integer(PympzObject *self)
{
#ifdef PY3
    return Pympz2PyLong(self);
#else
    if (mpz_fits_slong_p(self->z))
        return PyInt_FromLong(mpz_get_si(self->z));
    else
        return Pympz2PyLong(self);
#endif
}

static PyObject *
Pyxmpz_To_Integer(PyxmpzObject *self)
{
#ifdef PY3
    return Pyxmpz2PyLong(self);
#else
    if (mpz_fits_slong_p(self->z))
        return PyInt_FromLong(mpz_get_si(self->z));
    else
        return Pyxmpz2PyLong(self);
#endif
}

#ifdef PY2
static PyObject *
Pympq2PyInt(PympqObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympq2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz_To_Integer(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}
#endif

static PyObject *
Pympz2PyFloat(PympzObject *self)
{
    double res = mpz_get_d(self->z);

    return PyFloat_FromDouble(res);
}

static PyObject *
Pympq2PyFloat(PympqObject *self)
{
    double res = mpq_get_d(self->q);

    return PyFloat_FromDouble(res);
}

/*
 *  build binary representation of mpz (base-256 little-endian)
 *  Note: design limitation used to forbid binary repr of <0 mpz;
 *  this has now been remedied, but at the price of full compatibility
 *  with files saved in gmpy releases 0.6 and earlier.
 */
static PyObject *
mpz2binary(mpz_ptr z)
{
    size_t size, usize;
    int negative, needtrail;
    char *buffer;
    PyObject *s;

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z); /* Change the sign temporarily! */
    }
    else {
        negative = 0;
    }

    size = mpz_sizeinbase(z, 2);
    needtrail = (size%8)==0;
    usize = size = (size + 7) / 8;
    if (negative || needtrail)
        ++size;

    TEMP_ALLOC(buffer, size);
    buffer[0] = 0x00;
    mpz_export(buffer, NULL, -1, sizeof(char), 0, 0, z);
    if (usize < size) {
        buffer[usize] = negative?0xff:0x00;
    }
    if (negative) {
        mpz_neg(z, z);
    }
    s = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return s;
}

static PyObject *
Pympz2binary(PympzObject *self)
{
    return mpz2binary(self->z);
}

static PyObject *
Pyxmpz2binary(PyxmpzObject *self)
{
    return mpz2binary(self->z);
}

/*
 *  build binary representation of mpq (base-256 little-endian
 *  for num, then den; before those, 4 bytes for _length_ of
 *  numerator, which also encode sign as the single top bit).
 */

static PyObject *
Pympq2binary(PympqObject *self)
{
    size_t sizenum, sizeden, size, sizetemp;
    int negative = 0;
    char *buffer;
    int i;
    PyObject *result;

    if (mpq_sgn(self->q) < 0) {
        negative = 1;
        mpz_abs(mpq_numref(self->q), mpq_numref(self->q));
    }
    else {
        negative = 0;
    }
    assert(mpz_sgn(mpq_denref(self->q)) > 0);

    sizenum = (mpz_sizeinbase(mpq_numref(self->q), 2) + 7) / 8;
    sizeden = (mpz_sizeinbase(mpq_denref(self->q), 2) + 7) / 8;
    size = sizenum + sizeden + 4;

    TEMP_ALLOC(buffer, size);

    sizetemp = sizenum;
    for (i=0; i<4; i++) {
        buffer[i] = (char)(sizetemp & 0xff);
        sizetemp >>= 8;
    }
    if (negative)
        buffer[3] |= 0x80;
    buffer[4] = 0x00;

    mpz_export(buffer+4, NULL, -1, sizeof(char), 0, 0, mpq_numref(self->q));
    mpz_export(buffer+sizenum+4, NULL, -1, sizeof(char), 0, 0, mpq_denref(self->q));
    if (negative) {
        mpz_neg( mpq_numref(self->q), mpq_numref(self->q));
    }
    result = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return result;
}

/* Format an mpz into any base (2 to 62). Bits in the option parameter
 * control various behaviors:
 *   bit 0: if set, output is wrapped with 'mpz(...)'
 *   bit 1: if set, a '+' is included for positive numbers
 *   bit 2: if set, a ' ' is included for positive nubmers
 *   bit 3: if set, a '0b', '0o', or '0x' is included for binary, octal, hex
 *   bit 4: if set, no prefix is included for binary, octal, hex
 *
 * Note: if neither bit 3 or 4 is set, prefixes that match the platform default
 * are included.
 *
 * If base < 0, capital letters are used
 */
static char* ztag = "mpz(";
static PyObject *
mpz_ascii(mpz_t z, int base, int option)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if (!((base == 0) || ((base >= -36) && (base <= -2)) ||
        ((base >= 2) && (base <= 62)) )) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 62");
        return NULL;
    }

    /* Allocate extra space for:
     *
     * minus sign and trailing NULL byte (2)
     * 'mpz()' tag                       (5)
     * '0x' prefix                       (2)
     * 'L' suffix                        (1)
     *                                  -----
     *                                   10
     */
    size = mpz_sizeinbase(z, base) + 11;
    TEMP_ALLOC(buffer, size);

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z);
    }

    p = buffer;
    if (option & 1) {
       strcpy(p, ztag);
       p += strlen(p);
    }

    if (negative)
        *(p++) = '-';
    else if (option & 2)
        *(p++) = '+';
    else if (option & 4)
        *(p++) = ' ';

    if (option & 8) {
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }
    else if (!(option & 24)) {
    #ifdef PY2
        if (base == 8) {
            *(p++) = '0';
        }
    #else
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
    #endif
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }

    mpz_get_str(p, base, z);     /* Doesn't return number of characters */
    p = buffer + strlen(buffer); /* Address of NULL byte */
#ifdef PY2
    if ((option & 1) && !mpz_fits_slong_p(z))
        *(p++) = 'L';
#endif
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';

    result = Py_BuildValue("s", buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}

/* Format an xmpz into any base (2 to 62). If with_tag != 0, the the output
 * is wrapped with 'xmpz(...)'. If with_sign > 0, a plus or minus leading
 * sign is always included. If with_sign < 0, a space is included instead of
 * a plus sign.
 */
static char* xztag = "xmpz(";
static PyObject *
xmpz_ascii(mpz_t z, int base, int option)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if (!((base == 0) || ((base >= -36) && (base <= -2)) ||
        ((base >= 2) && (base <= 62)) )) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 62");
        return NULL;
    }

    /* Allocate extra space for:
     *
     * minus sign and trailing NULL byte (2)
     * 'gmpy2.xmpz()' tag                (6)
     * '0x' prefix                       (2)
     * 'L' suffix                        (1)
     *                                  -----
     *                                   11
     */
    size = mpz_sizeinbase(z, base) + 12;
    TEMP_ALLOC(buffer, size);

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z);
    }

    p = buffer;
    if (option & 1) {
       strcpy(p, xztag);
       p += strlen(p);
    }

    if (negative)
        *(p++) = '-';
    else if (option & 2)
        *(p++) = '+';
    else if (option & 4)
        *(p++) = ' ';

    if (option & 8) {
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }
    else if (!(option & 24)) {
    #ifdef PY2
        if (base == 8) {
            *(p++) = '0';
        }
    #else
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
    #endif
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }

    mpz_get_str(p, base, z);     /* Doesn't return number of characters */
    p = buffer + strlen(buffer); /* Address of NULL byte */
#ifdef PY2
    if ((option & 1) && !mpz_fits_slong_p(z))
        *(p++) = 'L';
#endif
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';

    result = Py_BuildValue("s", buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}

static PyObject *
Pympz_ascii(PympzObject *self, int base, int option)
{
    return mpz_ascii(self->z, base, option);
}

static PyObject *
Pyxmpz_ascii(PyxmpzObject *self, int base, int option)
{
    return xmpz_ascii(self->z, base, option);
}

static int qden_1(mpq_t q)
{
    return 0 == mpz_cmp_ui(mpq_denref(q),1);
}

static PyObject *
Pympq_ascii(PympqObject *self, int base, int option)
{
    PyObject *result = 0, *numstr = 0, *denstr = 0;
    char buffer[50], *p;

    numstr = mpz_ascii(mpq_numref(self->q), base, 0);
    if (!numstr)
        return NULL;

    /* Check if denominator is 1 and no tag is requested. If so, just
     * return the numerator.
     */
    if (!(option & 1) && qden_1(self->q))
        return numstr;

    denstr = mpz_ascii(mpq_denref(self->q), base, 0);
    if (!denstr) {
        Py_DECREF(numstr);
        return NULL;
    }

    /* Build the format string. */
    p = buffer;
    if (option & 1) {
        *(p++) = 'm';
        *(p++) = 'p';
        *(p++) = 'q';
        *(p++) = '(';
    }
#ifdef PY2
    *(p++) = '%';
    *(p++) = 's';
    if (!mpz_fits_slong_p(mpq_numref(self->q)))
        *(p++) = 'L';
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 's';
    if (!mpz_fits_slong_p(mpq_denref(self->q)))
        *(p++) = 'L';
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';
    result = PyString_FromFormat(buffer, PyString_AS_STRING(numstr),
                                 PyString_AS_STRING(denstr));
#else
    *(p++) = '%';
    *(p++) = 'U';
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 'U';
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';
    result = PyUnicode_FromFormat(buffer, numstr, denstr);
#endif
    Py_DECREF(numstr);
    Py_DECREF(denstr);
    return result;
}

static int isFraction(PyObject* obj)
{
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isDecimal(PyObject* obj)
{
    if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) return 1;

    return 0;
}

#ifdef WITHMPC
/* Classify an object as a type of number. If an object is recognized as a
 * number, it must be properly converted by the routines below.
 */

static int isComplex(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (Pympc_Check(obj))       return 1;
    if (PyFloat_Check(obj))     return 1;
    if (PyComplex_Check(obj))   return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}
#endif

#ifdef WITHMPFR
static int isReal(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (PyFloat_Check(obj))     return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}
#endif

static int isRational(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

static int isInteger(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pyxmpz_Check(obj))      return 1;

    return 0;
}

/* Number conversion routines
 *
 * The routines anynum2mpX will attempt to convert any number-like object into
 * into a gmpy object. These routines are intended for construction of mpXs.
 * The accepted number-like objects are:
 *      1) int (Python 2.x)
 *      2) long (Python 2.x and 3.x)
 *      3) float
 *      4) Decimal
 *      5) Fraction
 *      6) other gmpy objects
 *
 * The routine Pympz_From_Integer will only convert integer-like objects into to a
 * gmpy mpz. The accepted integer-like objects are:
 *      1) int
 *      2) long
 *      3) mpz
 *      4) xmpz
 *
 * The routine Pympq_From_Rational will convert an integer- and rational-like
 * object into a gmpy mpq. The accepted objects are:
 *      1) int
 *      2) long
 *      3) Fraction
 *      4) mpz
 *      5) mpq
 *      6) xmpz
 */

/* NOTE: Pympq_From_Decimal returns an invalid mpq object when attempting to
 *       convert a NaN or inifinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 */

static PympqObject*
Pympq_From_Decimal(PyObject* obj)
{
    PympqObject *result;
    PyObject *d_exp, *d_int, *d_sign, *d_is_special;
    long exp;
    mpz_t temp;
    const char *string;

    if (!(result = Pympq_new()))
        return NULL;
    mpq_set_si(result->q, 0, 1);

    d_exp = PyObject_GetAttrString(obj, "_exp");
    d_int = PyObject_GetAttrString(obj, "_int");
    d_sign = PyObject_GetAttrString(obj, "_sign");
    d_is_special = PyObject_GetAttrString(obj, "_is_special");
    if (!d_exp || !d_int || !d_sign || !d_is_special) {
        SYSTEM_ERROR("Object does not appear to be Decimal");
        Py_XDECREF(d_exp);
        Py_XDECREF(d_int);
        Py_XDECREF(d_sign);
        Py_XDECREF(d_is_special);
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    if (PyObject_IsTrue(d_is_special)) {
        string = Py2or3String_AsString(d_exp);
        if (string[0] == 'N' || string[0] == 'n') {
            mpz_set_si(mpq_denref(result->q), 0);
            Py_DECREF(d_exp);
            Py_DECREF(d_int);
            Py_DECREF(d_sign);
            Py_DECREF(d_is_special);
            return result;
        }
        if (string[0] == 'F') {
            if (PyObject_IsTrue(d_sign))
                mpq_set_si(result->q, -1, 0);
            else
                mpq_set_si(result->q, 1, 0);
            Py_DECREF(d_exp);
            Py_DECREF(d_int);
            Py_DECREF(d_sign);
            Py_DECREF(d_is_special);
            return result;
        }
        SYSTEM_ERROR("Cannot convert Decimal to mpq");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_set_PyStr(mpq_numref(result->q), d_int, 10);

    if (PyObject_IsTrue(d_sign))
        mpz_mul_si(mpq_numref(result->q), mpq_numref(result->q), -1);

    exp = PyIntOrLong_AsLong(d_exp);
    if (exp == -1 && PyErr_Occurred()) {
        SYSTEM_ERROR("Decimal _exp is not valid or overflow occurred");
        Py_DECREF((PyObject*)result);
        Py_DECREF(d_exp);
        Py_DECREF(d_int);
        Py_DECREF(d_sign);
        Py_DECREF(d_is_special);
        return NULL;
    }

    mpz_inoc(temp);
    if (exp <= 0)
        mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(-exp));
    else {
        mpz_inoc(temp);
        mpz_ui_pow_ui(temp, 10, (unsigned long)(exp));
        mpz_mul(mpq_numref(result->q), mpq_numref(result->q), temp);
        mpz_cloc(temp);
    }

    mpq_canonicalize(result->q);
    Py_DECREF(d_exp);
    Py_DECREF(d_int);
    Py_DECREF(d_sign);
    Py_DECREF(d_is_special);

    return result;
}

static PympqObject*
Pympq_From_Real(PyObject* obj)
{
    PympqObject* newob = 0;

    if (Pympq_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympqObject *) obj;
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympq(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympq(obj);
#endif
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pympq(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympq(obj);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympq(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympq(obj);
    }
    else if (isDecimal(obj)) {
        if ((newob = Pympq_From_Decimal(obj))) {
            if (!mpz_cmp_si(mpq_denref(newob->q), 0)) {
                if (mpz_get_si(mpq_numref(newob->q)) == 0)
                    VALUE_ERROR("'mpq' does not support NaN");
                else
                    VALUE_ERROR("'mpq' does not support Infinity");
                Py_DECREF((PyObject*)newob);
                newob = NULL;
            }
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }

    return newob;
}

/* Convert an integer or mpz to mpq. */

static PympqObject*
Pympq_From_Rational(PyObject* obj)
{
    PympqObject* newob = 0;

    if (Pympq_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympqObject *) obj;
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympq(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympq(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympq(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympq(obj);
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }

    return newob;
}

static PympzObject*
anynum2Pympz(PyObject* obj)
{
    PympzObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympz_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympzObject *) obj;
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympz(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympz(obj);
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pympz(obj);
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pympz(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympz(obj);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = PyLong2Pympz(s);
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pympz((PyObject *)temp);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }

    return newob;
}

static PyxmpzObject*
anynum2Pyxmpz(PyObject* obj)
{
    PyxmpzObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympz_Check(obj)) {
        newob = Pympz2Pyxmpz(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pyxmpz(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pyxmpz(obj);
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pyxmpz(obj);
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pyxmpz(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pyxmpz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pyxmpz(obj);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = PyLong2Pyxmpz(s);
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pyxmpz((PyObject *)temp);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }

    return newob;
}

/*
 * Convert an Integer-like object (as determined by isInteger) to
 * a Pympz. Returns NULL and raises a TypeError if obj is not an
 * Integer-like object.
 */

static PympzObject*
Pympz_From_Integer(PyObject* obj)
{
    PympzObject* newob = 0;

    if (Pympz_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympzObject*) obj;
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympz(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympz(obj);
    }
    if (!newob)
        TYPE_ERROR("conversion error in Pympz_From_Integer");
    return newob;
}

/*
 * Convert an Integer-like object (as determined by isInteger) to
 * a C long. Returns -1 and raises OverflowError if the the number is
 * too large. Returns -1 and raises TypeError if obj was not an
 * Integer-like object.
 */

#ifndef WIN64
static long
clong_From_Integer(PyObject *obj)
{
    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in clong_From_Integer");
            return -1;
        }
    }
    TYPE_ERROR("conversion error in clong_From_Integer");
    return -1;
}

static unsigned long
culong_From_Integer(PyObject *obj)
{
    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsUnsignedLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_ulong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_ui(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in culong_From_Integer");
            return (unsigned long)-1;
        }
    }
    TYPE_ERROR("conversion error in culong_From_Integer");
    return (unsigned long)-1;
}

#define MP_BITCNT_FROM_INTEGER(obj) culong_From_Integer(obj)

#else
static long long
clonglong_From_Integer(PyObject *obj)
{
    long long val;
    PyObject* temp;

    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsLongLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return (long long)mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            /* This section should only be called on Win64. */
            temp = mpz_get_PyLong(Pympz_AS_MPZ(obj));
            if (!temp) {
                TYPE_ERROR("conversion error in clonglong_From_Integer");
                return -1;
            }
            else {
                val = PyLong_AsLongLong(temp);
                Py_DECREF(temp);
                return val;
            }
        }
    }
    TYPE_ERROR("conversion error in clonglong_From_Integer");
    return -1;
}

static unsigned long long
culonglong_From_Integer(PyObject *obj)
{
    unsigned long long val;
    PyObject* temp;

    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsUnsignedLongLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return (unsigned long long)mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            /* This section should only be called on Win64. */
            temp = mpz_get_PyLong(Pympz_AS_MPZ(obj));
            if (!temp) {
                TYPE_ERROR("conversion error in culonglong_From_Integer");
                return -1;
            }
            else {
                val = PyLong_AsUnsignedLongLong(temp);
                Py_DECREF(temp);
                return val;
            }
        }
    }
    TYPE_ERROR("conversion error in culonglong_From_Integer");
    return -1;
}

#define MP_BITCNT_FROM_INTEGER(obj) culonglong_From_Integer(obj)
#endif

/*
 * Convert an Integer-like object (as determined by isInteger) to
 * a Py_ssize_t. Returns -1 and raises OverflowError if the the number is
 * too large. Returns -1 and raises TypeError if obj was not an
 * Integer-like object.
 */

static Py_ssize_t
ssize_t_From_Integer(PyObject *obj)
{
    Py_ssize_t val;
    PyObject* temp;

    if (PyLong_Check(obj)) {
        return PyLong_AsSsize_t(obj);
    }
#ifdef PY2
    else if (PyInt_Check(obj)) {
        return PyInt_AsSsize_t(obj);
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return (Py_ssize_t)mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            /* This section should only be called on Win64. */
            temp = mpz_get_PyLong(Pympz_AS_MPZ(obj));
            if (!temp) {
                TYPE_ERROR("conversion error in ssize_t_From_Integer");
                return -1;
            }
            else {
                val = PyLong_AsSsize_t(temp);
                Py_DECREF(temp);
                return val;
            }
        }
    }
    TYPE_ERROR("conversion error in ssize_t_From_Integer");
    return -1;
}

/*
 * coerce any number to a mpz
 */
int
Pympz_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympzObject* newob = Pympz_From_Integer(arg);
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpz'");
        return 0;
    }
}

/*
 * coerce any number to a mpq
 */

int
Pympq_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympqObject* newob = Pympq_From_Rational(arg);
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        if (!PyErr_Occurred()) {
            TYPE_ERROR("argument can not be converted to 'mpq'");
        }
        return 0;
    }
}

/* str and repr implementations for mpz */
static PyObject *
Pympz2str(PympzObject *self)
{
    /* base-10, no tag */
    return Pympz_ascii(self, 10, 0);
}

static PyObject *
Pympz2repr(PympzObject *self)
{
    /* base-10, with tag */
    return Pympz_ascii(self, 10, 1);
}

/* str and repr implementations for xmpz */
static PyObject *
Pyxmpz2str(PyxmpzObject *self)
{
    /* base-10, no tag */
    return Pyxmpz_ascii(self, 10, 0);
}

static PyObject *
Pyxmpz2repr(PyxmpzObject *self)
{
    /* base-10, with tag */
    return Pyxmpz_ascii(self, 10, 1);
}

/* str and repr implementations for mpq */
static PyObject *
Pympq2str(PympqObject *self)
{
    /* base-10, no tag */
    return Pympq_ascii(self, 10, 0);
}

static PyObject *
Pympq2repr(PympqObject *self)
{
    /* base-10, with tag */
    return Pympq_ascii(self, 10, 1);
}

/* CONSTRUCTORS */
static char doc_mpz[] = "\
mpz(n):\n\
      builds an 'mpz' object with a numeric value n (truncating n\n\
      to its integer part if it's a Fraction, 'mpq', Decimal, float\n\
      or 'mpfr')\n\
mpz(s,base=0):\n\
      builds an 'mpz' object from a string s made up of digits in the\n\
      given base.  If base=0, binary, octal, or hex Python strings\n\
      are recognized by leading 0b, 0o, or 0x characters, otherwise\n\
      the string is assumed to be decimal.\n\
";
static PyObject *
Pygmpy_mpz(PyObject *self, PyObject *args, PyObject *keywds)
{
    PympzObject *result = 0;
    PyObject *n = 0;
    long base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"n", "base", NULL };

    TRACE("Pygmpy_mpz() called...\n");

    /* Optimize the most common use case */
    argc = PyTuple_Size(args);
    if (argc == 1) {
        n = PyTuple_GetItem(args, 0);
#ifdef WITHMPFR
        if (isReal(n) && !keywds) {
#else
        if ((isRational(n) || PyFloat_Check(n)) && !keywds) {
#endif
            result = anynum2Pympz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("mpz() requires numeric or string argument");
            return (PyObject*)result;
        }
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist,
                                     &n, &base))
        return NULL;

    if ((base!=0) && ((base<2)||(base>62))) {
        VALUE_ERROR("base for mpz() must be 0 or in the "
                    "interval 2 ... 62");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        /* build-from-string (ascii or unicode) */
        result = PyStr2Pympz(n, base);
    }
    else {
        if (argc==2 || (argc == 1 && keywds))
            TYPE_ERROR("mpz() with non-string argument needs exactly "
                       "1 argument");
        else {
            result = anynum2Pympz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("mpz() requires numeric or string argument");
        }
    }
    return (PyObject*)result;
}

static char doc_xmpz[] = "\
xmpz(n):\n\
      builds an 'xmpz' object with a numeric value n (truncating n\n\
      to its integer part if it's a Fraction, 'mpq', Decimal, float\n\
      or 'mpfr')\n\
xmpz(s,base=0):\n\
      builds an 'xmpz' object from a string s made up of digits in the\n\
      given base.  If base=0, binary, octal, or hex Python strings\n\
      are recognized by leading 0b, 0o, or 0x characters, otherwise\n\
      the string is assumed to be decimal.\n\
";
static PyObject *
Pygmpy_xmpz(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyxmpzObject *result = 0;
    PyObject *n = 0;
    long base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"n", "base", NULL };

    TRACE("Pygmpy_xmpz() called...\n");

    /* Optimize the most common use case */
    argc = PyTuple_Size(args);
    if (argc == 1) {
        n = PyTuple_GetItem(args, 0);
#ifdef WITHMPFR
        if (isReal(n) && !keywds) {
#else
        if ((isRational(n) || PyFloat_Check(n)) && !keywds) {
#endif
            result = anynum2Pyxmpz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("xmpz() requires numeric or string argument");
            return (PyObject*)result;
        }
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist,
                                     &n, &base))
        return NULL;

    if ((base!=0) && ((base<2)||(base>62))) {
        VALUE_ERROR("base for xmpz() must be 0 or in the "
                    "interval 2 ... 62");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        /* build-from-string (ascii or unicode) */
        result = PyStr2Pyxmpz(n, base);
    }
    else {
        if (argc==2 || (argc == 1 && keywds))
            TYPE_ERROR("xmpz() with non-string argument needs exactly "
                       "1 argument");
        else {
            result = anynum2Pyxmpz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("xmpz() requires numeric or string argument");
        }
    }
    return (PyObject*)result;
}

static char doc_mpq[] = "\
mpq(n):\n\
      builds an 'mpq' object with a numeric value n. Decimal and\n\
      Fraction values are converted exactly.\n\
mpq(n,m):\n\
      builds an 'mpq' object with a numeric value n/m.\n\
mpq(s,base=10):\n\
      builds an 'mpq' object from a string s made up of digits in the\n\
      given base.  s may be made up of two numbers in the same base\n\
      separated by a '/' character.\n\
";
static PyObject *
Pygmpy_mpq(PyObject *self, PyObject *args, PyObject *keywds)
{
    PympqObject *result = 0, *temp;
    PyObject *n = 0, *m = 0;
    long base = 10;
    Py_ssize_t argc;
    static char *kwlist[] = {"s", "base", NULL };

    TRACE("Pygmpy_mpq() called...\n");

    argc = PyTuple_Size(args);
    if (argc < 1 || argc > 2) {
        TYPE_ERROR("mpq() requires 1 or 2 arguments");
        return NULL;
    }

    n = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(n)) {
        /* keyword base is legal */
        if (PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist, &n, &base)) {
            if ((base!=0) && ((base<2)||(base>62))) {
                VALUE_ERROR("base for mpq() must be 0 or in the "
                            "interval 2 ... 62");
            }
            else {
                result = PyStr2Pympq(n, base);
            }
        }
        return (PyObject*)result;
    }

    if (isDecimal(n)) {
        result = Pympq_From_Decimal(n);
        if (!mpz_cmp_si(mpq_denref(result->q), 0)) {
            if (mpz_get_si(mpq_numref(result->q)) == 0)
                VALUE_ERROR("'mpq' does not support NaN");
            else
                VALUE_ERROR("'mpq' does not support Infinity");
            Py_DECREF((PyObject*)result);
            result = NULL;
        }
        return (PyObject*)result;
    }

    if (argc == 2)
        m = PyTuple_GetItem(args, 1);

#ifdef WITHMPFR
    if (!isReal(n) || (m && !isReal(m))) {
#else
    if (!(isRational(n) || PyFloat_Check(n)) ||
        (m && !(isRational(m) || PyFloat_Check(m)))) {
#endif
        TYPE_ERROR("mpq() requires numeric or string argument");
        return NULL;
    }

    /* should now have one or two numeric values */
    result = Pympq_From_Real(n);
    if (!result && !PyErr_Occurred()) {
        TYPE_ERROR("mpq() requires numeric or string argument");
        return NULL;
    }
    if (m) {
        temp = Pympq_From_Real(m);
        if (!temp && !PyErr_Occurred()) {
            TYPE_ERROR("mpq() requires numeric or string argument");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpq_sgn(temp->q) == 0) {
            ZERO_ERROR("zero denominator in 'mpq'");
            Py_DECREF((PyObject*)result);
            Py_DECREF((PyObject*)temp);
            return NULL;
        }
        mpq_div(result->q, result->q, temp->q);
        Py_DECREF((PyObject*)temp);
    }
    return (PyObject*)result;
}

/* Helper function for hashing for Python < 3.2 */
/* Currently only used for mpq. Should refactor the mpq code and remove. */
#ifndef _PyHASH_MODULUS
static long
dohash(PyObject* tempPynum)
{
    long hash;
    if (!tempPynum) return -1;
    hash = PyObject_Hash(tempPynum);
    Py_DECREF(tempPynum);
    return hash;
}
#endif

/* ARITHMETIC */

#include "gmpy_mpz.c"
#include "gmpy_mpz_divmod2exp.c"
#include "gmpy_mpz_divmod.c"
#include "gmpy_mpz_inplace.c"
#include "gmpy_xmpz_inplace.c"
#include "gmpy_mpq.c"

#ifdef WITHMPFR
#include "gmpy_mpfr.c"
#endif

#ifdef WITHMPC
#include "gmpy_mpc.c"
#endif

#include "gmpy_basic.c"

/* COMPARING */

static PyObject *_cmp_to_object(int c, int op)
{
    PyObject *result;
    switch (op) {
    case Py_LT: c = c <  0; break;
    case Py_LE: c = c <= 0; break;
    case Py_EQ: c = c == 0; break;
    case Py_NE: c = c != 0; break;
    case Py_GT: c = c >  0; break;
    case Py_GE: c = c >= 0; break;
    }
    result = c ? Py_True : Py_False;
    Py_INCREF(result);
    return result;
}
static PyObject *
mpany_richcompare(PyObject *a, PyObject *b, int op)
{
    int c, overflow;
    long temp;
    mpz_t tempz;
    PyObject *tempa = 0, *tempb = 0;
    PyObject *result = 0;

    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            TRACE("compare (mpz,int)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                c = mpz_cmp(Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else {
                c = mpz_cmp_si(Pympz_AS_MPZ(a), temp);
            }
            return _cmp_to_object(c, op);
        }
        if (CHECK_MPZANY(b)) {
            TRACE("compare (mpz,mpz)\n");
            return _cmp_to_object(mpz_cmp(Pympz_AS_MPZ(a), Pympz_AS_MPZ(b)), op);
        }
        if (isInteger(b)) {
            TRACE("compare (mpz,integer)\n");
            tempb = (PyObject*)Pympz_From_Integer(b);
            if (!tempb)
                return NULL;
            c = mpz_cmp(Pympz_AS_MPZ(a), Pympz_AS_MPZ(tempb));
            Py_DECREF(tempb);
            return _cmp_to_object(c, op);
        }
        if (isRational(b)) {
            TRACE("compare (mpz,rational)\n");
            tempa = (PyObject*)Pympq_From_Rational(a);
            tempb = (PyObject*)Pympq_From_Rational(b);
            if (!tempa || !tempb) {
                Py_XDECREF(a);
                Py_XDECREF(b);
                return NULL;
            }
            c = mpq_cmp(Pympq_AS_MPQ(a), Pympq_AS_MPQ(tempb));
            Py_DECREF(tempa);
            Py_DECREF(tempb);
            return _cmp_to_object(c, op);
        }
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            if (Py_IS_NAN(d)) {
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else if (Py_IS_INFINITY(d)) {
                if (d < 0.0)
                    return _cmp_to_object(1, op);
                else
                    return _cmp_to_object(-1, op);
            }
            else {
                return _cmp_to_object(mpz_cmp_d(Pympz_AS_MPZ(a), d), op);
            }
        }
        if (isDecimal(b)) {
            tempa = (PyObject*)Pympq_From_Rational(a);
            tempb = (PyObject*)Pympq_From_Decimal(b);
            if (!tempa || !tempb) {
                Py_XDECREF(a);
                Py_XDECREF(b);
                return NULL;
            }
            if (!mpz_cmp_si(mpq_denref(Pympq_AS_MPQ(tempb)), 0)) {
                if (!mpz_cmp_si(mpq_numref(Pympq_AS_MPQ(tempb)), 0)) {
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_DECREF(tempa);
                    Py_DECREF(tempb);
                    Py_INCREF(result);
                    return result;
                }
                else if (mpz_cmp_si(mpq_numref(Pympq_AS_MPQ(tempb)), 0) < 0) {
                    Py_DECREF(tempa);
                    Py_DECREF(tempb);
                    return _cmp_to_object(1, op);
                }
                else {
                    Py_DECREF(tempa);
                    Py_DECREF(tempb);
                    return _cmp_to_object(-1, op);
                }
            }
            else {
                c = mpq_cmp(Pympq_AS_MPQ(tempa), Pympq_AS_MPQ(tempb));
                Py_DECREF(tempa);
                Py_DECREF(tempb);
                return _cmp_to_object(c, op);
            }
        }
    }
    if (Pympq_Check(a)) {
        if (Pympq_Check(b)) {
            TRACE("compare (mpq,mpq)\n");
            return _cmp_to_object(mpq_cmp(Pympq_AS_MPQ(a), Pympq_AS_MPQ(b)), op);
        }
        if (isRational(b)) {
            TRACE("compare (mpq,rational)\n");
            tempb = (PyObject*)Pympq_From_Rational(b);
            c = mpq_cmp(Pympq_AS_MPQ(a), Pympq_AS_MPQ(tempb));
            Py_DECREF(tempb);
            return _cmp_to_object(c, op);
        }
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            if (Py_IS_NAN(d)) {
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else if (Py_IS_INFINITY(d)) {
                if (d < 0.0)
                    return _cmp_to_object(1, op);
                else
                    return _cmp_to_object(-1, op);
            }
            else {
                tempb = (PyObject*)Pympq_new();
                if (!tempb)
                    return NULL;
                mpq_set_d(Pympq_AS_MPQ(tempb), d);
                c = mpq_cmp(Pympq_AS_MPQ(a), Pympq_AS_MPQ(tempb));
                Py_DECREF(tempb);
                return _cmp_to_object(c, op);
            }
        }
        if (isDecimal(b)) {
            if (!(tempb = (PyObject*)Pympq_From_Decimal(b)))
                return NULL;
            if (!mpz_cmp_si(mpq_denref(Pympq_AS_MPQ(tempb)), 0)) {
                if (!mpz_cmp_si(mpq_numref(Pympq_AS_MPQ(tempb)), 0)) {
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_DECREF(tempb);
                    Py_INCREF(result);
                    return result;
                }
                else if (mpz_cmp_si(mpq_numref(Pympq_AS_MPQ(tempb)), 0) < 0) {
                    Py_DECREF(tempb);
                    return _cmp_to_object(1, op);
                }
                else {
                    Py_DECREF(tempb);
                    return _cmp_to_object(-1, op);
                }
            }
            else {
                c = mpq_cmp(Pympq_AS_MPQ(a), Pympq_AS_MPQ(tempb));
                Py_DECREF(tempb);
                return _cmp_to_object(c, op);
            }
        }
    }
#ifdef WITHMPFR
    if (Pympfr_Check(a)) {
        if (Pympfr_Check(b)) {
            TRACE("compare (mpfr,mpfr)\n");
            mpfr_clear_flags();
            c = mpfr_cmp(Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(b));
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            TRACE("compare (mpfr,float)\n");
            mpfr_clear_flags();
            c = mpfr_cmp_d(Pympfr_AS_MPFR(a), d);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (isInteger(b)) {
            TRACE("compare (mpfr,integer)\n");
            tempb = (PyObject*)Pympz_From_Integer(b);
            if (!tempb)
                return NULL;
            mpfr_clear_flags();
            c = mpfr_cmp_z(Pympfr_AS_MPFR(a), Pympz_AS_MPZ(tempb));
            Py_DECREF(tempb);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (isRational(b)) {
            TRACE("compare (mpfr,rational)\n");
            tempb = (PyObject*)Pympq_From_Rational(b);
            if (!tempb)
                return NULL;
            mpfr_clear_flags();
            c = mpfr_cmp_q(Pympfr_AS_MPFR(a), Pympq_AS_MPQ(tempb));
            Py_DECREF(tempb);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (isDecimal(b)) {
            TRACE("compare (mpfr,decimal)\n");
            tempb = (PyObject*)Pympq_From_Decimal(b);
            if (!tempb)
                return NULL;
            if (!mpz_cmp_si(mpq_denref(Pympq_AS_MPQ(tempb)), 0)) {
                if (!mpz_cmp_si(mpq_numref(Pympq_AS_MPQ(tempb)), 0)) {
                    context->now.erange = 1;
                    if (context->now.trap_erange) {
                        GMPY_ERANGE("comparison with NaN");
                        return NULL;
                    }
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_DECREF(tempb);
                    Py_INCREF(result);
                    return result;
                }
                else if (mpz_cmp_si(mpq_numref(Pympq_AS_MPQ(tempb)), 0) < 0) {
                    Py_DECREF(tempb);
                    return _cmp_to_object(1, op);
                }
                else {
                    Py_DECREF(tempb);
                    return _cmp_to_object(-1, op);
                }
            }
            else {
                mpfr_clear_flags();
                c = mpfr_cmp_q(Pympfr_AS_MPFR(a), Pympq_AS_MPQ(tempb));
                Py_DECREF(tempb);
                if (mpfr_erangeflag_p()) {
                    /* Set erange and check if an exception should be raised. */
                    context->now.erange = 1;
                    if (context->now.trap_erange) {
                        GMPY_ERANGE("comparison with NaN");
                        return NULL;
                    }
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_INCREF(result);
                    return result;
                }
                else {
                    return _cmp_to_object(c, op);
                }
            }
        }
        if (isReal(b)) {
            TRACE("compare (mpfr,real)\n");
            tempb = (PyObject*)Pympfr_From_Real(b, 0);
            if (!tempb)
                return NULL;
            mpfr_clear_flags();
            c = mpfr_cmp(Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(tempb));
            Py_DECREF(tempb);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
    }
#endif

#ifdef WITHMPC
    if (Pympc_Check(a)) {
        if (!(op == Py_EQ || op == Py_NE)) {
            TYPE_ERROR("no ordering relation is defined for complex numbers");
            return NULL;
        }
        if (Pympc_Check(b)) {
            TRACE("compare (mpc,mpc)\n");
            mpfr_clear_flags();
            c = mpc_cmp(Pympc_AS_MPC(a), Pympc_AS_MPC(b));
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (PyComplex_Check(b)) {
            PympcObject *tempmpc;

            if (!(tempmpc = PyComplex2Pympc(b, 53, 53)))
                return NULL;
            mpfr_clear_flags();
            c = mpc_cmp(Pympc_AS_MPC(a), Pympc_AS_MPC(tempmpc));
            Py_DECREF((PyObject*)tempmpc);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        /* a.imag must be 0 or else all further comparisons will be NE */
        if (!mpfr_zero_p(mpc_imagref(Pympc_AS_MPC(a)))) {
            /* if a.real is NaN, possibly raise exception */
            if (mpfr_nan_p(mpc_realref(Pympc_AS_MPC(a)))) {
                context->now.erange = 1;
                if (context->now.trap_erange) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
            }
            result = (op == Py_NE) ? Py_True : Py_False;
            Py_INCREF(result);
            return result;
        }
        else {
            PympfrObject *tempmpfr;

            tempmpfr = Pympfr_new(mpfr_get_prec(mpc_realref(Pympc_AS_MPC(a))));
            if (!tempmpfr)
                return NULL;
            mpc_real(tempmpfr->f, Pympc_AS_MPC(a), context->now.mpfr_round);
            result = mpany_richcompare((PyObject*)tempmpfr, b, op);
            Py_DECREF((PyObject*)tempmpfr);
            return result;
        }
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}


/* Include the module-level methods that call the type-specific methods. */

#include "gmpy_mpany.c"

/* Include helper functions for mpmath. */

#include "gmpy_mpmath.c"

/* method-tables */

#ifdef PY3
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympz_neg,               /* nb_negative             */
    (unaryfunc) Pympz_pos,               /* nb_positive             */
    (unaryfunc) Pympz_abs,               /* nb_absolute             */
    (inquiry) Pympz_nonzero,             /* nb_bool                 */
    (unaryfunc) Pympz_com,               /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
    (unaryfunc) Pympz2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (binaryfunc) Pympz_inplace_add,      /* nb_inplace_add          */
    (binaryfunc) Pympz_inplace_sub,      /* nb_inplace_subtract     */
    (binaryfunc) Pympz_inplace_mul,      /* nb_inplace_multiply     */
    (binaryfunc) Pympz_inplace_rem,      /* nb_inplace_remainder    */
    (ternaryfunc) Pympz_inplace_pow,     /* nb_inplace_power        */
    (binaryfunc) Pympz_inplace_lshift,   /* nb_inplace_lshift       */
    (binaryfunc) Pympz_inplace_rshift,   /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv, /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc)  Pympz_To_Integer,       /* nb_index                */
};

#else
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympz_neg,               /* nb_negative             */
    (unaryfunc) Pympz_pos,               /* nb_positive             */
    (unaryfunc) Pympz_abs,               /* nb_absolute             */
    (inquiry) Pympz_nonzero,             /* nb_bool                 */
    (unaryfunc) Pympz_com,               /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympz_To_Integer,        /* nb_int                  */
    (unaryfunc) Pympz2PyLong,            /* nb_long                 */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (unaryfunc) Pympz_oct,               /* nb_oct                  */
    (unaryfunc) Pympz_hex,               /* nb_hex                  */
    (binaryfunc) Pympz_inplace_add,      /* nb_inplace_add          */
    (binaryfunc) Pympz_inplace_sub,      /* nb_inplace_subtract     */
    (binaryfunc) Pympz_inplace_mul,      /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
    (binaryfunc) Pympz_inplace_rem,      /* nb_inplace_remainder    */
    (ternaryfunc) Pympz_inplace_pow,     /* nb_inplace_power        */
    (binaryfunc) Pympz_inplace_lshift,   /* nb_inplace_lshift       */
    (binaryfunc) Pympz_inplace_rshift,   /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv, /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc) Pympz_To_Integer,        /* nb_index                */
};
#endif

static PyMappingMethods mpz_mapping_methods = {
    (lenfunc)Pympz_nbits,
    (binaryfunc)Pympz_subscript,
    NULL
};

#ifdef PY3
static PyNumberMethods xmpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pyxmpz_neg,              /* nb_negative             */
    (unaryfunc) Pyxmpz_pos,              /* nb_positive             */
    (unaryfunc) Pyxmpz_abs,              /* nb_absolute             */
    (inquiry) Pyxmpz_nonzero,            /* nb_bool                 */
    (unaryfunc) Pyxmpz_com,              /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
    (unaryfunc) Pympz2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (binaryfunc) Pyxmpz_inplace_add,     /* nb_inplace_add          */
    (binaryfunc) Pyxmpz_inplace_sub,     /* nb_inplace_subtract     */
    (binaryfunc) Pyxmpz_inplace_mul,     /* nb_inplace_multiply     */
    (binaryfunc) Pyxmpz_inplace_rem,     /* nb_inplace_remainder    */
    (ternaryfunc) Pyxmpz_inplace_pow,    /* nb_inplace_power        */
    (binaryfunc) Pyxmpz_inplace_lshift,  /* nb_inplace_lshift       */
    (binaryfunc) Pyxmpz_inplace_rshift,  /* nb_inplace_rshift       */
    (binaryfunc) Pyxmpz_inplace_and,     /* nb_inplace_and          */
    (binaryfunc) Pyxmpz_inplace_xor,     /* nb_inplace_xor          */
    (binaryfunc) Pyxmpz_inplace_ior,     /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pyxmpz_inplace_floordiv,/* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc)  Pyxmpz_To_Integer,      /* nb_index                */
};

#else
static PyNumberMethods xmpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pyxmpz_neg,              /* nb_negative             */
    (unaryfunc) Pyxmpz_pos,              /* nb_positive             */
    (unaryfunc) Pyxmpz_abs,              /* nb_absolute             */
    (inquiry) Pyxmpz_nonzero,            /* nb_bool                 */
    (unaryfunc) Pyxmpz_com,              /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympz_To_Integer,        /* nb_int                  */
    (unaryfunc) Pympz2PyLong,            /* nb_long                 */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (unaryfunc) Pyxmpz_oct,              /* nb_oct                  */
    (unaryfunc) Pyxmpz_hex,              /* nb_hex                  */
    (binaryfunc) Pyxmpz_inplace_add,     /* nb_inplace_add          */
    (binaryfunc) Pyxmpz_inplace_sub,     /* nb_inplace_subtract     */
    (binaryfunc) Pyxmpz_inplace_mul,     /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
    (binaryfunc) Pyxmpz_inplace_rem,     /* nb_inplace_remainder    */
    (ternaryfunc) Pyxmpz_inplace_pow,    /* nb_inplace_power        */
    (binaryfunc) Pyxmpz_inplace_lshift,  /* nb_inplace_lshift       */
    (binaryfunc) Pyxmpz_inplace_rshift,  /* nb_inplace_rshift       */
    (binaryfunc) Pyxmpz_inplace_and,     /* nb_inplace_and          */
    (binaryfunc) Pyxmpz_inplace_xor,     /* nb_inplace_xor          */
    (binaryfunc) Pyxmpz_inplace_ior,     /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pyxmpz_inplace_floordiv,/* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc) Pyxmpz_To_Integer,       /* nb_index                */
};
#endif

static PyMappingMethods xmpz_mapping_methods = {
    (lenfunc)Pyxmpz_nbits,
    (binaryfunc)Pyxmpz_subscript,
    (objobjargproc)Pyxmpz_assign_subscript
};

#ifdef PY3
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympq_neg,               /* nb_negative             */
    (unaryfunc) Pympq_pos,               /* nb_positive             */
    (unaryfunc) Pympq_abs,               /* nb_absolute             */
    (inquiry) Pympq_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympq2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympq2PyFloat,           /* nb_float                */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympq_neg,               /* nb_negative             */
    (unaryfunc) Pympq_pos,               /* nb_positive             */
    (unaryfunc) Pympq_abs,               /* nb_absolute             */
    (inquiry) Pympq_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympq2PyInt,             /* nb_int                  */
    (unaryfunc) Pympq2PyLong,            /* nb_long                 */
    (unaryfunc) Pympq2PyFloat,           /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add;         */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympq_getseters[] =
{
    { "numerator", (getter)Pympq_getnumer, NULL, "numerator", NULL },
    { "denominator", (getter)Pympq_getdenom, NULL, "denominator", NULL },
    {NULL}
};

#ifdef WITHMPFR
#ifdef PY3
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,              /* nb_negative             */
    (unaryfunc) Pympfr_pos,              /* nb_positive             */
    (unaryfunc) Pympfr_abs,              /* nb_absolute             */
    (inquiry) Pympfr_nonzero,            /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympfr2PyLong,           /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympfr2PyFloat,          /* nb_float                */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,              /* nb_negative             */
    (unaryfunc) Pympfr_pos,              /* nb_positive             */
    (unaryfunc) Pympfr_abs,              /* nb_absolute             */
    (inquiry) Pympfr_nonzero,            /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympfr2PyInt,            /* nb_int                  */
    (unaryfunc) Pympfr2PyLong,           /* nb_long                 */
    (unaryfunc) Pympfr2PyFloat,          /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympfr_getseters[] =
{
    {"precision", (getter)Pympfr_getprec_attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)Pympfr_getrc_attrib, NULL, "return code", NULL},
    {"imag", (getter)Pympfr_getimag_attrib, NULL, "imaginary component", NULL},
    {"real", (getter)Pympfr_getreal_attrib, NULL, "real component", NULL},
    {NULL}
};
#endif

#ifdef WITHMPC
#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympc_neg,               /* nb_negative             */
    (unaryfunc) Pympc_pos,               /* nb_positive             */
    (unaryfunc) Pympc_abs,               /* nb_absolute             */
    (inquiry) Pympc_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympc2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympc2PyFloat,           /* nb_float                */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympc_neg,               /* nb_negative             */
    (unaryfunc) Pympc_pos,               /* nb_positive             */
    (unaryfunc) Pympc_abs,               /* nb_absolute             */
    (inquiry) Pympc_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympc2PyInt,             /* nb_int                  */
    (unaryfunc) Pympc2PyLong,            /* nb_long                 */
    (unaryfunc) Pympc2PyFloat,           /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympc_getseters[] =
{
    {"precision", (getter)Pympc_getprec_attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)Pympc_getrc_attrib, NULL, "return code", NULL},
    {"imag", (getter)Pympc_getimag_attrib, NULL, "imaginary component", NULL},
    {"real", (getter)Pympc_getreal_attrib, NULL, "real component", NULL},
    {NULL}
};
#endif

static PyMethodDef Pygmpy_methods [] =
{
    { "_cvsid", Pygmpy_get_cvsid, METH_NOARGS, doc_cvsid },
    { "add", Pympany_add, METH_VARARGS, doc_mpany_add },
    { "bit_clear", Pygmpy_bit_clear, METH_VARARGS, doc_bit_clearg },
    { "bit_flip", Pygmpy_bit_flip, METH_VARARGS, doc_bit_flipg },
    { "bit_length", Pympz_bit_length, METH_O, doc_bit_lengthg },
    { "bit_mask", Pympz_bit_mask, METH_O, doc_bit_maskg },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0g },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1g },
    { "bit_set", Pygmpy_bit_set, METH_VARARGS, doc_bit_setg },
    { "bit_test", Pygmpy_bit_test, METH_VARARGS, doc_bit_testg },
    { "binary", Pympany_binary, METH_O, doc_binaryg },
    { "bincoef", Pympz_bincoef, METH_VARARGS, doc_bincoefg },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combg },
    { "c_div", Pygmpy_c_div, METH_VARARGS, doc_gmpy_c_div },
    { "c_div_2exp", Pygmpy_c_div_2exp, METH_VARARGS, doc_gmpy_c_div_2exp },
    { "c_divmod", Pygmpy_c_divmod, METH_VARARGS, doc_gmpy_c_divmod },
    { "c_divmod_2exp", Pygmpy_c_divmod_2exp, METH_VARARGS, doc_gmpy_c_divmod_2exp },
    { "c_mod", Pygmpy_c_mod, METH_VARARGS, doc_gmpy_c_mod },
    { "c_mod_2exp", Pygmpy_c_mod_2exp, METH_VARARGS, doc_gmpy_c_mod_2exp },
    { "denom", Pympq_denom, METH_VARARGS, doc_denomg },
    { "digits", Pympany_digits, METH_VARARGS, doc_g_mpany_digits },
    { "div", Pympany_div, METH_VARARGS, doc_mpany_div },
    { "divexact", Pygmpy_divexact, METH_VARARGS, doc_divexactg },
    { "divm", Pygmpy_divm, METH_VARARGS, doc_divm },
    { "fac", Pygmpy_fac, METH_O, doc_fac },
    { "fib", Pygmpy_fib, METH_O, doc_fib },
    { "fib2", Pygmpy_fib2, METH_O, doc_fib2 },
    { "f_div", Pygmpy_f_div, METH_VARARGS, doc_gmpy_f_div },
    { "f_div_2exp", Pygmpy_f_div_2exp, METH_VARARGS, doc_gmpy_f_div_2exp },
    { "f_divmod", Pygmpy_f_divmod, METH_VARARGS, doc_gmpy_f_divmod },
    { "f_divmod_2exp", Pygmpy_f_divmod_2exp, METH_VARARGS, doc_gmpy_f_divmod_2exp },
    { "f_mod", Pygmpy_f_mod, METH_VARARGS, doc_gmpy_f_mod },
    { "f_mod_2exp", Pygmpy_f_mod_2exp, METH_VARARGS, doc_gmpy_f_mod_2exp },
    { "gcd", Pygmpy_gcd, METH_VARARGS, doc_gcd },
    { "gcdext", Pygmpy_gcdext, METH_VARARGS, doc_gcdext },
    { "get_cache", Pygmpy_get_cache, METH_NOARGS, doc_get_cache },
    { "hamdist", Pympz_hamdist, METH_VARARGS, doc_hamdistg },
    { "invert", Pygmpy_invert, METH_VARARGS, doc_invertg },
    { "isqrt", Pympz_isqrt, METH_O, doc_mpz_isqrt },
    { "isqrt_rem", Pympz_isqrt_rem, METH_VARARGS, doc_mpz_isqrt_rem },
    { "is_even", Pympz_is_even, METH_O, doc_is_eveng },
    { "is_odd", Pympz_is_odd, METH_O, doc_is_oddg },
    { "is_power", Pympz_is_power, METH_O, doc_is_powerg },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primeg },
    { "is_square", Pympz_is_square, METH_O, doc_is_squareg },
    { "jacobi", Pympz_jacobi, METH_VARARGS, doc_jacobig },
    { "kronecker", Pympz_kronecker, METH_VARARGS, doc_kroneckerg },
    { "lcm", Pygmpy_lcm, METH_VARARGS, doc_lcm },
    { "legendre", Pympz_legendre, METH_VARARGS, doc_legendreg },
    { "license", Pygmpy_get_license, METH_NOARGS, doc_license },
    { "lucas", Pygmpy_lucas, METH_O, doc_lucas },
    { "lucas2", Pygmpy_lucas2, METH_O, doc_lucas2 },
    { "mp_version", Pygmpy_get_mp_version, METH_NOARGS, doc_mp_version },
    { "mp_limbsize", Pygmpy_get_mp_limbsize, METH_NOARGS, doc_mp_limbsize },
    { "mpc_version", Pygmpy_get_mpc_version, METH_NOARGS, doc_mpc_version },
    { "mpfr_version", Pygmpy_get_mpfr_version, METH_NOARGS, doc_mpfr_version },
    { "mpq", (PyCFunction)Pygmpy_mpq, METH_VARARGS | METH_KEYWORDS, doc_mpq },
    { "mpq_from_old_binary", Pympq_From_Old_Binary, METH_O, doc_g_mpq_from_old_binary },
    { "mpz", (PyCFunction)Pygmpy_mpz, METH_VARARGS | METH_KEYWORDS, doc_mpz },
    { "mpz_from_old_binary", Pympz_From_Old_Binary, METH_O, doc_g_mpz_from_old_binary },
    { "mul", Pympany_mul, METH_VARARGS, doc_mpany_mul },
    { "next_prime", Pympz_next_prime, METH_O, doc_next_primeg },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsg },
    { "numer", Pympq_numer, METH_VARARGS, doc_numerg },
    { "pack", Pygmpy_pack, METH_VARARGS, doc_gmpy_pack },
    { "popcount", Pympz_popcount, METH_O, doc_popcountg },
    { "printf", Pympany_printf, METH_VARARGS, doc_printf },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivg },
    { "remove", Pympz_remove, METH_VARARGS, doc_removeg },
    { "iroot", Pympz_iroot, METH_VARARGS, doc_mpz_iroot },
    { "iroot_rem", Pympz_iroot_rem, METH_VARARGS, doc_mpz_iroot_rem },
    { "set_cache", Pygmpy_set_cache, METH_VARARGS, doc_set_cache },
    { "set_debug", Pygmpy_set_debug, METH_VARARGS, doc_set_debug },
    { "sign", Pympany_sign, METH_O, doc_g_mpany_sign },
    { "square", Pympany_square, METH_O, doc_mpany_square },
    { "sub", Pympany_sub, METH_VARARGS, doc_mpany_sub },
    { "t_div", Pygmpy_t_div, METH_VARARGS, doc_gmpy_t_div },
    { "t_div_2exp", Pygmpy_t_div_2exp, METH_VARARGS, doc_gmpy_t_div_2exp },
    { "t_divmod", Pygmpy_t_divmod, METH_VARARGS, doc_gmpy_t_divmod },
    { "t_divmod_2exp", Pygmpy_t_divmod_2exp, METH_VARARGS, doc_gmpy_t_divmod_2exp },
    { "t_mod", Pygmpy_t_mod, METH_VARARGS, doc_gmpy_t_mod },
    { "t_mod_2exp", Pygmpy_t_mod_2exp, METH_VARARGS, doc_gmpy_t_mod_2exp },
    { "unpack", Pygmpy_unpack, METH_VARARGS, doc_gmpy_unpack },
    { "version", Pygmpy_get_version, METH_NOARGS, doc_version },
    { "xbit_mask", Pyxmpz_xbit_mask, METH_O, doc_xbit_maskg },
    { "xmpz", (PyCFunction)Pygmpy_xmpz, METH_VARARGS | METH_KEYWORDS, doc_xmpz },
    { "_mpmath_normalize", Pympz_mpmath_normalize, METH_VARARGS, doc_mpmath_normalizeg },
    { "_mpmath_create", Pympz_mpmath_create, METH_VARARGS, doc_mpmath_createg },
#ifdef WITHMPFR
    { "acos", Pympany_acos, METH_O, doc_mpany_acos },
    { "acosh", Pympany_acosh, METH_O, doc_mpany_acosh },
    { "ai", Pympfr_ai, METH_O, doc_g_mpfr_ai },
    { "agm", Pympfr_agm, METH_VARARGS, doc_g_mpfr_agm },
    { "asin", Pympany_asin, METH_O, doc_mpany_asin },
    { "asinh", Pympany_asinh, METH_O, doc_mpany_asinh },
    { "atan", Pympany_atan, METH_O, doc_mpany_atan },
    { "atanh", Pympany_atanh, METH_O, doc_mpany_atanh },
    { "atan2", Pympfr_atan2, METH_VARARGS, doc_g_mpfr_atan2 },
    { "cbrt", Pympfr_cbrt, METH_O, doc_g_mpfr_cbrt },
    { "ceil", Pympfr_ceil, METH_O, doc_g_mpfr_ceil },
    { "check_range", Pympfr_check_range, METH_O, doc_g_mpfr_check_range },
    { "const_catalan", Pympfr_const_catalan, METH_NOARGS, doc_mpfr_const_catalan },
    { "const_euler", Pympfr_const_euler, METH_NOARGS, doc_mpfr_const_euler },
    { "const_log2", Pympfr_const_log2, METH_NOARGS, doc_mpfr_const_log2 },
    { "const_pi", Pympfr_const_pi, METH_VARARGS, doc_mpfr_const_pi },
    { "context", (PyCFunction)Pygmpy_context, METH_VARARGS | METH_KEYWORDS, doc_context },
    { "copy_sign", Pympfr_copy_sign, METH_VARARGS, doc_g_mpfr_copy_sign },
    { "cos", Pympany_cos, METH_O, doc_mpany_cos },
    { "cosh", Pympany_cosh, METH_O, doc_mpany_cosh },
    { "cot", Pympfr_cot, METH_O, doc_g_mpfr_cot },
    { "coth", Pympfr_coth, METH_O, doc_g_mpfr_coth },
    { "csc", Pympfr_csc, METH_O, doc_g_mpfr_csc },
    { "csch", Pympfr_csch, METH_O, doc_g_mpfr_csch },
    { "digamma", Pympfr_digamma, METH_O, doc_g_mpfr_digamma },
    { "div_2exp", Pympany_div_2exp, METH_VARARGS, doc_mpany_div_2exp },
    { "eint", Pympfr_eint, METH_O, doc_g_mpfr_eint },
    { "erf", Pympfr_erf, METH_O, doc_g_mpfr_erf },
    { "erfc", Pympfr_erfc, METH_O, doc_g_mpfr_erfc },
    { "exp", Pympany_exp, METH_O, doc_mpany_exp },
    { "expm1", Pympfr_expm1, METH_O, doc_g_mpfr_expm1 },
    { "exp10", Pympfr_exp10, METH_O, doc_g_mpfr_exp10 },
    { "exp2", Pympfr_exp2, METH_O, doc_g_mpfr_exp2 },
    { "f2q", Pympfr_f2q, METH_VARARGS, doc_g_mpfr_f2q },
    { "factorial", Pympfr_factorial, METH_O, doc_g_mpfr_factorial },
    { "floor", Pympfr_floor, METH_O, doc_g_mpfr_floor},
    { "fma", Pympany_fma, METH_VARARGS, doc_mpany_fma },
    { "fms", Pympany_fms, METH_VARARGS, doc_mpany_fms },
    { "fmod", Pympfr_fmod, METH_VARARGS, doc_g_mpfr_fmod },
    { "frac", Pympfr_frac, METH_O, doc_g_mpfr_frac },
    { "frexp", Pympfr_frexp, METH_O, doc_g_mpfr_frexp },
    { "gamma", Pympfr_gamma, METH_O, doc_g_mpfr_gamma },
    { "get_emax_max", Pympfr_get_emax_max, METH_NOARGS, doc_g_mpfr_get_emax_max },
    { "get_emin_min", Pympfr_get_emin_min, METH_NOARGS, doc_g_mpfr_get_emin_min },
    { "get_exp", Pympfr_get_exp, METH_O, doc_g_mpfr_get_exp },
    { "get_max_precision", Pympfr_get_max_precision, METH_NOARGS, doc_g_mpfr_get_max_precision },
    { "hypot", Pympfr_hypot, METH_VARARGS, doc_g_mpfr_hypot },
    { "inf", Pympfr_set_inf, METH_VARARGS, doc_g_mpfr_set_inf },
    { "is_inf", Pympany_is_inf, METH_O, doc_mpany_is_inf },
    { "is_integer", Pympfr_is_integer, METH_O, doc_g_mpfr_is_integer },
    { "is_lessgreater", Pympfr_is_lessgreater, METH_VARARGS, doc_g_mpfr_is_lessgreater },
    { "is_nan", Pympany_is_nan, METH_O, doc_mpany_is_nan },
    { "is_number", Pympfr_is_number, METH_O, doc_g_mpfr_is_number },
    { "is_regular", Pympfr_is_regular, METH_O, doc_g_mpfr_is_regular },
    { "is_signed", Pympfr_is_signed, METH_O, doc_g_mpfr_is_signed },
    { "is_unordered", Pympfr_is_unordered, METH_VARARGS, doc_g_mpfr_is_unordered },
    { "is_zero", Pympany_is_zero, METH_O, doc_mpany_is_zero },
    { "jn", Pympfr_jn, METH_VARARGS, doc_g_mpfr_jn },
    { "j0", Pympfr_j0, METH_O, doc_g_mpfr_j0 },
    { "j1", Pympfr_j1, METH_O, doc_g_mpfr_j1 },
    { "lgamma", Pympfr_lgamma, METH_O, doc_g_mpfr_lgamma },
    { "li2", Pympfr_li2, METH_O, doc_g_mpfr_li2 },
    { "lngamma", Pympfr_lngamma, METH_O, doc_g_mpfr_lngamma },
    { "log", Pympany_log, METH_O, doc_mpany_log },
    { "log1p", Pympfr_log1p, METH_O, doc_g_mpfr_log1p },
    { "log10", Pympfr_log10, METH_O, doc_g_mpfr_log10 },
    { "log2", Pympfr_log2, METH_O, doc_g_mpfr_log2 },
    { "max", Pympfr_max, METH_VARARGS, doc_g_mpfr_max },
    { "min", Pympfr_min, METH_VARARGS, doc_g_mpfr_min },
    { "modf", Pympfr_modf, METH_O, doc_g_mpfr_modf },
    { "mpfr", (PyCFunction)Pygmpy_mpfr, METH_VARARGS | METH_KEYWORDS, doc_mpfr },
    { "mpfr_from_old_binary", Pympfr_From_Old_Binary, METH_O, doc_g_mpfr_from_old_binary },
    { "mul_2exp", Pympany_mul_2exp, METH_VARARGS, doc_mpany_mul_2exp },
    { "nan", Pympfr_set_nan, METH_NOARGS, doc_g_mpfr_set_nan },
    { "new_context", (PyCFunction)Pygmpy_new_context, METH_VARARGS | METH_KEYWORDS, doc_new_context },
    { "next_above", Pympfr_nextabove, METH_O, doc_g_mpfr_nextabove },
    { "next_below", Pympfr_nextbelow, METH_O, doc_g_mpfr_nextbelow },
    { "next_toward", Pympfr_nexttoward, METH_VARARGS, doc_g_mpfr_nexttoward },
    { "rec_sqrt", Pympfr_rec_sqrt, METH_O, doc_g_mpfr_rec_sqrt },
    { "reldiff", Pympfr_reldiff, METH_VARARGS, doc_g_mpfr_reldiff },
    { "remainder", Pympfr_remainder, METH_VARARGS, doc_g_mpfr_remainder },
    { "remquo", Pympfr_remquo, METH_VARARGS, doc_g_mpfr_remquo },
    { "rint", Pympfr_rint, METH_O, doc_g_mpfr_rint },
    { "rint_ceil", Pympfr_rint_ceil, METH_O, doc_g_mpfr_rint_ceil },
    { "rint_floor", Pympfr_rint_floor, METH_O, doc_g_mpfr_rint_floor },
    { "rint_round", Pympfr_rint_round, METH_O, doc_g_mpfr_rint_round },
    { "rint_trunc", Pympfr_rint_trunc, METH_O, doc_g_mpfr_rint_trunc },
    { "root", Pympfr_root, METH_VARARGS, doc_mpfr_root },
    { "round", Pympfr_round, METH_VARARGS, doc_g_mpfr_round },
    { "round2", Pympfr_round2, METH_O, doc_g_mpfr_round2 },
    { "sec", Pympfr_sec, METH_O, doc_g_mpfr_sec },
    { "sech", Pympfr_sech, METH_O, doc_g_mpfr_sech },
    { "set_context", Pygmpy_set_context, METH_O, doc_set_context },
    { "set_exp", Pympfr_set_exp, METH_VARARGS, doc_g_mpfr_set_exp },
    { "set_sign", Pympfr_set_sign, METH_VARARGS, doc_g_mpfr_set_sign },
    { "sin", Pympany_sin, METH_O, doc_mpany_sin },
    { "sinh", Pympany_sinh, METH_O, doc_mpany_sinh },
    { "sinh_cosh", Pympfr_sinh_cosh, METH_O, doc_g_mpfr_sinh_cosh },
    { "sin_cos", Pympany_sin_cos, METH_O, doc_mpany_sin_cos },
    { "sqrt", Pympany_sqrt, METH_O, doc_mpany_sqrt },
    { "tan", Pympany_tan, METH_O, doc_mpany_tan },
    { "tanh", Pympany_tanh, METH_O, doc_mpany_tanh },
    { "trunc", Pympfr_trunc, METH_O, doc_g_mpfr_trunc },
    { "yn", Pympfr_yn, METH_VARARGS, doc_g_mpfr_yn },
    { "y0", Pympfr_y0, METH_O, doc_g_mpfr_y0 },
    { "y1", Pympfr_y1, METH_O, doc_g_mpfr_y1 },
    { "zero", Pympfr_set_zero, METH_VARARGS, doc_g_mpfr_set_zero },
    { "zeta", Pympfr_zeta, METH_O, doc_g_mpfr_zeta },
#endif

#ifdef WITHMPC
    { "mpc", (PyCFunction)Pygmpy_mpc, METH_VARARGS | METH_KEYWORDS, doc_g_mpc },
    { "norm", Pympc_norm, METH_O, doc_mpc_norm },
    { "polar", Pympc_polar, METH_O, doc_mpc_polar },
    { "phase", Pympc_phase, METH_O, doc_mpc_phase },
    { "proj", Pympc_proj, METH_O, doc_mpc_proj },
    { "rect", Pympc_rect, METH_VARARGS, doc_mpc_rect },
#endif
    { NULL, NULL, 1}
};

static PyMethodDef Pympz_methods [] =
{
    { "__format__", Pympz_format, METH_VARARGS, doc_mpz_format },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "digits", Pympz_digits, METH_VARARGS, doc_mpz_digits },
    { "is_even", Pympz_is_even, METH_NOARGS, doc_is_evenm },
    { "is_odd", Pympz_is_odd, METH_NOARGS, doc_is_oddm },
    { "is_square", Pympz_is_square, METH_NOARGS, doc_is_squarem },
    { "is_power", Pympz_is_power, METH_NOARGS, doc_is_powerm },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primem },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pyxmpz_methods [] =
{
    { "__format__", Pympz_format, METH_VARARGS, doc_mpz_format },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "copy", Pyxmpz_copy, METH_NOARGS, doc_xmpz_copy },
    { "digits", Pyxmpz_digits, METH_VARARGS, doc_mpz_digits },
    { "is_even", Pympz_is_even, METH_NOARGS, doc_is_evenm },
    { "is_odd", Pympz_is_odd, METH_NOARGS, doc_is_oddm },
    { "is_square", Pympz_is_square, METH_VARARGS, doc_is_squarem },
    { "is_power", Pympz_is_power, METH_VARARGS, doc_is_powerm },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primem },
    { "make_mpz", Pyxmpz_make_mpz, METH_NOARGS, doc_make_mpzm },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pympq_methods [] =
{
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "digits", Pympq_digits, METH_VARARGS, doc_qdigitsm },
    { NULL, NULL, 1 }
};

#ifdef WITHMPFR
static PyMethodDef Pympfr_methods [] =
{
    { "__ceil__", Pympfr_ceil, METH_NOARGS, doc_mpfr_ceil },
    { "__floor__", Pympfr_floor, METH_NOARGS, doc_mpfr_floor },
    { "__format__", Pympfr_format, METH_VARARGS, doc_mpfr_format },
    { "__trunc__", Pympfr_trunc, METH_NOARGS, doc_mpfr_trunc },
    { "as_integer_ratio", Pympfr_integer_ratio, METH_NOARGS, doc_mpfr_integer_ratio },
    { "as_mantissa_exp", Pympfr_mantissa_exp, METH_NOARGS, doc_mpfr_mantissa_exp },
    { "as_simple_fraction", (PyCFunction)Pympfr_simple_fraction, METH_VARARGS | METH_KEYWORDS, doc_mpfr_simple_fraction },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "conjugate", Pympfr_conjugate, METH_NOARGS, doc_mpfr_conjugate },
    { "digits", Pympfr_digits, METH_VARARGS, doc_mpfr_digits },
    { "is_inf", Pympfr_is_inf, METH_NOARGS, doc_mpfr_is_inf },
    { "is_integer", Pympfr_is_integer, METH_NOARGS, doc_mpfr_is_integer },
    { "is_lessgreater", Pympfr_is_lessgreater, METH_VARARGS, doc_mpfr_is_lessgreater },
    { "is_nan", Pympfr_is_nan, METH_NOARGS, doc_mpfr_is_nan },
    { "is_number", Pympfr_is_number, METH_NOARGS, doc_mpfr_is_number },
    { "is_regular", Pympfr_is_regular, METH_NOARGS, doc_mpfr_is_regular },
    { "is_signed", Pympfr_is_signed, METH_NOARGS, doc_mpfr_is_signed },
    { "is_unordered", Pympfr_is_unordered, METH_VARARGS, doc_mpfr_is_unordered },
    { "is_zero", Pympfr_is_zero, METH_NOARGS, doc_mpfr_is_zero },
    { NULL, NULL, 1 }
};
#endif

#ifdef WITHMPC
static PyMethodDef Pympc_methods[] =
{
    { "__format__", Pympc_format, METH_VARARGS, doc_mpc_format },
    { "conjugate", Pympc_conjugate, METH_NOARGS, doc_mpc_conjugate },
    { "digits", Pympc_digits, METH_VARARGS, doc_mpc_digits },
    { "is_inf", Pympc_is_INF, METH_NOARGS, doc_mpc_is_inf },
    { "is_nan", Pympc_is_NAN, METH_NOARGS, doc_mpc_is_nan },
    { "is_zero", Pympc_is_ZERO, METH_NOARGS, doc_mpc_is_zero },
    { NULL, NULL, 1 }
};
#endif

static PyTypeObject Pympz_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpz",                                  /* tp_name          */
    sizeof(PympzObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympz_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympz2repr,                  /* tp_repr          */
    &mpz_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &mpz_mapping_methods,                   /* tp_as_mapping    */
    (hashfunc) Pympz_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympz2str,                   /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_INDEX|Py_TPFLAGS_HAVE_RICHCOMPARE| \
    Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_HAVE_CLASS| \
    Py_TPFLAGS_HAVE_INPLACEOPS,
#endif
    "GNU Multi Precision signed integer",   /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympz_methods,                          /* tp_methods       */
};

static PyTypeObject Pyxmpz_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "xmpz",                                 /* tp_name          */
    sizeof(PyxmpzObject),                   /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pyxmpz_dealloc,            /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pyxmpz2repr,                 /* tp_repr          */
    &xmpz_number_methods,                   /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &xmpz_mapping_methods,                  /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pyxmpz2str,                  /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_INDEX|Py_TPFLAGS_HAVE_RICHCOMPARE| \
    Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_HAVE_CLASS| \
    Py_TPFLAGS_HAVE_INPLACEOPS,
#endif
    "GNU Multi Precision signed integer",   /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pyxmpz_methods,                         /* tp_methods       */
};

static PyTypeObject Pympq_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpq",                                  /* tp_name          */
    sizeof(PympqObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympq_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympq2repr,                  /* tp_repr          */
    &mpq_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympq_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympq2str,                   /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE |
        Py_TPFLAGS_CHECKTYPES,              /* tp_flags         */
#endif
    "GNU Multi Precision rational number",  /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympq_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympq_getseters,                        /* tp_getset        */
};

#ifdef WITHMPFR
static PyTypeObject Pympfr_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpfr",                                 /* tp_name          */
    sizeof(PympfrObject),                   /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympfr_dealloc,            /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympfr2repr,                 /* tp_repr          */
    &mpfr_number_methods,                   /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympfr_hash,                 /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympfr2str,                  /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "GNU Multi Precision floating point",   /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympfr_methods,                         /* tp_methods       */
        0,                                  /* tp_members       */
    Pympfr_getseters,                       /* tp_getset        */
};
#endif

#ifdef WITHMPC
static PyTypeObject Pympc_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpc",                                  /* tp_name          */
    sizeof(PympcObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympc_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympc2repr,                  /* tp_repr          */
    &mpc_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympc_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympc2str,                   /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "MPC-based complex number",             /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympc_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympc_getseters,                        /* tp_getset        */
};
#endif

#ifdef USE_PYMEM
static void *
gmpy_allocate(size_t size)
{
    void *res;
    size_t usize = size;

    if (usize < GMPY_ALLOC_MIN)
        usize = GMPY_ALLOC_MIN;

    if (!(res = PyObject_Malloc(usize)))
        Py_FatalError("mp_allocate failure");

    return res;
}

static void *
gmpy_reallocate(void *ptr, size_t old_size, size_t new_size)
{
    void *res;
    size_t uold = old_size;
    size_t unew = new_size;

    if (uold < GMPY_ALLOC_MIN)
        uold = GMPY_ALLOC_MIN;
    if (unew < GMPY_ALLOC_MIN)
        unew = GMPY_ALLOC_MIN;

    if (uold == unew)
        return ptr;

    if (!(res = PyObject_Realloc(ptr, unew)))
        Py_FatalError("mp_reallocate failure");

    return res;
}

static void
gmpy_free( void *ptr, size_t size)
{
    size_t usize=size;

    if (usize < GMPY_ALLOC_MIN)
        usize = GMPY_ALLOC_MIN;

    PyObject_Free(ptr);
}
#endif /* USE_PYMEM */

static void
_PyInitGMP(void)
{
#ifdef WITHMPFR
    PyObject *temp = NULL;
#endif

#ifdef USE_PYMEM
    mp_set_memory_functions(gmpy_allocate, gmpy_reallocate, gmpy_free);
#endif
    set_zcache();
    set_pympzcache();
    set_pympqcache();
    set_pyxmpzcache();
#ifdef WITHMPFR
    set_pympfrcache();
    context = GMPyContext_new();
    GMPyExc_GmpyError = PyErr_NewException("gmpy2.GmpyError",
                                           PyExc_ArithmeticError, NULL);
    GMPyExc_Erange = PyErr_NewException("gmpy2.RangeError",
                                        GMPyExc_GmpyError, NULL);
    GMPyExc_Invalid = PyErr_NewException("gmpy2.InvalidOperationError",
                                         GMPyExc_GmpyError, NULL);
    GMPyExc_Inexact = PyErr_NewException("gmpy2.InexactError",
                                         GMPyExc_GmpyError, NULL);
    GMPyExc_Overflow = PyErr_NewException("gmpy2.OverflowError",
                                          GMPyExc_Inexact, NULL);
    GMPyExc_Underflow = PyErr_NewException("gmpy2.UnderflowError",
                                           GMPyExc_Inexact, NULL);
    GMPyExc_ExpBound = PyErr_NewException("gmpy2.ExponentOutOfBoundsError",
                                          GMPyExc_ExpBound, NULL);

    temp = PyTuple_Pack(2, GMPyExc_GmpyError, PyExc_ZeroDivisionError);
    GMPyExc_DivZero = PyErr_NewException("gmpy2.DivisionByZeroError",
                                         temp, NULL);
    Py_XDECREF(temp);

#endif
}

static char _gmpy_docs[] = "\
gmpy2 2.0.0a3 - General Multiple-precision arithmetic for Python\n\
\n\
Exposes functionality from the GMP or MPIR library to Python 2.6\n\
and later. If available, the MPFR and MPC libraries are used to\n\
support multiple-precision floating-point and complex numbers.\n\
\n\
Allows creation of multiple-precision integer (mpz), mutable\n\
integers (xmpz), rational (mpq), floating-point (mpfr), and complex\n\
(mpc) numbers. Supported operations include conversion between them\n\
and to/from Python numbers/strings, arithmetic, bitwise, and higher-\n\
level mathematical operations.\n\
\n\
'mpz' has comparable functionality to Python's builtin longs, but\n\
can be faster for some operations (particularly multiplication\n\
and raising-to-power) and has many further useful and speedy\n\
functions (prime testing and generation, factorial, fibonacci,\n\
binary-coefficients, gcd, lcm, square and other roots, ...).\n\
\n\
'mpfr' and 'mpc' provide multiple-precision real and complex numbers\n\
with user-definable precision, rounding, and exponent range. All the\n\
advanced functions from the MPFR and MPC libraries are available.\n\
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
PyMODINIT_FUNC PyInit_gmpy2(void)
#else
#define INITERROR return
PyMODINIT_FUNC initgmpy2(void)
#endif
{
    PyObject* copy_reg_module = NULL;

    /* Validate the sizes of the various typedef'ed integer types. */
    if (sizeof(mp_limb_t) != sizeof(size_t)) {
        SYSTEM_ERROR("Size of mp_limb_t and size_t not compatible");
        INITERROR;
    }
    if (sizeof(mp_bitcnt_t) != sizeof(size_t)) {
        SYSTEM_ERROR("Size of mp_bitcnt_t and size_t not compatible");
        INITERROR;
    }
    if (sizeof(mp_size_t) != sizeof(long)) {
        SYSTEM_ERROR("Size of mp_size_t and long not compatible");
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

    if (PyType_Ready(&Pympz_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pympq_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pyxmpz_Type) < 0)
        INITERROR;
#ifdef WITHMPFR
    if (PyType_Ready(&Pympfr_Type) < 0)
        INITERROR;
    if (PyType_Ready(&GMPyContext_Type) < 0)
        INITERROR;
#endif
#ifdef WITHMPC
    if (PyType_Ready(&Pympc_Type) < 0)
        INITERROR;
#endif

    _PyInitGMP();

#ifdef PY3
    gmpy_module = PyModule_Create(&moduledef);
#else
    gmpy_module = Py_InitModule3("gmpy2", Pygmpy_methods, _gmpy_docs);
#endif

    if (gmpy_module == NULL)
        INITERROR;

    /* Add the constants for defining rounding modes. */

#ifdef WITHMPFR
    PyModule_AddIntConstant(gmpy_module, "RoundToNearest", MPFR_RNDN);
    PyModule_AddIntConstant(gmpy_module, "RoundToZero", MPFR_RNDZ);
    PyModule_AddIntConstant(gmpy_module, "RoundUp", MPFR_RNDU);
    PyModule_AddIntConstant(gmpy_module, "RoundDown", MPFR_RNDD);
    PyModule_AddIntConstant(gmpy_module, "RoundAwayZero", MPFR_RNDA);
    PyModule_AddIntConstant(gmpy_module, "Default", GMPY_DEFAULT);
    Py_INCREF(GMPyExc_DivZero);
    PyModule_AddObject(gmpy_module, "DivisionByZeroError", GMPyExc_DivZero);
    Py_INCREF(GMPyExc_Inexact);
    PyModule_AddObject(gmpy_module, "InexactError", GMPyExc_Inexact);
    Py_INCREF(GMPyExc_Invalid);
    PyModule_AddObject(gmpy_module, "InvalidOperationError", GMPyExc_Invalid);
    Py_INCREF(GMPyExc_Overflow);
    PyModule_AddObject(gmpy_module, "OverflowError", GMPyExc_Overflow);
    Py_INCREF(GMPyExc_Underflow);
    PyModule_AddObject(gmpy_module, "UnderflowError", GMPyExc_Underflow);
    Py_INCREF(GMPyExc_Erange);
    PyModule_AddObject(gmpy_module, "RangeError", GMPyExc_Erange);
    Py_INCREF(GMPyExc_ExpBound);
    PyModule_AddObject(gmpy_module, "ExponentOutOfBoundsError", GMPyExc_Erange);
#endif

    /* Add support for pickling. */
#ifdef PY3
    copy_reg_module = PyImport_ImportModule("copyreg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def mpz_reducer(an_mpz): return (gmpy2.mpz_from_old_binary, (an_mpz.binary(),))\n"
            "def mpq_reducer(an_mpq): return (gmpy2.mpq_from_old_binary, (an_mpq.binary(),))\n"
            "copyreg.pickle(type(gmpy2.mpz(0)), mpz_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpq(0)), mpq_reducer)\n"
#ifdef WITHMPFR
            "def mpfr_reducer(an_mpfr): return (gmpy2.mpfr_from_old_binary, (an_mpfr.binary(),))\n"
            "copyreg.pickle(type(gmpy2.mpfr(0)), mpfr_reducer)\n"
#endif
        ;
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;
        PyDict_SetItemString(namespace, "copyreg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (result) {
        }
        else {
            PyErr_Clear();
        }
        Py_DECREF(namespace);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
    }
#else
    copy_reg_module = PyImport_ImportModule("copy_reg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def mpz_reducer(an_mpz): return (gmpy2.mpz_from_old_binary, (an_mpz.binary(),))\n"
            "def mpq_reducer(an_mpq): return (gmpy2.mpq_from_old_binary, (an_mpq.binary(),))\n"
            "copy_reg.pickle(type(gmpy2.mpz(0)), mpz_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpq(0)), mpq_reducer)\n"
#ifdef WITHMPFR
            "def mpfr_reducer(an_mpfr): return (gmpy2.mpfr_from_old_binary, (an_mpfr.binary(),))\n"
            "copy_reg.pickle(type(gmpy2.mpfr(0)), mpfr_reducer)\n"
#endif
        ;
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
