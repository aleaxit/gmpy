/* gmpy2.c
 *
 * Python interface to the GMP or MPIR multiple precision library,
 * Copyright (C) 2000 - 2009 Alex Martelli
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 *************************************************************************
 *
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

char gmpy_version[] = "2.0.0a1+";

char _gmpy_cvs[] = "$Id$";

/*
 * global data declarations
 */

static PyObject *gmpy_module = NULL;

static struct gmpy_global {
    int debug;               /* != 0 if debug messages desired on stderr */
    int cache_size;          /* size of cache, for all caches */
    int cache_obsize;        /* maximum size of the objects that are cached */
    int mpc_rc;              /* result code from MPC */
} global = {
    0,                       /* debug */
    100,                     /* cache_size */
    128,                     /* cache_obsize */
    0                        /* mpc_rc */
};

/* The context manager should really be thread-specific. Until then, gmpy2
 * is NOT thread-safe for mpfr and mpc calculations.
 */

static GMPyContextObject *context = NULL;

/* Define gmpy2 specific errors for mpfr and mpc data types. No change will
 * be made the exceptions raised by mpz, xmpz, and mpq.
 */

static PyObject *GMPyExc_DivZero = NULL;
static PyObject *GMPyExc_Inexact = NULL;
static PyObject *GMPyExc_Invalid = NULL;
static PyObject *GMPyExc_Overflow = NULL;
static PyObject *GMPyExc_Underflow = NULL;
static PyObject *GMPyExc_Erange = NULL;

/* forward declarations of type-objects and method-arrays for them */
#ifdef _MSC_VER
PyMethodDef Pympz_methods [];
PyMethodDef Pympq_methods [];
PyMethodDef Pympfr_methods [];
PyMethodDef Pympc_methods [];
PyMethodDef Pyxmpz_methods [];
#else
static PyMethodDef Pympz_methods [];
static PyMethodDef Pympq_methods [];
static PyMethodDef Pympfr_methods [];
static PyMethodDef Pympc_methods [];
static PyMethodDef Pyxmpz_methods [];
#endif

/* Include fast mpz to/from PyLong conversion from sage. */
#include "mpz_pylong.c"

/* The code for object creation, deletion, and caching is in gmpy_cache.c. */

#include "gmpy_cache.c"

/* Miscellaneous helper functions and simple methods are in gmpy_misc.c. */

#include "gmpy_misc.c"

/* Support for context manager */

#include "gmpy_context.c"

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

/*
  Make a copy of an mpfr object. If bits is 0, the new object will have
  the same precision as the original object. If the requested precision
  is less than the precision of the original object, the new object
  will be rounded to requested precision using the current rounding mode.
*/

static PympfrObject *
Pympfr2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob;

    assert(Pympfr_Check(self));
    if (bits == 0)
        bits = mpfr_get_prec(Pympfr_AS_MPFR(self));
    if ((newob = Pympfr_new(bits)))
        newob->rc = mpfr_set(newob->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    return newob;
}

static PympcObject *
Pympc2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *newob;

    assert(Pympc_Check(self));
    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, Pympc_AS_MPC(self));
    if ((newob = Pympc_new(rprec, iprec)))
        global.mpc_rc = mpc_set(newob->c, Pympc_AS_MPC(self), context->now.mpc_round);
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

static PympfrObject *
PyInt2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympfr_new(bits)))
        newob->rc = mpfr_set_si(newob->f, PyInt_AsLong(self), context->now.mpfr_round);
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
            VALUE_ERROR("gmpy2.mpz() does not handle nan");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpz() does not handle infinity");
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
            VALUE_ERROR("gmpy2.xmpz() does not handle nan");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.xmpz() does not handle infinity");
            return NULL;
        }
        mpz_set_d(newob->z, d);
    }
    return newob;
}

/* forward...: */
static PyObject *f2q_internal(PympfrObject* self, PympfrObject* err,
        unsigned int bits, int mayz);
static PyObject* Pympfr_f2q(PyObject *self, PyObject *args);

static PympqObject *
PyFloat2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pympq_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpq() does not handle nan");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpq() does not handle infinity");
            return NULL;
        }
        mpq_set_d(newob->q, d);
    }
    return newob;
}

/* forward */
static PympfrObject *PyStr2Pympfr(PyObject *s, long base, mpfr_prec_t bits);

static PympfrObject *
PyFloat2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob = 0;

    assert(PyFloat_Check(self));
    if (!bits)
        bits = DBL_MANT_DIG;
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "PyFloat2Pympfr(%p,%ld)\n", self, (long) bits);
#endif
    if ((newob = Pympfr_new(bits)))
        newob->rc = mpfr_set_d(newob->f, PyFloat_AS_DOUBLE(self), context->now.mpfr_round);
    return newob;
}

static PympfrObject *
Pympz2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pympfr_new(bits)))
        newob->rc = mpfr_set_z(newob->f, Pympz_AS_MPZ(self), context->now.mpfr_round);
    return newob;
}

static PympfrObject *
Pyxmpz2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob;

    assert(Pyxmpz_Check(self));
    if ((newob = Pympfr_new(bits)))
        newob->rc = mpfr_set_z(newob->f, Pympz_AS_MPZ(self), context->now.mpfr_round);
    return newob;
}

static PympzObject *
Pympfr2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pympfr_Check(self));
    if ((newob = Pympz_new())) {
        if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpz() does not handle nan");
            return NULL;
        }
        if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpz() does not handle infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(newob->z, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    }
    return newob;
}

static PyxmpzObject *
Pympfr2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pympfr_Check(self));
    if ((newob = Pyxmpz_new())) {
        if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.xmpz() does not handle nan");
            return NULL;
        }
        if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.xmpz() does not handle infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(newob->z, Pympfr_AS_MPFR(self), context->now.mpfr_round);
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

static PympqObject *
Pympfr2Pympq(PyObject *self)
{
    return (PympqObject*) Pympfr_f2q(self, 0);
}

static PympfrObject *
Pympq2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob;

    assert(Pympq_Check(self));
    if (!(newob = Pympfr_new(bits)))
        return NULL;
    newob->rc = mpfr_set_q(newob->f, Pympq_AS_MPQ(self), context->now.mpfr_round);
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

/*
 * long->mpfr delegates via long->mpz->mpfr to avoid duplicating
 * the above-seen dependencies; ditto long->mpq
 */
static PympfrObject *
PyLong2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *newob;
    PyObject *temp = (PyObject*)PyLong2Pympz(self);

    if (!temp)
        return NULL;
    newob = Pympz2Pympfr(temp, bits);
    Py_DECREF(temp);
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

/*
 * mpz conversion from string includes from-binary (base-256 LSB string
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

    if (256 == base) {
        /* Least significant octet first */
        int negative = 0;

        if (cp[len-1] == 0xFF) {
            negative = 1;
            --len;
        }
        mpz_set_si(z, 0);
        mpz_import(z, len, -1, sizeof(char), 0, 0, cp);
        if (negative)
            mpz_neg(z, z);
    }
    else {
        /* Don't allow NULL characters */
        for (i=0; i<len; i++) {
            if (cp[i] == '\0') {
                VALUE_ERROR("string without NULL characters expected");
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

    if (256 == base) {
        /* TODO: better factoring of str2mpz (for speed) */
        int topper, isnega, numlen;
        PyObject *s;
        PympzObject *numerator, *denominator;

        if (len < 6) {
            VALUE_ERROR("invalid mpq binary (too short)");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return 0;
        }
        topper = cp[3] & 0x7f;
        isnega = cp[3] & 0x80;
        numlen = cp[0] + 256 * (cp[1] + 256 * (cp[2] + 256 * topper));
        if (len < (4 + numlen + 1)) {
            VALUE_ERROR("invalid mpq binary (num len)");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return 0;
        }
        s = PyBytes_FromStringAndSize((char*)cp + 4, numlen);
        numerator = PyStr2Pympz(s, 256);
        Py_DECREF(s);
        if (!numerator) {
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return 0;
        }
        if (mpz_sgn(numerator->z) < 0) {
            VALUE_ERROR("invalid mpq binary (num sgn)");
            Py_DECREF((PyObject*)newob);
            Py_DECREF((PyObject*)numerator);
            Py_XDECREF(ascii_str);
            return 0;
        }
        if (isnega)
            mpz_neg(numerator->z, numerator->z);
        s = PyBytes_FromStringAndSize((char*)cp+4+numlen, len-4-numlen);
        denominator = PyStr2Pympz(s,256);
        Py_DECREF(s);
        if (!denominator) {
            Py_DECREF((PyObject*)newob);
            Py_DECREF((PyObject*)numerator);
            Py_XDECREF(ascii_str);
            return 0;
        }
        if (mpz_sgn(denominator->z) != 1) {
            VALUE_ERROR("invalid mpq binary (den sgn)");
            Py_DECREF((PyObject*)newob);
            Py_DECREF((PyObject*)numerator);
            Py_DECREF((PyObject*)denominator);
            Py_XDECREF(ascii_str);
            return 0;
        }
        mpq_set_num(newob->q, numerator->z);
        mpq_set_den(newob->q, denominator->z);
        mpq_canonicalize(newob->q);
        Py_DECREF((PyObject*)numerator);
        Py_DECREF((PyObject*)denominator);
    }
    else {
        /* Don't allow NULL characters */
        for (i=0; i<len; i++) {
            if (cp[i] == '\0') {
                VALUE_ERROR("string without NULL characters expected");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }
        }
        /* trickily delegate the rest to GMP avoiding allocations/copies */
        {
            char* whereslash = strchr((char*)cp,'/');
            char* wheredot = 0;
            if (whereslash) {
                *whereslash = 0;
            }
            else {
                wheredot = strchr((char*)cp, '.');
                if (wheredot) {
                    PympfrObject* temp = PyStr2Pympfr(stringarg, base, (mpfr_prec_t)len*4);
                    if (temp) {
                        newob = Pympfr2Pympq((PyObject*)temp);
                        Py_DECREF((PyObject*)temp);
                    }
                    return newob;
                }
            }
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
                    ZERO_ERROR("mpq: zero denominator");
                    return NULL;
                }
                mpq_canonicalize(newob->q);
            }
            else {
                mpz_set_ui(mpq_denref (newob->q), 1);
            }
        }
    }
    Py_XDECREF(ascii_str);
    return newob;
}

/*
 * mpfr conversion from string includes from-binary (base-256, format is
 * explained later) and 'true' from-string (bases 2 to 62), where exponent
 * if any is denoted by 'e' if base<=10, else by '@', and is always decimal.
 */
static PympfrObject *
PyStr2Pympfr(PyObject *s, long base, mpfr_prec_t bits)
{
    PympfrObject *newob;
    unsigned char *cp;
    mpfr_prec_t prec;
    Py_ssize_t i, len;
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
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }
    //~ if ((!bits) && (base ==256)) {
        //~ prec = 8 * (len - 5);
        //~ if ((len >= 5) && (cp[0] & 8)) {
            //~ bits = 0;
            //~ for(i=4; i>0; --i) {
                //~ bits = (bits<<8) | cp[i];
            //~ }
        //~ }
    //~ }

    if (bits > 0) {
        prec = bits;
    }
    else { /* precision to be defaulted or fetched */
        if (base == 256) {  /* it may be encoded for fetching */
            prec = (mpfr_prec_t)(8 * (len - 5));      /* default precision */
            if ((len>=5) && (cp[0]&8)) { /* precision must be fetched */
                prec = 0;
                for (i=4; i>0; --i) {
                    prec = (prec << 8) | cp[i];
                }
            }
        }
        else { /* true-string, never encoded, just default it */
            prec = context->now.mpfr_prec;
        }
    }
    if (prec < MPFR_PREC_MIN)
        prec = MPFR_PREC_MIN;

    if (!(newob = Pympfr_new(prec))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    if (256 == base) {
        /*
         * binary format for MP floats: first, a code-byte, then, a LSB
         * 4-byte unsigned int (exponent magnitude), then the "mantissa"
         * (actually, "significand", but "mantissa" is the usual term...)
         * in MSB form.
         *
         * The codebyte encodes both the signs, exponent and result, or
         * also the zeroness of the result (in which case, nothing more).
         */
        mpfr_t digit;
        int codebyte = cp[0];
        int resusign = codebyte & 1;
        int exposign = codebyte & 2;
        int resuzero = codebyte & 4;
        int precilen = (codebyte & 8)?4:0;
        unsigned int expomag = 0;

        /* mpfr zero has a very compact (1-byte) binary encoding!-) */
        if (resuzero) {
            newob->rc = mpfr_set_ui(newob->f, 0, context->now.mpfr_round);
            return newob;
        }

        /* all other numbers are 6+ bytes: codebyte, 4-byte exp, 1+
         * bytes for the mantissa; check this string is 6+ bytes
         */
        if (len < 6 + precilen) {
            VALUE_ERROR("string too short to be a gmpy2.mpfr binary encoding");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }
        /* reconstruct exponent */
        for (i = 4 + precilen; i > precilen; --i) {
            expomag = (expomag<<8) | cp[i];
        }

        /* reconstruct 'mantissa' (significand) */
        mpfr_set_si(newob->f, 0, context->now.mpfr_round);
        mpfr_init2(digit, prec);
        for (i = 5 + precilen; i<len; i++) {
            mpfr_set_ui(digit, cp[i], context->now.mpfr_round);
            mpfr_div_2ui(digit, digit, (unsigned long)((i-4-precilen) * 8),
                         context->now.mpfr_round);
            mpfr_add(newob->f, newob->f, digit, context->now.mpfr_round);
        }
        mpfr_clear(digit);
        /* apply exponent, with its appropriate sign */
        if (exposign)
            mpfr_div_2ui(newob->f, newob->f, 8*expomag, context->now.mpfr_round);
        else
            mpfr_mul_2ui(newob->f, newob->f, 8*expomag, context->now.mpfr_round);
        /* apply significand-sign (sign of the overall number) */
        if (resusign)
            mpfr_neg(newob->f, newob->f, context->now.mpfr_round);
    }
    else {
        /* Don't allow NULL characters */
        for (i=0; i<len; i++) {
            if (cp[i] == '\0') {
                VALUE_ERROR("string without NULL characters expected");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }
        }
        /* delegate the rest to MPFR */
        if (-1 == mpfr_set_str(newob->f, (char*)cp, base, context->now.mpfr_round)) {
            VALUE_ERROR("invalid digits");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
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

/*
 * mpfr->long delegates via mpfr->mpz->long to avoid duplicating
 * the above-seen thorny dependencies; ditto mpq->long
 */
static PyObject *
Pympfr2PyLong(PympfrObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympfr2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz2PyLong(temp);
    Py_DECREF((PyObject*)temp);
    return result;
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

/*
 * mpfr->int delegates via mpfr->mpz->int for convenience; ditto mpq->int
 */

#ifdef PY2
static PyObject *
Pympfr2PyInt(PympfrObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympfr2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz_To_Integer(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}

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
Pympfr2PyFloat(PympfrObject *self)
{
    double res = mpfr_get_d(self->f, context->now.mpfr_round);

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

/*
 * helper functions for mpfr->binary conversion
 * hof: maps a hex-digit character into 0..15
 * di256: maps two hex-digits chars into 0..255
 */

static int
hof(int hedi)
{
    static char table[] = "0123456789abcdef";
    char* p = strchr(table, tolower(hedi));

    assert(hedi && p);
    return (int)(p-table);
}

static char
di256(int di1, int di2)
{
    return (char)(hof(di2)+16*hof(di1));
}

/*
 * build binary representation of mpfr (see format description above)
 */

static PyObject *
Pympfr2binary(PympfrObject *self)
{
    size_t size, hexdigs, i, j;
    char *buffer, *aux;
    PyObject *result;
    int sign, codebyte;
    mpfr_exp_t the_exp;
    long lexp, lprec;
    int lexpodd, extrabyte;

    /* prepare codebyte */
    sign = mpfr_sgn(self->f);
    if (sign == 0) {
        /* 0 -> single codebyte with 'zerovalue' bit set */
#ifdef PY3
        return Py_BuildValue("y", "\004");
#else
        return Py_BuildValue("s", "\004");
#endif
        /* codebyte = 0; */
    }
    else if (sign < 0) {
        codebyte = 1;
        mpfr_neg(self->f, self->f, context->now.mpfr_round);
    }
    else {
        codebyte = 0;
    }

    /* get buffer of base-16 digits */
    buffer  = mpfr_get_str(0, &the_exp, 16, 0, self->f, context->now.mpfr_round);

    /* strip trailing zeros */
    hexdigs = strlen(buffer) - 1;
    while ((hexdigs >= 1) && (buffer[hexdigs] == '0'))
        buffer[hexdigs--] = 0x00;

    /* no need to worry about null-buffer as x->f==0.0 was
     * already handled above (see first test on 'sign').
     */
    /* restore correct sign to x->f if it was changed! */
    if (codebyte) {
        mpfr_neg(self->f, self->f, context->now.mpfr_round);
    }
    hexdigs = strlen(buffer);
    /* adjust exponent, & possibly set codebyte's expo-sign bit.
     * note the_exp is base-16 exp, while we need to have it in
     * base-256 -- so it's basically halved (but, with care...!).
     */
    if (the_exp<0) {
        codebyte |= 2;
        the_exp = -the_exp;
    }
    lexp = the_exp;
    lexpodd = lexp & 1;
    lexp = lexp/2 + lexpodd;
    if (lexpodd && (codebyte&2))
        --lexp;
    /* remember we also store precision explicitly */
    codebyte |= 8;

    /* allocate suitably-sized, uninitialized Python string */
    size = (hexdigs + 1) / 2;
    /* allocate an extra byte if lexpodd and hexdigs is even */
    extrabyte = lexpodd & ~hexdigs;
    result = PyBytes_FromStringAndSize(0, 1+4+size+4+extrabyte);
    if (!result)
        return NULL;
    /* set the data to the new Python string's buffer */
    aux = PyBytes_AS_STRING(result);
    /* codebyte first */
    aux[0] = (char)codebyte;
    /* then precision */
    lprec = mpfr_get_prec(self->f);
    for (i=0; i<4; ++i) {
        aux[i+1] = (char)(lprec & 0xFF);
        lprec >>= 8;
    }
    /* then exponent */
    for (i=0; i<4; ++i) {
        aux[4+i+1] = (char)(lexp & 0xFF);
        lexp >>= 8;
    }
    /* then mantissa, grouping 2 hex digits per base-256 digits;
     * with some care for the first & last ones...
     */
    j=0; i=0;
    if (lexpodd) {
        aux[i+9] = di256('0',buffer[0]);
        j=1; i=1;
    }
    for (; i<size+extrabyte; ++i) {
        int secdig = (j+1)<hexdigs? buffer[j+1] : '0';
        aux[i+9] = di256(buffer[j], secdig);
        j += 2;
    }

    PyMem_Free(buffer);
    return result;
}

/*
 * format mpz into any base (2 to 62)
 */
static char* ztag = "mpz(";
static PyObject *
mpz_ascii(mpz_t z, int base, int with_tag)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if ((base != 0) && ((base < 2) || (base > 62))) {
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
    if (with_tag) {
       strcpy(p, ztag);
       p += strlen(p);
    }
    if (negative)
        *(p++) = '-';
#ifdef PY2
    if (base == 8) {
        *(p++) = '0';
#else
    if (base == 2) {
        *(p++) = '0';
        *(p++) = 'b';
    }
    else if (base == 8) {
        *(p++) = '0';
        *(p++) = 'o';
#endif
    }
    else if (base == 16) {
        *(p++) = '0';
        *(p++) = 'x';
    }

    mpz_get_str(p, base, z);     /* Doesn't return number of characters */
    p = buffer + strlen(buffer); /* Address of NULL byte */
#ifdef PY2
    if (with_tag && !mpz_fits_slong_p(z))
        *(p++) = 'L';
#endif
    if (with_tag)
        *(p++) = ')';
    result = PyBytes_FromStringAndSize(buffer, p - buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}

/*
 * format xmpz into any base (2 to 62)
 */
static char* xztag = "xmpz(";
static PyObject *
xmpz_ascii(mpz_t z, int base, int with_tag)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if ((base != 0) && ((base < 2) || (base > 62))) {
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
    if (with_tag) {
       strcpy(p, xztag);
       p += strlen(p);
    }
    if (negative)
        *(p++) = '-';
#ifdef PY2
    if (base == 8) {
        *(p++) = '0';
#else
    if (base == 2) {
        *(p++) = '0';
        *(p++) = 'b';
    }
    else if (base == 8) {
        *(p++) = '0';
        *(p++) = 'o';
#endif
    }
    else if (base == 16) {
        *(p++) = '0';
        *(p++) = 'x';
    }

    mpz_get_str(p, base, z);     /* Doesn't return number of characters */
    p = buffer + strlen(buffer); /* Address of NULL byte */
#ifdef PY2
    if (with_tag && !mpz_fits_slong_p(z))
        *(p++) = 'L';
#endif
    if (with_tag)
        *(p++) = ')';
    result = PyBytes_FromStringAndSize(buffer, p - buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}

static PyObject *
Pympz_ascii(PympzObject *self, int base, int with_tag)
{
#ifdef PY3
    PyObject *result, *temp;

    assert(Pympz_Check((PyObject*)self));
    temp = mpz_ascii(self->z, base, with_tag);
    if (!temp)
        return NULL;
    result = PyUnicode_FromString(PyBytes_AS_STRING(temp));
    Py_DECREF(temp);
    return result;
#else
    assert(Pympz_Check((PyObject*)self));
    return mpz_ascii(self->z, base, with_tag);
#endif
}

static PyObject *
Pyxmpz_ascii(PyxmpzObject *self, int base, int with_tag)
{
#ifdef PY3
    PyObject *result, *temp;

    assert(Pyxmpz_Check((PyObject*)self));
    temp = xmpz_ascii(self->z, base, with_tag);
    if (!temp)
        return NULL;
    result = PyUnicode_FromString(PyBytes_AS_STRING(temp));
    Py_DECREF(temp);
    return result;
#else
    assert(Pyxmpz_Check((PyObject*)self));
    return xmpz_ascii(self->z, base, with_tag);
#endif
}

static char* qtag = "mpq(";
static int qden_1(mpq_t q)
{
    return 0 == mpz_cmp_ui(mpq_denref(q),1);
}

static PyObject *
Pympq_ascii(PympqObject *self, int base, int with_tag)
{
    PyObject *result = 0;
    PyObject *numstr = mpz_ascii(mpq_numref(self->q), base, 0);
    PyObject *denstr = 0;
    PyObject *temp = 0;

    if (!numstr)
        return NULL;

    if (with_tag || !qden_1(self->q)) {
        denstr = mpz_ascii(mpq_denref(self->q), base, 0);
        if (!denstr) {
            Py_DECREF(numstr);
            return NULL;
        }
    }

    if (with_tag) {
        result = PyBytes_FromString(qtag);
        if (!result) {
            Py_XDECREF(numstr);
            Py_XDECREF(denstr);
            return NULL;
        }
        PyBytes_ConcatAndDel(&result, numstr);
#ifdef PY2
        if (!mpz_fits_slong_p(mpq_numref(self->q))) {
            temp = PyBytes_FromString("L");
            if (!temp) {
                Py_XDECREF(denstr);
                Py_XDECREF(result);
                return NULL;
            }
            PyBytes_ConcatAndDel(&result, temp);
            if (!result) {
                Py_XDECREF(denstr);
                return NULL;
            }
        }
#endif
    }
    else {
        result = numstr;
        numstr = 0;
    }
    if (denstr) {
        char* separator = with_tag?",":"/";
        temp = PyBytes_FromString(separator);
        if (!temp) {
            Py_XDECREF(denstr);
            Py_XDECREF(result);
            return NULL;
        }
        PyBytes_ConcatAndDel(&result, temp);
        if (!result) {
            Py_DECREF(denstr);
            return NULL;
        }
        PyBytes_ConcatAndDel(&result, denstr);
#ifdef PY2
        if (with_tag && !mpz_fits_slong_p(mpq_denref(self->q))) {
            temp = PyBytes_FromString("L");
            if (!temp) {
                Py_XDECREF(result);
                return NULL;
            }
            PyBytes_ConcatAndDel(&result, temp);
            if (!result)
                return NULL;
        }
#endif
    }
    if (with_tag && result) {
        temp = PyBytes_FromString(")");
        if (!temp) {
            Py_XDECREF(result);
            return NULL;
        }
        PyBytes_ConcatAndDel(&result, temp);
        if (!result)
            return NULL;
    }
#ifdef PY3
    temp = PyUnicode_FromString(PyBytes_AS_STRING(result));
    Py_DECREF(result);
    return temp;
#else
    return result;
#endif
}

#define OP_TAG 1
#define OP_RAW 2
static char ftag[]="mpfr('";
/*
 * format mpfr into any base (2 to 62)
 * digits: number of digits to ask MPFR for (0=all of
 *     them) -- fewer will be given, if fewer significant
 * minexfi: format as mantissa-exponent if exp<minexfi
 * maxexfi: format as mantissa-exponent if exp>maxexfi
 *     note that, e.g., minexfi=0, maxexfi=-1, means
 *     "always format as mantissa-exponent".  If not
 *     mantissa-exponent, the number will be formatted
 *     as "fixed point" (FP).  Note the decimal point
 *     is _always_ explicitly inserted by this function
 *     (except when bit OP_RAW is set in optionflags).
 * optionflags: bitmap of option-values; currently:
 *     OP_TAG (1): add the mpfr('...') tag
 *     OP_RAW (2): ignore minexfi/maxexfi/OP_TAG
 *         and return a 3-element tuple digits/exponent/rprec
 *         (as GMP gives them) for Python formatting;
 *         'digits' may include a '-' sign, but no decimal
 *         point, nor tag, nor any exponent-indicator.
 *     other bits are currently ignored
 */
static PyObject *
Pympfr_ascii(PympfrObject *self, int base, int digits,
    int minexfi, int maxexfi, int optionflags)
{
    PyObject *result;
#ifdef PY3
    PyObject *temp;
#endif
    char *buffer;
    mpfr_exp_t the_exp;
    size_t buflen;

    /* check arguments are valid */
    assert(Pympfr_Check((PyObject*)self));
    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    /* Process special cases first */
    if (!(mpfr_regular_p(self->f))) {
        if (mpfr_nan_p(self->f)) {
            if (optionflags & OP_RAW)
                result = Py_BuildValue("(sii)", "nan", 0, 0);
            else
                if (optionflags & OP_TAG)
                    result = Py_BuildValue("s", "mpfr('nan')");
                else
                    result = Py_BuildValue("s", "nan");
        }
        else if (mpfr_inf_p(self->f) && !mpfr_signbit(self->f)) {
            if (optionflags & OP_RAW)
                result = Py_BuildValue("(sii)", "inf", 0, 0);
            else
                if (optionflags & OP_TAG)
                    result = Py_BuildValue("s", "mpfr('inf')");
                else
                    result = Py_BuildValue("s", "inf");
        }
        else if (mpfr_inf_p(self->f) && mpfr_signbit(self->f)) {
            if (optionflags & OP_RAW)
                result = Py_BuildValue("(sii)", "-inf", 0, 0);
            else
                if (optionflags & OP_TAG)
                    result = Py_BuildValue("s", "mpfr('-inf')");
                else
                    result = Py_BuildValue("s", "-inf");
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self->f)) {
            if (optionflags & OP_RAW)
                result = Py_BuildValue("(sii)", "-0", 0,
                                       mpfr_get_prec(self->f));
            else
                if (optionflags & OP_TAG)
                    result = Py_BuildValue("s", "mpfr('-0.0e0')");
                else
                    result = Py_BuildValue("s", "-0.0e0");
        }
        else {
            if (optionflags & OP_RAW)
                result = Py_BuildValue("(sii)", "0", 0,
                                       mpfr_get_prec(self->f));
            else
                if (optionflags & OP_TAG)
                    result = Py_BuildValue("s", "mpfr('0.0e0')");
                else
                    result = Py_BuildValue("s", "0.0e0");
        }
        return result;
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, context->now.mpfr_round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in Pympfr_ascii");
        return NULL;
    }

    if (optionflags & OP_RAW) {
        result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
        PyMem_Free(buffer);
        return result;
    }
    else {
        char expobuf[24];
        char auprebuf[24];
        int isfp = 1;   /* flag: fixed-point format (FP)? */
        int isnegative = 0;
        size_t size;

        /* check if it is negative */
        if (buffer[0]==0x2d)
            isnegative = 1;

        /* strip trailing zeros */
        buflen = strlen(buffer) - 1;
        while ((buflen >= (2 + isnegative)) && (buffer[buflen] == '0'))
            buffer[buflen--] = 0x00;

        /* insert formatting elements (decimal-point, leading or
           trailing 0's, other indication of exponent...) */
        buflen = strlen(buffer);

        /* account for the decimal point that is always inserted */
        size = buflen + 1;

        /* compute size of needed Python string */
        if (optionflags & OP_TAG) {
            size += strlen(ftag) + 2;
            if (mpfr_get_prec(self->f) != DBL_MANT_DIG) {
                sprintf(auprebuf, ",%ld", (long)mpfr_get_prec(self->f));
                size += strlen(auprebuf);
            }
        }

        /* exponential format */
        if (the_exp < minexfi || the_exp > maxexfi) {
            /* add exponent-length + 1 for '@' or 'e' marker */
            sprintf(expobuf, "%ld", the_exp - 1);
            size += strlen(expobuf) + 1;
            isfp = 0;
        }
        /* 'fixed-point' format */
        else {
            /* add number of leading or trailing 0's */
            if (the_exp <= 0) {
                /* add leading 0's */
                size += abs(the_exp) + 1;
            }
            else {
                /* add trailing 0's if needed */
                if (the_exp >= (buflen - isnegative))
                    size += (the_exp - (buflen - isnegative)) + 1;
            }
        }

        /* allocate the string itself (uninitialized, as yet) */
        result = PyBytes_FromStringAndSize(0, size);

        {
            /* proceed with building the string-buffer value */
            char* pd = PyBytes_AS_STRING(result);
            char* ps = buffer;

            /* insert leading tag if requested */
            if (optionflags & OP_TAG) {
                char* pt = ftag;
                while(*pt) *pd++ = *pt++;
            }

            /* copy sign if it's there */
            if (*ps == '-') {
                *pd++ = *ps++;
            }

            /* insert a leading-0 if needed for non-positive-exp FP,
             * else just copy the leading digit (goes before '.')
             */
            if (isfp && the_exp<=0)
                *pd++ = '0';
            else if (*ps)
                *pd++ = *ps++;
            else
                *pd++ = '0';

            /* insert what else goes before '.' for FP */
            if (isfp && the_exp > 1) {
                /* number of digits-to-copy before the '.' */
                int dtc = the_exp - 1;
                /* copy requested # of digits as long as there
                 * are still digits to copy in the buffer
                 */
                while (dtc && *ps) {
                    *pd++ = *ps++;
                    --dtc;
                }
                /* insert trailing 0's before the '.' if
                 * needed to make up the total # digits
                 * that go before the '.' in FP/large exp
                 */
                while (dtc > 0) {
                    *pd++ = '0';
                    --dtc;
                }
            }

            /* the decimal-point is _always_ explicitly there */
            *pd++ = '.';

            /* as is at least 1 trailing-digit after it, if FP,
             * so put a 0 if no more digits to copy
             */
            if (isfp && !*ps)
                *pd++ = '0';

            /* in FP with negative exp, we have more leading 0's
             * after the decimal-point before copying the digits
             * from the buffer
             */
            if (isfp && the_exp<0) {
                int dtc = abs(the_exp);
                while (dtc>0) {
                    *pd++ = '0';
                    --dtc;
                }
            }

            /* copy all remaining digits from buffer, if any */
            while(*ps) *pd++ = *ps++;

            /* insert marker-and-exponent if _not_ FP */
            if (!isfp) {
                char* pe = expobuf;
                *pd++ = (base<=10)?'e':'@';
                while(*pe) *pd++ = *pe++;
            }

            /* insert trailing-part of the tag if needed */
            if (optionflags & OP_TAG) {
                char* pe = auprebuf;
                *pd++ = '\'';
                if (mpfr_get_prec(self->f) != DBL_MANT_DIG)
                    while(*pe) *pd++ = *pe++;
                *pd++ = ')';
            }
        }
        PyMem_Free(buffer);
#ifdef PY3
        temp = PyUnicode_FromString(PyBytes_AS_STRING(result));
        Py_DECREF(result);
        return temp;
#else
        return result;
#endif
    }
}

/* Classify an object as a type of number. If an object is recognized as a
 * number, it must be properly converted by the routines below.
 */

static int isComplex(PyObject* obj)
{
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "isComplex: object type is %s\n", Py_TYPE(obj)->tp_name);
#endif
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (Pympc_Check(obj))       return 1;
    if (PyFloat_Check(obj))     return 1;
    if (PyComplex_Check(obj))   return 1;
    if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal"))  return 1;
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isReal(PyObject* obj)
{
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "isReal: object type is %s\n", Py_TYPE(obj)->tp_name);
#endif
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (PyFloat_Check(obj))     return 1;
    if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal"))  return 1;
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isRational(PyObject* obj)
{
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "isRational: object type is %s\n", Py_TYPE(obj)->tp_name);
#endif
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isInteger(PyObject* obj)
{
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "isInteger: object type is %s\n", Py_TYPE(obj)->tp_name);
#endif
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

static PympqObject*
anynum2Pympq(PyObject* obj)
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
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pympq(obj);
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympq(obj);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympq(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympq(obj);
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr,"anynum2Pympq(%p)->%p\n", obj, newob);
#endif

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
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr,"Pympq_From_Rational(%p)->%p\n", obj, newob);
#endif

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
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pympz(obj);
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympz(obj);
    }
    else if (PyNumber_Check(obj) && !strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = PyLong2Pympz(s);
            Py_DECREF(s);
        }
    }
    else if (PyNumber_Check(obj) && !strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pympz((PyObject *)temp);
            Py_DECREF(s); Py_DECREF((PyObject*)temp);
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr,"anynum2Pympz(%p)->%p\n", obj, newob);
#endif
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
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pyxmpz(obj);
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pyxmpz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pyxmpz(obj);
    }
    else if (PyNumber_Check(obj) && !strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = PyLong2Pyxmpz(s);
            Py_DECREF(s);
        }
    }
    else if (PyNumber_Check(obj) && !strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pyxmpz((PyObject *)temp);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr,"anynum2Pympz(%p)->%p\n", obj, newob);
#endif
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
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr,"Pympz_From_Integer(%p)->%p\n", obj, newob);
#endif
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

static long
clong_From_Integer(PyObject *obj)
{
    if (PyLong_Check(obj)) {
        return PyLong_AsLong(obj);
    }
#ifdef PY2
    else if (PyInt_Check(obj)) {
        return PyInt_AS_LONG(obj);
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
    }
    TYPE_ERROR("conversion error in clong_From_Integer");
    return -1;
}

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
 * If obj is a PyFloat and bits is 0, then the conversion is done exactly.
 * If obj is a Pympfr and bits is 0 or bits is the same as the precision of
 * obj, then a new reference is created.
 *
 * For all other numerical types with bits = 0, the conversion is rounded to
 * context->now.mpfr_prec.
 */

static PympfrObject*
Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits)
{
    PympfrObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympfr_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpfr is still
         * valid in the current context. */
        if (!bits || mpfr_get_prec(Pympfr_AS_MPFR(obj)) == bits) {
            newob = (PympfrObject*) obj;
            Py_INCREF(obj);
        }
        else {
            newob = Pympfr2Pympfr((PyObject*)obj, bits);
        }
    }
    else if (Pympfr_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer valid
         * and mpfr_check_range needs to be called. */
        if ((newob = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(obj))))) {
            mpfr_set(newob->f, Pympfr_AS_MPFR(obj), context->now.mpfr_round);
            newob->round_mode = ((PympfrObject*)obj)->round_mode;
            newob->rc = ((PympfrObject*)obj)->rc;
            newob->rc = mpfr_check_range(newob->f, newob->rc, newob->round_mode);
        }
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympfr(obj, bits);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympfr(obj, bits);
#endif
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pympfr(obj, bits);
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympfr(obj, bits);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympfr(obj, bits);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympfr(obj, bits);
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympfr(s, 10, bits);
            if (!newob) {
                Py_DECREF(s);
                return NULL;
            }
            Py_DECREF(s);
        }
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pympfr((PyObject *)temp, bits);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "Pympfr_From_Real(%p,%ld)->%p (%ld)\n", obj,
                (long)bits, newob, newob != 0 ? (long)mpfr_get_prec(newob->f) : -1);
#endif
    return newob;
}

/*
 * coerce any number to a mpz
 */
int
Pympz_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympzObject* newob = Pympz_From_Integer(arg);
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mpz_conv_arg(%p)->%p\n", arg, newob);
#endif
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to mpz");
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
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mpq_conv_arg(%p)->%p\n", arg, newob);
#endif
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        if (!PyErr_Occurred()) {
            TYPE_ERROR("argument can not be converted to mpq");
        }
        return 0;
    }
}

/*
 * coerce any number to a mpf
 */

int
Pympfr_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympfrObject* newob = Pympfr_From_Real(arg, 0);

#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mpfr_conv_arg(%p)->%p\n", arg, newob);
#endif
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to mpfr");
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

/* str and repr implementations for mpz */
static PyObject *
Pympfr2str(PympfrObject *self)
{
    /* base-10, FP for exp -2 to 8, no tag */
    return Pympfr_ascii(self, 10, 0, -2, 8, 0);
}

static PyObject *
Pympfr2repr(PympfrObject *self)
{
    /* base-10, always mantissa+exp, with tag */
    return Pympfr_ascii(self, 10, 0, 0, -1, OP_TAG);
}


/* CONSTRUCTORS */
static char doc_mpz[] = "\
mpz(n):\n\
      builds an mpz object with a numeric value n (truncating n\n\
      to its integer part if it's a float or mpfr)\n\
mpz(s,base=0):\n\
      builds an mpz object from a string s made up of digits in the\n\
      given base.  If base=0, binary, octal, or hex Python strings\n\
      are recognized by leading 0b, 0o, or 0x characters, otherwise\n\
      the string is assumed to be decimal. If base=256, s must be a\n\
      gmpy2.mpz portable binary representation as built by the function\n\
      gmpy2.binary (and the .binary method of mpz objects).\n\
";
static PyObject *
Pygmpy_mpz(PyObject *self, PyObject *args)
{
    PympzObject *newob;
    PyObject *obj;
    Py_ssize_t argc;

    TRACE("Pygmpy_mpz() called...\n");

    assert(PyTuple_Check(args));

    argc = PyTuple_Size(args);
    if ((argc < 1) || (argc > 2)) {
        TYPE_ERROR("gmpy2.mpz() requires 1 or 2 arguments");
        return NULL;
    }

    obj = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(obj)) {
        /* build-from-string (ascii or binary) */
        long base=0;
        if (argc == 2) {
            PyObject *pbase = PyTuple_GetItem(args, 1);
            base = clong_From_Integer(pbase);
            if (base == -1 && PyErr_Occurred()) {
                TYPE_ERROR("gmpy2.mpz(): base must be an integer");
                return NULL;
            }
            if ((base!=0) && (base!=256) && ((base<2)||(base>62))) {
                VALUE_ERROR("base for gmpy2.mpz() must be 0, 256, or in the "
                            "interval 2 ... 62");
                return NULL;
            }
        }
        newob = PyStr2Pympz(obj, base);
        if (!newob) {
            return NULL;
        }
    }
    else {
        if (argc==2) {
            TYPE_ERROR("gmpy2.mpz() with numeric argument needs exactly 1 argument");
            return NULL;
        }
        newob = anynum2Pympz(obj);
        if (!newob) {
            if (!PyErr_Occurred())
                TYPE_ERROR("gmpy2.mpz() requires numeric or string argument");
            return NULL;
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "Pygmpy_mpz: created mpz = %ld\n",
                mpz_get_si(newob->z));
#endif
    return (PyObject *) newob;
}

static char doc_xmpz[] = "\
xmpz(n):\n\
    builds an xmpz object from any number n (truncating n\n\
    to its integer part if it's a float or mpfr)\n\
xmpz(s, base=0):\n\
    builds an xmpz object from a string s made up of digits in the\n\
    given base.  If base=0, binary, octal, and hex Python strings\n\
    are recognized by leading 0b, 0o, or 0x characters, otherwise\n\
    the string is assumed to be decimal. If base=256, s must be a\n\
    gmpy2.xmpz portable binary representation as built by the function\n\
    gmpy2.binary (and the .binary method of xmpz objects).\n\
";
static PyObject *
Pygmpy_xmpz(PyObject *self, PyObject *args)
{
    PyxmpzObject *newob;
    PyObject *obj, *obj1;
    Py_ssize_t argc;
    long base = 0;

    TRACE("Pygmpy_xmpz() called...\n");

    assert(PyTuple_Check(args));

    argc = PyTuple_Size(args);
    if ((argc < 1) || (argc > 2)) {
        TYPE_ERROR("gmpy2.xmpz() requires 1 or 2 arguments");
        return NULL;
    }

    obj = PyTuple_GetItem(args, 0);

    if (argc == 2) {
        if (!PyStrOrUnicode_Check(obj)) {
            TYPE_ERROR("gmpy2.xmpz() with numeric argument accepts only 1 argument");
            return NULL;
        }
        obj1 = PyTuple_GetItem(args, 1);
        base = clong_From_Integer(obj1);
        if (base == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.xmpz(): base must be an integer");
            return NULL;
        }
        if ((base!=0) && (base!=256) && ((base<2)||(base>62))) {
            VALUE_ERROR("gmpy2.xmpz(): base must be 0, 256, or in the "
                        "interval 2 ... 62");
            return NULL;
        }
    }

    if (PyStrOrUnicode_Check(obj)) {
        /* build-from-string (ascii or binary) */
        newob = PyStr2Pyxmpz(obj, base);
        if (!newob) {
            if (!PyErr_Occurred()) {
                VALUE_ERROR("gmpy2.xmpz(): invalid string");
            }
            return NULL;
        }
    }
    else {
        newob = anynum2Pyxmpz(obj);
        if (!newob) {
            if (!PyErr_Occurred()) {
                TYPE_ERROR("gmpy2.xmpz() requires integer or string argument");
            }
            return NULL;
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "Pygmpy_xmpz: created xmpz = %ld\n",
                mpz_get_si(newob->z));
#endif
    return (PyObject *) newob;
}

static char doc_mpq[] = "\
mpq(n):\n\
    builds an mpq object with a numeric value n\n\
mpq(n,m):\n\
    builds an mpq object with a numeric value n/m\n\
mpq(s,base=10):\n\
    builds an mpq object from a string s made up of digits in the\n\
    given base.  s may be made up of two numbers in the same base\n\
    separated by a '/' character.  If base=256, s must be a\n\
    gmpy2.mpq portable binary representation as built by the\n\
    gmpy2.qbinary (or the .binary method of mpq objects).\n\
";
static PyObject *
Pygmpy_mpq(PyObject *self, PyObject *args)
{
    PympqObject *newob;
    PyObject *obj;
    int wasnumeric;
    Py_ssize_t argc;

    TRACE("Pygmpy_mpq() called...\n");

    assert(PyTuple_Check(args));

    argc = PyTuple_Size(args);
    if ((argc < 1) || (argc > 2)) {
        TYPE_ERROR("gmpy2.mpq() requires 1 or 2 arguments");
        return NULL;
    }

    obj = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(obj)) {
        /* build-from-string (ascii or binary) */
        long base=10;
        wasnumeric=0;
        if (argc == 2) {
            PyObject *pbase = PyTuple_GetItem(args, 1);
            base = clong_From_Integer(pbase);
            if (base == -1 && PyErr_Occurred()) {
                TYPE_ERROR("gmpy2.mpq(): base must be an integer");
                return NULL;
            }
            if ((base!=0) && (base!=256) && ((base<2)||(base>62))) {
                VALUE_ERROR("base for gmpy2.mpq() must be 0, 256, or in the "
                            "interval 2 ... 62");
                return NULL;
            }
        }
        newob = PyStr2Pympq(obj, base);
        if (!newob) {
            return NULL;
        }
    }
    else {
        wasnumeric=1;
        newob = anynum2Pympq(obj);
        if (!newob) {
            if (!PyErr_Occurred())
                TYPE_ERROR("gmpy2.mpq() requires numeric or string argument");
            return NULL;
        }
    }
#ifdef DEBUG
    if (global.debug) {
        fputs("Pygmpy_mpq: created mpq = ", stderr);
        mpq_out_str(stderr, 10, newob->q);
        putc('\n', stderr);
    }
#endif
    if (wasnumeric && argc==2) {
        PympqObject *denominator;
        denominator = anynum2Pympq(PyTuple_GET_ITEM(args, 1));
        if (!denominator) {
            TYPE_ERROR("argument can not be converted to mpq");
            Py_DECREF((PyObject*)newob);
            return NULL;
        }
        if (0==mpq_sgn(Pympq_AS_MPQ(denominator))) {
            ZERO_ERROR("mpq: zero denominator");
            Py_DECREF((PyObject*) newob);
            Py_DECREF((PyObject*)denominator);
            return NULL;
        }
        mpq_div(newob->q, newob->q, denominator->q);
        Py_DECREF((PyObject*)denominator);
    }

    return (PyObject *) newob;
}

static char doc_mpfr[] = "\
mpfr(n):\n\
    builds an mpfr object with a numeric value n (n may be any\n\
    Python number, or an mpz, mpq, or mpfr object) and a default\n\
    precision (in bits) depending on the nature of n\n\
mpfr(n,bits=0):\n\
    as above, but with the specified number of bits (0\n\
    means to use default precision, as above)\n\
mpfr(s,bits=0,base=10):\n\
    builds an mpfr object from a string s made up of\n\
    digits in the given base, possibly with fraction-part (with\n\
    period as a separator) and/or exponent-part (with exponent\n\
    marker 'e' for base<=10, else '@'). If base=256, s must be\n\
    a gmpy2.mpfr portable binary representation as built by the\n\
    function gmpy2.binary (and the .binary method of mpfr objects).\n\
    The resulting mpfr object is built with a default precision (in\n\
    bits) if bits is 0 or absent, else with the specified number\n\
    of bits.\n\
";
static PyObject *
Pygmpy_mpfr(PyObject *self, PyObject *args)
{
    PympfrObject *newob;
    PyObject *obj;
    Py_ssize_t argc;
    mpfr_prec_t bits=0;

    TRACE("Pygmpy_mpfr() called...\n");

    assert(PyTuple_Check(args));

    argc = PyTuple_Size(args);
    if ((argc < 1) || (argc > 3)) {
        TYPE_ERROR("gmpy2.mpfr() requires 1 to 3 arguments");
        return NULL;
    }

    obj = PyTuple_GetItem(args, 0);

    if (2 <= argc) {
        long sbits;
        PyObject *pbits = PyTuple_GetItem(args, 1);
        sbits = clong_From_Integer(pbits);
        if (sbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.mpfr(): bits must be an integer");
            return NULL;
        }
        if (sbits<0) {
            VALUE_ERROR("bits for gmpy2.mpfr() must be >= 0");
            return NULL;
        }
        bits = sbits;
    }

    if (PyStrOrUnicode_Check(obj) || PyUnicode_Check(obj)) {
        /* build-from-string (ascii or binary) */
        long base=10;
        if (3 == argc) {
            PyObject *pbase = PyTuple_GetItem(args, 2);
            base = clong_From_Integer(pbase);
            if (base == -1 && PyErr_Occurred()) {
                TYPE_ERROR("gmpy2.mpfr(): base must be an integer");
                return NULL;
            }
            if ((base!=0) && (base!=256) && ((base<2)||(base>62))) {
                VALUE_ERROR("base for gmpy2.mpfr() must be 0, 256, or in the "
                            "interval 2 ... 62");
                return NULL;
            }
        }
        newob = PyStr2Pympfr(obj, base, bits);
        if (!newob) {
            return NULL;
        }
    }
    else {
        if (argc==3) {
            TYPE_ERROR("gmpy2.mpfr() with numeric 1st argument needs 1 or 2 arguments");
            return NULL;
        }
        newob = Pympfr_From_Real(obj, bits);
        if (!newob) {
            if (!PyErr_Occurred())
                TYPE_ERROR("gmpy2.mpfr() requires numeric or string argument");
            return NULL;
        }
    }
#ifdef DEBUG
    if (global.debug) {
        fputs("Pygmpy_mpfr: created mpfr = ", stderr);
        mpfr_out_str(stderr, 10, 0, newob->f, context->now.mpfr_round);
        fprintf(stderr," bits=%ld (%ld)\n",
                (long)mpfr_get_prec(newob->f), (long)bits);
    }
#endif
    return (PyObject *) newob;
} /* Pygmpy_mpfr() */

PyDoc_STRVAR(doc_g_mpc,
"mpc(n, [prec=(i,i) | rnd=(i,i]) -> mpc object\n\n"
"Return an mpc object by converting a numeric value 'n' into a\n"
"complex number. If 'prec' is omitted, then get_mpc_precision() is\n"
"used. If 'rnd' is omitted, then get_mpc_round() is used.\n\n"
"mpc(s, [base = 10 | prec=(i,i) | rnd=(i,i)]) -> mpc object\n\n"
"Return an mpc object by converting a string 's' into a complex\n"
"number. If 'base' is omitted, then a base 10 representation is\n"
"assumed otherwise a base between 2 and 36 can be specified.\n"
"If 'prec' is omitted, then get_mpc_precision() is used. If 'rnd'\n"
"is omitted, then get_mpc_round() is used.");
static PyObject *
Pygmpy_mpc(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PympcObject *result;
    PyObject *arg0 = NULL, *prec_obj = NULL;
    Py_ssize_t dummy, rprec, iprec;
    int base, rmode;
    static char *kwlist[] = {"n", "base", "prec", "rnd", NULL};

    dummy = 0;
    base = 10;
    rprec = context->now.mpc_rprec;
    iprec = context->now.mpc_iprec;
    rmode = context->now.mpc_round;

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "|OiOi", kwlist,
                                      &arg0, &base, &prec_obj, &rmode)))
        return NULL;

    if (base < 2 || base > 36) {
        VALUE_ERROR("base for gmpy2.mpc() must be in the interval 2..36.");
        return NULL;
    }

    if (rmode != context->now.mpc_round && Pymisc_verify_mpc_round(rmode) == -1) {
        VALUE_ERROR("invalid rounding mode for complex arithmetic.");
        return NULL;
    }

    if (prec_obj && PyIntOrLong_Check(prec_obj)) {
        rprec = PyLong_AsSsize_t(prec_obj);
        iprec = rprec;

        /* Also catches return value of -1. */
        if (!Pymisc_verify_mpc_precision(rprec, iprec)) {
            VALUE_ERROR("invalid precision for complex arithmetic");
            return NULL;
        }
    }
    else if (prec_obj && PyTuple_Check(prec_obj)) {
        if (!(PyArg_ParseTuple(prec_obj,
                    "nn;invalid precision for complex arithmetic.",
                    &rprec, &iprec)))
            return NULL;

        if (!Pymisc_verify_mpc_precision(rprec, iprec)) {
            VALUE_ERROR("invalid precision for complex arithmetic");
            return NULL;
        }
    }
    else {
        VALUE_ERROR("invalid precision for complex arithmetic");
        return NULL;
    }

    if (arg0 == NULL) {
        if (!(result = Pympc_new(rprec, iprec)))
            return NULL;
        /* May want to return (NaN NaN) ?? */
        mpc_set_ui(result->c, 0, rmode);
        return (PyObject*)result;
    }

    Py_INCREF(arg0);
    return Py_BuildValue("OiOi", arg0, base,
                         Py_BuildValue("nn", rprec, iprec),
                         rmode);
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
#include "gmpy_mpq.c"
#include "gmpy_mpfr.c"
#include "gmpy_mpc.c"
#include "gmpy_basic.c"
#include "gmpy_mpz_inplace.c"
#include "gmpy_xmpz_inplace.c"

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
    PyObject *tempa = 0, *tempb = 0, *result = 0;

#ifdef DEBUG
    if (global.debug) {
        fprintf(stderr, "rich_compare: type(a) is %s\n", Py_TYPE(a)->tp_name);
        fprintf(stderr, "rich_compare: type(b) is %s\n", Py_TYPE(b)->tp_name);
    }
#endif
    if (CHECK_MPZANY(a) && PyIntOrLong_Check(b)) {
        TRACE("compare (mpz,small_int)\n");
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
    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        TRACE("compare (mpz,mpz)\n");
        return _cmp_to_object(mpz_cmp(Pympz_AS_MPZ(a), Pympz_AS_MPZ(b)), op);
    }
    if (Pympq_Check(a) && Pympq_Check(b)) {
        TRACE("compare (mpq,mpq)\n");
        return _cmp_to_object(mpq_cmp(Pympq_AS_MPQ(a), Pympq_AS_MPQ(b)), op);
    }
    if (Pympfr_Check(a) && Pympfr_Check(b)) {
        TRACE("compare (mpfr,mpfr)\n");
        return _cmp_to_object(mpfr_cmp(Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(b)), op);
    }
    if (isInteger(a) && isInteger(b)) {
        TRACE("compare (mpz,int)\n");
        tempa = (PyObject*)Pympz_From_Integer(a);
        tempb = (PyObject*)Pympz_From_Integer(b);
        c = mpz_cmp(Pympz_AS_MPZ(tempa), Pympz_AS_MPZ(tempb));
        Py_DECREF(tempa);
        Py_DECREF(tempb);
        return _cmp_to_object(c, op);
    }
    if (isRational(a) && isRational(b)) {
        TRACE("compare (mpq,rational)\n");
        tempa = (PyObject*)Pympq_From_Rational(a);
        tempb = (PyObject*)Pympq_From_Rational(b);
        c = mpq_cmp(Pympq_AS_MPQ(tempa), Pympq_AS_MPQ(tempb));
        Py_DECREF(tempa);
        Py_DECREF(tempb);
        return _cmp_to_object(c, op);
    }
    if (isReal(a) && isReal(b)) {
        TRACE("compare (mpfr,float)\n");
        /* Handle non-numbers separately. */
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            if (Py_IS_NAN(d)) {
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else if (Py_IS_INFINITY(d)) {
                if (d < 0.0) {
                    return _cmp_to_object(1, op);
                }
                else {
                    return _cmp_to_object(-1, op);
                }
            }
        }
        tempa = (PyObject*)Pympfr_From_Real(a, 0);
        tempb = (PyObject*)Pympfr_From_Real(b, 0);
        c = mpfr_cmp(Pympfr_AS_MPFR(tempa), Pympfr_AS_MPFR(tempb));
        Py_DECREF(tempa);
        Py_DECREF(tempb);
        return _cmp_to_object(c, op);
    }
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
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv, /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc)  Pympz_To_Integer,       /* nb_index                */
};

#else
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_div2,           /* nb_divide               */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
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
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
    (binaryfunc) Pyxmpz_inplace_floordiv,/* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc)  Pyxmpz_To_Integer,      /* nb_index                */
};

#else
static PyNumberMethods xmpz_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_div2,           /* nb_divide               */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
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
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_div2,           /* nb_divide               */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
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

#ifdef PY3
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_div2,           /* nb_divide               */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,              /* nb_negative             */
    (unaryfunc) Pympfr_pos,               /* nb_positive             */
    (unaryfunc) Pympfr_abs,              /* nb_absolute             */
    (inquiry) Pympfr_nonzero,             /* nb_bool                 */
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
    (binaryfunc) Pympany_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pympany_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympfr_getseters[] =
{
    {"precision", (getter)Pympfr_getprec_attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)Pympfr_getrc_attrib, NULL, "return code", NULL},
    {NULL}
};

#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
    0,            /* nb_add                  */
    0,            /* nb_subtract             */
    0,            /* nb_multiply             */
    0,            /* nb_remaider             */
    0,         /* nb_divmod               */
    0,           /* nb_power                */
    0,              /* nb_negative             */
    0,              /* nb_positive             */
    0,              /* nb_absolute             */
    0,            /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    0,           /* nb_int                  */
        0,                               /* nb_reserved             */
    0,          /* nb_float                */
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
    0,       /* nb_floor_divide         */
    0,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
    0,            /* nb_add                  */
    0,            /* nb_subtract             */
    0,            /* nb_multiply             */
    0,           /* nb_divide               */
    0,            /* nb_remaider             */
    0,         /* nb_divmod               */
    0,           /* nb_power                */
    0,              /* nb_negative             */
    0,               /* nb_positive             */
    0,              /* nb_absolute             */
    0,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    0,            /* nb_int                  */
    0,           /* nb_long                 */
    0,          /* nb_float                */
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
    0,       /* nb_floor_divide         */
    0,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyMethodDef Pygmpy_methods [] =
{
    { "_cvsid", Pygmpy_get_cvsid, METH_NOARGS, doc_cvsid },
    { "_copy", Pympany_copy, METH_O, doc_copyg },
    { "acos", Pympfr_acos, METH_O, doc_g_mpfr_acos },
    { "acosh", Pympfr_acosh, METH_O, doc_g_mpfr_acosh },
    { "add", Pympfr_add, METH_VARARGS, doc_g_mpfr_add },
    { "ai", Pympfr_ai, METH_O, doc_g_mpfr_ai },
    { "agm", Pympfr_agm, METH_VARARGS, doc_g_mpfr_agm },
    { "asin", Pympfr_asin, METH_O, doc_g_mpfr_asin },
    { "asinh", Pympfr_asinh, METH_O, doc_g_mpfr_asinh },
    { "atan", Pympfr_atan, METH_O, doc_g_mpfr_atan },
    { "atanh", Pympfr_atanh, METH_O, doc_g_mpfr_atanh },
    { "atan2", Pympfr_atan2, METH_VARARGS, doc_g_mpfr_atan2 },
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
    { "cbrt", Pympfr_cbrt, METH_O, doc_g_mpfr_cbrt },
    { "ceil", Pympfr_ceil, METH_O, doc_g_mpfr_ceil },
    { "check_range", Pympfr_check_range, METH_O, doc_g_mpfr_check_range },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combg },
    { "const_catalan", Pympfr_const_catalan, METH_NOARGS, doc_mpfr_const_catalan },
    { "const_euler", Pympfr_const_euler, METH_NOARGS, doc_mpfr_const_euler },
    { "const_log2", Pympfr_const_log2, METH_NOARGS, doc_mpfr_const_log2 },
    { "const_pi", Pympfr_const_pi, METH_VARARGS, doc_mpfr_const_pi },
    { "context", Pygmpy_context, METH_NOARGS, doc_context },
    { "cos", Pympfr_cos, METH_O, doc_g_mpfr_cos },
    { "cosh", Pympfr_cosh, METH_O, doc_g_mpfr_cosh },
    { "cot", Pympfr_cot, METH_O, doc_g_mpfr_cot },
    { "coth", Pympfr_coth, METH_O, doc_g_mpfr_coth },
    { "csc", Pympfr_csc, METH_O, doc_g_mpfr_csc },
    { "csch", Pympfr_csch, METH_O, doc_g_mpfr_csch },
    { "c_div", Pygmpy_c_div, METH_VARARGS, doc_gmpy_c_div },
    { "c_div_2exp", Pygmpy_c_div_2exp, METH_VARARGS, doc_gmpy_c_div_2exp },
    { "c_divmod", Pygmpy_c_divmod, METH_VARARGS, doc_gmpy_c_divmod },
    { "c_divmod_2exp", Pygmpy_c_divmod_2exp, METH_VARARGS, doc_gmpy_c_divmod_2exp },
    { "c_mod", Pygmpy_c_mod, METH_VARARGS, doc_gmpy_c_mod },
    { "c_mod_2exp", Pygmpy_c_mod_2exp, METH_VARARGS, doc_gmpy_c_mod_2exp },
    { "denom", Pympq_denom, METH_VARARGS, doc_denomg },
    { "digamma", Pympfr_digamma, METH_O, doc_g_mpfr_digamma },
    { "digits", Pympany_digits, METH_VARARGS, doc_g_mpany_digits },
    { "div", Pympfr_div, METH_VARARGS, doc_g_mpfr_div },
    { "divexact", Pygmpy_divexact, METH_VARARGS, doc_divexactg },
    { "divm", Pygmpy_divm, METH_VARARGS, doc_divm },
    { "eint", Pympfr_eint, METH_O, doc_g_mpfr_eint },
    { "erf", Pympfr_erf, METH_O, doc_g_mpfr_erf },
    { "erfc", Pympfr_erfc, METH_O, doc_g_mpfr_erfc },
    { "exp", Pympfr_exp, METH_O, doc_g_mpfr_exp },
    { "expm1", Pympfr_expm1, METH_O, doc_g_mpfr_expm1 },
    { "exp10", Pympfr_exp10, METH_O, doc_g_mpfr_exp10 },
    { "exp2", Pympfr_exp2, METH_O, doc_g_mpfr_exp2 },
    { "fac", Pygmpy_fac, METH_O, doc_fac },
    { "factorial", Pympfr_factorial, METH_O, doc_g_mpfr_factorial },
    { "fib", Pygmpy_fib, METH_O, doc_fib },
    { "fib2", Pygmpy_fib2, METH_O, doc_fib2 },
    { "floor", Pympfr_floor, METH_O, doc_g_mpfr_floor},
    { "fma", Pympfr_fma, METH_VARARGS, doc_g_mpfr_fma },
    { "fms", Pympfr_fms, METH_VARARGS, doc_g_mpfr_fms },
    { "fmod", Pympfr_fmod, METH_VARARGS, doc_g_mpfr_fmod },
    { "frac", Pympfr_frac, METH_O, doc_g_mpfr_frac },
    { "f2q", Pympfr_f2q, METH_VARARGS, doc_g_mpfr_f2q },
    { "f_div", Pygmpy_f_div, METH_VARARGS, doc_gmpy_f_div },
    { "f_div_2exp", Pygmpy_f_div_2exp, METH_VARARGS, doc_gmpy_f_div_2exp },
    { "f_divmod", Pygmpy_f_divmod, METH_VARARGS, doc_gmpy_f_divmod },
    { "f_divmod_2exp", Pygmpy_f_divmod_2exp, METH_VARARGS, doc_gmpy_f_divmod_2exp },
    { "f_mod", Pygmpy_f_mod, METH_VARARGS, doc_gmpy_f_mod },
    { "f_mod_2exp", Pygmpy_f_mod_2exp, METH_VARARGS, doc_gmpy_f_mod_2exp },
    { "gamma", Pympfr_gamma, METH_O, doc_g_mpfr_gamma },
    { "gcd", Pygmpy_gcd, METH_VARARGS, doc_gcd },
    { "gcdext", Pygmpy_gcdext, METH_VARARGS, doc_gcdext },
    { "get_cache", Pygmpy_get_cache, METH_NOARGS, doc_get_cache },
    { "get_emax_max", Pympfr_get_emax_max, METH_NOARGS, doc_g_mpfr_get_emax_max },
    { "get_emin_min", Pympfr_get_emin_min, METH_NOARGS, doc_g_mpfr_get_emin_min },
    { "get_max_precision", Pympfr_get_max_precision, METH_NOARGS, doc_g_mpfr_get_max_precision },
    { "hamdist", Pympz_hamdist, METH_VARARGS, doc_hamdistg },
    { "hypot", Pympfr_hypot, METH_VARARGS, doc_g_mpfr_hypot },
    { "inf", Pympfr_set_inf, METH_O, doc_g_mpfr_set_inf },
    { "invert", Pygmpy_invert, METH_VARARGS, doc_invertg },
    { "is_even", Pympz_is_even, METH_O, doc_is_eveng },
    { "is_inf", Pympfr_is_inf, METH_O, doc_g_mpfr_is_inf },
    { "is_integer", Pympfr_is_integer, METH_O, doc_g_mpfr_is_integer },
    { "is_lessgreater", Pympfr_is_lessgreater, METH_VARARGS, doc_g_mpfr_is_lessgreater },
    { "is_nan", Pympfr_is_nan, METH_O, doc_g_mpfr_is_nan },
    { "is_number", Pympfr_is_number, METH_O, doc_g_mpfr_is_number },
    { "is_odd", Pympz_is_odd, METH_O, doc_is_oddg },
    { "is_power", Pympz_is_power, METH_O, doc_is_powerg },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primeg },
    { "is_regular", Pympfr_is_regular, METH_O, doc_g_mpfr_is_regular },
    { "is_square", Pympz_is_square, METH_O, doc_is_squareg },
    { "is_unordered", Pympfr_is_unordered, METH_VARARGS, doc_g_mpfr_is_unordered },
    { "is_zero", Pympfr_is_zero, METH_O, doc_g_mpfr_is_zero },
    { "jacobi", Pympz_jacobi, METH_VARARGS, doc_jacobig },
    { "jn", Pympfr_jn, METH_VARARGS, doc_g_mpfr_jn },
    { "j0", Pympfr_j0, METH_O, doc_g_mpfr_j0 },
    { "j1", Pympfr_j1, METH_O, doc_g_mpfr_j1 },
    { "kronecker", Pympz_kronecker, METH_VARARGS, doc_kroneckerm },
    { "lcm", Pygmpy_lcm, METH_VARARGS, doc_lcm },
    { "legendre", Pympz_legendre, METH_VARARGS, doc_legendreg },
    { "lgamma", Pympfr_lgamma, METH_O, doc_g_mpfr_lgamma },
    { "license", Pygmpy_get_license, METH_NOARGS, doc_license },
    { "li2", Pympfr_li2, METH_O, doc_g_mpfr_li2 },
    { "lngamma", Pympfr_lngamma, METH_O, doc_g_mpfr_lngamma },
    { "log", Pympfr_log, METH_O, doc_g_mpfr_log },
    { "log1p", Pympfr_log1p, METH_O, doc_g_mpfr_log1p },
    { "log10", Pympfr_log10, METH_O, doc_g_mpfr_log10 },
    { "log2", Pympfr_log2, METH_O, doc_g_mpfr_log2 },
    { "lucas", Pygmpy_lucas, METH_O, doc_lucas },
    { "lucas2", Pygmpy_lucas2, METH_O, doc_lucas2 },
    { "max", Pympfr_max, METH_VARARGS, doc_g_mpfr_max },
    { "min", Pympfr_min, METH_VARARGS, doc_g_mpfr_min },
    { "modf", Pympfr_modf, METH_O, doc_g_mpfr_modf },
    { "mpfr", Pygmpy_mpfr, METH_VARARGS, doc_mpfr },
    { "mp_version", Pygmpy_get_mp_version, METH_NOARGS, doc_mp_version },
    { "mp_limbsize", Pygmpy_get_mp_limbsize, METH_NOARGS, doc_mp_limbsize },
    { "mpc", (PyCFunction)Pygmpy_mpc, METH_VARARGS | METH_KEYWORDS, doc_g_mpc },
    { "mpc_version", Pygmpy_get_mpc_version, METH_NOARGS, doc_mpc_version },
    { "mpfr_version", Pygmpy_get_mpfr_version, METH_NOARGS, doc_mpfr_version },
    { "mpq", Pygmpy_mpq, METH_VARARGS, doc_mpq },
    { "mpz", Pygmpy_mpz, METH_VARARGS, doc_mpz },
    { "mul", Pympfr_mul, METH_VARARGS, doc_g_mpfr_mul },
    { "nan", Pympfr_set_nan, METH_NOARGS, doc_g_mpfr_set_nan },
    { "new_context", (PyCFunction)Pygmpy_new_context, METH_VARARGS | METH_KEYWORDS, doc_new_context },
    { "next_above", Pympfr_nextabove, METH_O, doc_g_mpfr_nextabove },
    { "next_below", Pympfr_nextbelow, METH_O, doc_g_mpfr_nextbelow },
    { "next_prime", Pympz_next_prime, METH_O, doc_next_primeg },
    { "next_toward", Pympfr_nexttoward, METH_VARARGS, doc_g_mpfr_nexttoward },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsg },
    { "numer", Pympq_numer, METH_VARARGS, doc_numerg },
    { "pack", Pygmpy_pack, METH_VARARGS, doc_gmpy_pack },
    { "popcount", Pympz_popcount, METH_O, doc_popcountg },
    { "pow", Pympfr_pow, METH_VARARGS, doc_g_mpfr_pow },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivg },
    { "rec_sqrt", Pympfr_rec_sqrt, METH_O, doc_g_mpfr_rec_sqrt },
    { "reldiff", Pympfr_reldiff, METH_VARARGS, doc_g_mpfr_reldiff },
    { "remainder", Pympfr_remainder, METH_VARARGS, doc_g_mpfr_remainder },
    { "remove", Pympz_remove, METH_VARARGS, doc_removeg },
    { "remquo", Pympfr_remquo, METH_VARARGS, doc_g_mpfr_remquo },
    { "rint", Pympfr_rint, METH_O, doc_g_mpfr_rint },
    { "rint_ceil", Pympfr_rint_ceil, METH_O, doc_g_mpfr_rint_ceil },
    { "rint_floor", Pympfr_rint_floor, METH_O, doc_g_mpfr_rint_floor },
    { "rint_round", Pympfr_rint_round, METH_O, doc_g_mpfr_rint_round },
    { "rint_trunc", Pympfr_rint_trunc, METH_O, doc_g_mpfr_rint_trunc },
    { "root", Pympany_root, METH_VARARGS, doc_g_mpany_root },
    { "rootrem", Pympz_rootrem, METH_VARARGS, doc_rootremg },
    { "round", Pympfr_round, METH_VARARGS, doc_g_mpfr_round },
    { "round2", Pympfr_round2, METH_O, doc_g_mpfr_round2 },
    { "sec", Pympfr_sec, METH_O, doc_g_mpfr_sec },
    { "sech", Pympfr_sech, METH_O, doc_g_mpfr_sech },
    { "set_cache", Pygmpy_set_cache, METH_VARARGS, doc_set_cache },
    { "set_context", Pygmpy_set_context, METH_O, doc_set_context },
    { "set_debug", Pygmpy_set_debug, METH_VARARGS, doc_set_debug },
    { "sign", Pympany_sign, METH_O, doc_g_mpany_sign },
    { "sin", Pympfr_sin, METH_O, doc_g_mpfr_sin },
    { "sinh", Pympfr_sinh, METH_O, doc_g_mpfr_sinh },
    { "sinh_cosh", Pympfr_sinh_cosh, METH_O, doc_g_mpfr_sinh_cosh },
    { "sin_cos", Pympfr_sin_cos, METH_O, doc_g_mpfr_sin_cos },
    { "square", Pympany_square, METH_O, doc_g_mpany_square },
    { "sqrt", Pympany_sqrt, METH_O, doc_g_mpany_sqrt },
    { "sqrtrem", Pympz_sqrtrem, METH_VARARGS, doc_sqrtremg },
    { "sub", Pympfr_sub, METH_VARARGS, doc_g_mpfr_sub },
    { "tan", Pympfr_tan, METH_O, doc_g_mpfr_tan },
    { "tanh", Pympfr_tanh, METH_O, doc_g_mpfr_tanh },
    { "t_div", Pygmpy_t_div, METH_VARARGS, doc_gmpy_t_div },
    { "t_div_2exp", Pygmpy_t_div_2exp, METH_VARARGS, doc_gmpy_t_div_2exp },
    { "t_divmod", Pygmpy_t_divmod, METH_VARARGS, doc_gmpy_t_divmod },
    { "t_divmod_2exp", Pygmpy_t_divmod_2exp, METH_VARARGS, doc_gmpy_t_divmod_2exp },
    { "t_mod", Pygmpy_t_mod, METH_VARARGS, doc_gmpy_t_mod },
    { "t_mod_2exp", Pygmpy_t_mod_2exp, METH_VARARGS, doc_gmpy_t_mod_2exp },
    { "trunc", Pympfr_trunc, METH_O, doc_g_mpfr_trunc },
    { "unpack", Pygmpy_unpack, METH_VARARGS, doc_gmpy_unpack },
    { "version", Pygmpy_get_version, METH_NOARGS, doc_version },
    { "xbit_mask", Pyxmpz_xbit_mask, METH_O, doc_xbit_maskg },
    { "xmpz", Pygmpy_xmpz, METH_VARARGS, doc_xmpz },
    { "yn", Pympfr_yn, METH_VARARGS, doc_g_mpfr_yn },
    { "y0", Pympfr_y0, METH_O, doc_g_mpfr_y0 },
    { "y1", Pympfr_y1, METH_O, doc_g_mpfr_y1 },
    { "zero", Pympfr_set_zero, METH_O, doc_g_mpfr_set_zero },
    { "zeta", Pympfr_zeta, METH_O, doc_g_mpfr_zeta },
    { "_mpmath_normalize", Pympz_mpmath_normalize, METH_VARARGS, doc_mpmath_normalizeg },
    { "_mpmath_create", Pympz_mpmath_create, METH_VARARGS, doc_mpmath_createg },
    { NULL, NULL, 1}
};

static PyMethodDef Pympz_methods [] =
{
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "bincoef", Pympz_bincoef, METH_VARARGS, doc_bincoefm },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combm },
    { "c_div", Pympz_c_div, METH_O, doc_mpz_c_div },
    { "c_div_2exp", Pympz_c_div_2exp, METH_O, doc_mpz_c_div_2exp },
    { "c_mod", Pympz_c_mod, METH_O, doc_mpz_c_mod },
    { "c_mod_2exp", Pympz_c_mod_2exp, METH_O, doc_mpz_c_mod_2exp },
    { "digits", Pympz_digits, METH_VARARGS, doc_mpz_digits },
    { "divexact", Pympz_divexact, METH_O, doc_divexactm },
    { "f_div", Pympz_f_div, METH_O, doc_mpz_f_div },
    { "f_div_2exp", Pympz_f_div_2exp, METH_O, doc_mpz_f_div_2exp },
    { "f_mod", Pympz_f_mod, METH_O, doc_mpz_f_mod },
    { "f_mod_2exp", Pympz_f_mod_2exp, METH_O, doc_mpz_f_mod_2exp },
    { "hamdist", Pympz_hamdist, METH_VARARGS, doc_hamdistm },
    { "invert", Pympz_invert, METH_O, doc_invertm },
    { "is_even", Pympz_is_even, METH_NOARGS, doc_is_evenm },
    { "is_odd", Pympz_is_odd, METH_NOARGS, doc_is_oddm },
    { "is_square", Pympz_is_square, METH_NOARGS, doc_is_squarem },
    { "is_power", Pympz_is_power, METH_NOARGS, doc_is_powerm },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primem },
    { "jacobi", Pympz_jacobi, METH_VARARGS, doc_jacobim },
    { "kronecker", Pympz_kronecker, METH_VARARGS, doc_kroneckerg },
    { "legendre", Pympz_legendre, METH_VARARGS, doc_legendrem },
    { "remove", Pympz_remove, METH_VARARGS, doc_removem },
    { "next_prime", Pympz_next_prime, METH_NOARGS, doc_next_primem },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsm },
    { "popcount", Pympz_popcount, METH_NOARGS, doc_popcountm },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivm },
    { "root", Pympz_root, METH_VARARGS, doc_mpz_root },
    { "rootrem", Pympz_rootrem, METH_VARARGS, doc_rootremm },
    { "sign", Pympz_sign, METH_NOARGS, doc_signm },
    { "sqrt", Pympz_sqrt, METH_NOARGS, doc_sqrtm },
    { "sqrtrem", Pympz_sqrtrem, METH_VARARGS, doc_sqrtremm },
    { "t_div", Pympz_t_div, METH_O, doc_mpz_t_div },
    { "t_div_2exp", Pympz_t_div_2exp, METH_O, doc_mpz_t_div_2exp },
    { "t_mod", Pympz_t_mod, METH_O, doc_mpz_t_mod },
    { "t_mod_2exp", Pympz_t_mod_2exp, METH_O, doc_mpz_t_mod_2exp },
    { NULL, NULL, 1 }
};

static PyMethodDef Pyxmpz_methods [] =
{
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "bincoef", Pympz_bincoef, METH_VARARGS, doc_bincoefm },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combm },
    { "c_div", Pympz_c_div, METH_O, doc_mpz_c_div },
    { "c_div_2exp", Pympz_c_div_2exp, METH_O, doc_mpz_c_div_2exp },
    { "c_mod", Pympz_c_mod, METH_O, doc_mpz_c_mod },
    { "c_mod_2exp", Pympz_c_mod_2exp, METH_O, doc_mpz_c_mod_2exp },
    { "digits", Pyxmpz_digits, METH_VARARGS, doc_mpz_digits },
    { "divexact", Pympz_divexact, METH_O, doc_divexactm },
    { "f_div", Pympz_f_div, METH_O, doc_mpz_f_div },
    { "f_div_2exp", Pympz_f_div_2exp, METH_O, doc_mpz_f_div_2exp },
    { "f_mod", Pympz_f_mod, METH_O, doc_mpz_f_mod },
    { "f_mod_2exp", Pympz_f_mod_2exp, METH_O, doc_mpz_f_mod_2exp },
    { "hamdist", Pympz_hamdist, METH_VARARGS, doc_hamdistm },
    { "invert", Pympz_invert, METH_O, doc_invertm },
    { "is_even", Pympz_is_even, METH_NOARGS, doc_is_evenm },
    { "is_odd", Pympz_is_odd, METH_NOARGS, doc_is_oddm },
    { "is_square", Pympz_is_square, METH_VARARGS, doc_is_squarem },
    { "is_power", Pympz_is_power, METH_VARARGS, doc_is_powerm },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primem },
    { "jacobi", Pympz_jacobi, METH_VARARGS, doc_jacobim },
    { "kronecker", Pympz_kronecker, METH_VARARGS, doc_kroneckerg },
    { "legendre", Pympz_legendre, METH_VARARGS, doc_legendrem },
    { "make_mpz", Pyxmpz_make_mpz, METH_NOARGS, doc_make_mpzm },
    { "next_prime", Pympz_next_prime, METH_NOARGS, doc_next_primem },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsm },
    { "popcount", Pympz_popcount, METH_NOARGS, doc_popcountm },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivm },
    { "remove", Pympz_remove, METH_VARARGS, doc_removem },
    { "root", Pympz_root, METH_VARARGS, doc_mpz_root },
    { "rootrem", Pympz_rootrem, METH_VARARGS, doc_rootremm },
    { "sign", Pympz_sign, METH_NOARGS, doc_signm },
    { "sqrt", Pympz_sqrt, METH_NOARGS, doc_sqrtm },
    { "sqrtrem", Pympz_sqrtrem, METH_VARARGS, doc_sqrtremm },
    { "t_div", Pympz_t_div, METH_O, doc_mpz_t_div },
    { "t_div_2exp", Pympz_t_div_2exp, METH_O, doc_mpz_t_div_2exp },
    { "t_mod", Pympz_t_mod, METH_O, doc_mpz_t_mod },
    { "t_mod_2exp", Pympz_t_mod_2exp, METH_O, doc_mpz_t_mod_2exp },
    { NULL, NULL, 1 }
};

static PyMethodDef Pympq_methods [] =
{
    { "sign", Pympq_sign, METH_NOARGS, doc_qsignm },
    { "numer", Pympq_numer, METH_VARARGS, doc_numerm },
    { "denom", Pympq_denom, METH_VARARGS, doc_denomm },
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "digits", Pympq_digits, METH_VARARGS, doc_qdigitsm },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pympfr_methods [] =
{
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "acos", Pympfr_acos, METH_NOARGS, doc_mpfr_acos },
    { "acosh", Pympfr_acosh, METH_NOARGS, doc_mpfr_acosh },
    { "add", Pympfr_add, METH_VARARGS, doc_mpfr_add },
    { "agm", Pympfr_agm, METH_VARARGS, doc_mpfr_agm },
    { "ai", Pympfr_ai, METH_NOARGS, doc_mpfr_ai },
    { "asin", Pympfr_asin, METH_NOARGS, doc_mpfr_asin },
    { "asinh", Pympfr_asinh, METH_NOARGS, doc_mpfr_asinh },
    { "atan", Pympfr_atan, METH_NOARGS, doc_mpfr_atan },
    { "atanh", Pympfr_atanh, METH_NOARGS, doc_mpfr_atanh },
    { "atan2", Pympfr_atan2, METH_VARARGS, doc_mpfr_atan2 },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "cbrt", Pympfr_cbrt, METH_NOARGS, doc_mpfr_cbrt },
    { "ceil", Pympfr_ceil, METH_NOARGS, doc_mpfr_ceil },
    { "check_range", Pympfr_check_range, METH_NOARGS, doc_mpfr_check_range },
    { "cos", Pympfr_cos, METH_NOARGS, doc_mpfr_cos },
    { "cosh", Pympfr_cosh, METH_NOARGS, doc_mpfr_cosh },
    { "cot", Pympfr_cot, METH_NOARGS, doc_mpfr_cot },
    { "coth", Pympfr_coth, METH_NOARGS, doc_mpfr_coth },
    { "csc", Pympfr_csc, METH_NOARGS, doc_mpfr_csc },
    { "csch", Pympfr_csch, METH_NOARGS, doc_mpfr_csch },
    { "digamma", Pympfr_digamma, METH_NOARGS, doc_mpfr_digamma },
    { "digits", Pympfr_digits, METH_VARARGS, doc_mpfr_digits },
    { "div", Pympfr_div, METH_VARARGS, doc_mpfr_div },
    { "eint", Pympfr_eint, METH_NOARGS, doc_mpfr_eint },
    { "erf", Pympfr_erf, METH_NOARGS, doc_mpfr_erf },
    { "erfc", Pympfr_erfc, METH_NOARGS, doc_mpfr_erfc },
    { "exp", Pympfr_exp, METH_NOARGS, doc_mpfr_exp },
    { "expm1", Pympfr_expm1, METH_NOARGS, doc_mpfr_expm1 },
    { "exp10", Pympfr_exp10, METH_NOARGS, doc_mpfr_exp10 },
    { "exp2", Pympfr_exp2, METH_NOARGS, doc_mpfr_exp2 },
    { "floor", Pympfr_floor, METH_NOARGS, doc_mpfr_floor },
    { "fmod", Pympfr_fmod, METH_VARARGS, doc_mpfr_fmod },
    { "frac", Pympfr_frac, METH_NOARGS, doc_mpfr_frac },
    { "f2q", Pympfr_f2q, METH_VARARGS, doc_mpfr_f2q },
    { "gamma", Pympfr_gamma, METH_NOARGS, doc_mpfr_gamma },
    { "hypot", Pympfr_hypot, METH_VARARGS, doc_mpfr_hypot },
    { "is_inf", Pympfr_is_inf, METH_NOARGS, doc_mpfr_is_inf },
    { "is_integer", Pympfr_is_integer, METH_NOARGS, doc_mpfr_is_integer },
    { "is_lessgreater", Pympfr_is_lessgreater, METH_VARARGS, doc_mpfr_is_lessgreater },
    { "is_nan", Pympfr_is_nan, METH_NOARGS, doc_mpfr_is_nan },
    { "is_number", Pympfr_is_number, METH_NOARGS, doc_mpfr_is_number },
    { "is_regular", Pympfr_is_regular, METH_NOARGS, doc_mpfr_is_regular },
    { "is_unordered", Pympfr_is_unordered, METH_VARARGS, doc_mpfr_is_unordered },
    { "is_zero", Pympfr_is_zero, METH_NOARGS, doc_mpfr_is_zero },
    { "jn", Pympfr_jn, METH_VARARGS, doc_mpfr_jn },
    { "j0", Pympfr_j0, METH_NOARGS, doc_mpfr_j0 },
    { "j1", Pympfr_j1, METH_NOARGS, doc_mpfr_j1 },
    { "lgamma", Pympfr_lgamma, METH_NOARGS, doc_mpfr_lgamma },
    { "li2", Pympfr_li2, METH_NOARGS, doc_mpfr_li2 },
    { "lngamma", Pympfr_lngamma, METH_NOARGS, doc_mpfr_lngamma },
    { "log", Pympfr_log, METH_NOARGS, doc_mpfr_log },
    { "log1p", Pympfr_log1p, METH_NOARGS, doc_mpfr_log1p },
    { "log10", Pympfr_log10, METH_NOARGS, doc_mpfr_log10 },
    { "log2", Pympfr_log2, METH_NOARGS, doc_mpfr_log2 },
    { "max", Pympfr_max, METH_VARARGS, doc_mpfr_max },
    { "min", Pympfr_min, METH_VARARGS, doc_mpfr_min },
    { "modf", Pympfr_modf, METH_NOARGS, doc_mpfr_modf },
    { "mul", Pympfr_mul, METH_VARARGS, doc_mpfr_mul },
    { "next_above", Pympfr_nextabove, METH_NOARGS, doc_mpfr_nextabove },
    { "next_below", Pympfr_nextbelow, METH_NOARGS, doc_mpfr_nextbelow },
    { "next_toward", Pympfr_nexttoward, METH_VARARGS, doc_mpfr_nexttoward },
    { "pow", Pympfr_pow, METH_VARARGS, doc_mpfr_pow },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivm },
    { "rec_sqrt", Pympfr_rec_sqrt, METH_NOARGS, doc_mpfr_rec_sqrt },
    { "reldiff", Pympfr_reldiff, METH_VARARGS, doc_mpfr_reldiff },
    { "remainder", Pympfr_remainder, METH_VARARGS, doc_mpfr_remainder },
    { "remquo", Pympfr_remquo, METH_VARARGS, doc_mpfr_remquo },
    { "rint", Pympfr_rint, METH_NOARGS, doc_mpfr_rint },
    { "rint_ceil", Pympfr_rint_ceil, METH_NOARGS, doc_mpfr_rint_ceil },
    { "rint_floor", Pympfr_rint_floor, METH_NOARGS, doc_mpfr_rint_floor },
    { "rint_round", Pympfr_rint_round, METH_NOARGS, doc_mpfr_rint_round },
    { "rint_trunc", Pympfr_rint_trunc, METH_NOARGS, doc_mpfr_rint_trunc },
    { "root", Pympfr_root, METH_VARARGS, doc_mpfr_root },
    { "round", Pympfr_round, METH_VARARGS, doc_mpfr_round },
    { "round2", Pympfr_round2, METH_NOARGS, doc_mpfr_round2 },
    { "sec", Pympfr_sec, METH_NOARGS, doc_mpfr_sec },
    { "sech", Pympfr_sech, METH_NOARGS, doc_mpfr_sech },
    { "sign", Pympfr_sign, METH_NOARGS, doc_mpfr_sign },
    { "sin", Pympfr_sin, METH_NOARGS, doc_mpfr_sin },
    { "sin_cos", Pympfr_sin_cos, METH_NOARGS, doc_mpfr_sin_cos },
    { "sinh", Pympfr_sinh, METH_NOARGS, doc_mpfr_sinh },
    { "sinh_cosh", Pympfr_sinh_cosh, METH_NOARGS, doc_mpfr_sinh_cosh },
    { "square", Pympfr_sqr, METH_NOARGS, doc_mpfr_sqr },
    { "sqrt", Pympfr_sqrt, METH_NOARGS, doc_mpfr_sqrt },
    { "sub", Pympfr_sub, METH_VARARGS, doc_mpfr_sub },
    { "tan", Pympfr_tan, METH_NOARGS, doc_mpfr_tan },
    { "tanh", Pympfr_tanh, METH_NOARGS, doc_mpfr_tanh },
    { "trunc", Pympfr_trunc, METH_NOARGS, doc_mpfr_trunc },
    { "yn", Pympfr_yn, METH_VARARGS, doc_mpfr_yn },
    { "y0", Pympfr_y0, METH_NOARGS, doc_mpfr_y0 },
    { "y1", Pympfr_y1, METH_NOARGS, doc_mpfr_y1 },
    { "zeta", Pympfr_zeta, METH_NOARGS, doc_mpfr_zeta },
    { NULL, NULL, 1 }
};

static PyMethodDef Pympc_methods[] =
{
    { NULL, NULL, 1 }
};

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
    Pympfr_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympfr_getseters,                        /* tp_getset        */
};

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
        0,                                  /* tp_repr          */
    &mpc_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) 0,                           /* tp_str           */
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
        0,                                  /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympc_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
        0,                                  /* tp_getset        */
};


static void *
gmpy_allocate(size_t size)
{
    void *res;
    size_t usize = size;

    if (usize < GMPY_ALLOC_MIN)
        usize = GMPY_ALLOC_MIN;
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mp_allocate( %llu->%llu )\n",
            (unsigned long long)size, (unsigned long long)usize);
#endif
    if (!(res = PyMem_Malloc(usize))) {
        fprintf(stderr, "mp_allocate( %llu->%llu )\n",
            (unsigned long long)size, (unsigned long long)usize);
        Py_FatalError("mp_allocate failure");
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mp_allocate( %llu->%llu ) ->%8p\n",
            (unsigned long long)size, (unsigned long long)usize, res);
#endif
    return res;
} /* mp_allocate() */


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
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr,
            "mp_reallocate: old address %8p, old size %llu(%llu), new %llu(%llu)\n",
            ptr, (unsigned long long)old_size, (unsigned long long)uold,
            (unsigned long long)new_size, (unsigned long long)unew);
#endif
    if (uold==unew) {
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "mp_reallocate: avoided realloc for %llu\n",
                (unsigned long long)unew);
#endif
        return ptr;
    }

    if (!(res = PyMem_Realloc(ptr, unew))) {
        fprintf(stderr,
            "mp_reallocate: old address %8p, old size %llu(%llu), new %llu(%llu)\n",
            ptr, (unsigned long long)old_size, (unsigned long long)uold,
            (unsigned long long)new_size, (unsigned long long)unew);
        Py_FatalError("mp_reallocate failure");
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mp_reallocate: newob address %8p, newob size %llu(%llu)\n",
        res, (unsigned long long)new_size, (unsigned long long)unew);
#endif
    return res;
} /* mp_reallocate() */

static void
gmpy_free( void *ptr, size_t size)
{
    size_t usize=size;
    if (usize < GMPY_ALLOC_MIN)
        usize = GMPY_ALLOC_MIN;
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "mp_free      : old address %8p, old size %llu(%llu)\n",
            ptr, (unsigned long long)size, (unsigned long long)usize);
#endif
    PyMem_Free(ptr);
} /* mp_free() */

static void
_PyInitGMP(void)
{
    mp_set_memory_functions(gmpy_allocate, gmpy_reallocate, gmpy_free);
    set_zcache();
    set_pympzcache();
    set_pympqcache();
    set_pympfrcache();
    set_pyxmpzcache();
    context = GMPyContext_new();
    GMPyExc_DivZero = PyErr_NewException("gmpy2.DivisionByZeroError", NULL, NULL);
    GMPyExc_Inexact = PyErr_NewException("gmpy2.InexactError", NULL, NULL);
    GMPyExc_Invalid = PyErr_NewException("gmpy2.InvalidOperationError", NULL, NULL);
    GMPyExc_Overflow = PyErr_NewException("gmpy2.OverflowError", NULL, NULL);
    GMPyExc_Underflow = PyErr_NewException("gmpy2.UnderflowError", NULL, NULL);
    GMPyExc_Erange = PyErr_NewException("gmpy2.RangeError", NULL, NULL);
}

static char _gmpy_docs[] = "\
gmpy2 2.0.0a1 - General Multiprecision arithmetic for Python:\n\
exposes functionality from the GMP or MPIR library to Python 2.6\n\
and later.\n\
\n\
Allows creation of multiprecision integer (mpz), mutable integers\n\
(xmpz), float (mpfr), and rational (mpq) numbers, conversion between\n\
them and to/from Python numbers/strings, arithmetic, bitwise, and\n\
some other higher-level mathematical operations.\n\
\n\
mpz has comparable functionality to Python's builtin longs, but\n\
can be faster for some operations (particularly multiplication\n\
and raising-to-power) and has many further useful and speedy\n\
functions (prime testing and generation, factorial, fibonacci,\n\
binary-coefficients, gcd, lcm, square and other roots, ...).\n\
\n\
mpfr and mpq only offer basic arithmetic abilities, but they\n\
do add the ability to have floating-point numbers ensuring at\n\
least a predefined number of bits' worth of precision (and with\n\
potentially-huge or extremely-tiny magnitudes), as well as\n\
unlimited-precision rationals, with reasonably-fast operations,\n\
which are not built-in features of Python.\n\
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
    char *do_debug = getenv("GMPY_DEBUG");
    if (PyType_Ready(&Pympz_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pympq_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pympfr_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pyxmpz_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pympc_Type) < 0)
        INITERROR;
    if (PyType_Ready(&GMPyContext_Type) < 0)
        INITERROR;

    if (do_debug)
        sscanf(do_debug, "%d", &global.debug);
#ifdef DEBUG
    if (global.debug)
        fputs( "initgmpy2() called...\n", stderr );
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

#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "gmpy_module at %p\n", gmpy_module);
#endif
    /* Add support for pickling. */
#ifdef PY3
    copy_reg_module = PyImport_ImportModule("copyreg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def mpz_reducer(an_mpz): return (gmpy2.mpz, (an_mpz.binary(), 256))\n"
            "def mpq_reducer(an_mpq): return (gmpy2.mpq, (an_mpq.binary(), 256))\n"
            "def mpfr_reducer(an_mpfr): return (gmpy2.mpfr, (an_mpfr.binary(), 0, 256))\n"
            "copyreg.pickle(type(gmpy2.mpz(0)), mpz_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpq(0)), mpq_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpfr(0)), mpfr_reducer)\n"
        ;
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "gmpy_module imported copyreg OK\n");
#endif
        PyDict_SetItemString(namespace, "copyreg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (result) {
#ifdef DEBUG
            if (global.debug)
                fprintf(stderr, "gmpy_module enable pickle OK\n");
#endif
        }
        else {
#ifdef DEBUG
            if (global.debug)
                fprintf(stderr, "gmpy_module could not enable pickle\n");
#endif
            PyErr_Clear();
        }
        Py_DECREF(namespace);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "gmpy_module could not import copyreg\n");
#endif
    }
#else
    copy_reg_module = PyImport_ImportModule("copy_reg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def mpz_reducer(an_mpz): return (gmpy2.mpz, (an_mpz.binary(), 256))\n"
            "def mpq_reducer(an_mpq): return (gmpy2.mpq, (an_mpq.binary(), 256))\n"
            "def mpfr_reducer(an_mpfr): return (gmpy2.mpfr, (an_mpfr.binary(), 0, 256))\n"
            "copy_reg.pickle(type(gmpy2.mpz(0)), mpz_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpq(0)), mpq_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpfr(0)), mpfr_reducer)\n"
        ;
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "gmpy_module imported copy_reg OK\n");
#endif
        PyDict_SetItemString(namespace, "copy_reg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (result) {
#ifdef DEBUG
            if (global.debug)
                fprintf(stderr, "gmpy_module enable pickle OK\n");
#endif
        }
        else {
#ifdef DEBUG
            if (global.debug)
                fprintf(stderr, "gmpy_module could not enable pickle\n");
#endif
            PyErr_Clear();
        }
        Py_DECREF(namespace);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "gmpy_module could not import copy_reg\n");
#endif
    }
#endif

#ifdef PY3
    return gmpy_module;
#endif
}