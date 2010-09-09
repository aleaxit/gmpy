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
 *   added options.fcoform for optional use of intermediate string in
 *       float2mpf (used for any float->mpf conversion)
 *   added set_fcoform function for options.fcoform access
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
 *      conditioned by options.debug, & a couple of very obscure cases)
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

/* Include fast mpz to/from PyLong conversion from sage. */
#include "mpz_pylong.c"

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
#endif

char gmpy_version[] = "2.0.0a1";

char _gmpy_cvs[] = "$Id$";

/*
 * global data declarations
 */

static PyObject *gmpy_module = NULL;

static struct gmpy_options {
    int debug;               /* != 0 if debug messages desired on stderr */
    mpfr_prec_t precision;   /* current precision in bits */
    mpfr_rnd_t rounding;     /* current rounding mode */
    int cache_size;          /* size of cache, for all caches */
    int cache_obsize;        /* maximum size of the objects that are cached */
    int prefer_mutable;      /* != 0 if mixed mpz/xmpz operation results in xmpz */
    PyObject* fcoform;       /* if non-NULL, format for float->mpf (via string) */
} options = { 0, DBL_MANT_DIG, MPFR_RNDN, 100, 128, 0, 0 };

/* forward declarations of type-objects and method-arrays for them */
#ifdef _MSC_VER
PyMethodDef Pympz_methods [];
PyMethodDef Pympq_methods [];
PyMethodDef Pympf_methods [];
PyMethodDef Pyxmpz_methods [];
#else
static PyMethodDef Pympz_methods [];
static PyMethodDef Pympq_methods [];
static PyMethodDef Pympf_methods [];
static PyMethodDef Pyxmpz_methods [];
#endif

/* gmpy2 caches objects so they can be reused quickly without involving a new
 * memory allocation or object construction. There are two different types of
 * object caches used in gmpy2.
 *
 * "zcache" and "qcache" are used to cache mpz_t and mpq_t objects. The cache
 * is accessed via the functions mpz_inoc/mpz_cloc and mpq_inoc/mpq_cloc. The
 * functions set_zcache and set_qcache are used to change the size of the
 * array used to store the cached objects.
 *
 * "pympzcache" and "pympqcache" are used to cache Pympz and Pympq objects.
 * The cache is accessed via Pympz_new/Pympz_dealloc and Pympq_new/
 * Pympq_dealloc. The functions set_pympzcache and set_pympqcache are used
 * to change the size of the array used to the store the cached objects.
 */

static mpz_t* zcache;
static int in_zcache;

static void
set_zcache(void)
{
    TRACE("Entering set_zcache\n");
    if (in_zcache > options.cache_size) {
        int i;
        for(i = options.cache_size; i < in_zcache; ++i)
            mpz_clear(zcache[i]);
        in_zcache = options.cache_size;
    }
    zcache = PyMem_Realloc(zcache, sizeof(mpz_t) * options.cache_size);
}

static mpq_t* qcache;
static int in_qcache;

static void
set_qcache(void)
{
    TRACE("Entering set_qcache\n");
    if (in_qcache > options.cache_size) {
        int i;
        for (i = options.cache_size; i < in_qcache; ++i)
            mpq_clear(qcache[i]);
        in_qcache = options.cache_size;
    }
    qcache = PyMem_Realloc(qcache, sizeof(mpq_t) * options.cache_size);
}

static void
mpz_inoc(mpz_t newo)
{
    if (in_zcache) {
        if (options.debug)
            fprintf(stderr, "Getting %d from zcache\n", in_zcache);
        newo[0] = (zcache[--in_zcache])[0];
    }
    else {
        TRACE("Initing new not in zcache\n");
        mpz_init(newo);
    }
}

static void
mpz_cloc(mpz_t oldo)
{
    if (in_zcache<options.cache_size && oldo->_mp_alloc <= options.cache_obsize) {
        (zcache[in_zcache++])[0] = oldo[0];
        if (options.debug)
            fprintf(stderr, "Stashed %d to zcache\n", in_zcache);
    }
    else {
        if (options.debug)
            fprintf(stderr, "Not placing in full zcache(%d/%d)\n",
                    in_zcache, options.cache_size);
        mpz_clear(oldo);
    }
}

static void
mpq_inoc(mpq_t newo)
{
    if (in_qcache) {
        if (options.debug)
            fprintf(stderr, "Getting %d from qcache\n", in_qcache);
        newo[0] = (qcache[--in_qcache])[0];
    }
    else {
        TRACE("Initing new not in qcache\n");
        mpq_init(newo);
    }
}

static void
mpq_cloc(mpq_t oldo)
{
    if (in_qcache<options.cache_size &&
            mpq_numref(oldo)->_mp_alloc <= options.cache_obsize &&
            mpq_denref(oldo)->_mp_alloc <= options.cache_obsize) {
        (qcache[in_qcache++])[0] = oldo[0];
        if (options.debug)
            fprintf(stderr, "Stashed %d to qcache\n", in_qcache);
    }
    else {
        if (options.debug)
            fprintf(stderr, "Not placing in full qcache(%d/%d)\n",
                    in_qcache, options.cache_size);
        mpq_clear(oldo);
    }
}

/* Cache Pympz objects directly */

static PympzObject **pympzcache;
static int in_pympzcache;

static void
set_pympzcache(void)
{
    TRACE("Entering set_pympzcache\n");
    if (in_pympzcache > options.cache_size) {
        int i;
        for (i = options.cache_size; i < in_pympzcache; ++i) {
            mpz_cloc(pympzcache[i]->z);
            PyObject_Del(pympzcache[i]);
        }
        in_pympzcache = options.cache_size;
    }
    pympzcache = PyMem_Realloc(pympzcache, sizeof(PympzObject)*options.cache_size);
}

/* Cache Pyxmpz objects directly */

static PyxmpzObject **pyxmpzcache;
static int in_pyxmpzcache;

static void
set_pyxmpzcache(void)
{
    TRACE("Entering set_pyxmpzcache\n");
    if (in_pyxmpzcache > options.cache_size) {
        int i;
        for (i = options.cache_size; i < in_pyxmpzcache; ++i) {
            mpz_cloc(pyxmpzcache[i]->z);
            PyObject_Del(pyxmpzcache[i]);
        }
        in_pyxmpzcache = options.cache_size;
    }
    pyxmpzcache = PyMem_Realloc(pyxmpzcache, sizeof(PyxmpzObject)*options.cache_size);
}

/* Cache Pympq objects directly */

static PympqObject **pympqcache;
static int in_pympqcache;

static void
set_pympqcache(void)
{
    TRACE("Entering set_pympqcache\n");
    if (in_pympqcache > options.cache_size) {
        int i;
        for (i = options.cache_size; i < in_pympqcache; ++i) {
            mpq_cloc(pympqcache[i]->q);
            PyObject_Del(pympqcache[i]);
        }
        in_pympqcache = options.cache_size;
    }
    pympqcache = PyMem_Realloc(pympqcache, sizeof(PympqObject)*options.cache_size);
}

/* generation of new, uninitialized objects; deallocations */
static PympzObject *
Pympz_new(void)
{
    PympzObject * self;

    TRACE("Entering Pympz_new\n");
    if (in_pympzcache) {
        TRACE("Pympz_new is reusing an old object\n");
        self = (pympzcache[--in_pympzcache]);
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        TRACE("Pympz_new is creating a new object\n");
        if (!(self = PyObject_New(PympzObject, &Pympz_Type)))
            return NULL;
        mpz_inoc(self->z);
    }
    self->hash_cache = -1;
    return self;
}

static PyxmpzObject *
Pyxmpz_new(void)
{
    PyxmpzObject *self;

    TRACE("Entering Pyxmpz_new\n");
    if (in_pyxmpzcache) {
        TRACE("Pyxmpz_new is reusing an old object\n");
        self = (pyxmpzcache[--in_pyxmpzcache]);
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        TRACE("Pyxmpz_new is creating a new object\n");
        if (!(self = PyObject_New(PyxmpzObject, &Pyxmpz_Type)))
            return NULL;
        mpz_inoc(self->z);
    }
    return self;
}

static PympqObject *
Pympq_new(void)
{
    PympqObject *self;

    TRACE("Entering Pympq_new\n");
    if (in_pympqcache) {
        TRACE("Pympq_new is reusing an old object\n");
        self = (pympqcache[--in_pympqcache]);
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        TRACE("Pympq_new is creating a new object\n");
        if (!(self = PyObject_New(PympqObject, &Pympq_Type)))
            return NULL;
        mpq_inoc(self->q);
    }
    self->hash_cache = -1;
    return self;
}

static PympfObject *
Pympf_new(mpfr_prec_t bits)
{
    PympfObject *self;

    TRACE("Entering Pympf_new\n");
    if (!(self = PyObject_New(PympfObject, &Pympf_Type)))
        return NULL;
    if (bits == 0)
        bits = options.precision;
    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    mpfr_init2(self->f, bits);
    self->hash_cache = -1;
    return self;
}

static void
Pympz_dealloc(PympzObject *self)
{
    TRACE("Pympz_dealloc\n");
    if (in_pympzcache < options.cache_size &&
        self->z->_mp_alloc <= options.cache_obsize) {
        (pympzcache[in_pympzcache++]) = self;
    }
    else {
        mpz_cloc(self->z);
        PyObject_Del(self);
    }
}

static void
Pyxmpz_dealloc(PyxmpzObject *self)
{
    TRACE("Pyxmpz_dealloc\n");
    if (in_pyxmpzcache < options.cache_size &&
        self->z->_mp_alloc <= options.cache_obsize) {
        (pyxmpzcache[in_pyxmpzcache++]) = self;
    }
    else {
        mpz_cloc(self->z);
        PyObject_Del(self);
    }
}

static void
Pympq_dealloc(PympqObject *self)
{
    TRACE("Pympq_dealloc\n");
    if (in_pympqcache<options.cache_size &&
        mpq_numref(self->q)->_mp_alloc <= options.cache_obsize &&
        mpq_denref(self->q)->_mp_alloc <= options.cache_obsize) {
        (pympqcache[in_pympqcache++]) = self;
    }
    else {
        mpq_cloc(self->q);
        PyObject_Del(self);
    }
}

static void
Pympf_dealloc(PympfObject *self)
{
    TRACE("Pympf_dealloc\n");
    mpfr_clear(self->f);
    PyObject_Del(self);
}

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
  Make a copy of an mpf object. If bits is 0, the new object will have
  the same precision as the original object. If the requested precision
  is less than the precision of the original object, the new object
  will be rounded to requested precision using the current rounding mode.
*/

static PympfObject *
Pympf2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob;

    assert(Pympf_Check(self));
    if (bits == 0)
        bits = mpfr_get_prec(Pympf_AS_MPF(self));
    if ((newob = Pympf_new(bits)))
        mpfr_set(newob->f, Pympf_AS_MPF(self), options.rounding);
    return newob;
}

#ifdef PY2
static PympzObject *
PyInt2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympz_new()))
        mpz_set_si(newob->z, PyInt_AsLong(self));
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

static PympfObject *
PyInt2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympf_new(bits)))
        mpfr_set_si(newob->f, PyInt_AsLong(self), options.rounding);
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
            VALUE_ERROR("gmpy2.mpz does not handle nan");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpz does not handle infinity");
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
            VALUE_ERROR("gmpy2.xmpz does not handle nan");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.xmpz does not handle infinity");
            return NULL;
        }
        mpz_set_d(newob->z, d);
    }
    return newob;
}

/* forward...: */
static PyObject *f2q_internal(PympfObject* self, PympfObject* err,
        unsigned int bits, int mayz);
static PyObject* Pympf_f2q(PyObject *self, PyObject *args);
static PympfObject* anynum2Pympf(PyObject* obj, mpfr_prec_t bits);

static PympqObject *
PyFloat2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pympq_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpq does not handle nan");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpq does not handle infinity");
            return NULL;
        }
        mpq_set_d(newob->q, d);
    }
    return newob;
}

/* forward */
static PympfObject *PyStr2Pympf(PyObject *s, long base, Py_ssize_t bits);

static PympfObject *
PyFloat2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob = 0;

    assert(PyFloat_Check(self));
    if (!bits)
        bits = DBL_MANT_DIG;
    if (options.debug)
        fprintf(stderr, "PyFloat2Pympf(%p,%ld)\n", self, (long) bits);

    if (options.fcoform) {
        /* 2-step float->mpf conversion process: first, get a
         * Python string by formatting the Python float; then,
         * use str2mpf to build the mpf from the string.
         */
        PyObject *tuple, *s, *fmt;

        if (!(tuple = Py_BuildValue("(O)", self)))
            return NULL;
        if (!(fmt = (PyObject*)Py2or3String_FromString("%s"))) {
            Py_DECREF(tuple);
            return NULL;
        }
        s = Py2or3String_Format(fmt, tuple);
        Py_DECREF(tuple);
        Py_DECREF(fmt);
        if (!s)
            return NULL;
        newob = PyStr2Pympf(s, 10, bits);
        Py_DECREF(s);
        if (!newob)
            return NULL;
    }
    else { /* direct float->mpf conversion, faster but rougher */
        if ((newob = Pympf_new(bits)))
            mpfr_set_d(newob->f, PyFloat_AsDouble(self), bits);
    }
    return newob;
}

static PympfObject *
Pympz2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pympf_new(bits)))
        mpfr_set_z(newob->f, Pympz_AS_MPZ(self), options.rounding);
    return newob;
}

static PympfObject *
Pyxmpz2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob;

    assert(Pyxmpz_Check(obj));
    if ((newob = Pympf_new(bits)))
        mpfr_set_z(newob->f, Pympz_AS_MPZ(self), options.rounding);
    return newob;
}

static PympzObject *
Pympf2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pympf_Check(self));
    if ((newob = Pympz_new())) {
        if (mpfr_nan_p(Pympf_AS_MPF(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpz does not handle nan");
            return NULL;
        }
        if (mpfr_inf_p(Pympf_AS_MPF(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.mpz does not handle infinity");
            return NULL;
        }
        mpfr_get_z(newob->z, Pympf_AS_MPF(self), options.rounding);
    }
    return newob;
}

static PyxmpzObject *
Pympf2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pympf_Check(self));
    if ((newob = Pyxmpz_new())) {
        if (mpfr_nan_p(Pympf_AS_MPF(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.xmpz does not handle nan");
            return NULL;
        }
        if (mpfr_inf_p(Pympf_AS_MPF(self))) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("gmpy2.xmpz does not handle infinity");
            return NULL;
        }
        mpfr_get_z(newob->z, Pympf_AS_MPF(self), options.rounding);
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
Pympf2Pympq(PyObject *self)
{
    return (PympqObject*) Pympf_f2q(self, 0);
}

static PympfObject *
Pympq2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob;

    assert(Pympq_Check(self));
    if (!(newob = Pympf_new(bits)))
        return NULL;
    mpfr_set_q(newob->f, Pympq_AS_MPQ(self), options.rounding);
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
 * long->mpf delegates via long->mpz->mpf to avoid duplicating
 * the above-seen dependencies; ditto long->mpq
 */
static PympfObject *
PyLong2Pympf(PyObject *self, mpfr_prec_t bits)
{
    PympfObject *newob;
    PyObject *temp = (PyObject*)PyLong2Pympz(self);

    if (!temp)
        return NULL;
    newob = Pympz2Pympf(temp, bits);
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
 * of bytes) and 'true' from-string (bases 2 to 36; bases 8 and 16 are
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
 * of bytes) and 'true' from-string (bases 2 to 36; bases 8 and 16 are
 * special -- decorations of leading 0/0x are allowed (but not required).
 * For 'true-bases' 2..36 there is a '/' separator between numerator and
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
                    PympfObject* temp = PyStr2Pympf(stringarg, base, 4*len);
                    if (temp) {
                        newob = Pympf2Pympq((PyObject*)temp);
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
 * mpf conversion from string includes from-binary (base-256, format is
 * explained later) and 'true' from-string (bases 2 to 36), where exponent
 * if any is denoted by 'e' if base<=10, else by '@', and is always decimal.
 */
static PympfObject *
PyStr2Pympf(PyObject *s, long base, mpfr_prec_t bits)
{
    PympfObject *newob;
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
    if ((!bits) && (base ==256)) {
        prec = 8 * (len - 5);
        if ((len >= 5) && (cp[0] & 8)) {
            bits = 0;
            for(i=4; i>0; --i) {
                bits = (bits<<8) | cp[i];
            }
        }
    }

    if (!(newob = Pympf_new(bits))) {
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

        /* mpf zero has a very compact (1-byte) binary encoding!-) */
        if (resuzero) {
            mpfr_set_ui(newob->f, bits, options.rounding);
            return newob;
        }

        /* all other numbers are 6+ bytes: codebyte, 4-byte exp, 1+
         * bytes for the mantissa; check this string is 6+ bytes
         */
        if (len < 6 + precilen) {
            VALUE_ERROR("string too short to be a gmpy2.mpf binary encoding");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }
        /* reconstruct exponent */
        for (i = 4 + precilen; i > precilen; --i) {
            expomag = (expomag<<8) | cp[i];
        }

        /* reconstruct 'mantissa' (significand) */
        mpfr_set_si(newob->f, bits, options.rounding);
        mpfr_init2(digit, bits);
        for (i = 5 + precilen; i<len; i++) {
            mpfr_set_ui(digit, cp[i], options.rounding);
            mpfr_div_2ui(digit, digit, (unsigned long)((i-4-precilen) * 8),
                         options.rounding);
            mpfr_add(newob->f, newob->f, digit, options.rounding);
        }
        mpfr_clear(digit);
        /* apply exponent, with its appropriate sign */
        if (exposign)
            mpfr_div_2ui(newob->f, newob->f, 8*expomag, options.rounding);
        else
            mpfr_mul_2ui(newob->f, newob->f, 8*expomag, options.rounding);
        /* apply significand-sign (sign of the overall number) */
        if (resusign)
            mpfr_neg(newob->f, newob->f, options.rounding);
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
        /* delegate the rest to GMP */
        if (-1 == mpfr_set_str(newob->f, (char*)cp, base, options.rounding)) {
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
 * mpf->long delegates via mpf->mpz->long to avoid duplicating
 * the above-seen thorny dependencies; ditto mpq->long
 */
static PyObject *
Pympf2PyLong(PympfObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympf2Pympz((PyObject*)self);

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
 * mpf->int delegates via mpf->mpz->int for convenience; ditto mpq->int
 */

#ifdef PY2
static PyObject *
Pympf2PyInt(PympfObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympf2Pympz((PyObject*)self);

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
Pympf2PyFloat(PympfObject *self)
{
    double res = mpfr_get_d(self->f, options.rounding);

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
 * helper functions for mpf->binary conversion
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
 * build binary representation of mpf (see format description above)
 */

static PyObject *
Pympf2binary(PympfObject *self)
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
        mpfr_neg(self->f, self->f, options.rounding);
    }
    else {
        codebyte = 0;
    }

    /* get buffer of base-16 digits */
    buffer  = mpfr_get_str(0, &the_exp, 16, 0, self->f, options.rounding);
    /* no need to worry about null-buffer as x->f==0.0 was
     * already handled above (see first test on 'sign').
     */
    /* restore correct sign to x->f if it was changed! */
    if (codebyte) {
        mpfr_neg(self->f, self->f, options.rounding);
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
 * format mpz into any base (2 to 36), optionally with
 * a "gmpy2.mpz(...)" tag around it so it can be recovered
 * through a Python eval of the resulting string
 */
static char* ztag = "mpz(";
static PyObject *
mpz_ascii(mpz_t z, int base, int with_tag)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if ((base != 0) && ((base < 2) || (base > 36))) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 36");
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
 * format xmpz into any base (2 to 36), optionally with
 * a "gmpy2.xmpz(...)" tag around it so it can be recovered
 * through a Python eval of the resulting string
 */
static char* xztag = "xmpz(";
static PyObject *
xmpz_ascii(mpz_t z, int base, int with_tag)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if ((base != 0) && ((base < 2) || (base > 36))) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 36");
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

static char* qtag = "gmpy2.mpq(";
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

    if (!qden_1(self->q)) {
        denstr = mpz_ascii(mpq_denref(self->q), base, 0);
        if (!denstr) {
            Py_DECREF(numstr);
            return NULL;
        }
    }

    if (with_tag) {
        result = PyBytes_FromString(qtag);
        if (!result) {
            Py_XDECREF(denstr);
            return NULL;
        }
        PyBytes_ConcatAndDel(&result, numstr);
#ifdef PY2
        if (!mpz_fits_slong_p(mpq_numref(self->q))) {
            temp = PyBytes_FromString("L");
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
        PyBytes_ConcatAndDel(&result, temp);
        if (!result) {
            Py_DECREF(denstr);
            return NULL;
        }
        PyBytes_ConcatAndDel(&result, denstr);
#ifdef PY2
        if (with_tag && !mpz_fits_slong_p(mpq_denref(self->q))) {
            temp = PyBytes_FromString("L");
            PyBytes_ConcatAndDel(&result, temp);
            if (!result) {
                return NULL;
            }
        }
#endif
    }
    if (with_tag && result) {
        temp = PyBytes_FromString(")");
        PyBytes_ConcatAndDel(&result, temp);
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
static char ftag[]="mpf('";
/*
 * format mpf into any base (2 to 36), optionally with
 * a "gmpy2.mpf('...')" tag around it so it can be recovered
 * through a Python eval of the resulting string.
 * digits: number of digits to ask GMP for (0=all of
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
 *     OP_TAG (1): add the gmpy2.mpf('...') tag
 *     OP_RAW (2): ignore minexfi/maxexfi/OP_TAG
 *         and return a 3-element tuple digits/exponent/rprec
 *         (as GMP gives them) for Python formatting;
 *         'digits' may include a '-' sign, but no decimal
 *         point, nor tag, nor any exponent-indicator.
 *     other bits are currently ignored
 */
static PyObject *
Pympf_ascii(PympfObject *self, int base, int digits,
    int minexfi, int maxexfi, int optionflags)
{
    PyObject *res;
#ifdef PY3
    PyObject *temp;
#endif
    char *buffer;
    mp_exp_t the_exp;

    /* check arguments are valid */
    assert(Pympf_Check((PyObject*)self));
    if (!( (base==0) || ((base >= 2) && (base <= 36)))) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 36");
        return NULL;
    }
    if (digits < 0) {
        VALUE_ERROR("digits must be >= 0");
        return NULL;
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, options.rounding);
    if (!*buffer) {
        /* need to use malloc here for uniformity with mpfr_get_str */
        PyMem_Free(buffer);
        buffer = PyMem_Malloc(2);
        strcpy(buffer, "0");
        the_exp = 1;
    }

    if (optionflags & OP_RAW) {
        res = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
        PyMem_Free(buffer);
        return res;
    }
    else {
        /* insert formatting elements (decimal-point, leading or
         * trailing 0's, other indication of exponent...)
         */
        size_t buflen = strlen(buffer);
        /* account for the decimal point that is always inserted */
        size_t size = buflen+1;
        char expobuf[24];
        char auprebuf[24];
        int isfp=1;   /* flag: fixed-point format (FP)? */
        int isnegative=0;
        if (buffer[0]==0x2d)
            isnegative = 1;

        /* compute size of needed Python string */
        if (optionflags & OP_TAG) {
            size += strlen(ftag) + 2;
            if (mpfr_get_prec(self->f) != DBL_MANT_DIG) {
                sprintf(auprebuf, ",%ld", (long)mpfr_get_prec(self->f));
                size += strlen(auprebuf);
            }
        }
        if (the_exp < minexfi || the_exp > maxexfi) { /* exponential format */
            /* add exponent-length + 1 for '@' or 'e' marker */
            sprintf(expobuf, "%ld", the_exp - 1);
            size += strlen(expobuf) + 1;
            isfp = 0;
        }
        else { /* 'fixed-point' format */
            /* add number of leading or trailing 0's */
            if (the_exp <= 0) {
                /* add leading 0's */
                size += abs(the_exp) + 1;
            }
            else {
                /* add trailing 0's if needed */
                if (the_exp >= (buflen-isnegative))
                    size += (the_exp - (buflen - isnegative)) + 1;
            }
        }

        /* allocate the string itself (uninitialized, as yet) */
        res = PyBytes_FromStringAndSize(0, size);

        {
            /* proceed with building the string-buffer value */
            char* pd = PyBytes_AS_STRING(res);
            char* ps = buffer;

            /* insert leading tag if requested */
            if (optionflags & OP_TAG) {
                char* pt = ftag;
                while(*pt) *pd++ = *pt++;
            }

            /* copy sign if it's there */
            if (*ps=='-') {
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
                while (dtc>0) {
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
        temp = PyUnicode_FromString(PyBytes_AS_STRING(res));
        Py_DECREF(res);
        return temp;
#else
        return res;
#endif
    }
}

/* Classify an object as a type of number. If an object is recognized as a
 * number, it must be properly converted by the routines below.
 */

static int isNumber(PyObject* obj)
{
    if (options.debug)
        fprintf(stderr, "isNumber: object type is %s\n", Py_TYPE(obj)->tp_name);
    if (Pympz_Check(obj)) {
        return 1;
    }
    else if (PyIntOrLong_Check(obj)) {
        return 1;
    }
    else if (Pympq_Check(obj)) {
        return 1;
    }
    else if (Pympf_Check(obj)) {
        return 1;
    }
    else if (Pyxmpz_Check(obj)) {
        return 1;
    }
    else if (PyFloat_Check(obj)) {
        return 1;
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
        return 1;
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        return 1;
    }
    return 0;
}

static int isRational(PyObject* obj)
{
    if (options.debug)
        fprintf(stderr, "isRational: object type is %s\n", Py_TYPE(obj)->tp_name);
    if (Pympz_Check(obj)) {
        return 1;
    }
    else if (PyIntOrLong_Check(obj)) {
        return 1;
    }
    else if (Pympq_Check(obj)) {
        return 1;
    }
    else if (Pyxmpz_Check(obj)) {
        return 1;
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
        return 1;
    }
    return 0;
}

static int isInteger(PyObject* obj)
{
    if (options.debug)
        fprintf(stderr, "isInteger: object type is %s\n", Py_TYPE(obj)->tp_name);
    if (Pympz_Check(obj)) {
        return 1;
    }
    else if (PyIntOrLong_Check(obj)) {
        return 1;
    }
    else if (Pyxmpz_Check(obj)) {
        return 1;
    }
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
 * The routine anyrational2Pympq will convert an integer- and rational-like
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
    else if (Pympf_Check(obj)) {
        newob = Pympf2Pympq(obj);
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

    if (options.debug)
        fprintf(stderr,"anynum2Pympq(%p)->%p\n", obj, newob);

    return newob;
}

/* Convert an integer or mpz to mpq. */

static PympqObject*
anyrational2Pympq(PyObject* obj)
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

    if (options.debug)
        fprintf(stderr,"anyrational2Pympq(%p)->%p\n", obj, newob);

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
    else if (Pympf_Check(obj)) {
        newob = Pympf2Pympz(obj);
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
    if (options.debug)
        fprintf(stderr,"anynum2Pympz(%p)->%p\n", obj, newob);

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
    else if (Pympf_Check(obj)) {
        newob = Pympf2Pyxmpz(obj);
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
    if (options.debug)
        fprintf(stderr,"anynum2Pympz(%p)->%p\n", obj, newob);

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
    if (options.debug)
        fprintf(stderr,"Pympz_From_Integer(%p)->%p\n", obj, newob);
    if (!newob)
        TYPE_ERROR("conversion error in Pympz_From_Integer");
    return newob;
}

static PyxmpzObject*
Pyxmpz_From_Integer(PyObject* obj)
{
    PyxmpzObject* newob = 0;

    if (Pyxmpz_Check(obj)) {
        Py_INCREF(obj);
        newob = (PyxmpzObject*) obj;
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pyxmpz(obj);
#endif
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pyxmpz(obj);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pyxmpz(obj);
    }
    if (options.debug)
        fprintf(stderr,"Pyxmpz_From_Integer(%p)->%p\n", obj, newob);
    if (!newob)
        TYPE_ERROR("conversion error in Pyxmpz_From_Integer");
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
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        return PyInt_AS_LONG(obj);
#endif
    }
    else if (Pympz_Check(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
    }
    else if (Pyxmpz_Check(obj)) {
        if (mpz_fits_slong_p(Pyxmpz_AS_MPZ(obj))) {
            return mpz_get_si(Pyxmpz_AS_MPZ(obj));
        }
    }
    TYPE_ERROR("conversion error in clong_From_Integer");
    return -1;
}

static PympfObject*
anynum2Pympf(PyObject* obj, mpfr_prec_t bits)
{
    PympfObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympf_Check(obj)) {
        newob = (PympfObject*) obj;
        if (!bits || mpfr_get_prec(newob->f) == bits) {
            Py_INCREF(obj);
        }
        else {
            newob = Pympf2Pympf((PyObject*)newob, bits);
        }
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympf(obj, bits);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympf(obj, bits);
#endif
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pympf(obj, bits);
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympf(obj, bits);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympf(obj, bits);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympf(obj, bits);
    }
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympf(s, 10, bits);
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
            newob = Pympq2Pympf((PyObject *)temp, bits);
            Py_DECREF(s); Py_DECREF((PyObject*)temp);
        }
    }

    if (options.debug)
        fprintf(stderr, "anynum2Pympf(%p,%ld)->%p (%ld)\n", obj,
                (long)bits, newob, newob != 0 ? (long)mpfr_get_prec(newob->f) : -1);

    return newob;
}

/*
 * coerce any number to a mpz
 */
int
Pympz_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympzObject* newob = Pympz_From_Integer(arg);

    if (options.debug)
        fprintf(stderr, "mpz_conv_arg(%p)->%p\n", arg, newob);
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
    PympqObject* newob = anyrational2Pympq(arg);

    if (options.debug)
        fprintf(stderr, "mpq_conv_arg(%p)->%p\n", arg, newob);
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
Pympf_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympfObject* newob = anynum2Pympf(arg,0);

    if (options.debug)
        fprintf(stderr, "mpf_conv_arg(%p)->%p\n", arg, newob);
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to mpf");
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
Pympf2str(PympfObject *self)
{
    /* base-10, FP for exp -2 to 8, no tag */
    return Pympf_ascii(self, 10, 0, -2, 8, 0);
}

static PyObject *
Pympf2repr(PympfObject *self)
{
    /* base-10, always mantissa+exp, with tag */
    return Pympf_ascii(self, 10, 0, 0, -1, OP_TAG);
}

/*
 * The following functions are located in gmpy_misc.c:
 *
 *   Pygmpy_get_license(PyObject *self, PyObject *args)
 *   Pygmpy_get_version(PyObject *self, PyObject *args)
 *   Pygmpy_get_cvsid(PyObject *self, PyObject *args)
 *   Pygmpy_get_gmp_version(PyObject *self, PyObject *args)
 *   Pygmpy_get_mpir_version(PyObject *self, PyObject *args)
 *   Pygmpy_get_gmp_limbsize(PyObject *self, PyObject *args)
 *   Pygmpy_get_zcache(PyObject *self, PyObject *args)
 *   Pygmpy_get_qcache(PyObject *self, PyObject *args)
 *   Pygmpy_get_fcache(PyObject *self, PyObject *args)
 *   Pygmpy_set_zcache(PyObject *self, PyObject *args)
 *   Pygmpy_set_qcache(PyObject *self, PyObject *args)
 *   Pygmpy_set_fcache(PyObject *self, PyObject *args)
 *   Pygmpy_set_debug(PyObject *self, PyObject *args)
 *   Pygmpy_set_precison(PyObject *self, PyObject *args)
 *   Pygmpy_get_precison(PyObject *self, PyObject *args)
 *   Pygmpy_set_fcoform(PyObject *self, PyObject *args)
 *
 */

#include "gmpy_misc.c"

/*
 * The following functions are located in gmpy_mpz.c:
 *
 *   Pympz_copy(PyObject *self, PyObject *args)
 *   Pympz_binary(PyObject *self, PyObject *args)
 *   Pympz_digits(PyObject *self, PyObject *args)
 *   Pympz_numdigits(PyObject *self, PyObject *args)
 *   Pympz_bit_length(PyObject *self, PyObject *args)
 *   Pympz_scan0(PyObject *self, PyObject *args)
 *   Pympz_scan1(PyObject *self, PyObject *args)
 *   Pympz_popcount(PyObject *self, PyObject *args)
 *   Pympz_lowbits(PyObject *self, PyObject *args)
 *   Pympz_bit_test(PyObject *self, PyObject *args)
 *   Pympz_setbit(PyObject *self, PyObject *args)
 *   Pympz_root(PyObject *self, PyObject *args)
 *   Pympz_rootrem(PyObject *self, PyObject *args)
 *   Pympz_sign(PyObject *self, PyObject *args)
 *   Pympz_abs(PympzObject *x)
 *   Pympz_neg(PympzObject *x)
 *   Pympz_pos(PympzObject *x)
 *   Pympz_pow(PyObject *in_b, PyObject *in_e, PyObject *in_m)
 *   Pympz_nonzero(PympzObject *x)
 *   Pympz_com(PympzObject *x)
 *   Pympz_and(PyObject *a, PyObject *b)
 *   Pympz_ior(PyObject *a, PyObject *b)
 *   Pympz_xor(PyObject *a, PyObject *b)
 *   Pympz_rshift(PyObject *a, PyObject *b)
 *   Pympz_lshift(PyObject *a, PyObject *b)
 *   Pympz_oct(PympzObject *self)
 *   Pympz_hex(PympzObject *self)
 *   Pympz_hash(PympzObject *self)
 *   Pygmpy_gcd(PyObject *self, PyObject *args)
 *   Pygmpy_lcm(PyObject *self, PyObject *args)
 *   Pygmpy_gcdext(PyObject *self, PyObject *args)
 *   Pygmpy_divm(PyObject *self, PyObject *args)
 *   Pygmpy_fac(PyObject *self, PyObject *args)
 *   Pygmpy_fib(PyObject *self, PyObject *args)
 *   Pygmpy_fib2(PyObject *self, PyObject *args)
 *   Pygmpy_lucas(PyObject *self, PyObject *args)
 *   Pygmpy_lucas2(PyObject *self, PyObject *args)
 *   Pympz_bincoef(PyObject *self, PyObject *args)
 *   Pympz_sqrt(PyObject *self, PyObject *args)
 *   Pympz_sqrtrem(PyObject *self, PyObject *args)
 *   Pympz_remove(PyObject *self, PyObject *args)
 *   Pympz_invert(PyObject *self, PyObject *args)
 *   Pympz_hamdist(PyObject *self, PyObject *args)
 *   Pympz_divexact(PyObject *self, PyObject *args)
 *   Pygmpy_cdivmod(PyObject *self, PyObject *args)
 *   Pygmpy_fdivmod(PyObject *self, PyObject *args)
 *   Pygmpy_tdivmod(PyObject *self, PyObject *args)
 *   Pympz_is_square(PyObject *self, PyObject *args)
 *   Pympz_is_power(PyObject *self, PyObject *args)
 *   Pympz_is_prime(PyObject *self, PyObject *args)
 *   Pympz_next_prime(PyObject *self, PyObject *args)
 *   Pympz_jacobi(PyObject *self, PyObject *args)
 *   Pympz_legendre(PyObject *self, PyObject *args)
 *   Pympz_kronecker(PyObject *self, PyObject *args)
 *
 */

#include "gmpy_mpz.c"
#include "gmpy_mpz_divmod2exp.c"
#include "gmpy_mpz_divmod.c"
#include "gmpy_mpf.c"

/* produce digits for an mpq in requested base, default 10 */
static char doc_qdigitsm[]="\
x.digits([base]): returns Python string representing x in the\n\
given base (2 to 36, default 10 if omitted or 0); leading '-'\n\
is present if x<0, but no leading '+' if x>=0.\n\
";
static char doc_qdigitsg[]="\
qdigits(x[,base]): returns Python string representing x in the\n\
given base (2 to 36, default 10 if omitted or 0); leading '-'\n\
present if x<0, but no leading '+' if x>=0. x must be an mpq,\n\
or else gets coerced into one.\n\
";
static PyObject *
Pympq_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    PyObject *result;

    SELF_MPQ_ONE_ARG("|i", &base);
    assert(Pympq_Check(self));
    result = Pympq_ascii((PympqObject*)self, base, 0);
    Py_DECREF(self);
    return result;
}

/* produce string for an mpf with requested/defaulted parameters */
static char doc_fdigitsm[]="\
x.digits(base=10, digs=0, mine=0, maxe=-1, opts=0): formats x.\n\
\n\
Returns up to digs digits in the given base (if digs is 0, as many\n\
digits as are available), but no more than available given x's\n\
precision; the resulting string is formatted in fixed point\n\
if the exponent is >=mine and <=maxe, else in exponential (the\n\
exponent-separator is 'e' for base up to 10, else '@' -- the\n\
exponent is always output as a signed, base-10 integer). If opts\n\
has bit 1 set, the whole is wrapped in 'gmpy2.mpf(...)', to ease\n\
later approximate reconstruction via builtin function eval\n\
\n\
If opts has bit 2 set, then opts bit 1, mine, and maxe, are\n\
ignored; the result is then a 2-element tuple, first element\n\
the raw string of base-digits without formatting, second the\n\
exponent in base as a Python int.\n\
";

static char doc_fdigitsg[]="\
fdigits(x, base=10, digs=0, mine=0, maxe=-1, opts=0): formats x,\n\
which is an mpf or else gets coerced to one.\n\
\n\
Returns up to digs digits in the given base (if digs is 0, as many\n\
digits as are available), but no more than available given x's\n\
precision; the resulting string is formatted in fixed point\n\
if the exponent is >=mine and <=maxe, else in exponential (the\n\
exponent-separator is 'e' for base up to 10, else '@' -- the\n\
exponent is always output as a signed, base-10 integer). If opts\n\
has bit 1 set, the whole is wrapped in 'gmpy2.mpf(...)', to ease\n\
later approximate reconstruction via builtin function eval\n\
\n\
If opts has bit 2 set, then opts bit 1, mine, and maxe, are\n\
ignored; the result is then a 2-element tuple, first element\n\
the raw string of base-digits without formatting, second the\n\
exponent in base as a Python int.\n\
";

static PyObject *
Pympf_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    int digs = 0;
    int mine = 0;
    int maxe = -1;
    int opts = 0;
    PyObject *result;

    if (self && Pympf_Check(self)) {
        if (!PyArg_ParseTuple(args, "|iiiii", &base, &digs, &mine, &maxe, &opts))
            return NULL;
        Py_INCREF(self);
    }
    else {
        if (!PyArg_ParseTuple(args, "O&|iiiii", Pympf_convert_arg, &self, &base,
                &digs, &mine, &maxe, &opts))
            return NULL;
    }
    assert(Pympf_Check(self));
    result = Pympf_ascii( (PympfObject*)self, base, digs, mine, maxe, opts);
    Py_DECREF(self);
    return result;
}

static char doc_qsignm[]="\
x.sign(): returns -1, 0, or +1, if x is negative, 0, positive.\n\
";
static char doc_qsigng[]="\
qsign(x): returns -1, 0, or +1, if x is negative, 0, positive;\n\
x must be an mpq, or else gets coerced to one.\n\
";
static PyObject *
Pympq_sign(PyObject *self, PyObject *args)
{
    PyObject *result;

    SELF_MPQ_NO_ARG;
    assert(Pympq_Check(self));
    result = Py_BuildValue("i", mpq_sgn(Pympq_AS_MPQ(self)));
    Py_DECREF(self);
    return result;
}

static char doc_numerm[]="\
x.numer(): returns numerator of x.\n\
";
static char doc_numerg[]="\
numer(x): returns numerator of x;\n\
x must be an mpq, or else gets coerced to one.\n\
";
static PyObject *
Pympq_numer(PyObject *self, PyObject *args)
{
    PympzObject *result;

    if (!(result = Pympz_new()))
        return NULL;

    SELF_MPQ_NO_ARG;
    assert(Pympq_Check(self));
    mpz_set(result->z, mpq_numref(Pympq_AS_MPQ(self)));
    Py_DECREF(self);
    return (PyObject*)result;
}

static PyObject *
Pympq_getnumer(PympqObject *self, void *closure)
{
    PympzObject *result;

    if ((result = Pympz_new()))
        mpz_set(result->z, mpq_numref(Pympq_AS_MPQ(self)));
    return (PyObject*)result;
}

static char doc_denomm[]="\
x.denom(): returns denominator of x.\n\
";
static char doc_denomg[]="\
denom(x): returns denominator of x;\n\
x must be an mpq, or else gets coerced to one.\n\
";
static PyObject *
Pympq_denom(PyObject *self, PyObject *args)
{
    PympzObject *result;

    if (!(result = Pympz_new()))
        return NULL;

    SELF_MPQ_NO_ARG;
    assert(Pympq_Check(self));
    mpz_set(result->z, mpq_denref(Pympq_AS_MPQ(self)));
    Py_DECREF(self);
    return (PyObject*)result;
}

static PyObject *
Pympq_getdenom(PympqObject *self, void *closure)
{
    PympzObject *result;

    if ((result = Pympz_new()))
        mpz_set(result->z, mpq_denref(Pympq_AS_MPQ(self)));
    return (PyObject*)result;
}

static char doc_qdivm[]="\
x.qdiv(y=1): returns x/y as mpz if possible, or as mpq\n\
if x is not exactly divisible by y.\n\
";
static char doc_qdivg[]="\
qdiv(x,y=1): returns x/y as mpz if possible, or as mpq\n\
if x is not exactly divisible by y.\n\
";
static int isOne(PyObject* obj)
{
    if (!obj)
        return 1;

    if (Pympq_Check(obj)) {
        return (0==mpz_cmp_ui(mpq_denref(Pympq_AS_MPQ(obj)),1)) &&
               (0==mpz_cmp_ui(mpq_numref(Pympq_AS_MPQ(obj)),1));
    }
    else if (Pympz_Check(obj)) {
        return 0==mpz_cmp_ui(Pympz_AS_MPZ(obj),1);
    }
    else if (Pyxmpz_Check(obj)) {
        return 0==mpz_cmp_ui(Pyxmpz_AS_MPZ(obj),1);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        return PyInt_AS_LONG(obj)==1;
#endif
    }
    else if (Pympf_Check(obj)) {
        return mpfr_get_d(Pympf_AS_MPF(obj), options.rounding)==1.0;
    }
    else if (PyFloat_Check(obj)) {
        return PyFloat_AS_DOUBLE(obj)==1.0;
    }
    else if (PyLong_Check(obj)) {
        return PyLong_AsLong(obj)==1;
    }
    return 0;
}
static PyObject *
Pympq_qdiv(PyObject *self, PyObject *args)
{
    PyObject *other = 0;
    PyObject *s = 0;
    int wasone;

    if ( self && Pympq_Check(self)) {
        if (!PyArg_ParseTuple(args, "|O", &other))
            return NULL;
    }
    else {
        if (!PyArg_ParseTuple(args, "O|O", &self, &other))
            return NULL;
    }
    wasone = isOne(other);
    /* optimize if self must be returned unchanged */
    if (Pympq_Check(self) && wasone) {
        /* optimize if self is mpq and result must==self */
        if (mpz_cmp_ui(mpq_denref(Pympq_AS_MPQ(self)), 1) != 0) {
            Py_INCREF(self);
            return self;
        }
        else {
            /* denominator is 1, optimize returning an mpz */
            s = (PyObject*)Pympz_new();
            mpz_set(Pympz_AS_MPZ(s), mpq_numref(Pympq_AS_MPQ(self)));
            return s;
        }
    }
    else if (Pympz_Check(self) && wasone) {
        /* optimize if self is mpz and result must==self */
        Py_INCREF(self);
        return self;
    }
    /* normal, non-optimized case: must make new object as result */
    self = (PyObject*)anyrational2Pympq(self);
    if (!self) {
        if (!PyErr_Occurred())
            TYPE_ERROR("first argument can not be converted to mpq");
        return NULL;
    }
    if (wasone) { /* self was mpf, float, int, long... */
        s = self;
    }
    else {     /* other explicitly present and !=1... must compute */
        other = (PyObject*)anyrational2Pympq(other);
        if (!other) {
            Py_DECREF(self);
            if (!PyErr_Occurred())
                TYPE_ERROR("second argument can not be converted to mpq");
            return NULL;
        }
        if (mpq_sgn(Pympq_AS_MPQ(other))==0) {
            PyObject* result = 0;
            ZERO_ERROR("qdiv: zero divisor");
            Py_DECREF(self);
            Py_DECREF(other);
            return result;
        }
        s = (PyObject*)Pympq_new();
        mpq_div(Pympq_AS_MPQ(s), Pympq_AS_MPQ(self), Pympq_AS_MPQ(other));
        Py_DECREF(self);
        Py_DECREF(other);
    }
    if (mpz_cmp_ui(mpq_denref(Pympq_AS_MPQ(s)), 1) != 0) {
        return s;
    }
    else {
        /* denominator is 1, return an mpz */
        PyObject* ss = (PyObject*)Pympz_new();
        if (ss)
            mpz_set(Pympz_AS_MPZ(ss), mpq_numref(Pympq_AS_MPQ(s)));
        Py_DECREF(s);
        return ss;
    }
}

static char doc_f2qm[]="\
x.f2q([err]): returns the 'best' mpq approximating x to\n\
within relative error err (default, x's precision); 'best'\n\
rationals as per Stern-Brocot tree; mpz if denom is 1.\n\
If err<0, error sought is 2.0 ** err.\n\
";
static char doc_f2qg[]="\
f2q(x[,err]): returns the 'best' mpq approximating x to\n\
within relative error err (default, x's precision); 'best'\n\
rationals as per Stern-Brocot tree; mpz if denom is 1.\n\
If err<0, error sought is 2.0 ** err.\n\
";
static PyObject *
Pympf_f2q(PyObject *self, PyObject *args)
{
    PympfObject *err = 0;
    PympfObject *fself;

    if (options.debug)
        fprintf(stderr, "Pympf_f2q: %p, %p\n", self, args);

    SELF_MPF_ONE_ARG_CONVERTED_OPT(&err);
    assert(Pympf_Check(self));
    fself = (PympfObject*)self;

    return f2q_internal(fself, err, mpfr_get_prec(fself->f), args!=0);
}

static PyObject *
f2q_internal(PympfObject* self, PympfObject* err, unsigned int bits, int mayz)
{
    PympqObject *res = 0;
    int i, negative, errsign;
    mpfr_t f, al, a, r1[3], r2[3], minerr, curerr, newerr, temp;

    assert(!err || Pympf_Check(err));
    errsign = err ? mpfr_sgn(err->f) : 0;
    if (errsign == 0) {
        if (err)
            Py_DECREF((PyObject*)err);
        if (!(err = Pympf_new(20))) {
            Py_DECREF((PyObject*)self);
            return NULL;
        }
        mpfr_set_si(err->f, 1, options.rounding);
        mpfr_div_2ui(err->f, err->f, bits, options.rounding);
    }
    else if (errsign < 0) {
        int ubits;
        mpfr_floor(err->f, err->f);
        ubits = (int)mpfr_get_d(err->f, options.rounding);
        mpfr_set_si(err->f, 1, options.rounding);
        mpfr_div_2si(err->f, err->f, -ubits, options.rounding);
    }
    if (!(res = Pympq_new()))
        return NULL;
    mpfr_init2(minerr, 20);
    mpfr_set(minerr, err->f, options.rounding);
    Py_DECREF((PyObject*)err);

    mpfr_init2(f, bits);
    if (mpfr_sgn(self->f) < 0) {
        negative = 1;
        mpfr_abs(f, self->f, options.rounding);
    }
    else {
        negative = 0;
        mpfr_set(f, self->f, options.rounding);
    }
    Py_DECREF((PyObject*)self);
    mpfr_init2(al, bits);
    mpfr_set(al, f, options.rounding);
    mpfr_init2(a, bits);
    mpfr_floor(a, al);
    mpfr_init2(temp, bits);
    for (i=0; i<3; ++i) {
        mpfr_init2(r1[i], bits);
        mpfr_init2(r2[i], bits);
    }
    mpfr_set_si(r1[0], 0, options.rounding);
    mpfr_set_si(r1[1], 0, options.rounding);
    mpfr_set_si(r1[2], 1, options.rounding);
    mpfr_set_si(r2[0], 0, options.rounding);
    mpfr_set_si(r2[1], 1, options.rounding);
    mpfr_set(r2[2], a, options.rounding);
    mpfr_init2(curerr, 20);
    mpfr_init2(newerr, 20);
    mpfr_reldiff(curerr, f, a, options.rounding);
    while (mpfr_cmp(curerr, minerr) > 0) {
        mpfr_sub(temp, al, a, options.rounding);
        mpfr_ui_div(al, 1, temp, options.rounding);
        mpfr_floor(a, al);
        mpfr_swap(r1[0], r1[1]);
        mpfr_swap(r1[1], r1[2]);
        mpfr_mul(r1[2], r1[1], a, options.rounding);
        mpfr_add(r1[2], r1[2], r1[0], options.rounding);
        mpfr_swap(r2[0], r2[1]);
        mpfr_swap(r2[1], r2[2]);
        mpfr_mul(r2[2], r2[1], a, options.rounding);
        mpfr_add(r2[2], r2[2], r2[0], options.rounding);
        mpfr_div(temp, r2[2], r1[2], options.rounding);
        mpfr_reldiff(newerr, f, temp, options.rounding);
        if (mpfr_cmp(curerr, newerr) <= 0) {
            mpfr_swap(r1[1],r1[2]);
            mpfr_swap(r2[1],r2[2]);
            break;
        }
        mpfr_swap(curerr, newerr);
    }
    if (mayz && (mpfr_cmp_ui(r1[2],1)==0)) {
        Py_DECREF((PyObject*)res);
        res = (PympqObject*)Pympz_new();
        mpfr_get_z(Pympz_AS_MPZ(res), r2[2], options.rounding);
        if (negative)
            mpz_neg(Pympz_AS_MPZ(res),Pympz_AS_MPZ(res));
    }
    else {
        mpfr_get_z(mpq_numref(res->q), r2[2], options.rounding);
        mpfr_get_z(mpq_denref(res->q), r1[2], options.rounding);
        if (negative)
            mpz_neg(mpq_numref(res->q), mpq_numref(res->q));
    }
    mpfr_clear(minerr);
    mpfr_clear(al);
    mpfr_clear(a);
    mpfr_clear(f);
    for (i=0; i<3; ++i) {
        mpfr_clear(r1[i]);
        mpfr_clear(r2[i]);
    }
    mpfr_clear(curerr);
    mpfr_clear(newerr);
    mpfr_clear(temp);
    return (PyObject*)res;
}

/* CONSTRUCTORS */
static char doc_mpz[] = "\
mpz(n):\n\
      builds an mpz object with a numeric value n (truncating n\n\
      to its integer part if it's a float or mpf)\n\
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
            if ((base!=0) && (base!=256) && ((base<2)||(base>36))) {
                VALUE_ERROR("base for gmpy2.mpz must be 0, 256, or in the "
                            "interval 2 ... 36 .");
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
    if (options.debug)
        fprintf(stderr, "Pygmpy_mpz: created mpz = %ld\n",
                mpz_get_si(newob->z));

    return (PyObject *) newob;
}

static char doc_xmpz[] = "\
xmpz(n):\n\
    builds an xmpz object from any number n (truncating n\n\
    to its integer part if it's a float or mpf)\n\
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
        if ((base!=0) && (base!=256) && ((base<2)||(base>36))) {
            VALUE_ERROR("gmpy2.xmpz(): base must be 0, 256, or in the "
                        "interval 2 ... 36 .");
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
    if (options.debug)
        fprintf(stderr, "Pygmpy_xmpz: created xmpz = %ld\n",
                mpz_get_si(newob->z));

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
    int argc;

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
            if ((base!=0) && (base!=256) && ((base<2)||(base>36))) {
                VALUE_ERROR("base for gmpy2.mpq() must be 0, 256, or in the "
                            "interval 2 ... 36 .");
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
    if (options.debug) {
        fputs("Pygmpy_mpq: created mpq = ", stderr);
        mpq_out_str(stderr, 10, newob->q);
        putc('\n', stderr);
    }
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

static char doc_mpf[] = "\
mpf(n):\n\
    builds an mpf object with a numeric value n (n may be any\n\
    Python number, or an mpz, mpq, or mpf object) and a default\n\
    precision (in bits) depending on the nature of n\n\
mpf(n,bits=0):\n\
    as above, but with the specified number of bits (0\n\
    means to use default precision, as above)\n\
mpf(s,bits=0,base=10):\n\
    builds an mpf object from a string s made up of\n\
    digits in the given base, possibly with fraction-part (with\n\
    period as a separator) and/or exponent-part (with exponent\n\
    marker 'e' for base<=10, else '@'). If base=256, s must be\n\
    a gmpy2.mpf portable binary representation as built by the\n\
    function gmpy2.fbinary (and the .binary method of mpf objects).\n\
    The resulting mpf object is built with a default precision (in\n\
    bits) if bits is 0 or absent, else with the specified number\n\
    of bits.\n\
";
static PyObject *
Pygmpy_mpf(PyObject *self, PyObject *args)
{
    PympfObject *newob;
    PyObject *obj;
    int argc;
    unsigned int bits=0;

    TRACE("Pygmpy_mpf() called...\n");

    assert(PyTuple_Check(args));

    argc = PyTuple_Size(args);
    if ((argc < 1) || (argc > 3)) {
        TYPE_ERROR("gmpy2.mpf() requires 1 to 3 arguments");
        return NULL;
    }

    obj = PyTuple_GetItem(args, 0);

    if (2 <= argc) {
        long sbits;
        PyObject *pbits = PyTuple_GetItem(args, 1);
        sbits = clong_From_Integer(pbits);
        if (sbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.mpf(): bits must be an integer");
            return NULL;
        }
        if (sbits<0) {
            VALUE_ERROR("bits for gmpy2.mpf must be >= 0");
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
                TYPE_ERROR("gmpy2.mpf(): base must be an integer");
                return NULL;
            }
            if ((base!=0) && (base!=256) && ((base<2)||(base>36))) {
                VALUE_ERROR("base for gmpy2.mpf must be 0, 256, or in the "
                            "interval 2 ... 36 .");
                return NULL;
            }
        }
        newob = PyStr2Pympf(obj, base, bits);
        if (!newob) {
            return NULL;
        }
    }
    else {
        if (argc==3) {
            TYPE_ERROR("gmpy2.mpf() with numeric 1st argument needs 1 or 2 arguments");
            return NULL;
        }
        newob = anynum2Pympf(obj, bits);
        if (!newob) {
            if (!PyErr_Occurred())
                TYPE_ERROR("gmpy2.mpf() requires numeric or string argument");
            return NULL;
        }
    }
    if (options.debug) {
        fputs("Pygmpy_mpf: created mpf = ", stderr);
        mpfr_out_str(stderr, 10, 0, newob->f, options.rounding);
        fprintf(stderr," bits=%ld (%d)\n", (long)mpfr_get_prec(newob->f), bits);
    }

    return (PyObject *) newob;
} /* Pygmpy_mpf() */

/* ARITHMETIC */

/* Include the basic arithmetic operations from gmpy_basic.py.
 *
 * The following routines are contained in gmpy_basic.py:
 *   Pympany_add()
 *   Pympany_sub()
 *   Pympany_mul()
 *   Pympany_floordiv()
 *   Pympany_truediv()
 *   Pympany_div2()     -- Python 2.x only!
 *   Pympany_divmod()
 */

#include "gmpy_utility.c"
#include "gmpy_basic.c"
#include "gmpy_mpz_inplace.c"
#include "gmpy_xmpz_inplace.c"

#define MPF_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *a, PyObject *b) \
{ \
  unsigned int bits, bbits; \
  PympfObject *r; \
  PympfObject *pa = 0; \
  PympfObject *pb = 0; \
  if (Pympf_Check(a) && Pympf_Check(b)) { \
    bits = mpfr_get_prec(((PympfObject*)a)->f); \
    bbits = mpfr_get_prec(((PympfObject*)b)->f); \
    if (bits>bbits) bits=bbits; \
    if (!(r = Pympf_new(bits))) { \
      return NULL; \
    } \
    NAME(r->f, ((PympfObject*)a)->f, ((PympfObject*)b)->f, options.rounding); \
    if (options.debug) fprintf(stderr, "Py" #NAME "-> %p", r); \
    return (PyObject *) r; \
  } else { \
    if (Pympf_Check(a)) { \
      bits = mpfr_get_prec(((PympfObject*)a)->f); \
    } else { \
      bits = mpfr_get_prec(((PympfObject*)b)->f); \
    } \
    pa = anynum2Pympf(a, bits); \
    pb = anynum2Pympf(b, bits); \
    if (!pa || !pb) { \
      Py_XDECREF((PyObject*)pa); \
      Py_XDECREF((PyObject*)pb); \
      Py_RETURN_NOTIMPLEMENTED; \
    } \
    if (options.debug) fprintf(stderr, "Py" #NAME ": %p, %p", pa, pb); \
    if (!(r = Pympf_new(bits))) { \
      Py_DECREF((PyObject*)pa); \
      Py_DECREF((PyObject*)pb); \
      return NULL; \
    } \
    NAME(r->f, pa->f, pb->f, options.rounding); \
    Py_DECREF((PyObject*)pa); \
    Py_DECREF((PyObject*)pb); \
    if (options.debug) fprintf(stderr, "Py" #NAME "-> %p", r); \
    return (PyObject *) r; \
  } \
}

#define MPQ_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *a, PyObject *b) \
{ \
  PympqObject *r; \
  PympqObject *pa = 0; \
  PympqObject *pb = 0; \
  pa = anyrational2Pympq(a); \
  pb = anyrational2Pympq(b); \
  if (!pa || !pb) { \
    Py_XDECREF((PyObject*)pa); \
    Py_XDECREF((PyObject*)pb); \
    Py_RETURN_NOTIMPLEMENTED; \
  } \
  if (options.debug) fprintf(stderr, "Py" #NAME ": %p, %p", pa, pb); \
  if (!(r = Pympq_new())) { \
    Py_DECREF((PyObject*)pa); \
    Py_DECREF((PyObject*)pb); \
    return NULL; \
  } \
  NAME(r->q, pa->q, pb->q); \
  Py_DECREF((PyObject*)pa); \
  Py_DECREF((PyObject*)pb); \
  if (options.debug) fprintf(stderr, "Py" #NAME "-> %p", r); \
  return (PyObject *) r; \
}

MPF_BINOP(mpfr_reldiff)

#define MPF_MONOP(NAME) \
static PyObject * \
Py##NAME(PympfObject *x) \
{ \
  PympfObject *r; \
  if (options.debug) fprintf(stderr, "Py" #NAME ": %p\n", x); \
  if (!(r = Pympf_new(mpfr_get_prec(x->f)))) return NULL; \
  NAME(r->f, x->f, options.rounding); \
  if (options.debug) fprintf(stderr, "Py" #NAME "-> %p\n", r); \
  return (PyObject *) r; \
}

#define MPQ_MONOP(NAME) \
static PyObject * \
Py##NAME(PympqObject *x) \
{ \
  PympqObject *r; \
  if (options.debug) fprintf(stderr, "Py" #NAME ": %p\n", x); \
  if (!(r = Pympq_new())) return NULL; \
  NAME(r->q, x->q); \
  if (options.debug) fprintf(stderr, "Py" #NAME "-> %p\n", r); \
  return (PyObject *) r; \
}

/* MPQ_MONOP(mpq_inv) */

MPQ_MONOP(mpq_neg)

static PyObject *
Pympq_abs(PympqObject *x)
{
    PympqObject *r;
    if (!(r = Pympq_new()))
        return NULL;
    mpq_set(r->q, x->q);
    mpz_abs(mpq_numref(r->q),mpq_numref(r->q));
    return (PyObject *) r;
}

MPF_MONOP(mpfr_abs)
MPF_MONOP(mpfr_neg)

static PyObject *
Pympq_pos(PympqObject *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject *) x;
}

static PyObject *
Pympf_pos(PympfObject *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject *) x;
}

static PyObject *
Pympq_pow(PyObject *in_b, PyObject *in_e, PyObject *m)
{
    PympqObject *r;
    PympqObject *b = anyrational2Pympq(in_b);
    PympqObject *e = anyrational2Pympq(in_e);

    int esign;
    unsigned long ultem;

    assert(Pympq_Check(b));
    assert(Pympq_Check(e));

    if (!b || !e) {
        Py_XDECREF((PyObject*)b);
        Py_XDECREF((PyObject*)e);
        Py_RETURN_NOTIMPLEMENTED;
    }

    if (options.debug)
        fprintf(stderr, "Pympq_pow: %p, %p, %p\n", b, e, m);

    if ((PyObject*)m != Py_None) {
        PyErr_SetString(PyExc_ValueError, "mpq.pow no modulo allowed");
        Py_DECREF((PyObject*)b);
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    if (!mpz_fits_slong_p(mpq_numref(e->q))) {
        PyErr_SetString(PyExc_ValueError, "mpq.pow outrageous exp num");
        Py_DECREF((PyObject*)b);
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    if (!mpz_fits_slong_p(mpq_denref(e->q))) {
        PyErr_SetString(PyExc_ValueError, "mpq.pow outrageous exp den");
        Py_DECREF((PyObject*)b);
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    if (!(r = Pympq_new())) {
        Py_DECREF((PyObject*)b);
        Py_DECREF((PyObject*)e);
        return NULL;
    }

    esign = mpq_sgn(e->q);
    if (esign == 0) {
        if (options.debug)
            fprintf(stderr, "Pympq_pow (ui,0) -> %p\n", r);
        mpq_set_si(r->q, 1, 1);
        Py_DECREF((PyObject*)b);
        Py_DECREF((PyObject*)e);
        return (PyObject*)r;
    }
    else if (esign < 0) {
        int bsign = mpq_sgn(b->q);
        if (bsign == 0) {
            PyObject* result = 0;
            ZERO_ERROR("mpq.pow 0 base to <0 exponent");
            Py_DECREF((PyObject*)r);
            Py_DECREF((PyObject*)b);
            Py_DECREF((PyObject*)e);
            return result;
        }
        if (bsign < 0) {
            mpz_neg(mpq_numref(r->q), mpq_denref(b->q));
        }
        else {
            mpz_set(mpq_numref(r->q), mpq_denref(b->q));
        }
        mpz_abs(mpq_denref(r->q), mpq_numref(b->q));
        ultem = -mpz_get_si(mpq_numref(e->q));
    }
    else {
        mpq_set(r->q, b->q);
        ultem = mpz_get_ui(mpq_numref(e->q));
    }
    if (ultem>1) {
        mpz_pow_ui(mpq_numref(r->q), mpq_numref(r->q), ultem);
        mpz_pow_ui(mpq_denref(r->q), mpq_denref(r->q), ultem);
    }
    ultem = mpz_get_ui(mpq_denref(e->q));
    if (ultem>1) {
        static char* msgi = "mpq.pow fractional exponent, inexact-root";
        char* msg = msgi;
        int exact=0;
        if (mpq_sgn(r->q)<0) {
            static char* msgi = "mpq.pow fractional exponent, nonreal-root";
            msg = msgi;
        }
        else {
            mpz_t temp; /* workaround mpz_root bug in GMP 3.1.1 */
            mpz_inoc(temp);
            exact = mpz_root(temp, mpq_numref(r->q), ultem);
            if (exact) {
                mpz_set(mpq_numref(r->q), temp);
                exact = mpz_root(temp, mpq_denref(r->q), ultem);
                mpz_set(mpq_denref(r->q), temp);
            }
        }
        if (!exact) {
            Py_DECREF((PyObject*)r);
            VALUE_ERROR(msg);
            Py_DECREF((PyObject*)b);
            Py_DECREF((PyObject*)e);
            return NULL;
        }
    }
    if (options.debug)
        fprintf(stderr, "Pympq_pow (ui) -> %p\n", r);
    Py_DECREF((PyObject*)b);
    Py_DECREF((PyObject*)e);
    return (PyObject*)r;
}

static PyObject *
Pympf_pow(PyObject *xb, PyObject *xe, PyObject *m)
{
    PympqObject *qb, *qe;
    PyObject *r;
    unsigned int bits;
    int iexpo;
    PympfObject *b = 0, *e = 0;

    if ((PyObject*)m != Py_None) {
        PyErr_SetString(PyExc_ValueError, "mpf.pow no modulo allowed");
        return NULL;
    }

    if ((Pympf_Check(xb) && Pympf_Check(xe))) {
        b = anynum2Pympf(xb, 0);
        e = anynum2Pympf(xe, 0);
    }
    else {
        if (Pympf_Check(xb)) {
            b = anynum2Pympf(xb, 0);
            e = anynum2Pympf(xe, mpfr_get_prec(((PympfObject*)xb)->f));
        }
        if (Pympf_Check(xe)) {
            b = anynum2Pympf(xb, mpfr_get_prec(((PympfObject*)xe)->f));
            e = anynum2Pympf(xe, 0);
        }
    }

    if (!e || !b) {
        Py_XDECREF((PyObject*)e);
        Py_XDECREF((PyObject*)b);
        Py_RETURN_NOTIMPLEMENTED;
    }

    bits = mpfr_get_prec(b->f);
    if (bits > mpfr_get_prec(e->f))
        bits = mpfr_get_prec(e->f);
    if (options.debug)
        fprintf(stderr, "Pympf_pow(%d): %p, %p, %p\n", bits, b, e, m);

    iexpo = (int)mpfr_get_d(e->f, options.rounding);
    if (iexpo>0 && 0==mpfr_cmp_si(e->f, iexpo)) {
        r = (PyObject*)Pympf_new(mpfr_get_prec(b->f));
        if (!r) {
            Py_DECREF((PyObject*)e);
            Py_DECREF((PyObject*)b);
            return 0;
        }
        mpfr_pow_ui(Pympf_AS_MPF(r), b->f, iexpo, options.rounding);
    }
    else {
        qb = Pympf2Pympq((PyObject*)b);
        qe = Pympf2Pympq((PyObject*)e);
        r = Pympq_pow((PyObject*)qb, (PyObject*)qe, (PyObject*)m);
        Py_DECREF((PyObject*)qb); Py_DECREF((PyObject*)qe);
        if (!r || !Pympq_Check(r)) {
            Py_DECREF((PyObject*)e);
            Py_DECREF((PyObject*)b);
            return r;
        }
        qb = (PympqObject*)r;
        r = (PyObject*)Pympq2Pympf((PyObject*)qb, bits);
        Py_DECREF((PyObject*)qb);
    }
    Py_DECREF((PyObject*)e);
    Py_DECREF((PyObject*)b);
    return r;
}

static PyObject *
Pympany_pow(PyObject *in_b, PyObject *in_e, PyObject *in_m)
{
    PyObject *temp_b = 0, *temp_e = 0, *temp_r = 0, *res = 0;

    if (isInteger(in_b) && isInteger(in_e)) {
        return Pympz_pow(in_b, in_e, in_m);
    }
    else if ((PyFloat_Check(in_b) && Pympz_Check(in_e)) || \
            (PyFloat_Check(in_e) && Pympz_Check(in_b))) {
        if (in_m != Py_None) {
            TYPE_ERROR("3rd argument not allowed");
            return NULL;
        }
        else {
            if (Pympz_Check(in_b)) {
                temp_b = Pympz2PyFloat((PympzObject*)in_b);
            }
            else if (Pympq_Check(in_b)) {
                temp_b = Pympq2PyFloat((PympqObject*)in_b);
            }
            else if (Pympz_Check(in_b)) {
                temp_b = Pympf2PyFloat((PympfObject*)in_b);
            }
            else if (PyFloat_Check(in_b)) {
                temp_b = in_b;
                Py_INCREF(temp_b);
            }
            if (!temp_b) {
                Py_RETURN_NOTIMPLEMENTED;
            }
            if (Pympz_Check(in_e)) {
                temp_e = Pympz2PyFloat((PympzObject*)in_e);
            }
            else if (Pympq_Check(in_e)) {
                temp_e = Pympq2PyFloat((PympqObject*)in_e);
            }
            else if (Pympz_Check(in_e)) {
                temp_e = Pympf2PyFloat((PympfObject*)in_e);
            }
            else if (PyFloat_Check(in_e)) {
                temp_e = in_e;
                Py_INCREF(temp_e);
            }
            if (!temp_e) {
                Py_DECREF((PyObject*)temp_b);
                Py_RETURN_NOTIMPLEMENTED;
            }
            temp_r = PyNumber_Power(temp_b, temp_e, Py_None);
            Py_DECREF((PyObject*)temp_b);
            Py_DECREF((PyObject*)temp_e);
            if (!temp_r) {
                return NULL;
            }
            res = (PyObject *)PyFloat2Pympf(temp_r, 0);
            Py_DECREF((PyObject*)temp_r);
            return (PyObject *)res;
        }
    }
    else if (isRational(in_b) && isRational(in_e)) {
        return Pympq_pow(in_b, in_e, in_m);
    }
    else if (isNumber(in_b) && isNumber(in_e)) {
        return Pympf_pow(in_b, in_e, in_m);
    }
    Py_RETURN_NOTIMPLEMENTED;
}

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

    if (options.debug) {
        fprintf(stderr, "rich_compare: type(a) is %s\n", Py_TYPE(a)->tp_name);
        fprintf(stderr, "rich_compare: type(b) is %s\n", Py_TYPE(b)->tp_name);
    }

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
    if (Pympf_Check(a) && Pympf_Check(b)) {
        TRACE("compare (mpf,mpf)\n");
        return _cmp_to_object(mpfr_cmp(Pympf_AS_MPF(a), Pympf_AS_MPF(b)), op);
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
        tempa = (PyObject*)anyrational2Pympq(a);
        tempb = (PyObject*)anyrational2Pympq(b);
        c = mpq_cmp(Pympq_AS_MPQ(tempa), Pympq_AS_MPQ(tempb));
        Py_DECREF(tempa);
        Py_DECREF(tempb);
        return _cmp_to_object(c, op);
    }
    if (isNumber(a) && isNumber(b)) {
        TRACE("compare (mpf,float)\n");
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
        tempa = (PyObject*)anynum2Pympf(a,0);
        tempb = (PyObject*)anynum2Pympf(b,0);
        c = mpfr_cmp(Pympf_AS_MPF(tempa), Pympf_AS_MPF(tempb));
        Py_DECREF(tempa);
        Py_DECREF(tempb);
        return _cmp_to_object(c, op);
    }
    Py_RETURN_NOTIMPLEMENTED;
}

static int
Pympq_nonzero(PympqObject *x)
{
    return mpq_sgn(x->q) != 0;
}
static int
Pympf_nonzero(PympfObject *x)
{
    return mpfr_sgn(x->f) != 0;
}

/* float-truncations (return still a float!) */

#define MPF_UNIOP(NAME) \
static PyObject * \
Py##NAME(PyObject* self, PyObject *args) \
{ \
  PympfObject *r; \
  if (self && Pympf_Check(self)) { \
      if (args && !PyArg_ParseTuple(args, "")) return NULL; \
      Py_INCREF(self); \
  } else { \
      if (!PyArg_ParseTuple(args, "O&", Pympf_convert_arg, &self)) return NULL; \
  } \
  assert(Pympf_Check(self)); \
  if (options.debug) fprintf(stderr, "Py" #NAME ": %p\n", self); \
  if (!(r = Pympf_new(mpfr_get_prec(((PympfObject*)self)->f)))) return NULL; \
  NAME(r->f, Pympf_AS_MPF(self)); \
  if (options.debug) fprintf(stderr, "Py" #NAME "-> %p\n", r); \
  Py_DECREF(self); \
  return (PyObject *) r; \
}

static char doc_ceilm[]="\
x.ceil(): returns an mpf that is the smallest integer >= x\n\
";
static char doc_ceilg[]="\
ceil(x): returns an mpf that is the smallest integer >= x\n\
x must be an mpf, or else gets coerced to one.\n\
";
MPF_UNIOP(mpfr_ceil)

static char doc_floorm[]="\
x.floor(): returns an mpf that is the smallest integer <= x\n\
";
static char doc_floorg[]="\
floor(x): returns an mpf that is the smallest integer <= x\n\
x must be an mpf, or else gets coerced to one.\n\
";
MPF_UNIOP(mpfr_floor)

static char doc_truncm[]="\
x.trunc(): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
";
static char doc_truncg[]="\
trunc(x): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
x must be an mpf, or else gets coerced to one.\n\
";
MPF_UNIOP(mpfr_trunc)

/* hashing */
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
static long
Pympf_hash(PympfObject *self)
{
#ifdef _PyHASH_MODULUS
    unsigned long hash = 0;
    long exp = 0;
    size_t mbits = 0, fsize;
    double notneeded;
    mpz_t hack;
    int sign;

    if (self->hash_cache != -1)
        return self->hash_cache;
    fsize = (self->f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;

    /* Calculate the hash of the mantissa. */
    if (self->f->_mpfr_sign>0) {
        hash = mpn_mod_1(self->f->_mpfr_d, fsize, _PyHASH_MODULUS);
        sign = 1;
    }
    else if (self->f->_mpfr_sign<0) {
        hash = mpn_mod_1(self->f->_mpfr_d, fsize, _PyHASH_MODULUS);
        sign = -1;
    }
    else {
        return 0;
    }

    /* Get the number of bits in the mantissa. Ugly hack. */
    hack->_mp_size = fsize;
    hack->_mp_d = self->f->_mpfr_d;
    mbits = mpz_sizeinbase(hack, 2);

    /* Get the exponent as a power of 2. */
    notneeded = self->f->_mpfr_exp;

    /* Calculate the final hash. */
    exp -= (long)mbits;
    exp = exp >= 0 ? exp % _PyHASH_BITS : _PyHASH_BITS-1-((-1-exp) % _PyHASH_BITS);
    hash = ((hash << exp) & _PyHASH_MODULUS) | hash >> (_PyHASH_BITS - exp);

    hash *= sign;
    if (hash==(unsigned long)-1)
        hash = (unsigned long)-2;
    return (self->hash_cache = (long)hash);
#else
    double temp;
    if (self->hash_cache != -1)
        return self->hash_cache;
    temp = mpfr_get_d(self->f, options.rounding);
    return (self->hash_cache = _Py_HashDouble(temp));
#endif
}
static long
Pympq_hash(PympqObject *self)
{
#ifdef _PyHASH_MODULUS
    long hash = 0;
    mpz_t temp, mask;

    if (self->hash_cache != -1)
        return self->hash_cache;

    mpz_inoc(temp);
    mpz_inoc(mask);
    mpz_set_si(mask, _PyHASH_MODULUS);

    if (!mpz_invert(temp, mpq_denref(self->q), mask)) {
        mpz_cloc(temp);
        mpz_cloc(mask);
        hash = _PyHASH_INF;
        if (mpz_sgn(mpq_numref(self->q))<0)
            hash = -hash;
        self->hash_cache = hash;
        return hash;
    }
    mpz_powm_ui(temp, mpq_denref(self->q), _PyHASH_MODULUS - 2, mask);

    hash = (long)mpz_tdiv_ui(mpq_numref(self->q), _PyHASH_MODULUS);
    mpz_mul_si(temp, temp, hash);
    hash = (long)mpz_tdiv_ui(temp, _PyHASH_MODULUS);

    if (mpz_sgn(mpq_numref(self->q))<0)
        hash = -hash;
    if (hash==-1) hash = -2;
    mpz_cloc(temp);
    mpz_cloc(mask);
    self->hash_cache = hash;
    return hash;
#else
    if (self->hash_cache != -1)
        return self->hash_cache;
    return (self->hash_cache = dohash(Pympq2PyFloat(self)));
#endif
}

static char doc_pi[]="\
pi(n): returns pi with n bits of precision in an mpf object\n\
";
static PyObject *
Pygmpy_pi(PyObject *self, PyObject *args)
{
    PympfObject *pi;

    if ((pi = Pympf_new(0)))
        mpfr_set_d(pi->f, 3.2, options.rounding);
    return (PyObject*)pi;
}


static char doc_fsqrtm[]="\
x.fsqrt(): returns the square root of x.  x must be >= 0.\n\
";
static char doc_fsqrtg[]="\
fsqrt(x): returns the square root of x.  x must be an mpf, or\n\
else gets coerced to one; further, x must be >= 0.\n\
";
static PyObject *
Pympf_sqrt(PyObject *self, PyObject *args)
{
    PympfObject *root;

    SELF_MPF_NO_ARG;

    assert(Pympf_Check(self));

    if (mpfr_sgn(Pympf_AS_MPF(self)) < 0) {
        PyErr_SetString(PyExc_ValueError, "sqrt of negative number");
        Py_DECREF(self);
        return NULL;
    }

    if (!(root = Pympf_new(mpfr_get_prec(((PympfObject*)self)->f)))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_sqrt(root->f, Pympf_AS_MPF(self), options.rounding);
    Py_DECREF(self);
    return (PyObject *) root;
}

/* Include helper functions for mpmath. */

#include "gmpy_mpmath.c"

static char doc_getprecm[]="\
x.getprec(): returns the number of bits of precision in x.\n\
";
static char doc_getprecg[]="\
getprec(x): returns the number of bits of precision in x,\n\
which must be an mpf or else gets coerced to one.\n\
";
static PyObject *
Pympf_getprec(PyObject *self, PyObject *args)
{
    long precres;

    SELF_MPF_NO_ARG;
    assert(Pympf_Check(self));

    precres = (long) mpfr_get_prec(Pympf_AS_MPF(self));
    Py_DECREF(self);
    return PyIntOrLong_FromLong(precres);
}

static char doc_froundm[] = "\
x.round(n): returns x rounded to least n bits. Actual precision will\n\
be a multiple of gmp_limbsize().\n\
";
static char doc_froundg[] = "\
fround(x, n): returns x rounded to least n bits. Actual precision will\n\
be a multiple of gmp_limbsize(). x an mpf or coerced to an mpf.\n\
";
static PyObject *
Pympf_round(PyObject *self, PyObject *args)
{
    /* Should really get default precision. */
    long prec = DBL_MANT_DIG;
    PyObject *s;

    SELF_MPF_ONE_ARG("|l",&prec);
    assert(Pympf_Check(self));
    s = (PyObject*)Pympf2Pympf(self, prec);
    Py_DECREF(self);
    return s;
}

static char doc_reldiffm[] = "\
x.reldiff(y): returns the relative difference between x and y,\n\
where y can be any number and gets coerced to an mpf; result is\n\
a non-negative mpf roughly equal to abs(x-y)/((abs(x)+abs(y))/2).\n\
";
static char doc_reldiffg[] = "\
reldiff(x,y): returns the relative difference between x and y,\n\
where x and y can be any numbers and get coerced to mpf; result is\n\
a non-negative mpf roughly equal to abs(x-y)/((abs(x)+abs(y))/2).\n\
";
static PyObject *
Pympf_doreldiff(PyObject *self, PyObject *args)
{
    PympfObject *op;
    PyObject *res;

    SELF_MPF_ONE_ARG_CONVERTED(&op);
    assert(Pympf_Check(self));

    res = Pympfr_reldiff((PyObject*)self, (PyObject*)op);
    Py_DECREF(self); Py_DECREF((PyObject*)op);

    return res;
}

static char doc_fsignm[]="\
x.sign(): returns -1, 0, or +1, if x is negative, 0, positive.\n\
";
static char doc_fsigng[]="\
fsign(x): returns -1, 0, or +1, if x is negative, 0, positive;\n\
x must be an mpf, or else gets coerced to one.\n\
";
static PyObject *
Pympf_sign(PyObject *self, PyObject *args)
{
    long sign;

    SELF_MPF_NO_ARG;
    assert(Pympf_Check(self));

    sign = (long) mpfr_sgn(Pympf_AS_MPF(self));
    Py_DECREF(self);
    return PyIntOrLong_FromLong(sign);
}

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
static PyNumberMethods mpf_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,              /* nb_negative             */
    (unaryfunc) Pympf_pos,               /* nb_positive             */
    (unaryfunc) Pympfr_abs,              /* nb_absolute             */
    (inquiry) Pympf_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympf2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympf2PyFloat,           /* nb_float                */
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
static PyNumberMethods mpf_number_methods =
{
    (binaryfunc) Pympany_add,            /* nb_add                  */
    (binaryfunc) Pympany_sub,            /* nb_subtract             */
    (binaryfunc) Pympany_mul,            /* nb_multiply             */
    (binaryfunc) Pympany_div2,           /* nb_divide               */
    (binaryfunc) Pympany_rem,            /* nb_remaider             */
    (binaryfunc) Pympany_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,               /* nb_negative             */
    (unaryfunc) Pympf_pos,               /* nb_positive             */
    (unaryfunc) Pympfr_abs,               /* nb_absolute             */
    (inquiry) Pympf_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympf2PyInt,             /* nb_int                  */
    (unaryfunc) Pympf2PyLong,            /* nb_long                 */
    (unaryfunc) Pympf2PyFloat,           /* nb_float                */
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

static PyGetSetDef Pympf_getseters[] =
{
    { "precision", (getter)Pympf_getprec2, NULL, "precision in bits", NULL },
    {NULL}
};

static PyMethodDef Pygmpy_methods [] =
{
    { "_cvsid", Pygmpy_get_cvsid, METH_NOARGS, doc_cvsid },
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
    { "cdiv", Pygmpy_cdiv, METH_VARARGS, doc_cdivg },
    { "cdiv2exp", Pygmpy_cdiv2exp, METH_VARARGS, doc_cdiv2expg },
    { "cdivmod", Pygmpy_cdivmod, METH_VARARGS, doc_cdivmodg },
    { "cdivmod2exp", Pygmpy_cdivmod2exp, METH_VARARGS, doc_cdivmod2expg },
    { "ceil", Pympfr_ceil, METH_VARARGS, doc_ceilg },
    { "cmod", Pygmpy_cmod, METH_VARARGS, doc_cmodg },
    { "cmod2exp", Pygmpy_cmod2exp, METH_VARARGS, doc_cmod2expg },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combg },
    { "_copy", Pympany_copy, METH_O, doc_copyg },
    { "denom", Pympq_denom, METH_VARARGS, doc_denomg },
    { "digits", Pympz_digits, METH_VARARGS, doc_digitsg },
    { "divexact", Pygmpy_divexact, METH_VARARGS, doc_divexactg },
    { "divm", Pygmpy_divm, METH_VARARGS, doc_divm },
    { "fac", Pygmpy_fac, METH_O, doc_fac },
    { "fdigits", Pympf_digits, METH_VARARGS, doc_fdigitsg },
    { "fdiv", Pygmpy_fdiv, METH_VARARGS, doc_fdivg },
    { "fdiv2exp", Pygmpy_fdiv2exp, METH_VARARGS, doc_fdiv2expg },
    { "fdivmod", Pygmpy_fdivmod, METH_VARARGS, doc_fdivmodg },
    { "fdivmod2exp", Pygmpy_fdivmod2exp, METH_VARARGS, doc_fdivmod2expg },
    { "fib", Pygmpy_fib, METH_O, doc_fib },
    { "fib2", Pygmpy_fib2, METH_O, doc_fib2 },
    { "floor", Pympfr_floor, METH_VARARGS, doc_floorg },
    { "fmod", Pygmpy_fmod, METH_VARARGS, doc_fmodg },
    { "fmod2exp", Pygmpy_fmod2exp, METH_VARARGS, doc_fmod2expg },
    { "fround", Pympf_round, METH_VARARGS, doc_froundg },
    { "fsign", Pympf_sign, METH_VARARGS, doc_fsigng },
    { "fsqrt", Pympf_sqrt, METH_VARARGS, doc_fsqrtg },
    { "f2q", Pympf_f2q, METH_VARARGS, doc_f2qg },
    { "gcd", Pygmpy_gcd, METH_VARARGS, doc_gcd },
    { "gcdext", Pygmpy_gcdext, METH_VARARGS, doc_gcdext },
    { "get_cache", Pygmpy_get_cache, METH_NOARGS, doc_get_cache },
    { "get_precision", Pygmpy_get_precision, METH_VARARGS, doc_get_precision },
    { "get_prefer_mutable", Pygmpy_get_prefer_mutable, METH_NOARGS, doc_get_prefer_mutable },
    { "gmp_version", Pygmpy_get_gmp_version, METH_NOARGS, doc_gmp_version },
    { "gmp_limbsize", Pygmpy_get_gmp_limbsize, METH_NOARGS, doc_gmp_limbsize },
    { "getprec", Pympf_getprec, METH_VARARGS, doc_getprecg },
    { "hamdist", Pympz_hamdist, METH_VARARGS, doc_hamdistg },
    { "invert", Pygmpy_invert, METH_VARARGS, doc_invertg },
    { "is_even", Pympz_is_even, METH_O, doc_is_eveng },
    { "is_odd", Pympz_is_odd, METH_O, doc_is_oddg },
    { "is_power", Pympz_is_power, METH_O, doc_is_powerg },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primeg },
    { "is_square", Pympz_is_square, METH_O, doc_is_squareg },
    { "jacobi", Pympz_jacobi, METH_VARARGS, doc_jacobig },
    { "kronecker", Pympz_kronecker, METH_VARARGS, doc_kroneckerm },
    { "lcm", Pygmpy_lcm, METH_VARARGS, doc_lcm },
    { "legendre", Pympz_legendre, METH_VARARGS, doc_legendreg },
    { "license", Pygmpy_get_license, METH_NOARGS, doc_license },
    { "lucas", Pygmpy_lucas, METH_O, doc_lucas },
    { "lucas2", Pygmpy_lucas2, METH_O, doc_lucas2 },
    { "mpf", Pygmpy_mpf, METH_VARARGS, doc_mpf },
    { "mpir_version", Pygmpy_get_mpir_version, METH_NOARGS, doc_mpir_version },
    { "mpq", Pygmpy_mpq, METH_VARARGS, doc_mpq },
    { "mpz", Pygmpy_mpz, METH_VARARGS, doc_mpz },
    { "next_prime", Pympz_next_prime, METH_O, doc_next_primeg },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsg },
    { "numer", Pympq_numer, METH_VARARGS, doc_numerg },
    { "pack", Pygmpy_pack, METH_VARARGS, doc_packg },
    { "pi", Pygmpy_pi, METH_VARARGS, doc_pi },
    { "popcount", Pympz_popcount, METH_O, doc_popcountg },
    { "qdigits", Pympq_digits, METH_VARARGS, doc_qdigitsg },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivg },
    { "qsign", Pympq_sign, METH_VARARGS, doc_qsigng },
    { "reldiff", Pympf_doreldiff, METH_VARARGS, doc_reldiffg },
    { "remove", Pympz_remove, METH_VARARGS, doc_removeg },
    { "root", Pympz_root, METH_VARARGS, doc_rootg },
    { "rootrem", Pympz_rootrem, METH_VARARGS, doc_rootremg },
    { "set_cache", Pygmpy_set_cache, METH_VARARGS, doc_set_cache },
    { "set_debug", Pygmpy_set_debug, METH_VARARGS, doc_set_debug },
    { "set_fcoform", Pygmpy_set_fcoform, METH_VARARGS, doc_set_fcoform },
    { "set_precision", Pygmpy_set_precision, METH_VARARGS, doc_set_precision },
    { "set_prefer_mutable", Pygmpy_set_prefer_mutable, METH_VARARGS, doc_set_prefer_mutable },
    { "sign", Pympz_sign, METH_O, doc_signg },
    { "sqrt", Pygmpy_sqrt, METH_O, doc_sqrtg },
    { "sqrtrem", Pympz_sqrtrem, METH_VARARGS, doc_sqrtremg },
    { "tdiv", Pygmpy_tdiv, METH_VARARGS, doc_tdivg },
    { "tdiv2exp", Pygmpy_tdiv2exp, METH_VARARGS, doc_tdiv2expg },
    { "tdivmod", Pygmpy_tdivmod, METH_VARARGS, doc_tdivmodg },
    { "tdivmod2exp", Pygmpy_tdivmod2exp, METH_VARARGS, doc_tdivmod2expg },
    { "tmod", Pygmpy_tmod, METH_VARARGS, doc_tmodg },
    { "tmod2exp", Pygmpy_tmod2exp, METH_VARARGS, doc_tmod2expg },
    { "trunc", Pympfr_trunc, METH_VARARGS, doc_truncg },
    { "unpack", Pygmpy_unpack, METH_VARARGS, doc_unpackg },
    { "version", Pygmpy_get_version, METH_NOARGS, doc_version },
    { "xmpz", Pygmpy_xmpz, METH_VARARGS, doc_xmpz },
    { "_mpmath_normalize", Pympz_mpmath_normalize, METH_VARARGS, doc_mpmath_normalizeg },
    { "_mpmath_create", Pympz_mpmath_create, METH_VARARGS, doc_mpmath_createg },
    { "_mpmath_trim", Pympz_mpmath_trim, METH_VARARGS, doc_mpmath_trimg },
    { "_mpmath_add", Pympz_mpmath_add, METH_VARARGS, doc_mpmath_addg },
    { "_mpmath_mult", Pympz_mpmath_mult, METH_VARARGS, doc_mpmath_multg },
    { "_mpmath_div", Pympz_mpmath_div, METH_VARARGS, doc_mpmath_divg },
    { "_mpmath_sqrt", Pympz_mpmath_sqrt, METH_VARARGS, doc_mpmath_sqrtg },

    { NULL, NULL, 1}
};

static PyMethodDef Pympz_methods [] =
{
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "bincoef", Pympz_bincoef, METH_VARARGS, doc_bincoefm },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "cdiv", Pympz_cdiv, METH_O, doc_cdivm },
    { "cdiv2exp", Pympz_cdiv2exp, METH_O, doc_cdiv2expm },
    { "cmod", Pympz_cmod, METH_O, doc_cmodm },
    { "cmod2exp", Pympz_cmod2exp, METH_O, doc_cmod2expm },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combm },
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "digits", Pympz_digits, METH_VARARGS, doc_digitsm },
    { "divexact", Pympz_divexact, METH_O, doc_divexactm },
    { "fdiv", Pympz_fdiv, METH_O, doc_fdivm },
    { "fdiv2exp", Pympz_fdiv2exp, METH_O, doc_fdiv2expm },
    { "fmod", Pympz_fmod, METH_O, doc_fmodm },
    { "fmod2exp", Pympz_fmod2exp, METH_O, doc_fmod2expm },
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
    { "root", Pympz_root, METH_VARARGS, doc_rootm },
    { "rootrem", Pympz_rootrem, METH_VARARGS, doc_rootremm },
    { "sign", Pympz_sign, METH_NOARGS, doc_signm },
    { "sqrt", Pympz_sqrt, METH_NOARGS, doc_sqrtm },
    { "sqrtrem", Pympz_sqrtrem, METH_VARARGS, doc_sqrtremm },
    { "tdiv", Pympz_tdiv, METH_O, doc_tdivm },
    { "tdiv2exp", Pympz_tdiv2exp, METH_O, doc_tdiv2expm },
    { "tmod", Pympz_tmod, METH_O, doc_tmodm },
    { "tmod2exp", Pympz_tmod2exp, METH_O, doc_tmod2expm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pyxmpz_methods [] =
{
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "bincoef", Pympz_bincoef, METH_VARARGS, doc_bincoefm },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "cdiv", Pympz_cdiv, METH_O, doc_cdivm },
    { "cdiv2exp", Pympz_cdiv2exp, METH_O, doc_cdiv2expm },
    { "cmod", Pympz_cmod, METH_O, doc_cmodm },
    { "cmod2exp", Pympz_cmod2exp, METH_O, doc_cmod2expm },
    { "comb", Pympz_bincoef, METH_VARARGS, doc_combm },
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "digits", Pyxmpz_digits, METH_VARARGS, doc_digitsm },
    { "divexact", Pympz_divexact, METH_O, doc_divexactm },
    { "fdiv", Pympz_fdiv, METH_O, doc_fdivm },
    { "fdiv2exp", Pympz_fdiv2exp, METH_O, doc_fdiv2expm },
    { "fmod", Pympz_fmod, METH_O, doc_fmodm },
    { "fmod2exp", Pympz_fmod2exp, METH_O, doc_fmod2expm },
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
    { "root", Pympz_root, METH_VARARGS, doc_rootm },
    { "rootrem", Pympz_rootrem, METH_VARARGS, doc_rootremm },
    { "sign", Pympz_sign, METH_NOARGS, doc_signm },
    { "sqrt", Pyxmpz_sqrt, METH_NOARGS, doc_sqrtm },
    { "sqrtrem", Pympz_sqrtrem, METH_VARARGS, doc_sqrtremm },
    { "tdiv", Pympz_tdiv, METH_O, doc_tdivm },
    { "tdiv2exp", Pympz_tdiv2exp, METH_O, doc_tdiv2expm },
    { "tmod", Pympz_tmod, METH_O, doc_tmodm },
    { "tmod2exp", Pympz_tmod2exp, METH_O, doc_tmod2expm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pympq_methods [] =
{
    { "sign", Pympq_sign, METH_VARARGS, doc_qsignm },
    { "numer", Pympq_numer, METH_VARARGS, doc_numerm },
    { "denom", Pympq_denom, METH_VARARGS, doc_denomm },
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "digits", Pympq_digits, METH_VARARGS, doc_qdigitsm },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pympf_methods [] =
{
    { "reldiff", Pympf_doreldiff, METH_VARARGS, doc_reldiffm },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "digits", Pympf_digits, METH_VARARGS, doc_fdigitsm },
    { "round", Pympf_round, METH_VARARGS, doc_froundm },
    { "getprec", Pympf_getprec, METH_VARARGS, doc_getprecm },
    { "_copy", Pympany_copy, METH_NOARGS, doc_copym },
    { "sign", Pympf_sign, METH_VARARGS, doc_fsignm },
    { "sqrt", Pympf_sqrt, METH_VARARGS, doc_fsqrtm },
    { "qdiv", Pympq_qdiv, METH_VARARGS, doc_qdivm },
    { "f2q", Pympf_f2q, METH_VARARGS, doc_f2qm },
    { "ceil", Pympfr_ceil, METH_VARARGS, doc_ceilm },
    { "floor", Pympfr_floor, METH_VARARGS, doc_floorm },
    { "trunc", Pympfr_trunc, METH_VARARGS, doc_truncm },
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


static PyTypeObject Pympf_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpf",                                  /* tp_name          */
    sizeof(PympfObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympf_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympf2repr,                  /* tp_repr          */
    &mpf_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympf_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympf2str,                   /* tp_str           */
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
    Pympf_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympf_getseters,                        /* tp_getset        */
};


static void *
gmpy_allocate(size_t size)
{
    void *res;
    size_t usize = size;
    if (usize < GMPY_ALLOC_MIN) usize = GMPY_ALLOC_MIN;

    if (options.debug)
        fprintf(stderr, "mp_allocate( %llu->%llu )\n",
            (unsigned long long)size, (unsigned long long)usize);
    if (!(res = PyMem_Malloc(usize))) {
        fprintf(stderr, "mp_allocate( %llu->%llu )\n",
            (unsigned long long)size, (unsigned long long)usize);
        Py_FatalError("mp_allocate failure");
    }

    if (options.debug)
        fprintf(stderr, "mp_allocate( %llu->%llu ) ->%8p\n",
            (unsigned long long)size, (unsigned long long)usize, res);

    return res;
} /* mp_allocate() */


static void *
gmpy_reallocate(void *ptr, size_t old_size, size_t new_size)
{
    void *res;
    size_t uold = old_size;
    size_t unew = new_size;
    if (uold < GMPY_ALLOC_MIN) uold = GMPY_ALLOC_MIN;
    if (unew < GMPY_ALLOC_MIN) unew = GMPY_ALLOC_MIN;

    if (options.debug)
        fprintf(stderr,
            "mp_reallocate: old address %8p, old size %llu(%llu), new %llu(%llu)\n",
            ptr, (unsigned long long)old_size, (unsigned long long)uold,
            (unsigned long long)new_size, (unsigned long long)unew);

    if (uold==unew) {
        if (options.debug)
            fprintf(stderr, "mp_reallocate: avoided realloc for %llu\n",
                (unsigned long long)unew);
        return ptr;
    }

    if (!(res = PyMem_Realloc(ptr, unew))) {
        fprintf(stderr,
            "mp_reallocate: old address %8p, old size %llu(%llu), new %llu(%llu)\n",
            ptr, (unsigned long long)old_size, (unsigned long long)uold,
            (unsigned long long)new_size, (unsigned long long)unew);
        Py_FatalError("mp_reallocate failure");
    }

    if (options.debug)
        fprintf(stderr, "mp_reallocate: newob address %8p, newob size %llu(%llu)\n",
        res, (unsigned long long)new_size, (unsigned long long)unew);

    return res;
} /* mp_reallocate() */

static void
gmpy_free( void *ptr, size_t size)
{
    size_t usize=size;
    if (usize < GMPY_ALLOC_MIN) usize = GMPY_ALLOC_MIN;

    if (options.debug)
        fprintf(stderr, "mp_free      : old address %8p, old size %llu(%llu)\n",
            ptr, (unsigned long long)size, (unsigned long long)usize);

    PyMem_Free(ptr);
} /* mp_free() */

static void _PyInitGMP(void)
{
    mp_set_memory_functions(gmpy_allocate, gmpy_reallocate, gmpy_free);
    options.precision = DBL_MANT_DIG;
    set_zcache();
    set_qcache();
    set_pympzcache();
    set_pympqcache();
    set_pyxmpzcache();
}

static char _gmpy_docs[] = "\
gmpy2 2.0.0a1 - General Multiprecision arithmetic for Python:\n\
exposes functionality from the GMP or MPIR library to Python 2.6\n\
and later.\n\
\n\
Allows creation of multiprecision integer (mpz), mutable integers\n\
(xmpz), float (mpf), and rational (mpq) numbers, conversion between\n\
them and to/from Python numbers/strings, arithmetic, bitwise, and\n\
some other higher-level mathematical operations.\n\
\n\
mpz has comparable functionality to Python's builtin longs, but\n\
can be faster for some operations (particularly multiplication\n\
and raising-to-power) and has many further useful and speedy\n\
functions (prime testing and generation, factorial, fibonacci,\n\
binary-coefficients, gcd, lcm, square and other roots, ...).\n\
\n\
mpf and mpq only offer basic arithmetic abilities, but they\n\
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
    if (PyType_Ready(&Pympf_Type) < 0)
        INITERROR;
    if (PyType_Ready(&Pyxmpz_Type) < 0)
        INITERROR;

    if (do_debug)
        sscanf(do_debug, "%d", &options.debug);

    if (options.debug)
        fputs( "initgmpy2() called...\n", stderr );
    _PyInitGMP();

#ifdef PY3
    gmpy_module = PyModule_Create(&moduledef);
#else
    gmpy_module = Py_InitModule3("gmpy2", Pygmpy_methods, _gmpy_docs);
#endif

    /* Todo: Add error checking for status of gmpy_module returned above. */

    if (options.debug)
        fprintf(stderr, "gmpy_module at %p\n", gmpy_module);

    /* Add support for pickling. */
#ifdef PY3
    copy_reg_module = PyImport_ImportModule("copyreg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def mpz_reducer(an_mpz): return (gmpy2.mpz, (an_mpz.binary(), 256))\n"
            "def mpq_reducer(an_mpq): return (gmpy2.mpq, (an_mpq.binary(), 256))\n"
            "def mpf_reducer(an_mpf): return (gmpy2.mpf, (an_mpf.binary(), 0, 256))\n"
            "copyreg.pickle(type(gmpy2.mpz(0)), mpz_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpq(0)), mpq_reducer)\n"
            "copyreg.pickle(type(gmpy2.mpf(0)), mpf_reducer)\n"
        ;
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;
        if (options.debug)
            fprintf(stderr, "gmpy_module imported copyreg OK\n");
        PyDict_SetItemString(namespace, "copyreg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (result) {
            if (options.debug)
                fprintf(stderr, "gmpy_module enable pickle OK\n");
        }
        else {
            if (options.debug)
                fprintf(stderr, "gmpy_module could not enable pickle\n");
            PyErr_Clear();
        }
        Py_DECREF(namespace);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
        if (options.debug)
            fprintf(stderr, "gmpy_module could not import copyreg\n");
    }
#else
    copy_reg_module = PyImport_ImportModule("copy_reg");
    if (copy_reg_module) {
        char* enable_pickle =
            "def mpz_reducer(an_mpz): return (gmpy2.mpz, (an_mpz.binary(), 256))\n"
            "def mpq_reducer(an_mpq): return (gmpy2.mpq, (an_mpq.binary(), 256))\n"
            "def mpf_reducer(an_mpf): return (gmpy2.mpf, (an_mpf.binary(), 0, 256))\n"
            "copy_reg.pickle(type(gmpy2.mpz(0)), mpz_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpq(0)), mpq_reducer)\n"
            "copy_reg.pickle(type(gmpy2.mpf(0)), mpf_reducer)\n"
        ;
        PyObject* namespace = PyDict_New();
        PyObject* result = NULL;
        if (options.debug)
            fprintf(stderr, "gmpy_module imported copy_reg OK\n");
        PyDict_SetItemString(namespace, "copy_reg", copy_reg_module);
        PyDict_SetItemString(namespace, "gmpy2", gmpy_module);
        PyDict_SetItemString(namespace, "type", (PyObject*)&PyType_Type);
        result = PyRun_String(enable_pickle, Py_file_input,
                              namespace, namespace);
        if (result) {
            if (options.debug)
                fprintf(stderr, "gmpy_module enable pickle OK\n");
        }
        else {
            if (options.debug)
                fprintf(stderr, "gmpy_module could not enable pickle\n");
            PyErr_Clear();
        }
        Py_DECREF(namespace);
        Py_XDECREF(result);
    }
    else {
        PyErr_Clear();
        if (options.debug)
            fprintf(stderr, "gmpy_module could not import copy_reg\n");
    }
#endif


#ifdef PY3
    return gmpy_module;
#endif
}
