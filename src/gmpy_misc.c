/* gmpy_misc.c
 *
 * Miscellaneous module-level functions.
 *
 * This file should be considered part of gmpy2.c
 */

/* Return license information. */

#ifdef __MPIR_VERSION
#define MPIR_VER \
__MPIR_VERSION * 10000 + \
__MPIR_VERSION_MINOR * 100 + \
__MPIR_VERSION_PATCHLEVEL

#if ((MPIR_VER < 10300) && (MPIR_VER > 10399))
static char gmpy_license[] = "\
The GMPY2 source code is licensed under LGPL 2.1 or later. This version \n\
of the MPIR library is licensed under LGPL 3 or later. Therefore, this \n\
combined module is licensed under LGPL 3 or later.";
#else
static char gmpy_license[] = "\
The GMPY2 source code is licensed under LGPL 2.1 or later. This version \n\
of the MPIR library is licensed under LGPL 2.1 or later. Therefore, this \n\
combined module is licensed under LGPL 2.1 or later.";
#endif
#else
char gmpy_license[] = "\
The GMPY2 source code is licensed under LGPL 2.1 or later. The GMP \n\
library is licensed under LGPL 3 or later. Therefore, this combined \n\
module is licensed under LGPL 3 or later.";
#endif
#undef MPIR_VER

PyDoc_STRVAR(doc_license,
"license() -> string\n\n"
"Return string giving license information.");

static PyObject *
Pygmpy_get_license(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_license);
}

PyDoc_STRVAR(doc_version,
"version() -> string\n\n"
"Return string giving current GMPY2 version.");

static PyObject *
Pygmpy_get_version(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_version);
}

PyDoc_STRVAR(doc_cvsid,
"_cvsid() -> string\n\n"
"Return string giving current GMPY2 cvs Id.");

static PyObject *
Pygmpy_get_cvsid(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", _gmpy_cvs);
}

PyDoc_STRVAR(doc_mp_version,
"mp_version() -> string\n\n"
"Return string giving the name and version of the multiple precision\n"
"library used.");

static PyObject *
Pygmpy_get_mp_version(PyObject *self, PyObject *args)
{
#ifndef __MPIR_VERSION
    return PyUnicode_FromFormat("%s %s", "GMP", gmp_version);
#else
    return PyUnicode_FromFormat("%s %s", "MPIR", mpir_version);
#endif
}

PyDoc_STRVAR(doc_mpfr_version,
"mpfr_version() -> string\n\n"
"Return string giving current MPFR version.");

static PyObject *
Pygmpy_get_mpfr_version(PyObject *self, PyObject *args)
{
    return PyUnicode_FromFormat("%s %s", "MPFR",
                                MPFR_VERSION_STRING);
}

PyDoc_STRVAR(doc_mp_limbsize,
"mp_limbsize() -> integer\n\n\
Return the number of bits per limb.");

static PyObject *
Pygmpy_get_mp_limbsize(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", mp_bits_per_limb);
}

/*
 * access cache options
 */

PyDoc_STRVAR(doc_get_cache,
"get_cache() -> (cache_size, object_size)\n\n\
Return the current cache size (number of objects) and maximum size\n\
per object (number of limbs) for all GMPY2 objects.");

static PyObject *
Pygmpy_get_cache(PyObject *self, PyObject *args)
{
    return Py_BuildValue("ii", options.cache_size, options.cache_obsize);
}

PyDoc_STRVAR(doc_set_cache,
"set_cache(cache_size, object_size)\n\n\
Set the current cache size (number of objects) and the maximum size\n\
per object (number of limbs). Raises ValueError if cache size exceeds\n\
1000 or object size exceeds 16384.");

static PyObject *
Pygmpy_set_cache(PyObject *self, PyObject *args)
{
    int newcache, newsize;

    if (!PyArg_ParseTuple(args, "ii", &newcache, &newsize))
        return NULL;
    if (newcache<0 || newcache>MAX_CACHE) {
        VALUE_ERROR("cache size must between 0 and 1000");
        return NULL;
    }
    if (newsize<0 || newsize>MAX_CACHE_LIMBS) {
        VALUE_ERROR("object size must between 0 and 16384");
        return NULL;
    }

    options.cache_size = newcache;
    options.cache_obsize = newsize;
    set_zcache();
    set_pympzcache();
    set_pympqcache();
    set_pympfcache();
    set_pympccache();
    set_pyxmpzcache();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_set_debug,
"set_debug(boolean) -> boolean\n\n\
Set (if True) or clear (if False) the module level 'debug' setting\n\
and returns the previous value. If set, diagnostic information is\n\
sent to stderr. Note: only useful to debug GMPY2's own internals!");

static PyObject *
Pygmpy_set_debug(PyObject *self, PyObject *args)
{
#ifdef DEBUG
    int old = options.debug;
    if (!PyArg_ParseTuple(args, "i", &options.debug))
        return NULL;
    return Py_BuildValue("i", old);
#else
    PyErr_SetString(PyExc_NotImplementedError,
                    "gmpy2 was compiled without debug support.");
    return NULL;
#endif
}

PyDoc_STRVAR(doc_clear_underflow,
"clear_underflow()\n\n"
"Clear the underflow flag.");
static PyObject *
Pygmpy_clear_underflow(PyObject *self, PyObject *args)
{
    mpfr_clear_underflow();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_clear_overflow,
"clear_overflow()\n\n"
"Clear the overflow flag.");
static PyObject *
Pygmpy_clear_overflow(PyObject *self, PyObject *args)
{
    mpfr_clear_overflow();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_clear_nanflag,
"clear_nanflag()\n\n"
"Clear the Not A Number (nan) flag.");
static PyObject *
Pygmpy_clear_nanflag(PyObject *self, PyObject *args)
{
    mpfr_clear_nanflag();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_clear_inexflag,
"clear_inexactflag()\n\n"
"Clear the inexact flag.");
static PyObject *
Pygmpy_clear_inexflag(PyObject *self, PyObject *args)
{
    mpfr_clear_inexflag();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_clear_erangeflag,
"clear_erangeflag()\n\n"
"Clear the range error flag.");
static PyObject *
Pygmpy_clear_erangeflag(PyObject *self, PyObject *args)
{
    mpfr_clear_erangeflag();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_clear_flags,
"clear_flags()\n\n"
"Clear all exception flags.");
static PyObject *
Pygmpy_clear_flags(PyObject *self, PyObject *args)
{
    mpfr_clear_flags();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_set_underflow,
"set_underflow()\n\n"
"Set the underflow flag.");
static PyObject *
Pygmpy_set_underflow(PyObject *self, PyObject *args)
{
    mpfr_set_underflow();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_set_overflow,
"set_overflow()\n\n"
"Set the overflow flag.");
static PyObject *
Pygmpy_set_overflow(PyObject *self, PyObject *args)
{
    mpfr_set_overflow();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_set_nanflag,
"set_nanflag()\n\n"
"Set the Not A Number (nan) flag.");
static PyObject *
Pygmpy_set_nanflag(PyObject *self, PyObject *args)
{
    mpfr_set_nanflag();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_set_inexflag,
"set_inexactflag()\n\n"
"Set the inexact flag.");
static PyObject *
Pygmpy_set_inexflag(PyObject *self, PyObject *args)
{
    mpfr_set_inexflag();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_set_erangeflag,
"set_erangeflag()\n\n"
"Set the range error flag.");
static PyObject *
Pygmpy_set_erangeflag(PyObject *self, PyObject *args)
{
    mpfr_set_erangeflag();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_is_underflow,
"is_underflow() -> boolean\n\n"
"Return True if underflow flag is set.");
static PyObject *
Pygmpy_is_underflow(PyObject *self, PyObject *args)
{
    if (mpfr_underflow_p())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_overflow,
"is_overflow() -> boolean\n\n"
"Return True if overflow flag is set.");
static PyObject *
Pygmpy_is_overflow(PyObject *self, PyObject *args)
{
    if (mpfr_overflow_p())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_nanflag,
"is_nanflag() -> boolean\n\n"
"Return True if Not A Number (nan) flag is set.");
static PyObject *
Pygmpy_is_nanflag(PyObject *self, PyObject *args)
{
    if (mpfr_nanflag_p())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_inexflag,
"is_inexactflag() -> boolean\n\n"
"Return True if inexact flag is set.");
static PyObject *
Pygmpy_is_inexflag(PyObject *self, PyObject *args)
{
    if (mpfr_inexflag_p())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_erangeflag,
"is_erangeflag() -> boolean\n\n"
"Return True if range error flag is set.");
static PyObject *
Pygmpy_is_erangeflag(PyObject *self, PyObject *args)
{
    if (mpfr_erangeflag_p())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_get_mpfr_rc,
"get_mpfr_rc() -> integer\n\n\
Return the ternary result code from the most recent MPFR operation.\n"
"If the ternary value is 0, the result of the operation is exact.\n"
"If the ternary value is > 0, the result of the operation is greater\n"
"than the exact result. If the ternary value < 0, then the result\n"
"of the operation is less than the exact result.");
static PyObject *
Pygmpy_get_mpfr_rc(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", mpfr_rc);
}

PyDoc_STRVAR(doc_get_emin,
"get_emin() -> integer\n\n"
"Return the minimum exponent currently allowed.");
static PyObject *
Pygmpy_get_emin(PyObject *self, PyObject *args)
{
    return Py_BuildValue("n", (Py_ssize_t)mpfr_get_emin());
}

PyDoc_STRVAR(doc_get_emin_min,
"get_emin_min() -> integer\n\n"
"Return the minimum possible exponent that can be set.");
static PyObject *
Pygmpy_get_emin_min(PyObject *self, PyObject *args)
{
    return Py_BuildValue("n", (Py_ssize_t)mpfr_get_emin_min());
}

PyDoc_STRVAR(doc_set_emin,
"set_min(n)\n\n"
"Set the minimum allowed exponent.");
static PyObject *
Pygmpy_set_emin(PyObject *self, PyObject *args)
{
    Py_ssize_t exp;

    if (!PyArg_ParseTuple(args, "n", &exp))
        return NULL;
    if (mpfr_set_emin((mpfr_prec_t)exp)) {
        VALUE_ERROR("requested minimum exponent is invalid");
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_get_emax,
"get_emax() -> integer\n\n"
"Return the maximum exponent currently allowed.");
static PyObject *
Pygmpy_get_emax(PyObject *self, PyObject *args)
{
    return Py_BuildValue("n", (Py_ssize_t)mpfr_get_emax());
}

PyDoc_STRVAR(doc_get_emax_max,
"get_emax_max() -> integer\n\n"
"Return the maximum possible exponent that can be set.");
static PyObject *
Pygmpy_get_emax_max(PyObject *self, PyObject *args)
{
    return Py_BuildValue("n", (Py_ssize_t)mpfr_get_emax_max());
}

PyDoc_STRVAR(doc_set_emax,
"set_max(n)\n\n"
"Set the maximum allowed exponent.");
static PyObject *
Pygmpy_set_emax(PyObject *self, PyObject *args)
{
    Py_ssize_t exp;

    if (!PyArg_ParseTuple(args, "n", &exp))
        return NULL;
    if (mpfr_set_emax((mpfr_prec_t)exp)) {
        VALUE_ERROR("requested maximum exponent is invalid");
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_get_mode,
"get_mode() -> integer\n\n"
"Return the active mode for handling errors: ModePython raises\n"
"exception, ModeMPFR returns 'nan'.");
static PyObject *
Pygmpy_get_mode(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", options.raise);
}

PyDoc_STRVAR(doc_set_mode,
"set_mode(n)\n\n"
"Set the active mode for handling errors: ModePython raises\n"
"exception, ModeMPFR returns 'nan'.");
static PyObject *
Pygmpy_set_mode(PyObject *self, PyObject *args)
{
    int mode;

    if (!PyArg_ParseTuple(args, "i", &mode))
        return NULL;
    if (mode == GMPY_MODE_PYTHON)
        options.raise = GMPY_MODE_PYTHON;
    else if (mode == GMPY_MODE_MPFR)
        options.raise = GMPY_MODE_MPFR;
    else {
        VALUE_ERROR("invalid value for error handling mode");
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_get_precision,
"get_precision() -> integer\n\n"
"Return the number of bits of precision used for calculations.");
static PyObject *
Pygmpy_get_precision(PyObject *self, PyObject *args)
{
    return Py_BuildValue("n", options.precision);
}

PyDoc_STRVAR(doc_get_max_precision,
"get_max_precision() -> integer\n\n"
"Return the maximum bits of precision that can be used for calculations.\n"
"Note: to allow extra precision for intermediate calculations, avoid\n"
"setting precision close the maximum precisicon");
static PyObject *
Pygmpy_get_max_precision(PyObject *self, PyObject *args)
{
    return Py_BuildValue("n", MPFR_PREC_MAX);
}

PyDoc_STRVAR(doc_set_precision,
"set_precision(n)\n\n"
"Set the number of bits of precision to use for calculations.");
static PyObject *
Pygmpy_set_precision(PyObject *self, PyObject *args)
{
    int bits;

    if(!PyArg_ParseTuple(args, "i", &bits))
        return NULL;
    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    options.precision = bits;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_get_mpfr_round,
"get_mpfr_round() -> integer\n\n"
"Return the rounding mode for 'mpf' arithmetic. Rounding mode can"
"be one of RoundToNearest, RoundToZero, RoundUp, RoundDown, or"
"RoundAwayZero.");

static PyObject *
Pygmpy_get_mpfr_round(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", options.mpfr_round);
}

PyDoc_STRVAR(doc_set_mpfr_round,
"set_mpfr_round(n)\n\n"
"Set the rounding mode for 'mpf' arithmetic. Valid rounding modes"
"are RoundToNearest, RoundToZero, RoundUp, RoundDown, or"
"RoundAwayZero.");

static PyObject *
Pygmpy_set_mpfr_round(PyObject *self, PyObject *args)
{
    int mode;

    if(!PyArg_ParseTuple(args, "i", &mode))
        return NULL;

    if (mode == MPFR_RNDN)
        options.mpfr_round = mode;
    else if (mode == MPFR_RNDZ)
        options.mpfr_round = mode;
    else if (mode == MPFR_RNDU)
        options.mpfr_round = mode;
    else if (mode == MPFR_RNDD)
        options.mpfr_round = mode;
    else if (mode == MPFR_RNDA)
        options.mpfr_round = mode;
    else {
        VALUE_ERROR("invalid rounding mode");
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_get_mpc_round,
"get_mpc_round() -> (integer, integer)\n\n"
"Return the rounding mode for 'mpc' arithmetic. The first value is"
"the rounding mode for the real portion. The second value is the"
"rounding mode for the imaginary portion. Rounding mode can be"
"RoundToNearest, RoundToZero, RoundUp, or RoundDown.");

static PyObject *
Pygmpy_get_mpc_round(PyObject *self, PyObject *args)
{
    return Py_BuildValue("ii",
                         MPC_RND_RE(options.mpc_round),
                         MPC_RND_IM(options.mpc_round));
}

PyDoc_STRVAR(doc_set_mpc_round,
"set_mpc_round(n,n)\n\n"
"Set the rounding mode for complex arithmetic. The first value is"
"the rounding more for the real portion. The second value is the"
"rounding mode for the imaginary portion. Valid rounding modes are"
"RoundToNearest, RoundToZero, RoundUp, or RoundDown.");

static PyObject *
Pygmpy_set_mpc_round(PyObject *self, PyObject *args)
{
    int rmode, imode;

    if(!PyArg_ParseTuple(args, "ii", &rmode, &imode))
        return NULL;

    if (rmode == MPFR_RNDN && imode == MPFR_RNDN)
        options.mpc_round = MPC_RNDNN;
    else if (rmode == MPFR_RNDN && imode == MPFR_RNDZ)
        options.mpc_round = MPC_RNDNZ;
    else if (rmode == MPFR_RNDN && imode == MPFR_RNDU)
        options.mpc_round = MPC_RNDNU;
    else if (rmode == MPFR_RNDN && imode == MPFR_RNDD)
        options.mpc_round = MPC_RNDND;

    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDN)
        options.mpc_round = MPC_RNDZN;
    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDZ)
        options.mpc_round = MPC_RNDZZ;
    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDU)
        options.mpc_round = MPC_RNDZU;
    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDD)
        options.mpc_round = MPC_RNDZD;

    else if (rmode == MPFR_RNDU && imode == MPFR_RNDN)
        options.mpc_round = MPC_RNDUN;
    else if (rmode == MPFR_RNDU && imode == MPFR_RNDZ)
        options.mpc_round = MPC_RNDUZ;
    else if (rmode == MPFR_RNDU && imode == MPFR_RNDU)
        options.mpc_round = MPC_RNDUU;
    else if (rmode == MPFR_RNDU && imode == MPFR_RNDD)
        options.mpc_round = MPC_RNDUD;

    else if (rmode == MPFR_RNDD && imode == MPFR_RNDN)
        options.mpc_round = MPC_RNDDN;
    else if (rmode == MPFR_RNDD && imode == MPFR_RNDZ)
        options.mpc_round = MPC_RNDDZ;
    else if (rmode == MPFR_RNDD && imode == MPFR_RNDU)
        options.mpc_round = MPC_RNDDU;
    else if (rmode == MPFR_RNDD && imode == MPFR_RNDD)
        options.mpc_round = MPC_RNDDD;

    else {
        VALUE_ERROR("invalid rounding mode");
        return NULL;
    }
    Py_RETURN_NONE;
}

/* create a copy of a gmpy2 object */
PyDoc_STRVAR(doc_copym,
"x._copy() -> gmpy2_object\n\n"
"Return a copy of x.\n");
PyDoc_STRVAR(doc_copyg,
"_copy(x): -> gmpy2_object\n\n"
"Return a copy of x. Raises TypeError if x is not a gmpy2 object.");
static PyObject *
Pympany_copy(PyObject *self, PyObject *other)
{
    if (self && Pympz_Check(self))
        return (PyObject*)Pympz2Pympz(self);
    else if (self && Pyxmpz_Check(self))
        return (PyObject*)Pyxmpz2Pyxmpz(self);
    else if (self && Pympq_Check(self))
        return (PyObject*)Pympq2Pympq(self);
    else if (self && Pympf_Check(self))
        return (PyObject*)Pympf2Pympf(self, 0);
    else if (Pympz_Check(other))
        return (PyObject*)Pympz2Pympz(other);
    else if (Pyxmpz_Check(other))
        return (PyObject*)Pyxmpz2Pyxmpz(other);
    else if (Pympq_Check(other))
        return (PyObject*)Pympq2Pympq(other);
    else if (Pympf_Check(other))
        return (PyObject*)Pympf2Pympf(other, 0);
    TYPE_ERROR("_copy() requires a gmpy2 object as argument");
    return NULL;
}

PyDoc_STRVAR(doc_binarym,
"x.binary() -> binary string\n\n"
"Return a Python string (or bytes for Python 3+) that is a portable\n"
"binary representation of a gmpy2 object x. The binary string can\n"
"later be passed to the appropriate constructor function to obtain\n"
"an exact copy of x's value.");
PyDoc_STRVAR(doc_binaryg,
"binary(x) -> binary string\n\n"
"Return a Python string (or bytes for Python 3+) that is a portable\n"
"binary representation of a gmpy2 object x. The binary string can\n"
"later be passed to the appropriate constructor function to obtain\n"
"an exact copy of x's value. Raises TypeError if x is not a gmpy2\n"
"object.");

static PyObject *
Pympany_binary(PyObject *self, PyObject *other)
{
    if(self && Pympz_Check(self))
        return Pympz2binary((PympzObject*)self);
    else if(self && Pyxmpz_Check(self))
        return Pyxmpz2binary((PyxmpzObject*)self);
    else if(self && Pympq_Check(self))
        return Pympq2binary((PympqObject*)self);
    else if(self && Pympf_Check(self))
        return Pympf2binary((PympfObject*)self);
    else if(Pympz_Check(other))
        return Pympz2binary((PympzObject*)other);
    else if(Pyxmpz_Check(other))
        return Pyxmpz2binary((PyxmpzObject*)other);
    else if(Pympq_Check(other))
        return Pympq2binary((PympqObject*)other);
    else if(Pympf_Check(other))
        return Pympf2binary((PympfObject*)other);
    TYPE_ERROR("binary() requires a gmpy2 object as argument");
    return NULL;
}
