/* gmpy_misc.c
 *
 * Miscellaneous module-level functions and helper functions.
 *
 * This file should be considered part of gmpy2.c
 */

/* Return license information. */

char gmpy_license[] = "\
The GMPY2 source code is licensed under LGPL 2.1 or later. The GMP/MPIR, \n\
MPFR, and MPC libraries are licensed under LGPL 3 or later. Therefore, this \n\
combined module is licensed under LGPL 3 or later.";

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

PyDoc_STRVAR(doc_mpc_version,
"mpc_version() -> string\n\n"
"Return string giving current MPC version.");

static PyObject *
Pygmpy_get_mpc_version(PyObject *self, PyObject *args)
{
    return PyUnicode_FromFormat("%s %s", "MPC",
                                MPC_VERSION_STRING);
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
    return Py_BuildValue("ii", global.cache_size, global.cache_obsize);
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

    global.cache_size = newcache;
    global.cache_obsize = newsize;
    set_zcache();
    set_pympzcache();
    set_pympqcache();
    set_pympfrcache();
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
    int old = global.debug;
    if (!PyArg_ParseTuple(args, "i", &global.debug))
        return NULL;
    return Py_BuildValue("i", old);
#else
    PyErr_SetString(PyExc_NotImplementedError,
                    "gmpy2 was compiled without debug support.");
    return NULL;
#endif
}

PyDoc_STRVAR(doc_get_mode,
"get_mode() -> integer\n\n"
"Return the active mode for handling errors: ModePython raises\n"
"exception, ModeMPFR returns 'nan' or 'inf'.");

static PyObject *
Pygmpy_get_mode(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", context.raise);
}

PyDoc_STRVAR(doc_set_mode,
"set_mode(n)\n\n"
"Set the active mode for handling errors: ModePython raises\n"
"exception, ModeMPFR returns 'nan' or 'inf'.");

static PyObject *
Pygmpy_set_mode(PyObject *self, PyObject *args)
{
    int mode;

    if (!PyArg_ParseTuple(args, "i", &mode))
        return NULL;
    if (mode == GMPY_MODE_RAISE)
        context.raise = GMPY_MODE_RAISE;
    else if (mode == GMPY_MODE_NONSTOP)
        context.raise = GMPY_MODE_NONSTOP;
    else {
        VALUE_ERROR("invalid value for error handling mode");
        return NULL;
    }
    Py_RETURN_NONE;
}

/*
 * Helper functions.
 */

/* Verify that a valid rounding mode is specified for complex arithmetic.
 * Returns 0 (false) if the rounding mode is not valid else returns 1 (true).
 */

static int
Pymisc_verify_mpc_round(int rmode)
{
    if ( rmode == MPC_RNDNN || rmode == MPC_RNDNZ ||
         rmode == MPC_RNDNU || rmode == MPC_RNDND ||
         rmode == MPC_RNDZN || rmode == MPC_RNDZZ ||
         rmode == MPC_RNDZU || rmode == MPC_RNDZD ||
         rmode == MPC_RNDUN || rmode == MPC_RNDUZ ||
         rmode == MPC_RNDUU || rmode == MPC_RNDUD ||
         rmode == MPC_RNDDN || rmode == MPC_RNDDZ ||
         rmode == MPC_RNDDU || rmode == MPC_RNDDD )
        return 1;
    else
        return 0;
}

/* Verify that valid precisions are requested for complex arithmetic.
 * Returns 0 if the precisions are not valid else returns 1.
 */

static int
Pymisc_verify_mpc_precision(Py_ssize_t rprec, Py_ssize_t iprec)
{
    if ( rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
         iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX )
        return 0;
    else
        return 1;
}

