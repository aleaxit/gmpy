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
    return Py_BuildValue("i", global.raise);
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
    if (mode == GMPY_MODE_PYTHON)
        global.raise = GMPY_MODE_PYTHON;
    else if (mode == GMPY_MODE_MPFR)
        global.raise = GMPY_MODE_MPFR;
    else {
        VALUE_ERROR("invalid value for error handling mode");
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
    else if (self && Pympfr_Check(self))
        return (PyObject*)Pympfr2Pympfr(self, 0);
    else if (Pympz_Check(other))
        return (PyObject*)Pympz2Pympz(other);
    else if (Pyxmpz_Check(other))
        return (PyObject*)Pyxmpz2Pyxmpz(other);
    else if (Pympq_Check(other))
        return (PyObject*)Pympq2Pympq(other);
    else if (Pympfr_Check(other))
        return (PyObject*)Pympfr2Pympfr(other, 0);
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
    else if(self && Pympfr_Check(self))
        return Pympfr2binary((PympfrObject*)self);
    else if(Pympz_Check(other))
        return Pympz2binary((PympzObject*)other);
    else if(Pyxmpz_Check(other))
        return Pyxmpz2binary((PyxmpzObject*)other);
    else if(Pympq_Check(other))
        return Pympq2binary((PympqObject*)other);
    else if(Pympfr_Check(other))
        return Pympfr2binary((PympfrObject*)other);
    TYPE_ERROR("binary() requires a gmpy2 object as argument");
    return NULL;
}
