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
"license() -> string\n\n\
Return string giving license information.");

static PyObject *
Pygmpy_get_license(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_license);
}

PyDoc_STRVAR(doc_version,
"version() -> string\n\n\
Return string giving current GMPY2 version.");

static PyObject *
Pygmpy_get_version(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_version);
}

PyDoc_STRVAR(doc_cvsid,
"_cvsid() -> string\n\n\
Return string giving current GMPY2 cvs Id.");

static PyObject *
Pygmpy_get_cvsid(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", _gmpy_cvs);
}

PyDoc_STRVAR(doc_gmp_version,
"gmp_version() -> string\n\n\
Return string giving current GMP version. Empty string returned if\n\
MPIR is used.");

static PyObject *
Pygmpy_get_gmp_version(PyObject *self, PyObject *args)
{
#ifndef __MPIR_VERSION
    return Py_BuildValue("s", gmp_version);
#else
    return Py_BuildValue("s", "");
#endif
}

PyDoc_STRVAR(doc_mpir_version,
"mpir_version() -> string\n\n\
Return string giving current MPIR version. Empty string is returned if\n\
GMP was used.");

static PyObject *
Pygmpy_get_mpir_version(PyObject *self, PyObject *args)
{
#ifdef __MPIR_VERSION
    return Py_BuildValue("s", mpir_version);
#else
    return Py_BuildValue("s", "");
#endif
}

PyDoc_STRVAR(doc_gmp_limbsize,
"gmp_limbsize() -> integer\n\n\
Return the number of bits per limb.");

static PyObject *
Pygmpy_get_gmp_limbsize(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", GMP_NUMB_BITS);
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
    set_qcache();
    set_pympzcache();
    set_pympqcache();
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
    long old = options.debug;

    if (!PyArg_ParseTuple(args, "l", &options.debug))
        return NULL;
    return Py_BuildValue("l", old);
}

PyDoc_STRVAR(doc_set_prefer_mutable,
"set_prefer_mutable(boolean) -> boolean\n\n\
If set to True, the result type of an ambiguous operation (i.e.\n\
adding an 'mpz' and 'xmpz', or creating a GMPY2 type from a Python\n\
number) will be an 'xmpz'. Returns the previous value. Default is\n\
is False.");

static PyObject *
Pygmpy_set_prefer_mutable(PyObject *self, PyObject *args)
{
    long old = options.prefer_mutable;

    if(!PyArg_ParseTuple(args, "l", &options.prefer_mutable))
        return NULL;
    return Py_BuildValue("l", old);
}

static char doc_set_tagoff[]="\
set_tagoff(n): resets (if n==0) or sets (if n!=0) the module\n\
level 'tagoff' setting, removing the 'gmpy2.' prefix of the tag\n\
strings used by repr and (optionally) digits/fdigits/qdigits;\n\
also returns the previous value of this module-level setting.\n\
";
static PyObject *
Pygmpy_set_tagoff(PyObject *self, PyObject *args)
{
    int old = options.tagoff;

    if(!PyArg_ParseTuple(args, "i", &options.tagoff))
        return NULL;
    if(options.tagoff)
        options.tagoff = GMPY2_TAGOFF;
    return Py_BuildValue("i", old!=0);
}

static char doc_set_minprec[]="\
set_minprec(n): sets number of bits of precision to be at\n\
least n for all mpf objects generated from now on; also\n\
returns the previous value of this module-level setting.\n\
";
static PyObject *
Pygmpy_set_minprec(PyObject *self, PyObject *args)
{
    long old = options.minprec;
    int i;

    if(!PyArg_ParseTuple(args, "i", &i))
        return NULL;
    if(i<0) {
        VALUE_ERROR("minimum precision must be >= 0");
        return NULL;
    }
    options.minprec = i;
    return Py_BuildValue("l", old);
}

static char doc_set_fcoform[]="\
set_fcoform(s=None): resets (if s is None) or sets the module level\n\
'fcoform' setting, the format in which to build an intermediate string\n\
to be used in float->mpf conversion (direct, if no fcoform); also\n\
returns the previous value of this module-level setting.  Note that\n\
s must be a string usable for s%f formatting; or, s may be a Python\n\
int, 0<s<=30, in which case format is set to '%.<s>e'.\n\
";
static PyObject *
Pygmpy_set_fcoform(PyObject *self, PyObject *args)
{
    PyObject *old = options.fcoform;
    PyObject *new = 0;
    long inew;

    if(!PyArg_ParseTuple(args, "|O", &new))
        return NULL;
    if(new == Py_None) { /* none == missing-argument (reset string use) */
        new = 0;
    } else if(new) {
        char buf[20];
        if(isInteger(new)) {
            /* int arg (1 to 30) used as # of digits for intermediate string */
            inew = clong_From_Integer(new);
            if(inew==-1 && PyErr_Occurred()) {
                VALUE_ERROR("number of digits n must be 0<n<=30");
                return NULL;
            }
            /* check range for number-of-digits setting */
            if(inew<=0 || inew>30) {
                VALUE_ERROR("number of digits n must be 0<n<=30");
                return NULL;
            }
            /* prepare Python format-string '%.12e' or whatever */
            sprintf(buf,"%%.%lde",inew);
#ifdef PY3
            new = PyUnicode_FromString(buf);
#else
            new = PyString_FromString(buf);
#endif
        } else { /* else arg must be string directly usable in formatting */
            if(!Py2or3String_Check(new)) {
                TYPE_ERROR("set_fcoform argument must be int, string, or None");
                return NULL;
            }
            Py_INCREF(new);
        }
    }
    /* set new 'float conversion format' and return old one if any */
    options.fcoform = new;
    if(old)
        return old;
    else
        return Py_BuildValue("");
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
