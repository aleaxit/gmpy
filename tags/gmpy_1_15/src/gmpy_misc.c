/* gmpy_misc.c
 *
 * Miscellaneous module-level functions.
 *
 * This file should be considered part of gmpy.c.
 */

/* Return license information. */
static char doc_license[]="\
license(): returns string giving license information\n\
";
static PyObject *
Pygmpy_get_license(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("license expects 0 arguments");
    return Py_BuildValue("s", gmpy_license);
}

/* return GMPY, resp. GMP, versions, or CVS Id, as strings */
static char doc_version[]="\
version(): returns string giving current GMPY version\n\
";
static PyObject *
Pygmpy_get_version(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("version expects 0 arguments");
    return Py_BuildValue("s", gmpy_version);
}

static char doc_cvsid[]="\
_cvsid(): returns string giving current GMPY cvs Id\n\
";
static PyObject *
Pygmpy_get_cvsid(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("get_cvsid expects 0 arguments");
    return Py_BuildValue("s", _gmpy_cvs);
}

static char doc_gmp_version[]="\
gmp_version(): returns string giving current GMP version. Empty string\n\
returned if MPIR was used.\n\
";
static PyObject *
Pygmpy_get_gmp_version(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("gmp_version expects 0 arguments");
#ifndef __MPIR_VERSION
    return Py_BuildValue("s", gmp_version);
#else
    return Py_BuildValue("s", "");
#endif
}

static char doc_mpir_version[]="\
mpir_version(): returns string giving current MPIR version. Empty string\n\
returned if GMP was used.\n\
";
static PyObject *
Pygmpy_get_mpir_version(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("mpir_version expects 0 arguments");
#ifdef __MPIR_VERSION
    return Py_BuildValue("s", mpir_version);
#else
    return Py_BuildValue("s", "");
#endif
}

static char doc_gmp_limbsize[]="\
gmp_limbsize(): returns the number of bits per limb\n\
";
static PyObject *
Pygmpy_get_gmp_limbsize(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("gmp_limbsize expects 0 arguments");
    return Py_BuildValue("i", GMP_NUMB_BITS);
}

/*
 * access cache options
 */
static char doc_get_cache[]="\
get_cache(): returns the current cache-size (number of objects)\n\
and maximum size per object (number of limbs) for all objects.\n\
";
static PyObject *
Pygmpy_get_cache(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("get_cache expects 0 arguments");
    return Py_BuildValue("ii", options.cache_size, options.cache_obsize);
}

static char doc_set_cache[]="\
set_cache(n,size): sets the current cache-size (number of objects) to\n\
'n' and the maximum size per object (number of limbs) to 'size'.\n\
Note: cache size 'n' must be between 0 and 1000, included. Object size\n\
'size' must be between 0 and 16384, included.\n\
";
static PyObject *
Pygmpy_set_cache(PyObject *self, PyObject *args)
{
    int newcache, newsize;
    if(!PyArg_ParseTuple(args, "ii", &newcache, &newsize))
        return NULL;
    if(newcache<0 || newcache>MAX_CACHE) {
        PyErr_SetString(PyExc_ValueError, "cache must between 0 and 1000");
        return NULL;
    }
    if(newsize<0 || newsize>MAX_CACHE_LIMBS) {
        PyErr_SetString(PyExc_ValueError, "object size must between 0 and 16384");
        return NULL;
    }
    options.cache_size=newcache;
    options.cache_obsize=newsize;
    set_zcache();
    set_qcache();
    set_fcache();
    set_pympzcache();
    return Py_BuildValue("");
}

/* set a module-global flag, return previously-set value */
static char doc_set_debug[]="\
set_debug(n): resets (if n==0) or sets (if n!=0) the module\n\
level 'debug' setting, giving detailed info to stderr; also\n\
returns the previous value of this module-level setting.\n\
Note: only useful to debug gmpy's own internals!\n\
";
static PyObject *
Pygmpy_set_debug(PyObject *self, PyObject *args)
{
    int old = options.debug;

    ONE_ARG("set_debug", "i", &options.debug);
    return Py_BuildValue("i", old);
}

static char doc_set_tagoff[]="\
set_tagoff(n): resets (if n==0) or sets (if n!=0) the module\n\
level 'tagoff' setting, removing the 'gmpy.' prefix of the tag\n\
strings used by repr and (optionally) digits/fdigits/qdigits;\n\
also returns the previous value of this module-level setting.\n\
";
static PyObject *
Pygmpy_set_tagoff(PyObject *self, PyObject *args)
{
    int old = options.tagoff;

    ONE_ARG("set_tagoff", "i", &options.tagoff);
    if(options.tagoff) options.tagoff=5;
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
    long i;

    ONE_ARG("set_minprec", "l", &i);
    if(i<0) {
        PyErr_SetString(PyExc_ValueError,
            "minimum precision must be >= 0");
        return 0;
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

    ONE_ARG("set_fcoform", "|O", &new);
    if(new == Py_None) { /* none == missing-argument (reset string use) */
        new = 0;
    } else if(new) {
        char buf[20];
        if(isInteger(new)) {
            /* int arg (1 to 30) used as # of digits for intermediate string */
            inew = clong_From_Integer(new);
            if(inew==-1 && PyErr_Occurred()) {
                PyErr_SetString(PyExc_ValueError,
                    "number of digits n must be 0<n<=30");
                return 0;
            }
            /* check range for number-of-digits setting */
            if(inew<=0 || inew>30) {
                PyErr_SetString(PyExc_ValueError,
                    "number of digits n must be 0<n<=30");
                return 0;
            }
            /* prepare Python format-string '%.12e' or whatever */
            sprintf(buf,"%%.%lde",inew);
            new = Py2or3String_FromString(buf);
        } else { /* else arg must be string directly usable in formatting */
            if(!Py2or3String_Check(new)) {
                PyErr_SetString(PyExc_TypeError,
                    "set_fcoform argument must be int, string, or None");
                return 0;
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
