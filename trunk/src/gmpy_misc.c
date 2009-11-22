/* gmpy_misc.c
 *
 * Miscellaneous module-level functions.
 *
 * This file should be considered part of gmpy.c
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
static char doc_get_zcache[]="\
get_zcache(): returns the current cache-size (number of objects)\n\
for mpz objects.\n\
";
static PyObject *
Pygmpy_get_zcache(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("get_zcache expects 0 arguments");
    return Py_BuildValue("i", options.zcache);
}

static char doc_get_qcache[]="\
get_qcache(): returns the current cache-size (number of objects)\n\
for mpq objects.\n\
";
static PyObject *
Pygmpy_get_qcache(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("get_qcache expects 0 arguments");
    return Py_BuildValue("i", options.qcache);
}

static char doc_get_fcache[]="\
get_fcache(): returns the current cache-size (number of objects)\n\
for mpf objects.\n\
";
static PyObject *
Pygmpy_get_fcache(PyObject *self, PyObject *args)
{
    PARSE_NO_ARGS("get_fcache expects 0 arguments");
    return Py_BuildValue("i", options.fcache);
}

static char doc_set_zcache[]="\
set_zcache(n): sets the current cache-size (number of objects)\n\
for mpz objects to n (does not immediately flush or enlarge the\n\
cache, but rather lets it grow/shrink during later normal use).\n\
Note: cache size n must be between 0 and 1000, included.\n\
";
static PyObject *
Pygmpy_set_zcache(PyObject *self, PyObject *args)
{
    int newval;
    ONE_ARG("set_zcache", "i", &newval);
    if(newval<0 || newval>MAX_CACHE) {
        PyErr_SetString(PyExc_ValueError, "cache must between 0 and 1000");
        return 0;
    }
    set_zcache(newval);
    return Py_BuildValue("");
}

static char doc_set_qcache[]="\
set_qcache(n): sets the current cache-size (number of objects)\n\
for mpq objects to n (does not immediately flush or enlarge the\n\
cache, but rather lets it grow/shrink during later normal use).\n\
Note: cache size n must be between 0 and 1000, included.\n\
";
static PyObject *
Pygmpy_set_qcache(PyObject *self, PyObject *args)
{
    int newval;
    ONE_ARG("set_qcache", "i", &newval);
    if(newval<0 || newval>MAX_CACHE) {
        PyErr_SetString(PyExc_ValueError, "cache must between 0 and 1000");
        return 0;
    }
    set_qcache(newval);
    return Py_BuildValue("");
}

static char doc_set_fcache[]="\
set_fcache(n): sets the current cache-size (number of objects)\n\
for mpf objects to n (does not immediately flush or enlarge the\n\
cache, but rather lets it grow/shrink during later normal use).\n\
Note: cache size n must be between 0 and 1000, included.\n\
";
static PyObject *
Pygmpy_set_fcache(PyObject *self, PyObject *args)
{
    int newval;
    ONE_ARG("set_fcache", "i", &newval);
    if(newval<0 || newval>MAX_CACHE) {
        PyErr_SetString(PyExc_ValueError, "cache must between 0 and 1000");
        return 0;
    }
    set_fcache(newval);
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
    long old = options.debug;

    ONE_ARG("set_debug", "l", &options.debug);
    return Py_BuildValue("l", old);
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
    int i;

    ONE_ARG("set_minprec", "i", &i);
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
