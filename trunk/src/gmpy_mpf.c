/* gmpy_mpf.c
 *
 * Functions that operate strictly on mpf.
 *
 * This file should be considered part of gmpy2.c
 */

static PyObject *
Pympf_getprec2(PympfObject *self, void *closure)
{
    return PyIntOrLong_FromSize_t((size_t)mpfr_get_prec(self->f));
}

static int
Pympf_nonzero(PympfObject *x)
{
    return mpfr_sgn(x->f) != 0;
}

static long
Pympf_hash(PympfObject *self)
{
#ifdef _PyHASH_MODULUS
    unsigned long hash = 0;
    long exp;
    size_t msize;
    int sign;

    if (self->hash_cache != -1)
        return self->hash_cache;

    /* Handle special cases first */
    if (!mpfr_number_p(self->f)) {
        if (mpfr_inf_p(self->f))
            if (mpfr_sgn(self->f) > 0)
                return (self->hash_cache = _PyHASH_INF);
            else
                return (self->hash_cache = -_PyHASH_INF);
        else
            return (self->hash_cache = _PyHASH_NAN);
    }

    /* Calculate the number of limbs in the mantissa. */
    msize = (self->f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;

    /* Calculate the hash of the mantissa. */
    if (mpfr_sgn(self->f) > 0) {
        hash = mpn_mod_1(self->f->_mpfr_d, msize, _PyHASH_MODULUS);
        sign = 1;
    }
    else if (mpfr_sgn(self->f) < 0) {
        hash = mpn_mod_1(self->f->_mpfr_d, msize, _PyHASH_MODULUS);
        sign = -1;
    }
    else {
        return (self->hash_cache = 0);
    }

    /* Calculate the final hash. */
    exp = self->f->_mpfr_exp - (msize * mp_bits_per_limb);
    exp = exp >= 0 ? exp % _PyHASH_BITS : _PyHASH_BITS-1-((-1-exp) % _PyHASH_BITS);
    hash = ((hash << exp) & _PyHASH_MODULUS) | hash >> (_PyHASH_BITS - exp);

    hash *= sign;
    if (hash == (unsigned long)-1)
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

/* float-truncations (return still a float!) */

static char doc_ceilm[]="\
x.ceil(): returns an mpf that is the smallest integer >= x\n\
";
static char doc_ceilg[]="\
ceil(x): returns an mpf that is the smallest integer >= x\n\
x must be an mpf, or else gets coerced to one.\n\
";
static PyObject *
Pympf_ceil(PyObject *self, PyObject *args)
{
    PympfObject *result;

    PARSE_ONE_MPF("ceil() requires 'mpf' argument");

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_ceil(result->f, Pympf_AS_MPF(self));
    Py_DECREF(self);
    return (PyObject*)result;
}

static char doc_floorm[]="\
x.floor(): returns an mpf that is the smallest integer <= x\n\
";
static char doc_floorg[]="\
floor(x): returns an mpf that is the smallest integer <= x\n\
x must be an mpf, or else gets coerced to one.\n\
";
static PyObject *
Pympf_floor(PyObject *self, PyObject *args)
{
    PympfObject *result;

    PARSE_ONE_MPF("ceil() requires 'mpf' argument");

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_floor(result->f, Pympf_AS_MPF(self));
    Py_DECREF(self);
    return (PyObject*)result;
}

static char doc_truncm[]="\
x.trunc(): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
";
static char doc_truncg[]="\
trunc(x): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
x must be an mpf, or else gets coerced to one.\n\
";
static PyObject *
Pympf_trunc(PyObject *self, PyObject *args)
{
    PympfObject *result;

    PARSE_ONE_MPF("ceil() requires 'mpf' argument");

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_trunc(result->f, Pympf_AS_MPF(self));
    Py_DECREF(self);
    return (PyObject*)result;
}

static char doc_pi[]="\
pi(n): returns pi with n bits of precision in an mpf object\n\
";
static PyObject *
Pygmpy_pi(PyObject *self, PyObject *args)
{
    PympfObject *pi;
    long bits;

    if (PyTuple_GET_SIZE(args) == 0)
        bits = options.precision;
    else if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("pi() requires 0 or 1 arguments");
        return NULL;
    }
    else {
        bits = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (bits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("pi() requires 'int' argument");
            return NULL;
        }
        if (bits < 0) {
            VALUE_ERROR("precision must be >= 0");
            return NULL;
        }
    }

    if (!(pi = Pympf_new(bits)))
        return NULL;
    mpfr_const_pi(pi->f, options.rounding);
    return (PyObject*)pi;
}

static char doc_const_pi[]="\
const_pi(): returns the constant pi using default precision\n\
";
static PyObject *
Pygmpy_const_pi(PyObject *self, PyObject *args)
{
    PympfObject *result;

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_const_pi(result->f, options.rounding);
    return (PyObject*)result;
}

static char doc_const_euler[]="\
const_euler(): returns the euler constant using default precision\n\
";
static PyObject *
Pygmpy_const_euler(PyObject *self, PyObject *args)
{
    PympfObject *result;

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_const_euler(result->f, options.rounding);
    return (PyObject*)result;
}

static char doc_const_log2[]="\
const_log2(): returns the log2 constant using default precision\n\
";
static PyObject *
Pygmpy_const_log2(PyObject *self, PyObject *args)
{
    PympfObject *result;

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_const_log2(result->f, options.rounding);
    return (PyObject*)result;
}

static char doc_const_catalan[]="\
const_catalan(): returns the catalan constant using default precision\n\
";
static PyObject *
Pygmpy_const_catalan(PyObject *self, PyObject *args)
{
    PympfObject *result;

    if (!(result = Pympf_new(0)))
        return NULL;
    mpfr_const_catalan(result->f, options.rounding);
    return (PyObject*)result;
}

static char doc_fsqrtm[]="\
x.sqrt(): returns the square root of x.  x must be >= 0.\n\
";
static char doc_fsqrtg[]="\
fsqrt(x): returns the square root of x.  x must be an mpf, or\n\
else gets coerced to one; further, x must be >= 0.\n\
";
static PyObject *
Pympf_sqrt(PyObject *self, PyObject *args)
{
    PympfObject *result;

    PARSE_ONE_MPF("sqrt() requires 'mpf' argument");

    if (mpfr_sgn(Pympf_AS_MPF(self)) < 0) {
        VALUE_ERROR("sqrt() of negative number");
        Py_DECREF(self);
        return NULL;
    }

    if (!(result = Pympf_new(0))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_sqrt(result->f, Pympf_AS_MPF(self), options.rounding);
    Py_DECREF(self);
    return (PyObject*)result;
}

static char doc_froundm[] = "\
x.round(n): returns x rounded to n bits. Uses default precision if\n\
n is not specified.\n\
";
static char doc_froundg[] = "\
fround(x, n): returns x rounded to n bits. Uses default precision\n\
if n is not specified.\n\
";
static PyObject *
Pympf_round(PyObject *self, PyObject *args)
{
    mpfr_prec_t prec = options.precision;
    PympfObject *result;

    PARSE_ONE_MPF_OPT_CLONG(&prec,
            "round() requires 'mpf',['int'] arguments");
    if (!(result = Pympf_new(prec))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_set(result->f, Pympf_AS_MPF(self), options.rounding);
    Py_DECREF(self);
    return (PyObject*)result;
}

static char doc_reldiffm[] = "\
x.reldiff(y): returns the relative difference between x and y,\n\
where y can be any number and gets coerced to an mpf; result is\n\
an mpf roughly equal to abs(x-y)/x.\n\
";
static char doc_reldiffg[] = "\
reldiff(x,y): returns the relative difference between x and y,\n\
where x and y can be any numbers and get coerced to mpf; result is\n\
an mpf roughly equal to abs(x-y)/x.\n\
";
static PyObject *
Pympf_doreldiff(PyObject *self, PyObject *args)
{
    PympfObject *result;
    PyObject *other;

    PARSE_TWO_MPF(other, "reldiff() requires 'mpf,'mpf' arguments");

    if (!(result = Pympf_new(0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_reldiff(result->f, Pympf_AS_MPF(self), Pympf_AS_MPF(other),
                options.rounding);
    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
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

