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


static char doc_pi[]="\
pi(n): returns pi with n bits of precision in an mpf object\n\
";
static PyObject *
Pygmpy_pi(PyObject *self, PyObject *args)
{
    PympfObject *pi;
    mpfr_prec_t bits;

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
Pympf_sqrt(PyObject *self, PyObject *other)
{
    PympfObject *result, *tempx;

    if (!(result = Pympf_new(0)))
        return NULL;
    if(self && Pympf_Check(self)) {
        if (mpfr_sgn(Pympf_AS_MPF(self)) < 0) {
            VALUE_ERROR("sqrt() of negative number");
            return NULL;
        }
        mpfr_sqrt(result->f, Pympf_AS_MPF(self), options.rounding);
    }
    else if (Pympf_Check(other)) {
        if (mpfr_sgn(Pympf_AS_MPF(other)) < 0) {
            VALUE_ERROR("sqrt() of negative number");
            return NULL;
        }
        mpfr_sqrt(result->f, Pympf_AS_MPF(other), options.rounding);
    }
    else {
        if (!(tempx = Pympf_From_Float(other, 0))) {
            TYPE_ERROR("sqrt() requires 'mpf' argument");
            return NULL;
        }
        else {
            if (mpfr_sgn(Pympf_AS_MPF(tempx)) < 0) {
                VALUE_ERROR("sqrt() of negative number");
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpfr_sqrt(result->f, tempx->f, options.rounding);
            Py_DECREF((PyObject*)tempx);
        }
    }
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
Pympf_sign(PyObject *self, PyObject *other)
{
    long sign;

    PympfObject *tempx;

    if (self && (Pympf_Check(self))) {
        sign = mpfr_sgn(Pympf_AS_MPF(self));
    }
    else if (Pympf_Check(other)) {
        sign = mpfr_sgn(Pympf_AS_MPF(other));
    }
    else {
        if (!(tempx = Pympf_From_Float(other, 0))) {
            TYPE_ERROR("sign() requires 'mpf' argument");
            return NULL;
        }
        else {
            sign = mpfr_sgn(tempx->f);
            Py_DECREF((PyObject*)tempx);
        }
    }
    return PyIntOrLong_FromLong(sign);
}

#define MPF_UNIOP_NOROUND(NAME) \
static PyObject * \
Pympf_##NAME(PyObject* self, PyObject *other) \
{ \
    PympfObject *result, *tempx; \
    if (!(result = Pympf_new(0))) return NULL; \
    if(self && Pympf_Check(self)) { \
        mpfr_##NAME(result->f, Pympf_AS_MPF(self)); \
    } \
    else if (Pympf_Check(other)) { \
        mpfr_##NAME(result->f, Pympf_AS_MPF(other)); \
    } \
    else { \
        if (!(tempx = Pympf_From_Float(other, 0))) { \
            TYPE_ERROR(#NAME "() requires 'mpf' argument"); \
            return NULL; \
        } \
        else { \
            mpfr_##NAME(result->f, tempx->f); \
            Py_DECREF((PyObject*)tempx); \
        } \
    } \
    return (PyObject*)result; \
}

static char doc_ceilm[]="\
x.ceil(): returns an mpf that is the smallest integer >= x\n\
";
static char doc_ceilg[]="\
ceil(x): returns an mpf that is the smallest integer >= x\n\
x must be an mpf, or else gets coerced to one.\n\
";

MPF_UNIOP_NOROUND(ceil)

static char doc_floorm[]="\
x.floor(): returns an mpf that is the smallest integer <= x\n\
";
static char doc_floorg[]="\
floor(x): returns an mpf that is the smallest integer <= x\n\
x must be an mpf, or else gets coerced to one.\n\
";

MPF_UNIOP_NOROUND(floor);

static char doc_truncm[]="\
x.trunc(): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
";
static char doc_truncg[]="\
trunc(x): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
x must be an mpf, or else gets coerced to one.\n\
";

MPF_UNIOP_NOROUND(trunc)

#define MPF_UNIOP(NAME) \
static PyObject * \
Pympf_##NAME(PyObject* self, PyObject *other) \
{ \
    PympfObject *result, *tempx; \
    if (!(result = Pympf_new(0))) return NULL; \
    if(self && Pympf_Check(self)) { \
        mpfr_##NAME(result->f, Pympf_AS_MPF(self), options.rounding); \
    } \
    else if (Pympf_Check(other)) { \
        mpfr_##NAME(result->f, Pympf_AS_MPF(other), options.rounding); \
    } \
    else { \
        if (!(tempx = Pympf_From_Float(other, 0))) { \
            TYPE_ERROR(#NAME "() requires 'mpf' argument"); \
            return NULL; \
        } \
        else { \
            mpfr_##NAME(result->f, tempx->f, options.rounding); \
            Py_DECREF((PyObject*)tempx); \
        } \
    } \
    return (PyObject*)result; \
}

PyDoc_STRVAR(doc_mpf_sqr,
"x.square() ->: mpf\n\n"
"Return x * x.");

MPF_UNIOP(sqr)

static char doc_flogm[]="\
x.log(): returns natural logarithm of x.\n\
";
static char doc_flogg[]="\
log(x): returns natural logarithm of x.\n\
";

MPF_UNIOP(log)

static char doc_flog2m[]="\
x.log2(): returns base-2 logarithm of x.\n\
";
static char doc_flog2g[]="\
log2(x): returns base-2 logarithm of x.\n\
";

MPF_UNIOP(log2)

static char doc_flog10m[]="\
x.log10(): returns base-10 logarithm of x.\n\
";
static char doc_flog10g[]="\
log10(x): returns base-10 logarithm of x.\n\
";

MPF_UNIOP(log10)

static char doc_fexpm[]="\
x.exp(): returns exponential x.\n\
";
static char doc_fexpg[]="\
exp(x): returns exponential of x.\n\
";

MPF_UNIOP(exp)

static char doc_fexp2m[]="\
x.exp2(): returns 2**x.\n\
";
static char doc_fexp2g[]="\
exp2(x): returns 2**x.\n\
";

MPF_UNIOP(exp2)

static char doc_fexp10m[]="\
x.exp10(): returns 10**x.\n\
";
static char doc_fexp10g[]="\
exp10(x): returns 10**x.\n\
";

MPF_UNIOP(exp10)

static char doc_fsinm[]="\
x.sin(): returns sine of x; x in radians.\n\
";
static char doc_fsing[]="\
sin(x): returns sine of x; x in radians.\n\
";

MPF_UNIOP(sin)

static char doc_fcosm[]="\
x.cos(): returns cosine of x; x in radians.\n\
";
static char doc_fcosg[]="\
cos(x): returns cosine of x; x in radians.\n\
";

MPF_UNIOP(cos)

static char doc_ftanm[]="\
x.tan(): returns tangent of x; x in radians.\n\
";
static char doc_ftang[]="\
tan(x): returns tangent of x; x in radians.\n\
";

MPF_UNIOP(tan)

static char doc_fsecm[]="\
x.sec(): returns secant of x; x in radians.\n\
";
static char doc_fsecg[]="\
sec(x): returns secant of x; x in radians.\n\
";

MPF_UNIOP(sec)

static char doc_fcscm[]="\
x.csc(): returns cosecant of x; x in radians.\n\
";
static char doc_fcscg[]="\
csc(x): returns cosecant of x; x in radians.\n\
";

MPF_UNIOP(csc)

static char doc_fcotm[]="\
x.cot(): returns cotangent of x; x in radians.\n\
";
static char doc_fcotg[]="\
cot(x): returns cotangent of x; x in radians.\n\
";

MPF_UNIOP(cot)

static char doc_facosm[]="\
x.acos(): returns arc-cosine of x; x in radians.\n\
";
static char doc_facosg[]="\
acos(x): returns arc-cosine of x; x in radians.\n\
";

MPF_UNIOP(acos)

static char doc_fasinm[]="\
x.asin(): returns arc-sine of x; x in radians.\n\
";
static char doc_fasing[]="\
asin(x): returns arc-sine of x; x in radians.\n\
";

MPF_UNIOP(asin)

static char doc_fatanm[]="\
x.atan(): returns arc-tangent of x; x in radians.\n\
";
static char doc_fatang[]="\
atan(x): returns arc-tangent of x; x in radians.\n\
";

MPF_UNIOP(atan)

static char doc_fcoshm[]="\
x.cosh(): returns hyperbolic cosine of x.\n\
";
static char doc_fcoshg[]="\
cosh(x): returns hyperbolic cosine of x.\n\
";

MPF_UNIOP(cosh)

static char doc_fsinhm[]="\
x.sinh(): returns hyperbolic sine of x.\n\
";
static char doc_fsinhg[]="\
sinh(x): returns hyperbolic sine of x.\n\
";

MPF_UNIOP(sinh)

static char doc_ftanhm[]="\
x.tanh(): returns hyperbolic tangent of x.\n\
";
static char doc_ftanhg[]="\
tanh(x): returns hyperbolic tangent of x.\n\
";

MPF_UNIOP(tanh)

static char doc_fsechm[]="\
x.sech(): returns hyperbolic secant of x.\n\
";
static char doc_fsechg[]="\
sech(x): returns hyperbolic secant of x.\n\
";

MPF_UNIOP(sech)

static char doc_fcschm[]="\
x.csch(): returns hyperbolic cosecant of x.\n\
";
static char doc_fcschg[]="\
csch(x): returns hyperbolic cosecant of x.\n\
";

MPF_UNIOP(csch)

static char doc_fcothm[]="\
x.coth(): returns hyperbolic cotangent of x.\n\
";
static char doc_fcothg[]="\
coth(x): returns hyperbolic cotangent of x.\n\
";

MPF_UNIOP(coth)

static char doc_facoshm[]="\
x.acosh(): returns inverse hyperbolic cosine of x.\n\
";
static char doc_facoshg[]="\
acosh(x): returns inverse hyperbolic cosine of x.\n\
";

MPF_UNIOP(acosh)

static char doc_fasinhm[]="\
x.asinh(): returns inverse hyperbolic sine of x.\n\
";
static char doc_fasinhg[]="\
asinh(x): returns inverse hyperbolic sine of x.\n\
";

MPF_UNIOP(asinh)

static char doc_fatanhm[]="\
x.atanh(): returns inverse hyperbolic tangent of x.\n\
";
static char doc_fatanhg[]="\
atanh(x): returns inverse hyperbolic tangent of x.\n\
";

MPF_UNIOP(atanh)

static char doc_flog1pm[]="\
x.log1p(): returns logarithm of (1+x).\n\
";
static char doc_flog1pg[]="\
log1p(x): returns logarithm of (1+x).\n\
";

MPF_UNIOP(log1p)

static char doc_fexpm1m[]="\
x.expm1(): returns exponential(x) - 1.\n\
";
static char doc_fexpm1g[]="\
expm1(x): returns exponential(x) - 1.\n\
";

MPF_UNIOP(expm1)

static char doc_feintm[]="\
x.eint(): returns exponential integral of x.\n\
";
static char doc_feintg[]="\
eint(x): returns exponential integral of x.\n\
";

MPF_UNIOP(eint)

static char doc_fli2m[]="\
x.li2(): returns real part of dilogarithm of x.\n\
";
static char doc_fli2g[]="\
li2(x): returns real part of dilogarithm of x.\n\
";

MPF_UNIOP(li2)

static char doc_fgammam[]="\
x.gamma(): returns gamma of x.\n\
";
static char doc_fgammag[]="\
gamma(x): returns gamma of x.\n\
";

MPF_UNIOP(gamma)

static char doc_flngammam[]="\
x.lngamma(): returns logarithm of gamma(x).\n\
";
static char doc_flngammag[]="\
lngamma(x): returns logarithm of gamma(x).\n\
";

MPF_UNIOP(lngamma)

static char doc_fdigammam[]="\
x.digamma(): returns digamma of x.\n\
";
static char doc_fdigammag[]="\
digamma(x): returns digamma of x.\n\
";

MPF_UNIOP(digamma)

static char doc_fzetam[]="\
x.zeta(): returns Riemann zeta of x.\n\
";
static char doc_fzetag[]="\
zeta(x): returns Riemann zeta of x.\n\
";

MPF_UNIOP(zeta)

static char doc_ferfm[]="\
x.erf(): returns error function of x.\n\
";
static char doc_ferfg[]="\
erf(x): returns error function of x.\n\
";

MPF_UNIOP(erf)

static char doc_ferfcm[]="\
x.erfc(): returns complementary error function of x.\n\
";
static char doc_ferfcg[]="\
erfc(x): returns complementary error function of x.\n\
";

MPF_UNIOP(erfc)

static char doc_fj0m[]="\
x.j0(): returns first kind Bessel function of order 0 of x.\n\
";
static char doc_fj0g[]="\
j0(x): returns first kind Bessel function of order 0 of x.\n\
";

MPF_UNIOP(j0)

static char doc_fj1m[]="\
x.j1(): returns first kind Bessel function of order 1 of x.\n\
";
static char doc_fj1g[]="\
j1(x): returns first kind Bessel function of order 1 of x.\n\
";

MPF_UNIOP(j1)

static char doc_fy0m[]="\
x.y0(): returns second kind Bessel function of order 0 of x.\n\
";
static char doc_fy0g[]="\
y0(x): returns second kind Bessel function of order 0 of x.\n\
";

MPF_UNIOP(y0)

static char doc_fy1m[]="\
x.y1(): returns second kind Bessel function of order 1 of x.\n\
";
static char doc_fy1g[]="\
y1(x): returns second kind Bessel function of order 1 of x.\n\
";

MPF_UNIOP(y1)

static char doc_faim[]="\
x.ai(): returns Airy function of x.\n\
";
static char doc_faig[]="\
ai(x): returns Airy function of x.\n\
";

MPF_UNIOP(ai)

