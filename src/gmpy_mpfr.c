/* gmpy_mpfr.c
 *
 * Functions that operate strictly on mpfr.
 *
 * This file should be considered part of gmpy2.c
 */

/* Implement the .precision attribute of an mpfr. */

static PyObject *
Pympfr_getprec_attrib(PympfrObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_prec(self->f));
}

/* Implement the .rc attribute of an mpfr. */

static PyObject *
Pympfr_getrc_attrib(PympfrObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)self->rc);
}

/* Implement the nb_bool slot. */

static int
Pympfr_nonzero(PympfrObject *self)
{
    return !mpfr_zero_p(self->f);
}

/* Implement the nb_positive slot. */

static PyObject *
Pympfr_pos(PympfrObject *self)
{
    PympfrObject *result;

    if (!(result = Pympfr_new(mpfr_get_prec(self->f))))
        return NULL;

    mpfr_clear_flags();

    /* Since result has the same precision as self, no rounding occurs. */
    mpfr_set(result->f, self->f, context->now.mpfr_round);
    result->round_mode = self->round_mode;
    result->rc = self->rc;
    /* Force the exponents to be valid. */
    result->rc = mpfr_check_range(result->f, result->rc, result->round_mode);
    /* Now round result to the current precision. */
    result->rc = mpfr_prec_round(result->f, context->now.mpfr_prec,
                                 context->now.mpfr_round);

    SUBNORMALIZE(result);
    MERGE_FLAGS;
    CHECK_FLAGS("pos");
  done:
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_get_emin_min,
"get_emin_min() -> integer\n\n"
"Return the minimum possible exponent that can be set for 'mpfr'.");

static PyObject *
Pympfr_get_emin_min(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_emin_min());
}

PyDoc_STRVAR(doc_g_mpfr_get_emax_max,
"get_emax_max() -> integer\n\n"
"Return the maximum possible exponent that can be set for 'mpfr'.");

static PyObject *
Pympfr_get_emax_max(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)mpfr_get_emax_max());
}

PyDoc_STRVAR(doc_g_mpfr_get_max_precision,
"get_max_precision() -> integer\n\n"
"Return the maximum bits of precision that can be used for calculations.\n"
"Note: to allow extra precision for intermediate calculations, avoid\n"
"setting precision close the maximum precisicon.");

static PyObject *
Pympfr_get_max_precision(PyObject *self, PyObject *args)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)MPFR_PREC_MAX);
}

PyDoc_STRVAR(doc_g_mpfr_set_nan,
"nan() -> mpfr\n\n"
"Return an 'mpfr' inialized to 'nan'.");

static PyObject *
Pympfr_set_nan(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    if ((result = Pympfr_new(0)))
        mpfr_set_nan(result->f);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_inf,
"inf(n) -> mpfr\n\n"
"Return an 'mpfr' inialized to 'inf' with the same sign as n.");

static PyObject *
Pympfr_set_inf(PyObject *self, PyObject *other)
{
    PympfrObject *result;
    long s = clong_From_Integer(other);

    if (s == -1 && PyErr_Occurred()) {
        TYPE_ERROR("inf() requires 'int' argument");
        return NULL;
    }

    if ((result = Pympfr_new(0)))
        mpfr_set_inf(result->f, s<0?-1:1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_zero,
"zero() -> mpfr\n\n"
"Return an 'mpfr' inialized to 0.0 with the same sign as n.");

static PyObject *
Pympfr_set_zero(PyObject *self, PyObject *other)
{
    PympfrObject *result;
    long s = clong_From_Integer(other);

    if (s == -1 && PyErr_Occurred()) {
        TYPE_ERROR("zero() requires 'int' argument");
        return NULL;
    }

    if ((result = Pympfr_new(0)))
        mpfr_set_zero(result->f, s<0?-1:1);
    return (PyObject*)result;
}

#define MPFR_TEST_OTHER(NAME, msg) \
static PyObject * \
Pympfr_is_##NAME(PyObject *self, PyObject *other)\
{\
    int res;\
    if(self && Pympfr_Check(self)) {\
        Py_INCREF(self);\
    }\
    else if(Pympfr_Check(other)) {\
        self = other;\
        Py_INCREF((PyObject*)self);\
    }\
    else if (!(self = (PyObject*)Pympfr_From_Real(other, 0))) {\
        PyErr_SetString(PyExc_TypeError, msg);\
        return NULL;\
    }\
    res = mpfr_##NAME##_p(Pympfr_AS_MPFR(self));\
    Py_DECREF(self);\
    if (res)\
        Py_RETURN_TRUE;\
    else\
        Py_RETURN_FALSE;\
}

PyDoc_STRVAR(doc_mpfr_is_nan,
"x.is_nan() -> boolean\n\n"
"Return True if x is 'nan' (Not-A-Number), False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_nan,
"is_nan(x) -> boolean\n\n"
"Return True if x is 'nan' (Not-A-Number), False otherwise.");

MPFR_TEST_OTHER(nan, "is_nan() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_inf,
"x.is_inf() -> boolean\n\n"
"Return True if x is +Infinity or -Infinity, False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_inf,
"is_inf(x) -> boolean\n\n"
"Return True if x is +Infinity or -Infinity, False otherwise.");

MPFR_TEST_OTHER(inf, "is_inf() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_number,
"x.is_number() -> boolean\n\n"
"Return True if x is an actual number (i.e. not NaN or Infinity),\n"
"False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_number,
"is_number(x) -> boolean\n\n"
"Return True if x is an actual number (i.e. not NaN or Infinity),\n"
"False otherwise.");

MPFR_TEST_OTHER(number, "is_number() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_zero,
"x.is_zero() -> boolean\n\n"
"Return True if x is zero, False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_zero,
"is_zero(x) -> boolean\n\n"
"Return True if x is zero, False otherwise.");

MPFR_TEST_OTHER(zero, "is_zero() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_regular,
"x.is_regular() -> boolean\n\n"
"Return True if x is not zero, NaN, or Infinity, False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_regular,
"is_regular(x) -> boolean\n\n"
"Return True if x is not zero, NaN, or Infinity, False otherwise.");

MPFR_TEST_OTHER(regular, "is_regular() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_integer,
"x.is_integer() -> boolean\n\n"
"Return True if x is an integer, False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_integer,
"is_integer(x) -> boolean\n\n"
"Return True if x is an integer, False otherwise.");

MPFR_TEST_OTHER(integer, "is_integer() requires 'mpfr' argument");

/* produce string for an mpfr with requested/defaulted parameters */

PyDoc_STRVAR(doc_mpfr_digits,
"x.digits(base=10, prec=0) -> (mantissa, exponent, bits)\n\n"
"Returns up to 'prec' digits in the given base. If 'prec' is 0, as many\n"
"digits that are available are returned. No more digits than available\n"
"given x's precision are returned. 'base' must be between 2 and 62,\n"
"inclusive. The result is a three element tuple containig the mantissa,\n"
"the exponent, and the number of bits of precision.");

/* TODO: support keyword arguments. */

static PyObject *
Pympfr_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    int prec = 0;
    PyObject *result;

    if (self && Pympfr_Check(self)) {
        if (!PyArg_ParseTuple(args, "|ii", &base, &prec))
            return NULL;
        Py_INCREF(self);
    }
    else {
        if(!PyArg_ParseTuple(args, "O&|ii", Pympfr_convert_arg, &self,
                            &base, &prec))
        return NULL;
    }
    result = Pympfr_ascii((PympfrObject*)self, base, prec, 0, 0, 2);
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_mpfr_f2q,
"x.f2q([err]) -> mpq\n\n"
"Return the 'best' mpq approximating x to within relative error 'err'.\n"
"Default is the precision of x. Uses Stern-Brocot tree to find the\n"
"'best' approximation. An 'mpz' is returned if the the denominator\n"
"is 1. If 'err'<0, error sought is 2.0 ** err.");

PyDoc_STRVAR(doc_g_mpfr_f2q,
"f2q(x,[err]) -> mpq\n\n"
"Return the 'best' mpq approximating x to within relative error 'err'.\n"
"Default is the precision of x. Uses Stern-Brocot tree to find the\n"
"'best' approximation. An 'mpz' is returned if the the denominator\n"
"is 1. If 'err'<0, error sought is 2.0 ** err.");

/* TODO: Redo f2q. Ref-counting looks strange. */

static PyObject *
Pympfr_f2q(PyObject *self, PyObject *args)
{
    PympfrObject *err = 0;
    PympfrObject *fself;
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "Pympfr_f2q: %p, %p\n", self, args);
#endif
    SELF_MPFR_ONE_ARG_CONVERTED_OPT(&err);

    fself = (PympfrObject*)self;
    return f2q_internal(fself, err, mpfr_get_prec(fself->f), args!=0);
}

static PyObject *
f2q_internal(PympfrObject* self, PympfrObject* err, unsigned int bits, int mayz)
{
    PympqObject *res = 0;
    int i, negative, errsign;
    mpfr_t f, al, a, r1[3], r2[3], minerr, curerr, newerr, temp;

    assert(!err || Pympfr_Check(err));
    errsign = err ? mpfr_sgn(err->f) : 0;
    if (errsign == 0) {
        if (err) {
            Py_DECREF((PyObject*)err);
        }
        if (!(err = Pympfr_new(20))) {
            Py_DECREF((PyObject*)self);
            return NULL;
        }
        mpfr_set_si(err->f, 1, context->now.mpfr_round);
        mpfr_div_2ui(err->f, err->f, bits, context->now.mpfr_round);
    }
    else if (errsign < 0) {
        int ubits;
        mpfr_floor(err->f, err->f);
        ubits = (int)mpfr_get_d(err->f, context->now.mpfr_round);
        mpfr_set_si(err->f, 1, context->now.mpfr_round);
        mpfr_div_2si(err->f, err->f, -ubits, context->now.mpfr_round);
    }
    if (!(res = Pympq_new()))
        return NULL;
    mpfr_init2(minerr, 20);
    mpfr_set(minerr, err->f, context->now.mpfr_round);
    Py_DECREF((PyObject*)err);

    mpfr_init2(f, bits);
    if (mpfr_sgn(self->f) < 0) {
        negative = 1;
        mpfr_abs(f, self->f, context->now.mpfr_round);
    }
    else {
        negative = 0;
        mpfr_set(f, self->f, context->now.mpfr_round);
    }
    Py_DECREF((PyObject*)self);
    mpfr_init2(al, bits);
    mpfr_set(al, f, context->now.mpfr_round);
    mpfr_init2(a, bits);
    mpfr_floor(a, al);
    mpfr_init2(temp, bits);
    for (i=0; i<3; ++i) {
        mpfr_init2(r1[i], bits);
        mpfr_init2(r2[i], bits);
    }
    mpfr_set_si(r1[0], 0, context->now.mpfr_round);
    mpfr_set_si(r1[1], 0, context->now.mpfr_round);
    mpfr_set_si(r1[2], 1, context->now.mpfr_round);
    mpfr_set_si(r2[0], 0, context->now.mpfr_round);
    mpfr_set_si(r2[1], 1, context->now.mpfr_round);
    mpfr_set(r2[2], a, context->now.mpfr_round);
    mpfr_init2(curerr, 20);
    mpfr_init2(newerr, 20);
    mpfr_reldiff(curerr, f, a, context->now.mpfr_round);
    while (mpfr_cmp(curerr, minerr) > 0) {
        mpfr_sub(temp, al, a, context->now.mpfr_round);
        mpfr_ui_div(al, 1, temp, context->now.mpfr_round);
        mpfr_floor(a, al);
        mpfr_swap(r1[0], r1[1]);
        mpfr_swap(r1[1], r1[2]);
        mpfr_mul(r1[2], r1[1], a, context->now.mpfr_round);
        mpfr_add(r1[2], r1[2], r1[0], context->now.mpfr_round);
        mpfr_swap(r2[0], r2[1]);
        mpfr_swap(r2[1], r2[2]);
        mpfr_mul(r2[2], r2[1], a, context->now.mpfr_round);
        mpfr_add(r2[2], r2[2], r2[0], context->now.mpfr_round);
        mpfr_div(temp, r2[2], r1[2], context->now.mpfr_round);
        mpfr_reldiff(newerr, f, temp, context->now.mpfr_round);
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
        mpfr_get_z(Pympz_AS_MPZ(res), r2[2], context->now.mpfr_round);
        if (negative)
            mpz_neg(Pympz_AS_MPZ(res),Pympz_AS_MPZ(res));
    }
    else {
        mpfr_get_z(mpq_numref(res->q), r2[2], context->now.mpfr_round);
        mpfr_get_z(mpq_denref(res->q), r1[2], context->now.mpfr_round);
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

static Py_hash_t
Pympfr_hash(PympfrObject *self)
{
#ifdef _PyHASH_MODULUS
    Py_uhash_t hash = 0;
    Py_ssize_t exp;
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
    if (hash == (Py_uhash_t)-1)
        hash = (Py_uhash_t)-2;
    return (self->hash_cache = (Py_hash_t)hash);
#else
    double temp;
    if (self->hash_cache != -1)
        return self->hash_cache;
    temp = mpfr_get_d(self->f, context->now.mpfr_round);
    return (self->hash_cache = _Py_HashDouble(temp));
#endif
}

/* This function is used in gmpy_basic. */

static PyObject *
Pympfr2_pow(PyObject *base, PyObject *exp, PyObject *m)
{
    PympfrObject *tempb, *tempe, *result;

    if ((PyObject*)m != Py_None) {
        TYPE_ERROR("mpfr.pow() no modulo allowed");
        return NULL;
    }

    tempb = Pympfr_From_Real(base, 0);
    tempe = Pympfr_From_Real(exp, 0);
    result = Pympfr_new(0);

    if (!tempe || !tempb || !result) {
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempb);
        Py_XDECREF((PyObject*)result);
        Py_RETURN_NOTIMPLEMENTED;
    }

    if (mpfr_zero_p(tempb->f) && (mpfr_sgn(tempe->f) < 0)) {
        context->now.divzero = 1;
        if (context->now.trap_divzero) {
            GMPY_DIVZERO("zero cannot be raised to a negative power");
            goto done;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_pow(result->f, tempb->f, tempe->f,
                          context->now.mpfr_round);
    SUBNORMALIZE(result)
    MERGE_FLAGS
    CHECK_FLAGS("pow()")
  done:
    Py_DECREF((PyObject*)tempe);
    Py_DECREF((PyObject*)tempb);
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

#define MPFR_CONST(NAME) \
static PyObject * \
Pympfr_##NAME(PyObject *self, PyObject *args) \
{ \
    PympfrObject *result; \
    if ((result = Pympfr_new(0))) { \
        mpfr_clear_flags(); \
        result->rc = mpfr_##NAME(result->f, context->now.mpfr_round); \
        MERGE_FLAGS \
        CHECK_FLAGS(#NAME "()") \
    } \
  done: \
    return (PyObject*)result; \
}

PyDoc_STRVAR(doc_mpfr_const_pi,
"const_pi() -> mpfr\n\n"
"Return the constant pi using the default precision.");

MPFR_CONST(const_pi)

PyDoc_STRVAR(doc_mpfr_const_euler,
"const_euler() -> mpfr\n\n"
"Return the euler constant using the default precision.");

MPFR_CONST(const_euler)

PyDoc_STRVAR(doc_mpfr_const_log2,
"const_log2() -> mpfr\n\n"
"Return the log2 constant using the default precision.");

MPFR_CONST(const_log2)

PyDoc_STRVAR(doc_mpfr_const_catalan,
"const_catalan() -> mpfr\n\n"
"Return the catalan constant using the default precision.");

MPFR_CONST(const_catalan)

PyDoc_STRVAR(doc_mpfr_sqrt,
"x.sqrt() -> mpfr\n\n"
"Return the square root of x.  x must be >= 0.");

static PyObject *
Pympfr_sqrt(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("sqrt() requires 'mpfr' argument");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_sqrt(result->f, Pympfr_AS_MPFR(self),
                           context->now.mpfr_round);

    MPFR_CLEANUP_SELF("sqrt()");
}

PyDoc_STRVAR(doc_mpfr_rec_sqrt,
"x.rec_sqrt() -> mpfr\n\n"
"Return the reciprocal of the square root of x.");

PyDoc_STRVAR(doc_g_mpfr_rec_sqrt,
"rec_sqrt(x) -> mpfr\n\n"
"Return the reciprocal of the square root of x.");

static PyObject *
Pympfr_rec_sqrt(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("rec_sqrt() requires 'mpfr' argument");

    if (!(result = Pympfr_new(0)))
        goto done;

    if (mpfr_zero_p(Pympfr_AS_MPFR(self))) {
        context->now.divzero = 1;
        if (context->now.trap_divzero) {
            GMPY_DIVZERO("division by zero in 'mpfr' reciprocal square root");
            goto done;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_rec_sqrt(result->f, Pympfr_AS_MPFR(self),
                               context->now.mpfr_round);

    MPFR_CLEANUP_SELF("rec_sqrt()");
}

PyDoc_STRVAR(doc_mpfr_root,
"x.root(n) -> mpfr\n\n"
"Return the n-th root of x.");
static PyObject *
Pympfr_root(PyObject *self, PyObject *args)
{
    long n;
    PympfrObject *result;

    PARSE_ONE_MPFR_REQ_CLONG(&n, "root() requires 'mpfr','int' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    if (n <= 0) {
        VALUE_ERROR("n must be > 0");
        goto done;
    }

    mpfr_clear_flags();
    result->rc = mpfr_root(result->f, Pympfr_AS_MPFR(self), n,
                           context->now.mpfr_round);

    MPFR_CLEANUP_SELF("root()");
}


static char doc_mpfr_round[] = "\
x.round(n): returns x rounded to n bits. Uses default precision if\n\
n is not specified. See round2() to access the mpfr_round() function.";
static char doc_g_mpfr_round[] = "\
round(x, n): returns x rounded to n bits. Uses default precision\n\
if n is not specified. See round2() to access the mpfr_round()\n\
function.";
static PyObject *
Pympfr_round(PyObject *self, PyObject *args)
{
    mpfr_prec_t prec = context->now.mpfr_prec;
    PympfrObject *result;

    PARSE_ONE_MPFR_OPT_CLONG(&prec,
            "round() requires 'mpfr',['int'] arguments");

    if (!(result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(self))))) {
        goto done;
    }

    mpfr_clear_flags();
    /* Duplicate the code from Pympfr_pos. */
    mpfr_set(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    result->round_mode = ((PympfrObject*)self)->round_mode;
    result->rc = ((PympfrObject*)self)->rc;
    result->rc = mpfr_check_range(result->f, result->rc, result->round_mode);
    result->rc = mpfr_prec_round(result->f, prec, context->now.mpfr_round);

    MPFR_CLEANUP_SELF("round()");
}

static char doc_mpfr_reldiff[] = "\
x.reldiff(y): returns the relative difference between x and y,\n\
where y can be any number and gets coerced to an mpfr; result is\n\
an mpfr roughly equal to abs(x-y)/x.\n\
";
static char doc_g_mpfr_reldiff[] = "\
reldiff(x,y): returns the relative difference between x and y,\n\
where x and y can be any numbers and get coerced to mpfr; result is\n\
an mpfr roughly equal to abs(x-y)/x.\n\
";
static PyObject *
Pympfr_reldiff(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "reldiff() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    /* mpfr_reldiff doesn't guarantee correct rounding and doesn't appear
     * to set any exceptions.
     */
    mpfr_reldiff(result->f, Pympfr_AS_MPFR(self), Pympfr_AS_MPFR(other),
                 context->now.mpfr_round);
    result->rc = 0;
    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static char doc_mpfr_sign[]="\
x.sign(): returns -1, 0, or +1, if x is negative, 0, or positive.\n\
";
static PyObject *
Pympfr_sign(PyObject *self, PyObject *other)
{
    long sign;

    PARSE_ONE_MPFR_OTHER("sign() requires 'mpfr' argument");

    mpfr_clear_flags();
    sign = mpfr_sgn(Pympfr_AS_MPFR(self));

    MERGE_FLAGS;
    CHECK_ERANGE("range error in 'mpfr' sign(), NaN argument");

  done:
    Py_DECREF((PyObject*)self);
    if (PyErr_Occurred())
        return NULL;
    else
        return PyIntOrLong_FromLong(sign);
}

#define MPFR_MONOP(NAME) \
static PyObject * \
Py##NAME(PympfrObject *x) \
{ \
    PympfrObject *r; \
    if (!(r = Pympfr_new(mpfr_get_prec(x->f)))) \
        return NULL; \
    if (Pympfr_CheckAndExp(x)) { \
        r->rc = NAME(r->f, x->f, context->now.mpfr_round); \
    } \
    else { \
        mpfr_set(r->f, x->f, context->now.mpfr_round); \
        r->round_mode = x->round_mode; \
        r->rc = x->rc; \
        mpfr_clear_flags(); \
        mpfr_check_range(r->f, r->rc, r->round_mode); \
        r->rc = NAME(r->f, r->f, context->now.mpfr_round); \
        MERGE_FLAGS; \
        CHECK_FLAGS(#NAME "()"); \
    } \
  done: \
    return (PyObject *) r; \
}

MPFR_MONOP(mpfr_abs)
MPFR_MONOP(mpfr_neg)

#define MPFR_UNIOP_NOROUND(NAME) \
static PyObject * \
Pympfr_##NAME(PyObject* self, PyObject *other) \
{ \
    PympfrObject *result; \
    PARSE_ONE_MPFR_OTHER(#NAME "() requires 'mpfr' argument"); \
    if (!(result = Pympfr_new(0))) goto done; \
    mpfr_clear_flags(); \
    result->rc = mpfr_##NAME(result->f, Pympfr_AS_MPFR(self)); \
    MPFR_CLEANUP_SELF(#NAME "()"); \
}

static char doc_mpfr_ceil[]="\
x.ceil(): returns an mpfr that is the smallest integer >= x\n\
";
static char doc_g_mpfr_ceil[]="\
ceil(x): returns an mpfr that is the smallest integer >= x\n\
x must be an mpfr, or else gets coerced to one.\n\
";

MPFR_UNIOP_NOROUND(ceil)

static char doc_mpfr_floor[]="\
x.floor(): returns an mpfr that is the smallest integer <= x\n\
";
static char doc_g_mpfr_floor[]="\
floor(x): returns an mpfr that is the smallest integer <= x\n\
x must be an mpfr, or else gets coerced to one.\n\
";

MPFR_UNIOP_NOROUND(floor);

static char doc_mpfr_trunc[]="\
x.trunc(): returns an mpfr that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
";
static char doc_g_mpfr_trunc[]="\
trunc(x): returns an mpfr that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
x must be an mpfr, or else gets coerced to one.\n\
";

MPFR_UNIOP_NOROUND(trunc)

static char doc_mpfr_round2[]="\
x.round2(): returns an mpfr that is x rounded to the nearest\n\
integer, with ties rounded away from 0.";
static char doc_g_mpfr_round2[]="\
round2(x): returns an mpfr that is x rounded to the nearest integer,\n\
with ties rounded away from 0.";

static PyObject *
Pympfr_round2(PyObject* self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("round2() requires 'mpfr' argument");
    if (!(result = Pympfr_new(0)))
        goto done;
    mpfr_clear_flags();
    result->rc = mpfr_round(result->f, Pympfr_AS_MPFR(self));
    MPFR_CLEANUP_SELF("round2()");
}

#define MPFR_UNIOP(NAME) \
static PyObject * \
Pympfr_##NAME(PyObject* self, PyObject *other) \
{ \
    PympfrObject *result; \
    PARSE_ONE_MPFR_OTHER(#NAME "() requires 'mpfr' argument"); \
    if (!(result = Pympfr_new(0))) goto done; \
    mpfr_clear_flags(); \
    result->rc = mpfr_##NAME(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round); \
    MPFR_CLEANUP_SELF(#NAME "()"); \
}

PyDoc_STRVAR(doc_mpfr_rint,
"x.rint() -> mpfr\n\n"
"Return x rounded to the nearest integer using the current rounding\n"
"mode.");
PyDoc_STRVAR(doc_g_mpfr_rint,
"rint(x) -> mpfr\n\n"
"Return x rounded to the nearest integer using the current rounding\n"
"mode.");

MPFR_UNIOP(rint)

PyDoc_STRVAR(doc_mpfr_rint_ceil,
"x.rint_ceil() -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the next\n"
"higher or equal integer and then, if needed, using the current\n"
"rounding mode.");
PyDoc_STRVAR(doc_g_mpfr_rint_ceil,
"rint_ceil(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the next\n"
"higher or equal integer and then, if needed, using the current\n"
"rounding mode.");

MPFR_UNIOP(rint_ceil)

PyDoc_STRVAR(doc_mpfr_rint_floor,
"x.rint_floor() -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the next\n"
"lower or equal integer and then, if needed, using the current\n"
"rounding mode.");
PyDoc_STRVAR(doc_g_mpfr_rint_floor,
"rint_floor(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the next\n"
"lower or equal integer and then, if needed, using the current\n"
"rounding mode.");

MPFR_UNIOP(rint_floor)

PyDoc_STRVAR(doc_mpfr_rint_round,
"x.rint_round() -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the current rounding mode.");
PyDoc_STRVAR(doc_g_mpfr_rint_round,
"rint_round(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the current rounding mode.");

MPFR_UNIOP(rint_round)

PyDoc_STRVAR(doc_mpfr_rint_trunc,
"x.rint_trunc() -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the current rounding mode.");
PyDoc_STRVAR(doc_g_mpfr_rint_trunc,
"rint_trunc(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the current rounding mode.");

MPFR_UNIOP(rint_trunc)

PyDoc_STRVAR(doc_mpfr_frac,
"x.frac() -> mpfr\n\n"
"Return fractional part of x.");
PyDoc_STRVAR(doc_g_mpfr_frac,
"frac(x) -> mpfr\n\n"
"Return fractional part of x.");

MPFR_UNIOP(frac)

PyDoc_STRVAR(doc_mpfr_modf,
"x.modf() -> (mpfr, mpfr)\n\n"
"Return a tuple containing the integer and fractional portions\n"
"of x.");

PyDoc_STRVAR(doc_g_mpfr_modf,
"modf(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the integer and fractional portions\n"
"of x.");

static PyObject *
Pympfr_modf(PyObject *self, PyObject *other)
{
    PympfrObject *s, *c;
    PyObject *result;
    int code;

    PARSE_ONE_MPFR_OTHER("modf() requires 'mpfr' argument");

    s = Pympfr_new(0);
    c = Pympfr_new(0);
    result = PyTuple_New(2);
    if (!s || !c || !result)
        goto done;

    mpfr_clear_flags();
    code = mpfr_modf(s->f, c->f, Pympfr_AS_MPFR(self),
                     context->now.mpfr_round);
    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;
    SUBNORMALIZE(s);
    SUBNORMALIZE(c);
    MERGE_FLAGS;
    CHECK_FLAGS("modf()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)s);
        PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    }
    return result;
}

PyDoc_STRVAR(doc_mpfr_sqr,
"x.square() ->: mpfr\n\n"
"Return x * x.");

MPFR_UNIOP(sqr)

PyDoc_STRVAR(doc_mpfr_cbrt,
"x.cbrt() -> mpfr\n\n"
"Return the cube root of x.");
PyDoc_STRVAR(doc_g_mpfr_cbrt,
"cbrt(x) -> mpfr\n\n"
"Return the cube root of x.");

MPFR_UNIOP(cbrt)

static char doc_mpfr_log[]="\
x.log(): returns natural logarithm of x.\n\
";
static char doc_g_mpfr_log[]="\
log(x): returns natural logarithm of x.\n\
";

MPFR_UNIOP(log)

static char doc_mpfr_log2[]="\
x.log2(): returns base-2 logarithm of x.\n\
";
static char doc_g_mpfr_log2[]="\
log2(x): returns base-2 logarithm of x.\n\
";

MPFR_UNIOP(log2)

static char doc_mpfr_log10[]="\
x.log10(): returns base-10 logarithm of x.\n\
";
static char doc_g_mpfr_log10[]="\
log10(x): returns base-10 logarithm of x.\n\
";

MPFR_UNIOP(log10)

static char doc_mpfr_exp[]="\
x.exp(): returns exponential of x.\n\
";
static char doc_g_mpfr_exp[]="\
exp(x): returns exponential of x.\n\
";

MPFR_UNIOP(exp)

static char doc_mpfr_exp2[]="\
x.exp2(): returns 2**x.\n\
";
static char doc_g_mpfr_exp2[]="\
exp2(x): returns 2**x.\n\
";

MPFR_UNIOP(exp2)

static char doc_mpfr_exp10[]="\
x.exp10(): returns 10**x.\n\
";
static char doc_g_mpfr_exp10[]="\
exp10(x): returns 10**x.\n\
";

MPFR_UNIOP(exp10)

static char doc_mpfr_sin[]="\
x.sin(): returns sine of x; x in radians.\n\
";
static char doc_g_mpfr_sin[]="\
sin(x): returns sine of x; x in radians.\n\
";

MPFR_UNIOP(sin)

static char doc_mpfr_cos[]="\
x.cos(): returns cosine of x; x in radians.\n\
";
static char doc_g_mpfr_cos[]="\
cos(x): returns cosine of x; x in radians.\n\
";

MPFR_UNIOP(cos)

static char doc_mpfr_tan[]="\
x.tan(): returns tangent of x; x in radians.\n\
";
static char doc_g_mpfr_tan[]="\
tan(x): returns tangent of x; x in radians.\n\
";

MPFR_UNIOP(tan)

static char doc_mpfr_sec[]="\
x.sec(): returns secant of x; x in radians.\n\
";
static char doc_g_mpfr_sec[]="\
sec(x): returns secant of x; x in radians.\n\
";

MPFR_UNIOP(sec)

static char doc_mpfr_csc[]="\
x.csc(): returns cosecant of x; x in radians.\n\
";
static char doc_g_mpfr_csc[]="\
csc(x): returns cosecant of x; x in radians.\n\
";

MPFR_UNIOP(csc)

static char doc_mpfr_cot[]="\
x.cot(): returns cotangent of x; x in radians.\n\
";
static char doc_g_mpfr_cot[]="\
cot(x): returns cotangent of x; x in radians.\n\
";

MPFR_UNIOP(cot)

static char doc_mpfr_acos[]="\
x.acos(): returns arc-cosine of x; x in radians.\n\
";
static char doc_g_mpfr_acos[]="\
acos(x): returns arc-cosine of x; x in radians.\n\
";

MPFR_UNIOP(acos)

static char doc_mpfr_asin[]="\
x.asin(): returns arc-sine of x; x in radians.\n\
";
static char doc_g_mpfr_asin[]="\
asin(x): returns arc-sine of x; x in radians.\n\
";

MPFR_UNIOP(asin)

static char doc_mpfr_atan[]="\
x.atan(): returns arc-tangent of x; x in radians.\n\
";
static char doc_g_mpfr_atan[]="\
atan(x): returns arc-tangent of x; x in radians.\n\
";

MPFR_UNIOP(atan)

static char doc_mpfr_cosh[]="\
x.cosh(): returns hyperbolic cosine of x.\n\
";
static char doc_g_mpfr_cosh[]="\
cosh(x): returns hyperbolic cosine of x.\n\
";

MPFR_UNIOP(cosh)

static char doc_mpfr_sinh[]="\
x.sinh(): returns hyperbolic sine of x.\n\
";
static char doc_g_mpfr_sinh[]="\
sinh(x): returns hyperbolic sine of x.\n\
";

MPFR_UNIOP(sinh)

static char doc_mpfr_tanh[]="\
x.tanh(): returns hyperbolic tangent of x.\n\
";
static char doc_g_mpfr_tanh[]="\
tanh(x): returns hyperbolic tangent of x.\n\
";

MPFR_UNIOP(tanh)

static char doc_mpfr_sech[]="\
x.sech(): returns hyperbolic secant of x.\n\
";
static char doc_g_mpfr_sech[]="\
sech(x): returns hyperbolic secant of x.\n\
";

MPFR_UNIOP(sech)

static char doc_mpfr_csch[]="\
x.csch(): returns hyperbolic cosecant of x.\n\
";
static char doc_g_mpfr_csch[]="\
csch(x): returns hyperbolic cosecant of x.\n\
";

MPFR_UNIOP(csch)

static char doc_mpfr_coth[]="\
x.coth(): returns hyperbolic cotangent of x.\n\
";
static char doc_g_mpfr_coth[]="\
coth(x): returns hyperbolic cotangent of x.\n\
";

MPFR_UNIOP(coth)

static char doc_mpfr_acosh[]="\
x.acosh(): returns inverse hyperbolic cosine of x.\n\
";
static char doc_g_mpfr_acosh[]="\
acosh(x): returns inverse hyperbolic cosine of x.\n\
";

MPFR_UNIOP(acosh)

static char doc_mpfr_asinh[]="\
x.asinh(): returns inverse hyperbolic sine of x.\n\
";
static char doc_g_mpfr_asinh[]="\
asinh(x): returns inverse hyperbolic sine of x.\n\
";

MPFR_UNIOP(asinh)

static char doc_mpfr_atanh[]="\
x.atanh(): returns inverse hyperbolic tangent of x.\n\
";
static char doc_g_mpfr_atanh[]="\
atanh(x): returns inverse hyperbolic tangent of x.\n\
";

MPFR_UNIOP(atanh)

static char doc_mpfr_log1p[]="\
x.log1p(): returns logarithm of (1+x).\n\
";
static char doc_g_mpfr_log1p[]="\
log1p(x): returns logarithm of (1+x).\n\
";

MPFR_UNIOP(log1p)

static char doc_mpfr_expm1[]="\
x.expm1(): returns exponential(x) - 1.\n\
";
static char doc_g_mpfr_expm1[]="\
expm1(x): returns exponential(x) - 1.\n\
";

MPFR_UNIOP(expm1)

static char doc_mpfr_eint[]="\
x.eint(): returns exponential integral of x.\n\
";
static char doc_g_mpfr_eint[]="\
eint(x): returns exponential integral of x.\n\
";

MPFR_UNIOP(eint)

static char doc_mpfr_li2[]="\
x.li2(): returns real part of dilogarithm of x.\n\
";
static char doc_g_mpfr_li2[]="\
li2(x): returns real part of dilogarithm of x.\n\
";

MPFR_UNIOP(li2)

static char doc_mpfr_gamma[]="\
x.gamma(): returns gamma of x.\n\
";
static char doc_g_mpfr_gamma[]="\
gamma(x): returns gamma of x.\n\
";

MPFR_UNIOP(gamma)

static char doc_mpfr_lngamma[]="\
x.lngamma(): returns logarithm of gamma(x).\n\
";
static char doc_g_mpfr_lngamma[]="\
lngamma(x): returns logarithm of gamma(x).\n\
";

MPFR_UNIOP(lngamma)

PyDoc_STRVAR(doc_mpfr_lgamma,
"x.lgamma() -> (mpfr, int)\n\n"
"Return a 2-tuple containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

PyDoc_STRVAR(doc_g_mpfr_lgamma,
"lgamma(x) -> (mpfr, int)\n\n"
"Return a 2-tuple containing the logarithm of the absolute value of\n"
"gamma(x) and the sign of gamma(x)");

static PyObject *
Pympfr_lgamma(PyObject* self, PyObject *other)
{
    PyObject *result;
    PympfrObject *value;
    int signp = 0;

    PARSE_ONE_MPFR_OTHER("lgamma() requires 'mpfr' argument");

    value = Pympfr_new(0);
    result = PyTuple_New(2);
    if (!value || !result)
        goto done;

    mpfr_clear_flags();
    value->rc = mpfr_lgamma(value->f, &signp, Pympfr_AS_MPFR(self),
                            context->now.mpfr_round);
    SUBNORMALIZE(value);
    MERGE_FLAGS;
    CHECK_FLAGS("lgamma()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)value);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)value);
        PyTuple_SET_ITEM(result, 1, PyIntOrLong_FromLong((long)signp));
    }
    return result;
}

static char doc_mpfr_digamma[]="\
x.digamma(): returns digamma of x.\n\
";
static char doc_g_mpfr_digamma[]="\
digamma(x): returns digamma of x.\n\
";

MPFR_UNIOP(digamma)

static char doc_mpfr_zeta[]="\
x.zeta(): returns Riemann zeta of x.\n\
";
static char doc_g_mpfr_zeta[]="\
zeta(x): returns Riemann zeta of x.\n\
";

MPFR_UNIOP(zeta)

static char doc_mpfr_erf[]="\
x.erf(): returns error function of x.\n\
";
static char doc_g_mpfr_erf[]="\
erf(x): returns error function of x.\n\
";

MPFR_UNIOP(erf)

static char doc_mpfr_erfc[]="\
x.erfc(): returns complementary error function of x.\n\
";
static char doc_g_mpfr_erfc[]="\
erfc(x): returns complementary error function of x.\n\
";

MPFR_UNIOP(erfc)

static char doc_mpfr_j0[]="\
x.j0(): returns first kind Bessel function of order 0 of x.\n\
";
static char doc_g_mpfr_j0[]="\
j0(x): returns first kind Bessel function of order 0 of x.\n\
";

MPFR_UNIOP(j0)

static char doc_mpfr_j1[]="\
x.j1(): returns first kind Bessel function of order 1 of x.\n\
";
static char doc_g_mpfr_j1[]="\
j1(x): returns first kind Bessel function of order 1 of x.\n\
";

MPFR_UNIOP(j1)

PyDoc_STRVAR(doc_mpfr_jn,
"x.jn(n) -> mpfr\n\n"
"Return the first kind Bessel function of order n of x.");

PyDoc_STRVAR(doc_g_mpfr_jn,
"jn(x,n) -> mpfr\n\n"
"Return the first kind Bessel function of order n of x.");

static PyObject *
Pympfr_jn(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    long n = 0;

    PARSE_ONE_MPFR_REQ_CLONG(&n, "jn() requires 'mpfr','int' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_jn(result->f, n, Pympfr_AS_MPFR(self),
                         context->now.mpfr_round);
    MPFR_CLEANUP_SELF("jn()");
}

static char doc_mpfr_y0[]="\
x.y0(): returns second kind Bessel function of order 0 of x.\n\
";
static char doc_g_mpfr_y0[]="\
y0(x): returns second kind Bessel function of order 0 of x.\n\
";

MPFR_UNIOP(y0)

PyDoc_STRVAR(doc_mpfr_y1,
"x.y1() -> mpfr\n\n"
"Return second kind Bessel function of order 1 of x.");

PyDoc_STRVAR(doc_g_mpfr_y1,
"y1(x) -> mpfr\n\n"
"Return second kind Bessel function of order 1 of x.");

MPFR_UNIOP(y1)

PyDoc_STRVAR(doc_mpfr_yn,
"x.yn(n) -> mpfr\n\n"
"Return the second kind Bessel function of order n of x.");

PyDoc_STRVAR(doc_g_mpfr_yn,
"yn(x,n) -> mpfr\n\n"
"Return the second kind Bessel function of order n of x.");

static PyObject *
Pympfr_yn(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    long n = 0;

    PARSE_ONE_MPFR_REQ_CLONG(&n, "yn() requires 'mpfr','int' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_yn(result->f, n, Pympfr_AS_MPFR(self),
                         context->now.mpfr_round);
    MPFR_CLEANUP_SELF("yn()");
}

PyDoc_STRVAR(doc_mpfr_ai,
"x.ai() -> mpfr\n\n"
"Return Airy function of x.");

PyDoc_STRVAR(doc_g_mpfr_ai,
"ai(x) -> mpfr\n\n"
"Return Airy function of x.");

MPFR_UNIOP(ai)

PyDoc_STRVAR(doc_mpfr_add,
"x.add(y) -> mpfr\n\n"
"Return x + y.");

PyDoc_STRVAR(doc_g_mpfr_add,
"add(x, y) -> mpfr\n\n"
"Return x + y.");

static PyObject *
Pympfr_add(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "add() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_add(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("add()");
}

PyDoc_STRVAR(doc_mpfr_sub,
"x.sub(y) -> mpfr\n\n"
"Return x - y.");

PyDoc_STRVAR(doc_g_mpfr_sub,
"sub(x, y) -> mpfr\n\n"
"Return x - y.");

static PyObject *
Pympfr_sub(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "sub() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_sub(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("sub()");
}

PyDoc_STRVAR(doc_mpfr_mul,
"x.mul(y) -> mpfr\n\n"
"Return x * y.");

PyDoc_STRVAR(doc_g_mpfr_mul,
"mul(x, y) -> mpfr\n\n"
"Return x * y.");

static PyObject *
Pympfr_mul(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "mul() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_mul(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("mul()");
}

PyDoc_STRVAR(doc_mpfr_div,
"x.div(y) -> mpfr\n\n"
"Return x / y.");

PyDoc_STRVAR(doc_g_mpfr_div,
"div(x, y) -> mpfr\n\n"
"Return x / y.");

static PyObject *
Pympfr_div(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "div() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    if (mpfr_zero_p(Pympfr_AS_MPFR(other))) {
        context->now.divzero = 1;
        if (context->now.trap_divzero) {
            GMPY_DIVZERO("'mpfr' division by zero");
            goto done;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_div(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("div()");
}

PyDoc_STRVAR(doc_mpfr_fmod,
"x.fmod(y) -> mpfr\n\n"
"Return x - n * y where n is the integer quotient of x/y, rounded to 0.");

PyDoc_STRVAR(doc_g_mpfr_fmod,
"fmod(x, y) -> mpfr\n\n"
"Return x - n * y where n is the integer quotient of x/y, rounded to 0.");

static PyObject *
Pympfr_fmod(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "fmod() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    if (mpfr_zero_p(Pympfr_AS_MPFR(other))) {
        context->now.divzero = 1;
        if (context->now.trap_divzero) {
            GMPY_DIVZERO("'mpfr' division by zero");
            goto done;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_fmod(result->f, Pympfr_AS_MPFR(self),
                           Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("fmod()");
}

PyDoc_STRVAR(doc_mpfr_remainder,
"x.remainder(y) -> mpfr\n\n"
"Return x - n * y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

PyDoc_STRVAR(doc_g_mpfr_remainder,
"remainder(x, y) -> mpfr\n\n"
"Return x - n * y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

static PyObject *
Pympfr_remainder(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "remainder() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    if (mpfr_zero_p(Pympfr_AS_MPFR(other))) {
        context->now.divzero = 1;
        if (context->now.trap_divzero) {
            GMPY_DIVZERO("'mpfr' remainder by zero");
            goto done;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_remainder(result->f, Pympfr_AS_MPFR(self),
                                Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("remainder()");
}

PyDoc_STRVAR(doc_mpfr_remquo,
"x.remquo() -> (mpfr, int)\n\n"
"Return a 2-tuple containing the remainder(x,y) and the low bits of the\n"
"quotient.");

PyDoc_STRVAR(doc_g_mpfr_remquo,
"remquo(x) -> (mpfr, int)\n\n"
"Return a 2-tuple containing the remainder(x,y) and the low bits of the\n"
"quotient.");

static PyObject *
Pympfr_remquo(PyObject* self, PyObject *args)
{
    PyObject *result, *other;
    PympfrObject *value;
    long quobits = 0;

    PARSE_TWO_MPFR(other, "remquo() requires 'mpfr', 'mpfr' argument");

    value = Pympfr_new(0);
    result = PyTuple_New(2);
    if (!value || !result)
        goto done;

    mpfr_clear_flags();
    value->rc = mpfr_remquo(value->f, &quobits, Pympfr_AS_MPFR(self),
                            Pympfr_AS_MPFR(other), context->now.mpfr_round);
    SUBNORMALIZE(value);
    MERGE_FLAGS;
    CHECK_FLAGS("remquo()");

  done:
    Py_DECREF(self);
    Py_DECREF(other);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)value);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)value);
        PyTuple_SET_ITEM(result, 1, PyIntOrLong_FromLong(quobits));
    }
    return result;
}

PyDoc_STRVAR(doc_mpfr_pow,
"x.pow(y) -> mpfr\n\n"
"Return x ** y.");

PyDoc_STRVAR(doc_g_mpfr_pow,
"pow(x, y) -> mpfr\n\n"
"Return x ** y.");

static PyObject *
Pympfr_pow(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "pow() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    if ((mpfr_zero_p(Pympfr_AS_MPFR(self))) &&
        (mpfr_sgn(Pympfr_AS_MPFR(other)) < 0)) {
        context->now.divzero = 1;
        if (context->now.trap_divzero) {
            GMPY_DIVZERO("zero cannot be raised to a negative power");
            goto done;
        }
    }

    mpfr_clear_flags();
    result->rc = mpfr_pow(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("pow()");
}

PyDoc_STRVAR(doc_mpfr_atan2,
"y.atan2(x) -> mpfr\n\n"
"Return arc-tangent of (y/x).");

PyDoc_STRVAR(doc_g_mpfr_atan2,
"atan2(y, x) -> mpfr\n\n"
"Return arc-tangent of (y/x).");

static PyObject *
Pympfr_atan2(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "atan2() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_atan2(result->f, Pympfr_AS_MPFR(self),
                            Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("atan2()");
}

PyDoc_STRVAR(doc_mpfr_agm,
"x.agm(y) -> mpfr\n\n"
"Return arithmetic-geometric mean of x and y.");

PyDoc_STRVAR(doc_g_mpfr_agm,
"agm(x, y) -> mpfr\n\n"
"Return arithmetic-geometric mean of x and y.");

static PyObject *
Pympfr_agm(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "agm() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_agm(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("agm()");
}

PyDoc_STRVAR(doc_mpfr_hypot,
"y.hypot(x) -> mpfr\n\n"
"Return square root of (x**2 + y**2).");

PyDoc_STRVAR(doc_g_mpfr_hypot,
"hypot(y, x) -> mpfr\n\n"
"Return square root of (x**2 + y**2).");

static PyObject *
Pympfr_hypot(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "hypot() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_hypot(result->f, Pympfr_AS_MPFR(self),
                            Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("hypot()");
}

PyDoc_STRVAR(doc_mpfr_max,
"y.max(x) -> mpfr\n\n"
"Return maximum of x and y.");

PyDoc_STRVAR(doc_g_mpfr_max,
"max(y, x) -> mpfr\n\n"
"Return maximum of x and y.");

static PyObject *
Pympfr_max(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "max() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_max(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("max()");
}

PyDoc_STRVAR(doc_mpfr_min,
"y.min(x) -> mpfr\n\n"
"Return minimum of x and y.");

PyDoc_STRVAR(doc_g_mpfr_min,
"min(y, x) -> mpfr\n\n"
"Return minimum of x and y.");

static PyObject *
Pympfr_min(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "min() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_min(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("min()");
}

PyDoc_STRVAR(doc_mpfr_nexttoward,
"x.next_toward(y) -> mpfr\n\n"
"Return the next mpfr from x in the direction of y.");

PyDoc_STRVAR(doc_g_mpfr_nexttoward,
"next_toward(y, x) -> mpfr\n\n"
"Return the next mpfr from x in the direction of y.");

static PyObject *
Pympfr_nexttoward(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR(other, "next_toward() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(self)))))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    mpfr_nexttoward(result->f, Pympfr_AS_MPFR(other));
    result->rc = 0;
    MPFR_CLEANUP_SELF_OTHER("next_toward()");
}

PyDoc_STRVAR(doc_mpfr_nextabove,
"x.next_above() -> mpfr\n\n"
"Return the next mpfr from x toward +Infinity.");

PyDoc_STRVAR(doc_g_mpfr_nextabove,
"next_above(x) -> mpfr\n\n"
"Return the next mpfr from x toward +Infinity.");

static PyObject *
Pympfr_nextabove(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("next_above() requires 'mpfr' argument");

    if (!(result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(self)))))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    mpfr_nextabove(result->f);
    result->rc = 0;
    MPFR_CLEANUP_SELF("next_above()");
}

PyDoc_STRVAR(doc_mpfr_nextbelow,
"x.next_below() -> mpfr\n\n"
"Return the next mpfr from x toward -Infinity.");

PyDoc_STRVAR(doc_g_mpfr_nextbelow,
"next_below(x) -> mpfr\n\n"
"Return the next mpfr from x toward -Infinity.");

static PyObject *
Pympfr_nextbelow(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("next_below() requires 'mpfr' argument");

    if (!(result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(self)))))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    mpfr_nextbelow(result->f);
    result->rc = 0;
    MPFR_CLEANUP_SELF("next_below()");
}

PyDoc_STRVAR(doc_mpfr_sin_cos,
"x.sin_cos() -> (mpfr, mpfr)\n\n"
"Return a tuple containing the sine and cosine of x.");

PyDoc_STRVAR(doc_g_mpfr_sin_cos,
"sin_cos(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the sine and cosine of x.");

static PyObject *
Pympfr_sin_cos(PyObject *self, PyObject *other)
{
    PympfrObject *s, *c;
    PyObject *result;
    int code;

    PARSE_ONE_MPFR_OTHER("sin_cos() requires 'mpfr' argument");

    s = Pympfr_new(0);
    c = Pympfr_new(0);
    result = PyTuple_New(2);
    if (!s || !c || !result)
        goto done;

    mpfr_clear_flags();
    code = mpfr_sin_cos(s->f, c->f, Pympfr_AS_MPFR(self),
                        context->now.mpfr_round);
    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;
    SUBNORMALIZE(s);
    SUBNORMALIZE(c);
    MERGE_FLAGS;
    CHECK_FLAGS("sin_cos()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)s);
        PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    }
    return result;
}

PyDoc_STRVAR(doc_mpfr_sinh_cosh,
"x.sinh_cosh() -> (mpfr, mpfr)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

PyDoc_STRVAR(doc_g_mpfr_sinh_cosh,
"sinh_cosh(x) -> (mpfr, mpfr)\n\n"
"Return a tuple containing the hyperbolic sine and cosine of x.");

static PyObject *
Pympfr_sinh_cosh(PyObject *self, PyObject *other)
{
    PympfrObject *s, *c;
    PyObject *result;
    int code;

    PARSE_ONE_MPFR_OTHER("sinh_cosh() requires 'mpfr' argument");

    s = Pympfr_new(0);
    c = Pympfr_new(0);
    result = PyTuple_New(2);
    if (!s || !c || !result)
        goto done;

    mpfr_clear_flags();
    code = mpfr_sinh_cosh(s->f, c->f, Pympfr_AS_MPFR(self),
                          context->now.mpfr_round);
    s->rc = code & 0x03;
    c->rc = code >> 2;
    if (s->rc == 2) s->rc = -1;
    if (c->rc == 2) c->rc = -1;
    SUBNORMALIZE(s);
    SUBNORMALIZE(c);
    MERGE_FLAGS;
    CHECK_FLAGS("sin_cos()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)s);
        PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    }
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_fma,
"fma(x,y,z) -> mpfr\n\n"
"Return correctly rounded result of (x * y) + z.");

static PyObject *
Pympfr_fma(PyObject *self, PyObject *args)
{
    PympfrObject *result, *x, *y, *z;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fma() requires 'mpfr','mpfr','mpfr' arguments.");
        return NULL;
    }

    result = Pympfr_new(0);
    x = Pympfr_From_Real(PyTuple_GET_ITEM(args, 0), 0);
    y = Pympfr_From_Real(PyTuple_GET_ITEM(args, 1), 0);
    z = Pympfr_From_Real(PyTuple_GET_ITEM(args, 2), 0);
    if (!result || !x || !y || !z) {
        TYPE_ERROR("fma() requires 'mpfr','mpfr','mpfr' arguments.");
        goto done;
    }

    mpfr_clear_flags();
    result->rc = mpfr_fma(result->f, x->f, y->f, z->f,
                          context->now.mpfr_round);
    SUBNORMALIZE(result);
    MERGE_FLAGS;
    CHECK_FLAGS("fma()");

  done:
    Py_XDECREF((PyObject*)x);
    Py_XDECREF((PyObject*)y);
    Py_XDECREF((PyObject*)z);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_fms,
"fms(x,y,z) -> mpfr\n\n"
"Return correctly rounded result of (x * y) - z.");

static PyObject *
Pympfr_fms(PyObject *self, PyObject *args)
{
    PympfrObject *result, *x, *y, *z;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fms() requires 'mpfr','mpfr','mpfr' arguments.");
        return NULL;
    }

    result = Pympfr_new(0);
    x = Pympfr_From_Real(PyTuple_GET_ITEM(args, 0), 0);
    y = Pympfr_From_Real(PyTuple_GET_ITEM(args, 1), 0);
    z = Pympfr_From_Real(PyTuple_GET_ITEM(args, 2), 0);
    if (!result || !x || !y || !z) {
        TYPE_ERROR("fms() requires 'mpfr','mpfr','mpfr' arguments.");
        goto done;
    }

    mpfr_clear_flags();
    result->rc = mpfr_fms(result->f, x->f, y->f, z->f,
                          context->now.mpfr_round);
    SUBNORMALIZE(result);
    MERGE_FLAGS;
    CHECK_FLAGS("fms()");

  done:
    Py_XDECREF((PyObject*)x);
    Py_XDECREF((PyObject*)y);
    Py_XDECREF((PyObject*)z);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_factorial,
"factorial(n) -> mpfr\n\n"
"Return the floating-point approximation to the factorial of n.\n"
"See fac(n) to get the exact integer result.");

static PyObject *
Pympfr_factorial(PyObject *self, PyObject *other)
{
    PympfrObject *result;
    long n;

    n = clong_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("factorial() requires 'int' argument");
        return NULL;
    }

    if (n < 0) {
        VALUE_ERROR("factorial() of negative number");
        return NULL;
    }

    if (!(result = Pympfr_new(0)))
        return NULL;

    mpfr_clear_flags();
    mpfr_fac_ui(result->f, n, context->now.mpfr_round);

    MERGE_FLAGS;
    CHECK_FLAGS("factorial()");
  done:
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpfr_is_lessgreater,
"x.is_lessgreater(y) -> boolean\n\n"
"Return True if x > y or x < y. Return False if x == y or either x\n"
"and/or y is NaN.");

PyDoc_STRVAR(doc_g_mpfr_is_lessgreater,
"is_lessgreater(x,y) -> boolean\n\n"
"Return True if x > y or x < y. Return False if x == y or either x\n"
"and/or y is NaN.");

static PyObject *
Pympfr_is_lessgreater(PyObject *self, PyObject *args)
{
    PyObject *other;
    int temp;

    PARSE_TWO_MPFR(other, "is_lessgreater() requires 'mpfr','mpfr' arguments");

    temp = mpfr_lessgreater_p(Pympfr_AS_MPFR(self), Pympfr_AS_MPFR(other));
    Py_DECREF(self);
    Py_DECREF(other);
    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_mpfr_is_unordered,
"x.is_unordered(y) -> boolean\n\n"
"Return True if either x and/or y is NaN.");

PyDoc_STRVAR(doc_g_mpfr_is_unordered,
"is_unordered(x,y) -> boolean\n\n"
"Return True if either x and/or y is NaN.");

static PyObject *
Pympfr_is_unordered(PyObject *self, PyObject *args)
{
    PyObject *other;
    int temp;

    PARSE_TWO_MPFR(other, "unordered() requires 'mpfr','mpfr' arguments");

    temp = mpfr_unordered_p(Pympfr_AS_MPFR(self), Pympfr_AS_MPFR(other));
    Py_DECREF(self);
    Py_DECREF(other);
    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_mpfr_check_range,
"x.check_range() -> mpfr\n\n"
"Return a new 'mpfr' with exponent that lies within the current range\n"
"of emin and emax.");

PyDoc_STRVAR(doc_g_mpfr_check_range,
"check_range(x) -> mpfr\n\n"
"Return a new 'mpfr' with exponent that lies within the current range\n"
"of emin and emax.");

static PyObject *
Pympfr_check_range(PyObject *self, PyObject *other)
{
    PympfrObject *result = NULL;

    if (self && Pympfr_Check(self)) {
        if ((result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(self))))) {
            mpfr_set(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
            result->round_mode = ((PympfrObject*)self)->round_mode;
            result->rc = ((PympfrObject*)self)->rc;
            mpfr_clear_flags();
            result->rc = mpfr_check_range(result->f, result->rc,
                                          result->round_mode);
        }
    }
    else if (Pympfr_Check(other)) {
        if ((result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(other))))) {
            mpfr_set(result->f, Pympfr_AS_MPFR(other), context->now.mpfr_round);
            result->round_mode = ((PympfrObject*)other)->round_mode;
            result->rc = ((PympfrObject*)other)->rc;
            mpfr_clear_flags();
            result->rc = mpfr_check_range(result->f, result->rc,
                                          result->round_mode);
        }
    }
    else {
        TYPE_ERROR("check_range() requires 'mpfr' argument");
    }
    MERGE_FLAGS;
    CHECK_FLAGS("check_range()");
  done:
    return (PyObject*)result;
}

