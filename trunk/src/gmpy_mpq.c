/* gmpy_mpq.c
 *
 * Functions that operate strictly on mpq.
 *
 * This file should be considered part of gmpy2.c
 */

/* produce digits for an mpq in requested base, default 10 */
static char doc_qdigitsm[]="\
x.digits([base]): returns Python string representing x in the\n\
given base (2 to 36, default 10 if omitted or 0); leading '-'\n\
is present if x<0, but no leading '+' if x>=0.\n\
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

static char doc_qsignm[]="\
x.sign(): returns -1, 0, or +1, if x is negative, 0, positive.\n\
";
static PyObject *
Pympq_sign(PyObject *self, PyObject *other)
{
    long res;
    PympqObject* tempx;

    if (self && (Pympq_Check(self))) {
        res = mpq_sgn(Pympq_AS_MPQ(self));
    }
    else if (Pympq_Check(other)) {
        res = mpq_sgn(Pympq_AS_MPQ(other));
    }
    else {
        if (!(tempx = Pympq_From_Rational(other))) {
            TYPE_ERROR("sign() requires 'mpq' argument");
            return NULL;
        }
        else {
            res = mpq_sgn(tempx->q);
            Py_DECREF((PyObject*)tempx);
        }
    }
    return PyIntOrLong_FromLong(res);
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
    self = (PyObject*)Pympq_From_Rational(self);
    if (!self) {
        if (!PyErr_Occurred())
            TYPE_ERROR("first argument can not be converted to mpq");
        return NULL;
    }
    if (wasone) { /* self was mpf, float, int, long... */
        s = self;
    }
    else {     /* other explicitly present and !=1... must compute */
        other = (PyObject*)Pympq_From_Rational(other);
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

#define MPQ_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *a, PyObject *b) \
{ \
  PympqObject *r; \
  PympqObject *pa = 0; \
  PympqObject *pb = 0; \
  pa = Pympq_From_Rational(a); \
  pb = Pympq_From_Rational(b); \
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
    PympqObject *b = Pympq_From_Rational(in_b);
    PympqObject *e = Pympq_From_Rational(in_e);

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

static int
Pympq_nonzero(PympqObject *x)
{
    return mpq_sgn(x->q) != 0;
}

static Py_hash_t
Pympq_hash(PympqObject *self)
{
#ifdef _PyHASH_MODULUS
    Py_hash_t hash = 0;
    mpz_t temp, temp1, mask;

    if (self->hash_cache != -1)
        return self->hash_cache;

    mpz_inoc(temp);
    mpz_inoc(temp1);
    mpz_inoc(mask);
    mpz_set_si(mask, 1);
    mpz_mul_2exp(mask, mask, _PyHASH_BITS);
    mpz_sub_ui(mask, mask, 1);

    if (!mpz_invert(temp, mpq_denref(self->q), mask)) {
        mpz_cloc(temp);
        mpz_cloc(temp1);
        mpz_cloc(mask);
        hash = _PyHASH_INF;
        if (mpz_sgn(mpq_numref(self->q))<0)
            hash = -hash;
        self->hash_cache = hash;
        return hash;
    }
    mpz_set(temp1, mask);
    mpz_sub_ui(temp1, temp1, 2);
    mpz_powm(temp, mpq_denref(self->q), temp1, mask);
    
    mpz_tdiv_r(temp1, mpq_numref(self->q), mask);
    mpz_mul(temp, temp, temp1);
    hash = (Py_hash_t)mpn_mod_1(temp->_mp_d, mpz_size(temp), _PyHASH_MODULUS);

    if (mpz_sgn(mpq_numref(self->q))<0)
        hash = -hash;
    if (hash==-1) hash = -2;
    mpz_cloc(temp);
    mpz_cloc(temp1);
    mpz_cloc(mask);
    self->hash_cache = hash;
    return hash;
#else
    if (self->hash_cache != -1)
        return self->hash_cache;
    return (self->hash_cache = dohash(Pympq2PyFloat(self)));
#endif
}

