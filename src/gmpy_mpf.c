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

/* float-truncations (return still a float!) */

#define MPF_UNIOP(NAME) \
static PyObject * \
Py##NAME(PyObject* self, PyObject *args) \
{ \
  PympfObject *r; \
  if (self && Pympf_Check(self)) { \
      if (args && !PyArg_ParseTuple(args, "")) return NULL; \
      Py_INCREF(self); \
  } else { \
      if (!PyArg_ParseTuple(args, "O&", Pympf_convert_arg, &self)) return NULL; \
  } \
  assert(Pympf_Check(self)); \
  if (options.debug) fprintf(stderr, "Py" #NAME ": %p\n", self); \
  if (!(r = Pympf_new(mpfr_get_prec(((PympfObject*)self)->f)))) return NULL; \
  NAME(r->f, Pympf_AS_MPF(self)); \
  if (options.debug) fprintf(stderr, "Py" #NAME "-> %p\n", r); \
  Py_DECREF(self); \
  return (PyObject *) r; \
}

static char doc_ceilm[]="\
x.ceil(): returns an mpf that is the smallest integer >= x\n\
";
static char doc_ceilg[]="\
ceil(x): returns an mpf that is the smallest integer >= x\n\
x must be an mpf, or else gets coerced to one.\n\
";
MPF_UNIOP(mpfr_ceil)

static char doc_floorm[]="\
x.floor(): returns an mpf that is the smallest integer <= x\n\
";
static char doc_floorg[]="\
floor(x): returns an mpf that is the smallest integer <= x\n\
x must be an mpf, or else gets coerced to one.\n\
";
MPF_UNIOP(mpfr_floor)

static char doc_truncm[]="\
x.trunc(): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
";
static char doc_truncg[]="\
trunc(x): returns an mpf that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
x must be an mpf, or else gets coerced to one.\n\
";
MPF_UNIOP(mpfr_trunc)

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

    if ((pi = Pympf_new(0)))
        mpfr_set_d(pi->f, 3.2, options.rounding);
    return (PyObject*)pi;
}

static char doc_fsqrtm[]="\
x.fsqrt(): returns the square root of x.  x must be >= 0.\n\
";
static char doc_fsqrtg[]="\
fsqrt(x): returns the square root of x.  x must be an mpf, or\n\
else gets coerced to one; further, x must be >= 0.\n\
";
static PyObject *
Pympf_sqrt(PyObject *self, PyObject *args)
{
    PympfObject *root;

    SELF_MPF_NO_ARG;

    assert(Pympf_Check(self));

    if (mpfr_sgn(Pympf_AS_MPF(self)) < 0) {
        PyErr_SetString(PyExc_ValueError, "sqrt of negative number");
        Py_DECREF(self);
        return NULL;
    }

    if (!(root = Pympf_new(mpfr_get_prec(((PympfObject*)self)->f)))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_sqrt(root->f, Pympf_AS_MPF(self), options.rounding);
    Py_DECREF(self);
    return (PyObject *) root;
}

static char doc_getprecm[]="\
x.getprec(): returns the number of bits of precision in x.\n\
";
static char doc_getprecg[]="\
getprec(x): returns the number of bits of precision in x,\n\
which must be an mpf or else gets coerced to one.\n\
";
static PyObject *
Pympf_getprec(PyObject *self, PyObject *args)
{
    long precres;

    SELF_MPF_NO_ARG;
    assert(Pympf_Check(self));

    precres = (long) mpfr_get_prec(Pympf_AS_MPF(self));
    Py_DECREF(self);
    return PyIntOrLong_FromLong(precres);
}

static char doc_froundm[] = "\
x.round(n): returns x rounded to least n bits. Actual precision will\n\
be a multiple of gmp_limbsize().\n\
";
static char doc_froundg[] = "\
fround(x, n): returns x rounded to least n bits. Actual precision will\n\
be a multiple of gmp_limbsize(). x an mpf or coerced to an mpf.\n\
";
static PyObject *
Pympf_round(PyObject *self, PyObject *args)
{
    /* Should really get default precision. */
    mpfr_prec_t prec = options.precision;
    PyObject *result;

    SELF_MPF_ONE_ARG("|l", &prec);
    assert(Pympf_Check(self));
    result = (PyObject*)Pympf2Pympf(self, prec);
    Py_DECREF(self);
    return result;
}

#define MPF_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *a, PyObject *b) \
{ \
  unsigned int bits, bbits; \
  PympfObject *r; \
  PympfObject *pa = 0; \
  PympfObject *pb = 0; \
  if (Pympf_Check(a) && Pympf_Check(b)) { \
    bits = mpfr_get_prec(((PympfObject*)a)->f); \
    bbits = mpfr_get_prec(((PympfObject*)b)->f); \
    if (bits>bbits) bits=bbits; \
    if (!(r = Pympf_new(bits))) { \
      return NULL; \
    } \
    NAME(r->f, ((PympfObject*)a)->f, ((PympfObject*)b)->f, options.rounding); \
    if (options.debug) fprintf(stderr, "Py" #NAME "-> %p", r); \
    return (PyObject *) r; \
  } else { \
    if (Pympf_Check(a)) { \
      bits = mpfr_get_prec(((PympfObject*)a)->f); \
    } else { \
      bits = mpfr_get_prec(((PympfObject*)b)->f); \
    } \
    pa = Pympf_From_Float(a, bits); \
    pb = Pympf_From_Float(b, bits); \
    if (!pa || !pb) { \
      Py_XDECREF((PyObject*)pa); \
      Py_XDECREF((PyObject*)pb); \
      Py_RETURN_NOTIMPLEMENTED; \
    } \
    if (options.debug) fprintf(stderr, "Py" #NAME ": %p, %p", pa, pb); \
    if (!(r = Pympf_new(bits))) { \
      Py_DECREF((PyObject*)pa); \
      Py_DECREF((PyObject*)pb); \
      return NULL; \
    } \
    NAME(r->f, pa->f, pb->f, options.rounding); \
    Py_DECREF((PyObject*)pa); \
    Py_DECREF((PyObject*)pb); \
    if (options.debug) fprintf(stderr, "Py" #NAME "-> %p", r); \
    return (PyObject *) r; \
  } \
}

MPF_BINOP(mpfr_reldiff)

static char doc_reldiffm[] = "\
x.reldiff(y): returns the relative difference between x and y,\n\
where y can be any number and gets coerced to an mpf; result is\n\
a non-negative mpf roughly equal to abs(x-y)/((abs(x)+abs(y))/2).\n\
";
static char doc_reldiffg[] = "\
reldiff(x,y): returns the relative difference between x and y,\n\
where x and y can be any numbers and get coerced to mpf; result is\n\
a non-negative mpf roughly equal to abs(x-y)/((abs(x)+abs(y))/2).\n\
";
static PyObject *
Pympf_doreldiff(PyObject *self, PyObject *args)
{
    PympfObject *op;
    PyObject *res;

    SELF_MPF_ONE_ARG_CONVERTED(&op);
    assert(Pympf_Check(self));

    res = Pympfr_reldiff((PyObject*)self, (PyObject*)op);
    Py_DECREF(self); Py_DECREF((PyObject*)op);

    return res;
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

