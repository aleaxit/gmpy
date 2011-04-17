/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008 - 2011 Case Van Horsen                              *
 *                                                                         *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2.1 of the License, or  *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This library is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * Lesser General Public License for more details.                         *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with this library; if not, write to the Free Software     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA           *
 * 02110-1301  USA                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/* Functions that operate strictly on mpq. */

/* produce digits for an mpq in requested base, default 10 */
static char doc_qdigitsm[]="\
x.digits([base]): returns Python string representing x in the\n\
given base (2 to 62, default 10 if omitted or 0); leading '-'\n\
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
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        return mpfr_get_d(Pympfr_AS_MPFR(obj), context->now.mpfr_round)==1.0;
    }
#endif
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
  if (!(r = Pympq_new())) { \
    Py_DECREF((PyObject*)pa); \
    Py_DECREF((PyObject*)pb); \
    return NULL; \
  } \
  NAME(r->q, pa->q, pb->q); \
  Py_DECREF((PyObject*)pa); \
  Py_DECREF((PyObject*)pb); \
  return (PyObject *) r; \
}

#define MPQ_MONOP(NAME) \
static PyObject * \
Py##NAME(PympqObject *x) \
{ \
  PympqObject *r; \
  if (!(r = Pympq_new())) return NULL; \
  NAME(r->q, x->q); \
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

static PyObject *
Pympq_pos(PympqObject *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject *) x;
}

static PyObject *
Pympq_pow(PyObject *base, PyObject *exp, PyObject *m)
{
    PympqObject *rq, *tempbq;
    PympzObject *tempez;
#ifdef WITHMPFR
    PympfrObject *rf, *tempbf, *tempef;
#endif
    int esign, bsign;
    long tempexp;

    if ((PyObject*)m != Py_None) {
        TYPE_ERROR("mpq.pow() no modulo allowed");
        return NULL;
    }

    /* Only support mpq**int. Everything else gets converted to mpf. */
    if (isRational(base) && isInteger(exp)) {
        tempbq = Pympq_From_Rational(base);
        tempez = Pympz_From_Integer(exp);
        if (!tempbq || !tempez) {
            Py_XDECREF((PyObject*)tempbq);
            Py_XDECREF((PyObject*)tempez);
            return NULL;
        }
        if (!mpz_fits_slong_p(tempez->z)) {
            VALUE_ERROR("mpq.pow() outrageous exponent");
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return NULL;
        }
        esign = mpz_sgn(tempez->z);
        if (esign == 0) {
            mpq_set_si(rq->q, 1, 1);
            Py_DECREF((PyObject*)tempbq);
            Py_DECREF((PyObject*)tempez);
            return (PyObject*)rq;
        }
        bsign = mpq_sgn(tempbq->q);
        if (esign < 0) {
            if (bsign == 0) {
                ZERO_ERROR("mpq.pow() 0 base to negative exponent");
                Py_DECREF((PyObject*)rq);
                Py_DECREF((PyObject*)tempbq);
                Py_DECREF((PyObject*)tempez);
                return NULL;
            }
            if (bsign < 0) {
                mpz_neg(mpq_numref(rq->q), mpq_denref(tempbq->q));
            }
            else {
                mpz_set(mpq_numref(rq->q), mpq_denref(tempbq->q));
            }
            mpz_abs(mpq_denref(rq->q), mpq_numref(tempbq->q));
            tempexp = -mpz_get_si(tempez->z);
        }
        else {
            mpq_set(rq->q, tempbq->q);
            tempexp = mpz_get_si(tempez->z);
        }
        if (tempexp>1) {
            mpz_pow_ui(mpq_numref(rq->q), mpq_numref(rq->q), tempexp);
            mpz_pow_ui(mpq_denref(rq->q), mpq_denref(rq->q), tempexp);
        }
        Py_DECREF((PyObject*)tempbq);
        Py_DECREF((PyObject*)tempez);
        return (PyObject*)rq;
    }
    else {
#ifdef WITHMPFR
        tempbf = Pympfr_From_Real(base, 0);
        tempef = Pympfr_From_Real(exp, 0);
        rf = Pympfr_new(0);
        if (!tempbf || !tempef || !rf) {
            TYPE_ERROR("mpq.pow() unsupported operands");
            Py_XDECREF((PyObject*)tempbf);
            Py_XDECREF((PyObject*)tempef);
            Py_XDECREF((PyObject*)rf);
            return NULL;
        }
        rf->rc = mpfr_pow(rf->f, tempbf->f, tempef->f, context->now.mpfr_round);
        Py_DECREF((PyObject*)tempbf);
        Py_DECREF((PyObject*)tempef);
        return (PyObject*)rf;
#else
        TYPE_ERROR("mpq.pow() unsupported operands");
        return NULL;
#endif
    }
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

