/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen                  *
 *                                                                         *
 * This file is part of GMPY2.                                             *
 *                                                                         *
 * GMPY2 is free software: you can redistribute it and/or modify it under  *
 * the terms of the GNU Lesser General Public License as published by the  *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2 is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    *
 * License for more details.                                               *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with GMPY2; if not, see <http://www.gnu.org/licenses/>    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

PyDoc_STRVAR(doc_mpq,
"mpq(n) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n. Decimal and\n"
"     Fraction values are converted exactly.\n\n"
"mpq(n,m) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n/m.\n\n"
"mpq(s[, base=10]) -> mpq\n\n"
"     Return an 'mpq' object from a string s made up of digits in\n"
"     the given base. s may be made up of two numbers in the same\n"
"     base separated by a '/' character.\n");

static PyObject *
Pygmpy_mpq(PyObject *self, PyObject *args, PyObject *keywds)
{
    PympqObject *result = 0, *temp;
    PyObject *n = 0, *m = 0;
    long base = 10;
    Py_ssize_t argc;
    static char *kwlist[] = {"s", "base", NULL };

    TRACE("Pygmpy_mpq() called...\n");

    argc = PyTuple_Size(args);
    if (argc < 1 || argc > 2) {
        TYPE_ERROR("mpq() requires 1 or 2 arguments");
        return NULL;
    }

    n = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(n)) {
        /* keyword base is legal */
        if (PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist, &n, &base)) {
            if ((base!=0) && ((base<2)||(base>62))) {
                VALUE_ERROR("base for mpq() must be 0 or in the "
                            "interval 2 ... 62");
            }
            else {
                result = PyStr2Pympq(n, base);
            }
        }
        return (PyObject*)result;
    }

    if (isDecimal(n)) {
        result = Pympq_From_Decimal(n);
        if (!mpz_cmp_si(mpq_denref(result->q), 0)) {
            if (mpz_get_si(mpq_numref(result->q)) == 0)
                VALUE_ERROR("'mpq' does not support NaN");
            else
                VALUE_ERROR("'mpq' does not support Infinity");
            Py_DECREF((PyObject*)result);
            result = NULL;
        }
        return (PyObject*)result;
    }

    if (argc == 2)
        m = PyTuple_GetItem(args, 1);

#ifdef WITHMPFR
    if (!isReal(n) || (m && !isReal(m))) {
#else
    if (!(isRational(n) || PyFloat_Check(n)) ||
        (m && !(isRational(m) || PyFloat_Check(m)))) {
#endif
        TYPE_ERROR("mpq() requires numeric or string argument");
        return NULL;
    }

    /* should now have one or two numeric values */
    result = Pympq_From_Real(n);
    if (!result && !PyErr_Occurred()) {
        TYPE_ERROR("mpq() requires numeric or string argument");
        return NULL;
    }
    if (m) {
        temp = Pympq_From_Real(m);
        if (!temp && !PyErr_Occurred()) {
            TYPE_ERROR("mpq() requires numeric or string argument");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpq_sgn(temp->q) == 0) {
            ZERO_ERROR("zero denominator in 'mpq'");
            Py_DECREF((PyObject*)result);
            Py_DECREF((PyObject*)temp);
            return NULL;
        }
        mpq_div(result->q, result->q, temp->q);
        Py_DECREF((PyObject*)temp);
    }
    return (PyObject*)result;
}

/* Functions that operate strictly on mpq. */

/* produce digits for an mpq in requested base, default 10 */
PyDoc_STRVAR(doc_qdigitsm,
"x.digits([base=10]) -> string\n\n"
"Return a Python string representing x in the given base (2 to 62,\n"
"default is 10). A leading '-' is present if x<0, but no leading '+'\n"
"is present if x>=0.\n");

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

PyDoc_STRVAR(doc_numerg,
"numer(x) -> mpz\n\n"
"Return the numerator of x.");

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

PyDoc_STRVAR(doc_denomg,
"denom(x) -> mpz\n\n"
"Return the denominator of x.");

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

PyDoc_STRVAR(doc_qdivg,
"qdiv(x[, y=1]) -> number\n\n"
"Return x/y as 'mpz' if possible, or as 'mpq' if x is not exactly\n"
"divisible by y.");

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
            TYPE_ERROR("first argument cannot be converted to 'mpq'");
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
                TYPE_ERROR("second argument cannot be converted to 'mpq'");
            return NULL;
        }
        if (mpq_sgn(Pympq_AS_MPQ(other))==0) {
            PyObject* result = 0;
            ZERO_ERROR("division or modulo by zero in qdiv");
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
Pympq_square(PyObject *self, PyObject *other)
{
    PympqObject *tempx, *result;

    if (!(result = Pympq_new()))
        return NULL;

    if (self && (Pympq_Check(self))) {
        mpq_mul(result->q, Pympq_AS_MPQ(self), Pympq_AS_MPQ(self));
    }
    else if (Pympq_Check(other)) {
        mpq_mul(result->q, Pympq_AS_MPQ(other), Pympq_AS_MPQ(other));
    }
    else {
        if (!(tempx = Pympq_From_Rational(other))) {
            TYPE_ERROR("square() requires 'mpq' argument");
            return NULL;
        }
        else {
            mpq_mul(result->q, Pympq_AS_MPQ(tempx), Pympq_AS_MPQ(tempx));
            Py_DECREF((PyObject*)tempx);
        }
    }
    return (PyObject*)result;
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
    mpir_si tempexp;

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
    PyObject *temp;

    if (self->hash_cache != -1)
        return self->hash_cache;

    if (!(temp = Pympq2PyFloat(self))) {
        SYSTEM_ERROR("Could not convert 'mpq' to float.");
        return -1;
    }
    self->hash_cache = PyObject_Hash(temp);
    Py_DECREF(temp);
    return self->hash_cache;
#endif
}

static PyObject *
Pympq_add(PyObject *self, PyObject *args)
{
    PympqObject *result;
    PyObject *other;

    PARSE_TWO_MPQ(other, "add() requires 'mpq','mpq' arguments");

    if ((result = Pympq_new()))
        mpq_add(result->q, Pympq_AS_MPQ(self), Pympq_AS_MPQ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympq_sub(PyObject *self, PyObject *args)
{
    PympqObject *result;
    PyObject *other;

    PARSE_TWO_MPQ(other, "sub() requires 'mpq','mpq' arguments");

    if ((result = Pympq_new()))
        mpq_sub(result->q, Pympq_AS_MPQ(self), Pympq_AS_MPQ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympq_mul(PyObject *self, PyObject *args)
{
    PympqObject *result;
    PyObject *other;

    PARSE_TWO_MPQ(other, "mul() requires 'mpq','mpq' arguments");

    if ((result = Pympq_new()))
        mpq_mul(result->q, Pympq_AS_MPQ(self), Pympq_AS_MPQ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympq_div(PyObject *self, PyObject *args)
{
    PympqObject *result;
    PyObject *other;

    PARSE_TWO_MPQ(other, "div() requires 'mpq','mpq' arguments");

    if ((result = Pympq_new())) {
        if (mpq_sgn(Pympq_AS_MPQ(other)) == 0) {
            ZERO_ERROR("'mpq' division by zero");
            Py_DECREF((PyObject*)result);
            result = 0;
        }
        else {
            mpq_div(result->q, Pympq_AS_MPQ(self), Pympq_AS_MPQ(other));
        }
    }

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

#ifdef PY3
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympq_neg,               /* nb_negative             */
    (unaryfunc) Pympq_pos,               /* nb_positive             */
    (unaryfunc) Pympq_abs,               /* nb_absolute             */
    (inquiry) Pympq_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympq2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympq2PyFloat,           /* nb_float                */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympq_neg,               /* nb_negative             */
    (unaryfunc) Pympq_pos,               /* nb_positive             */
    (unaryfunc) Pympq_abs,               /* nb_absolute             */
    (inquiry) Pympq_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympq2PyInt,             /* nb_int                  */
    (unaryfunc) Pympq2PyLong,            /* nb_long                 */
    (unaryfunc) Pympq2PyFloat,           /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add;         */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympq_getseters[] =
{
    { "numerator", (getter)Pympq_getnumer, NULL, "numerator", NULL },
    { "denominator", (getter)Pympq_getdenom, NULL, "denominator", NULL },
    {NULL}
};

static PyMethodDef Pympq_methods [] =
{
    { "digits", Pympq_digits, METH_VARARGS, doc_qdigitsm },
    { NULL, NULL, 1 }
};

static PyTypeObject Pympq_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpq",                                  /* tp_name          */
    sizeof(PympqObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympq_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympq2repr,                  /* tp_repr          */
    &mpq_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympq_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympq2str,                   /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE |
        Py_TPFLAGS_CHECKTYPES,              /* tp_flags         */
#endif
    "Multiple precision rational",          /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympq_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympq_getseters,                        /* tp_getset        */
};

