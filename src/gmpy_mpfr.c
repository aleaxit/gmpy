/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpfr.c                                                             *
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

static int isReal(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (PyFloat_Check(obj))     return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

/* Functions that operate strictly on mpfr. */

/* Make a copy of an mpfr object. If bits is 0, the new object will have
 * the same precision as the original object. If the requested precision
 * is less than the precision of the original object, the new object
 * will be rounded to requested precision using the current rounding mode.
 */

static PympfrObject *
Pympfr2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if (bits == 0)
        bits = mpfr_get_prec(Pympfr_AS_MPFR(self));
    if ((result = Pympfr_new(bits)))
        result->rc = mpfr_set(result->f, Pympfr_AS_MPFR(self),
                              context->now.mpfr_round);
    return result;
}

static PympfrObject *
PyFloat2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if (!bits)
        bits = context->now.mpfr_prec;
    if ((result = Pympfr_new(bits)))
        result->rc = mpfr_set_d(result->f, PyFloat_AS_DOUBLE(self),
                                context->now.mpfr_round);
    return result;
}

static PympfrObject *
Pympz2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = Pympfr_new(bits)))
        result->rc = mpfr_set_z(result->f, Pympz_AS_MPZ(self),
                                context->now.mpfr_round);
    return result;
}

#define Pyxmpz2Pympfr Pympz2Pympfr

static PympzObject *
Pympfr2Pympz(PyObject *self)
{
    PympzObject *result;

    if ((result = Pympz_new())) {
        if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    }
    return result;
}

static PyxmpzObject *
Pympfr2Pyxmpz(PyObject *self)
{
    PyxmpzObject *result;

    if ((result = Pyxmpz_new())) {
        if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    }
    return result;
}

/* Return the simpliest rational number that approximates 'self' to the
 * requested precision 'err'. If 'err' is negative, then the requested
 * precision is -2**abs(int(err)). If 'err' is NULL, then the requested
 * precision is -2**prec. If 'prec' is 0, then the requested precision is
 * the precision of 'self'.
 */

static PympqObject *
stern_brocot(PympfrObject* self, PympfrObject *err, mpfr_prec_t prec, int mayz)
{
    PympqObject *result = 0;
    int i, negative, errsign;
    mpfr_t f, al, a, r1[3], r2[3], minerr, curerr, newerr, temp;
#define F2Q_PREC 20

    if (mpfr_nan_p(self->f)) {
        VALUE_ERROR("Cannot convert NaN to a number.");
        return NULL;
    }

    if (mpfr_inf_p(self->f)) {
        OVERFLOW_ERROR("Cannot convert Infinity to a number.");
        return NULL;
    }

    if (prec == 0)
        prec = mpfr_get_prec(self->f);

    errsign = err ? mpfr_sgn(err->f) : 0;
    if (errsign < 0)
        prec = (mpfr_prec_t)(-mpfr_get_si(err->f, context->now.mpfr_round));

    if (errsign <= 0 && (prec < 2 || prec > mpfr_get_prec(self->f))) {
        VALUE_ERROR("Requested precision out-of-bounds.");
        return NULL;
    }

    if (!(result = Pympq_new()))
        return NULL;

    mpfr_init2(minerr, F2Q_PREC);
    if (errsign <= 0) {
        mpfr_set_ui(minerr, 1, context->now.mpfr_round);
        mpfr_div_2si(minerr, minerr, prec, context->now.mpfr_round);
    }
    else {
        mpfr_set(minerr, err->f, context->now.mpfr_round);
    }

    mpfr_init2(f, prec);
    if (mpfr_sgn(self->f) < 0) {
        negative = 1;
        mpfr_abs(f, self->f, context->now.mpfr_round);
    }
    else {
        negative = 0;
        mpfr_set(f, self->f, context->now.mpfr_round);
    }

    mpfr_init2(al, prec);
    mpfr_set(al, f, context->now.mpfr_round);
    mpfr_init2(a, prec);
    mpfr_floor(a, al);
    mpfr_init2(temp, prec);
    for (i=0; i<3; ++i) {
        mpfr_init2(r1[i], prec);
        mpfr_init2(r2[i], prec);
    }
    mpfr_set_si(r1[0], 0, context->now.mpfr_round);
    mpfr_set_si(r1[1], 0, context->now.mpfr_round);
    mpfr_set_si(r1[2], 1, context->now.mpfr_round);
    mpfr_set_si(r2[0], 0, context->now.mpfr_round);
    mpfr_set_si(r2[1], 1, context->now.mpfr_round);
    mpfr_set(r2[2], a, context->now.mpfr_round);
    mpfr_init2(curerr, F2Q_PREC);
    mpfr_init2(newerr, F2Q_PREC);
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

    if (mayz && (mpfr_cmp_ui(r1[2],1) == 0)) {
        Py_DECREF((PyObject*)result);
        result = (PympqObject*)Pympz_new();
        mpfr_get_z(Pympz_AS_MPZ(result), r2[2], context->now.mpfr_round);
        if (negative)
            mpz_neg(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result));
    }
    else {
        mpfr_get_z(mpq_numref(result->q), r2[2], context->now.mpfr_round);
        mpfr_get_z(mpq_denref(result->q), r1[2], context->now.mpfr_round);
        if (negative)
            mpz_neg(mpq_numref(result->q), mpq_numref(result->q));
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
    return result;
}

static PympqObject *
Pympfr2Pympq(PyObject *self)
{
    return stern_brocot((PympfrObject*)self, 0, 0, 0);
}

static PympfrObject *
Pympq2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = Pympfr_new(bits)))
        result->rc = mpfr_set_q(result->f, Pympq_AS_MPQ(self),
                                context->now.mpfr_round);
    return result;
}

static PympfrObject *
PyLong2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;
    PyObject *temp = (PyObject*)PyLong2Pympz(self);

    if (!temp)
        return NULL;
    result = Pympz2Pympfr(temp, bits);
    Py_DECREF(temp);
    return result;
}

#ifdef PY2
static PympfrObject *
PyInt2Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = Pympfr_new(bits)))
        result->rc = mpfr_set_si(result->f, PyInt_AsLong(self),
                                 context->now.mpfr_round);
    return result;
}

static PyObject *
Pympfr2PyInt(PympfrObject *self)
{
    PyObject *result;
    PympzObject *temp = Pympfr2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz_To_Integer(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}
#endif

static PympfrObject *
PyStr2Pympfr(PyObject *s, long base, mpfr_prec_t bits)
{
    PympfrObject *result;
    unsigned char *cp;
    mpfr_prec_t prec;
    Py_ssize_t i, len;
    PyObject *ascii_str = NULL;

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = (unsigned char*)PyBytes_AsString(s);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }

    if (bits > 0)
        prec = bits;
    else
        prec = context->now.mpfr_prec;

    if (!(result = Pympfr_new(prec))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Don't allow NULL characters */
    for (i=0; i<len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            Py_DECREF((PyObject*)result);
            Py_XDECREF(ascii_str);
            return NULL;
        }
    }
    /* delegate the rest to MPFR */
    if (-1 == mpfr_set_str(result->f, (char*)cp, base,
                           context->now.mpfr_round)) {
        VALUE_ERROR("invalid digits");
        Py_DECREF((PyObject*)result);
        Py_XDECREF(ascii_str);
        return NULL;
    }
    Py_XDECREF(ascii_str);
    return result;
}

static PyObject *
Pympfr2PyLong(PympfrObject *self)
{
    PyObject *result;
    PympzObject *temp = Pympfr2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz2PyLong(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}

static PyObject *
Pympfr2PyFloat(PympfrObject *self)
{
    double res = mpfr_get_d(self->f, context->now.mpfr_round);

    return PyFloat_FromDouble(res);
}

/*
 * helper functions for mpfr->binary conversion
 * hof: maps a hex-digit character into 0..15
 * di256: maps two hex-digits chars into 0..255
 */

static int
hof(int hedi)
{
    static char table[] = "0123456789abcdef";
    char* p = strchr(table, tolower(hedi));

    assert(hedi && p);
    return (int)(p-table);
}

static char
di256(int di1, int di2)
{
    return (char)(hof(di2)+16*hof(di1));
}

/*
 * Build old-version binary string. This will be removed when the new
 * binary format is implemented.
 */

static PyObject *
Pympfr2binary(PympfrObject *self)
{
    size_t size, hexdigs, i, j;
    char *buffer, *aux;
    PyObject *result;
    int sign, codebyte;
    mpfr_exp_t the_exp;
    long lexp, lprec;
    int lexpodd, extrabyte;

    /* prepare codebyte */
    sign = mpfr_sgn(self->f);
    if (sign == 0) {
        /* 0 -> single codebyte with 'zerovalue' bit set */
#ifdef PY3
        return Py_BuildValue("y", "\004");
#else
        return Py_BuildValue("s", "\004");
#endif
        /* codebyte = 0; */
    }
    else if (sign < 0) {
        codebyte = 1;
        mpfr_neg(self->f, self->f, context->now.mpfr_round);
    }
    else {
        codebyte = 0;
    }

    /* get buffer of base-16 digits */
    buffer = mpfr_get_str(0, &the_exp, 16, 0, self->f,
                          context->now.mpfr_round);

    /* strip trailing zeros */
    hexdigs = strlen(buffer) - 1;
    while ((hexdigs >= 1) && (buffer[hexdigs] == '0'))
        buffer[hexdigs--] = 0x00;

    /* no need to worry about null-buffer as x->f==0.0 was
     * already handled above (see first test on 'sign').
     */
    /* restore correct sign to x->f if it was changed! */
    if (codebyte) {
        mpfr_neg(self->f, self->f, context->now.mpfr_round);
    }
    hexdigs = strlen(buffer);
    /* adjust exponent, & possibly set codebyte's expo-sign bit.
     * note the_exp is base-16 exp, while we need to have it in
     * base-256 -- so it's basically halved (but, with care...!).
     */
    if (the_exp<0) {
        codebyte |= 2;
        the_exp = -the_exp;
    }
    lexp = the_exp;
    lexpodd = lexp & 1;
    lexp = lexp/2 + lexpodd;
    if (lexpodd && (codebyte&2))
        --lexp;
    /* remember we also store precision explicitly */
    codebyte |= 8;

    /* allocate suitably-sized, uninitialized Python string */
    size = (hexdigs + 1) / 2;
    /* allocate an extra byte if lexpodd and hexdigs is even */
    extrabyte = lexpodd & ~hexdigs;
    result = PyBytes_FromStringAndSize(0, 1+4+size+4+extrabyte);
    if (!result)
        return NULL;
    /* set the data to the new Python string's buffer */
    aux = PyBytes_AS_STRING(result);
    /* codebyte first */
    aux[0] = (char)codebyte;
    /* then precision */
    lprec = mpfr_get_prec(self->f);
    for (i=0; i<4; ++i) {
        aux[i+1] = (char)(lprec & 0xFF);
        lprec >>= 8;
    }
    /* then exponent */
    for (i=0; i<4; ++i) {
        aux[4+i+1] = (char)(lexp & 0xFF);
        lexp >>= 8;
    }
    /* then mantissa, grouping 2 hex digits per base-256 digits;
     * with some care for the first & last ones...
     */
    j=0; i=0;
    if (lexpodd) {
        aux[i+9] = di256('0',buffer[0]);
        j=1; i=1;
    }
    for (; i<size+extrabyte; ++i) {
        int secdig = (j+1)<hexdigs? buffer[j+1] : '0';
        aux[i+9] = di256(buffer[j], secdig);
        j += 2;
    }

    mpfr_free_str(buffer);
    return result;
}

static PyObject*
Pympfr_ascii(PympfrObject *self, int base, int digits)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;

    /* check arguments are valid */
    assert(Pympfr_Check((PyObject*)self));
    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    /* Process special cases first */
    if (!(mpfr_regular_p(self->f))) {
        if (mpfr_nan_p(self->f)) {
            return Py_BuildValue("(sii)", "nan", 0, 0);
        }
        else if (mpfr_inf_p(self->f) && !mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "inf", 0, 0);
        }
        else if (mpfr_inf_p(self->f) && mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "-inf", 0, 0);
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "-0", 0, mpfr_get_prec(self->f));
        }
        else {
            return Py_BuildValue("(sii)", "0", 0, mpfr_get_prec(self->f));
        }
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, context->now.mpfr_round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in Pympfr_ascii");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
    mpfr_free_str(buffer);
    return result;
}

/*
 * If obj is a Pympfr and bits is 0 or bits is the same as the precision of
 * obj, then a new reference is created.
 *
 * For all other numerical types with bits = 0, the conversion is rounded to
 * context->now.mpfr_prec.
 */

static PympfrObject *
Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits)
{
    PympfrObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympfr_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpfr is still
         * valid in the current context. */
        if (!bits || mpfr_get_prec(Pympfr_AS_MPFR(obj)) == bits) {
            newob = (PympfrObject*) obj;
            Py_INCREF(obj);
        }
        else {
            newob = Pympfr2Pympfr((PyObject*)obj, bits);
        }
    }
    else if (Pympfr_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer valid
         * and mpfr_check_range needs to be called. */
        if (context->now.trap_expbound) {
            GMPY_EXPBOUND("exponent of existing 'mpfr' incompatible with current context");
            return NULL;
        }
        if ((newob = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(obj))))) {
            mpfr_set(newob->f, Pympfr_AS_MPFR(obj), context->now.mpfr_round);
            newob->round_mode = ((PympfrObject*)obj)->round_mode;
            newob->rc = ((PympfrObject*)obj)->rc;
            newob->rc = mpfr_check_range(newob->f, newob->rc, newob->round_mode);
        }
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympfr(obj, bits);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympfr(obj, bits);
#endif
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pympfr(obj, bits);
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympfr(obj, bits);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympfr(obj, bits);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympfr(obj, bits);
    }
    else if (isDecimal(obj)) {
        PympqObject *temp = Pympq_From_Decimal(obj);
        if (temp) {
            if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
                if (mpz_get_si(mpq_numref(temp->q)) == 0) {
                    if ((newob = Pympfr_new(0)))
                        mpfr_set_nan(newob->f);
                }
                else if (mpz_get_si(mpq_numref(temp->q)) < 0) {
                    if ((newob = Pympfr_new(0)))
                        mpfr_set_inf(newob->f, -1);
                }
                else {
                    if ((newob = Pympfr_new(0)))
                        mpfr_set_inf(newob->f, 1);
                }
            }
            else
                newob = Pympq2Pympfr((PyObject*)temp, bits);
            Py_DECREF((PyObject*)temp);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pympfr((PyObject*)temp, bits);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }
#ifdef DEBUG
    if (global.debug)
        fprintf(stderr, "Pympfr_From_Real(%p,%ld)->%p (%ld)\n", (void *)obj,
                (long)bits, (void *)newob, 
                newob != 0 ? (long)mpfr_get_prec(newob->f) : -1);
#endif
    if (!newob)
        TYPE_ERROR("object could not be converted to 'mpfr'");
    return newob;
}

/*
 * coerce any number to a mpf
 */

int
Pympfr_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympfrObject* newob = Pympfr_From_Real(arg, 0);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpfr'");
        return 0;
    }
}

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
    PyObject *result;

    if (!PyArg_ParseTuple(args, "O&|O&", Pympfr_convert_arg, &self,
                          Pympfr_convert_arg, &err)) {
        TYPE_ERROR("f2q() requires 'mpfr', ['mpfr'] arguments");
        return NULL;
    }

    result = (PyObject*)stern_brocot((PympfrObject*)self, err, 0, 1);
    Py_DECREF(self);
    Py_XDECREF((PyObject*)err);
    return result;
}

/* str and repr implementations for mpfr */
static PyObject *
Pympfr2str(PympfrObject *self)
{
    PyObject *result, *temp;
    long precision;
    char fmtstr[30];

    precision = (long)(log10(2) * (double)mpfr_get_prec(Pympfr_AS_MPFR(self))) + 2;

    sprintf(fmtstr, "{0:.%ldg}", precision);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympfr2repr(PympfrObject *self)
{
    PyObject *result, *temp;
    long precision, bits;
    char fmtstr[30];

    bits = mpfr_get_prec(Pympfr_AS_MPFR(self));
    precision = (long)(log10(2) * (double)bits) + 2;

    if (mpfr_number_p(Pympfr_AS_MPFR(self)) && bits != DBL_MANT_DIG)
        sprintf(fmtstr, "mpfr('{0:.%ldg}',%ld)", precision, bits);
    else
        sprintf(fmtstr, "mpfr('{0:.%ldg}')", precision);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static char doc_mpfr[] = "\
mpfr(n,[precison=0]):\n\
      Return an 'mpfr' object after converting a numeric value. If\n\
      no precision, or a precision of 0, is specified; the precison\n\
      is taken from the current context.\n\
mpfr(s,[precision=0],[base=0]):\n\
      Return 'mpfr' object after converting a string 's' made up of\n\
      digits in the given base, possibly with fraction-part (with\n\
      period as a separator) and/or exponent-part (with exponent\n\
      marker 'e' for base<=10, else '@'). If no precision, or a\n\
      precision of 0, is specified; the precison is taken from the\n\
      current context. The base of the string representation must\n\
      be 0 or in the interval 2 ... 62. If the base is 0, the leading\n\
      digits of the string are used to identify the base: 0b implies\n\
      base=2, 0x implies base=16, otherwise base=10 is assumed.\n\
";
static PyObject *
Pygmpy_mpfr(PyObject *self, PyObject *args, PyObject *keywds)
{
    PympfrObject *result = 0;
    PyObject *arg0;
    long base = 0;
    Py_ssize_t argc;
    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t bits = 0;
    static char *kwlist_s[] = {"s", "precision", "base", NULL};
    static char *kwlist_n[] = {"n", "precision", NULL};

    TRACE("Pygmpy_mpfr() called...\n");

    argc = PyTuple_Size(args);
    if ((argc < 1) || (argc > 3)) {
        TYPE_ERROR("mpfr() requires 1 to 3 arguments");
        return NULL;
    }

    arg0 = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(arg0)) {
        /* Can have both precision and/or base as keyword arguments. */
        if (PyArg_ParseTupleAndKeywords(args, keywds, "O|ll", kwlist_s,
                                        &arg0, &bits, &base)) {
            if ((base!=0) && ((base<2)||(base>62))) {
                VALUE_ERROR("base for mpfr() must be 0 or in the "
                            "interval 2 ... 62");
            }
            else if (bits < 0) {
                VALUE_ERROR("precision for mpfr() must be >= 0");
            }
            else {
                result = PyStr2Pympfr(arg0, base, bits);
            }
        }
        SUBNORMALIZE(result);
        return (PyObject*)result;
    }

    /* Optimize the common case */
    if (isReal(arg0) && argc == 1 && !keywds) {
        result = Pympfr_From_Real(arg0, bits);
        SUBNORMALIZE(result);
        return (PyObject*)result;
    }

    /* Can only have precision as keyword argument. */
    if (PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist_n, &arg0, &bits)) {
        if (bits < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
        }
        else {
            result = Pympfr_From_Real(arg0, bits);
            if (!result)
                TYPE_ERROR("mpfr() requires numeric or string argument");
        }
    }
    return (PyObject*)result;
}

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

/* Implement the .imag attribute of an mpfr. */

static PyObject *
Pympfr_getimag_attrib(PympfrObject *self, void *closure)
{
    PympfrObject *result;

    if ((result = Pympfr_new(0)))
        mpfr_set_zero(result->f, 1);
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpfr. */

static PyObject *
Pympfr_getreal_attrib(PympfrObject *self, void *closure)
{
    return (PyObject*)Pympfr2Pympfr((PyObject*)self, 0);
}

/* Implement the nb_bool slot. */

static int
Pympfr_nonzero(PympfrObject *self)
{
    return !mpfr_zero_p(self->f);
}

/* Implement the conjugate() method. */

PyDoc_STRVAR(doc_mpfr_conjugate,
"x.conjugate() -> mpfr\n\n"
"Return the conjugate of x (which is just a copy of x since x is\n"
"not a complex number).");

static PyObject *
Pympfr_conjugate(PyObject *self, PyObject *args)
{
    return (PyObject*)Pympfr2Pympfr(self, 0);
}

/* Implement the nb_positive slot. */

/* TODO: can probably just call Pympfr_From_Real. */

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
    CHECK_FLAGS("__pos__");
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

PyDoc_STRVAR(doc_g_mpfr_get_exp,
"get_exp(mpfr) -> integer\n\n"
"Return the exponent of an mpfr. Returns 0 for NaN or Infinity, sets\n"
"the erange flag and will raise an exception if trap_erange is set.");

static PyObject *
Pympfr_get_exp(PyObject *self, PyObject *other)
{
    PyObject *result = 0;
    Py_ssize_t exp;

    PARSE_ONE_MPFR_OTHER("get_exp() requires 'mpfr. argument");

    if (mpfr_regular_p(Pympfr_AS_MPFR(self))) {
        exp = (Py_ssize_t)mpfr_get_exp(Pympfr_AS_MPFR(self));
        result = PyIntOrLong_FromSsize_t((Py_ssize_t)exp);
    }
    else if (mpfr_zero_p(Pympfr_AS_MPFR(self))) {
        Py_DECREF(self);
        result = PyIntOrLong_FromSsize_t(0);
    }
    else {
        context->now.erange = 1;
        if (context->now.trap_erange) {
            GMPY_ERANGE("Can not get exponent from NaN or Infinity.");
        }
        else {
            result = PyIntOrLong_FromSsize_t(0);
        }
    }
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_set_exp,
"set_exp(mpfr, n) -> mpfr\n\n"
"Set the exponent of an mpfr to n. If n is outside the range of\n"
"valid exponents, set_exp() will set the erange flag and either\n"
"return the original value or raise an exception if trap_erange\n"
"is set.");

static PyObject *
Pympfr_set_exp(PyObject *self, PyObject *args)
{
    PympfrObject *result = 0;
    long exp = 0;

    if (!PyArg_ParseTuple(args, "O&l", Pympfr_convert_arg, &self, &exp)) {
        TYPE_ERROR("set_exp() requires 'mpfr', 'integer' arguments");
        return NULL;
    }

    if (!(result = Pympfr2Pympfr(self, 0)))
        return NULL;
    Py_DECREF(self);

    result->rc = mpfr_set_exp(Pympfr_AS_MPFR(result), exp);

    if (result->rc) {
        context->now.erange = 1;
        if (context->now.trap_erange) {
            GMPY_ERANGE("New exponent is out-of-bounds.");
            Py_DECREF(result);
            return NULL;
        }
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_sign,
"set_sign(mpfr, bool) -> mpfr\n\n"
"If 'bool' is True, then return an 'mpfr' with the sign bit set.");

static PyObject *
Pympfr_set_sign(PyObject *self, PyObject *args)
{
    PympfrObject *result = 0;
    PyObject *boolean = 0;
    int s;

    if (!PyArg_ParseTuple(args, "O&O", Pympfr_convert_arg, &self, &boolean)) {
        TYPE_ERROR("set_sign() requires 'mpfr', 'boolean' arguments");
        return NULL;
    }

    if (!(result = Pympfr_new(0)))
        return NULL;

    s = PyObject_IsTrue(boolean);
    if (s == -1) {
        TYPE_ERROR("set_sign() requires 'mpfr', 'boolean' arguments");
        Py_DECREF(self);
        Py_DECREF(boolean);
        Py_DECREF(result);
        return NULL;
    }

    result->rc = mpfr_setsign(Pympfr_AS_MPFR(result), Pympfr_AS_MPFR(self),
                              s, context->now.mpfr_round);

    Py_DECREF(self);
    Py_DECREF(boolean);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_copy_sign,
"copy_sign(mpfr, mpfr) -> mpfr\n\n"
"Return an 'mpfr' composed of the first argument with the sign of the\n"
"second argument.");

static PyObject *
Pympfr_copy_sign(PyObject *self, PyObject *args)
{
    PympfrObject *result = 0;
    PyObject *other = 0;

    if (!PyArg_ParseTuple(args, "O&O&", Pympfr_convert_arg, &self,
                          Pympfr_convert_arg, &other)) {
        TYPE_ERROR("copy_sign() requires 'mpfr', 'mpfr' arguments");
        return NULL;
    }

    if (!(result = Pympfr_new(0)))
        return NULL;

    result->rc = mpfr_copysign(Pympfr_AS_MPFR(result), Pympfr_AS_MPFR(self),
                              Pympfr_AS_MPFR(other), context->now.mpfr_round);

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympfr_div_2exp(PyObject *self, PyObject *args)
{
    PympfrObject *result = 0;
    unsigned long exp = 0;

    if (!PyArg_ParseTuple(args, "O&k", Pympfr_convert_arg, &self, &exp)) {
        TYPE_ERROR("div_2exp() requires 'mpfr', 'integer' arguments");
        return NULL;
    }

    if (!(result = Pympfr_new(0)))
        return NULL;

    mpfr_clear_flags();

    result->rc = mpfr_div_2ui(Pympfr_AS_MPFR(result), Pympfr_AS_MPFR(self),
                              exp, context->now.mpfr_round);

    MPFR_CLEANUP_SELF("div_2exp()");
}

static PyObject *
Pympfr_mul_2exp(PyObject *self, PyObject *args)
{
    PympfrObject *result = 0;
    unsigned long exp = 0;

    if (!PyArg_ParseTuple(args, "O&k", Pympfr_convert_arg, &self, &exp)) {
        TYPE_ERROR("mul_2exp() requires 'mpfr', 'integer' arguments");
        return NULL;
    }

    if (!(result = Pympfr_new(0)))
        return NULL;

    mpfr_clear_flags();

    result->rc = mpfr_mul_2ui(Pympfr_AS_MPFR(result), Pympfr_AS_MPFR(self),
                              exp, context->now.mpfr_round);

    MPFR_CLEANUP_SELF("mul_2exp()");
}

PyDoc_STRVAR(doc_g_mpfr_set_nan,
"nan() -> mpfr\n\n"
"Return an 'mpfr' initialized to NaN (Not-A-Number).");

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
"Return an 'mpfr' initialized to Infinity with the same sign as n. If n is\n"
"not given, +Infinity is returned.");

static PyObject *
Pympfr_set_inf(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    long s = 1;

    if (PyTuple_Size(args) == 1) {
        s = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (s == -1 && PyErr_Occurred()) {
            TYPE_ERROR("inf() requires 'int' argument");
            return NULL;
        }
    }

    if ((result = Pympfr_new(0)))
        mpfr_set_inf(result->f, s<0?-1:1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_set_zero,
"zero(n) -> mpfr\n\n"
"Return an 'mpfr' inialized to 0.0 with the same sign as n. If n is not\n"
"given, +0.0 is returned.");

static PyObject *
Pympfr_set_zero(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    long s = 1;

    if (PyTuple_Size(args) == 1) {
        s = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (s == -1 && PyErr_Occurred()) {
            TYPE_ERROR("zero() requires 'int' argument");
            return NULL;
        }
    }

    if ((result = Pympfr_new(0)))
        mpfr_set_zero(result->f, s<0?-1:1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpfr_is_signed,
"x.is_signed() -> boolean\n\n"
"Return True if the sign bit of x is set.");

PyDoc_STRVAR(doc_g_mpfr_is_signed,
"is_signed(x) -> boolean\n\n"
"Return True if the sign bit of x is set.");

static PyObject *
Pympfr_is_signed(PyObject *self, PyObject *other)
{
    int res;
    if(self && Pympfr_Check(self)) {
        Py_INCREF(self);
    }
    else if(Pympfr_Check(other)) {
        self = other;
        Py_INCREF((PyObject*)self);
    }
    else if (!(self = (PyObject*)Pympfr_From_Real(other, 0))) {
        TYPE_ERROR("is_signed() requires 'mpfr' argument");
        return NULL;
    }
    res = mpfr_signbit(Pympfr_AS_MPFR(self));
    Py_DECREF(self);
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
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
"Return True if x is 'nan' (Not-A-Number); False otherwise.");

MPFR_TEST_OTHER(nan, "is_nan() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_inf,
"x.is_inf() -> boolean\n\n"
"Return True if x is +Infinity or -Infinity; False otherwise.");

MPFR_TEST_OTHER(inf, "is_inf() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_number,
"x.is_number() -> boolean\n\n"
"Return True if x is an actual number (i.e. not NaN or Infinity);\n"
"False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_number,
"is_number(x) -> boolean\n\n"
"Return True if x is an actual number (i.e. not NaN or Infinity);\n"
"False otherwise.");

MPFR_TEST_OTHER(number, "is_number() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_zero,
"x.is_zero() -> boolean\n\n"
"Return True if x is zero, False otherwise.");

MPFR_TEST_OTHER(zero, "is_zero() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_regular,
"x.is_regular() -> boolean\n\n"
"Return True if x is not zero, NaN, or Infinity; False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_regular,
"is_regular(x) -> boolean\n\n"
"Return True if x is not zero, NaN, or Infinity; False otherwise.");

MPFR_TEST_OTHER(regular, "is_regular() requires 'mpfr' argument");

PyDoc_STRVAR(doc_mpfr_is_integer,
"x.is_integer() -> boolean\n\n"
"Return True if x is an integer; False otherwise.");

PyDoc_STRVAR(doc_g_mpfr_is_integer,
"is_integer(x) -> boolean\n\n"
"Return True if x is an integer; False otherwise.");

MPFR_TEST_OTHER(integer, "is_integer() requires 'mpfr' argument");

/* produce string for an mpfr with requested/defaulted parameters */

PyDoc_STRVAR(doc_mpfr_digits,
"x.digits(base=10, prec=0) -> (mantissa, exponent, bits)\n\n"
"Returns up to 'prec' digits in the given base. If 'prec' is 0, as many\n"
"digits that are available are returned. No more digits than available\n"
"given x's precision are returned. 'base' must be between 2 and 62,\n"
"inclusive. The result is a three element tuple containing the mantissa,\n"
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
    result = Pympfr_ascii((PympfrObject*)self, base, prec);
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_mpfr_integer_ratio,
"x.as_integer_ratio() -> (num,den)\n\n"
"Return the exact rational equivalent of an mpfr. Value is a tuple\n"
"for compatibility with Python's float.as_integer_ratio().");

static PyObject *
Pympfr_integer_ratio(PyObject *self, PyObject *args)
{
    PympzObject *num = 0, *den = 0;
    mpfr_exp_t temp, twocount;

    if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
        VALUE_ERROR("Cannot pass NaN to mpfr.as_integer_ratio.");
        return NULL;
    }

    if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
        OVERFLOW_ERROR("Cannot pass Infinity to mpfr.as_integer_ratio.");
        return NULL;
    }

    num = Pympz_new();
    den = Pympz_new();
    if (!num || !den) {
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        return NULL;
    }

    if (mpfr_zero_p(Pympfr_AS_MPFR(self))) {
        mpz_set_ui(num->z, 0);
        mpz_set_ui(den->z, 1);
    }
    else {
        temp = mpfr_get_z_2exp(num->z, Pympfr_AS_MPFR(self));
        twocount = (mpfr_exp_t)mpz_scan1(num->z, 0);
        if (twocount) {
            temp += twocount;
            mpz_div_2exp(num->z, num->z, twocount);
        }
        mpz_set_ui(den->z, 1);
        if (temp > 0)
            mpz_mul_2exp(num->z, num->z, temp);
        else if (temp < 0)
            mpz_mul_2exp(den->z, den->z, -temp);
    }
    return Py_BuildValue("(NN)", (PyObject*)num, (PyObject*)den);
}

PyDoc_STRVAR(doc_mpfr_mantissa_exp,
"x.as_mantissa_exp() -> (mantissa,exponent)\n\n"
"Return the mantissa and exponent of an mpfr.");

static PyObject *
Pympfr_mantissa_exp(PyObject *self, PyObject *args)
{
    PympzObject *mantissa = 0, *exponent = 0;
    mpfr_exp_t temp;

    if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
        VALUE_ERROR("Cannot pass NaN to mpfr.as_mantissa_exp.");
        return NULL;
    }

    if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
        OVERFLOW_ERROR("Cannot pass Infinity to mpfr.as_mantissa_exp.");
        return NULL;
    }

    mantissa = Pympz_new();
    exponent = Pympz_new();
    if (!mantissa || !exponent) {
        Py_XDECREF((PyObject*)mantissa);
        Py_XDECREF((PyObject*)exponent);
        return NULL;
    }

    if (mpfr_zero_p(Pympfr_AS_MPFR(self))) {
        mpz_set_ui(mantissa->z, 0);
        mpz_set_ui(exponent->z, 1);
    }
    else {
        temp = mpfr_get_z_2exp(mantissa->z, Pympfr_AS_MPFR(self));
        mpz_set_si(exponent->z, temp);
    }
    return Py_BuildValue("(NN)", (PyObject*)mantissa, (PyObject*)exponent);
}

PyDoc_STRVAR(doc_mpfr_simple_fraction,
"x.as_simple_fraction(precision=0) -> mpq\n\n"
"Return a simple rational approximation to x. The result will be\n"
"accurate to 'precision' bits. If 'precision' is 0, the precision\n"
"of 'x' will be used.");

static PyObject *
Pympfr_simple_fraction(PyObject *self, PyObject *args, PyObject *keywds)
{
    mpfr_prec_t prec = 0;
    static char *kwlist[] = {"precision", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "|l", kwlist, &prec))
        return NULL;

    return (PyObject*)stern_brocot((PympfrObject*)self, 0, prec, 0);
}

static Py_hash_t
_mpfr_hash(mpfr_t f)
{
#ifdef _PyHASH_MODULUS
    Py_uhash_t hash = 0;
    Py_ssize_t exp;
    size_t msize;
    int sign;

    /* Handle special cases first */
    if (!mpfr_number_p(f)) {
        if (mpfr_inf_p(f))
            if (mpfr_sgn(f) > 0)
                return _PyHASH_INF;
            else
                return -_PyHASH_INF;
        else
            return _PyHASH_NAN;
    }

    /* Calculate the number of limbs in the mantissa. */
    msize = (f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;

    /* Calculate the hash of the mantissa. */
    if (mpfr_sgn(f) > 0) {
        hash = mpn_mod_1(f->_mpfr_d, msize, _PyHASH_MODULUS);
        sign = 1;
    }
    else if (mpfr_sgn(f) < 0) {
        hash = mpn_mod_1(f->_mpfr_d, msize, _PyHASH_MODULUS);
        sign = -1;
    }
    else {
        return 0;
    }

    /* Calculate the final hash. */
    exp = f->_mpfr_exp - (msize * mp_bits_per_limb);
    exp = exp >= 0 ? exp % _PyHASH_BITS : _PyHASH_BITS-1-((-1-exp) % _PyHASH_BITS);
    hash = ((hash << exp) & _PyHASH_MODULUS) | hash >> (_PyHASH_BITS - exp);

    hash *= sign;
    if (hash == (Py_uhash_t)-1)
        hash = (Py_uhash_t)-2;
    return (Py_hash_t)hash;
#else
    double temp;
    temp = mpfr_get_d(f, context->now.mpfr_round);
    return _Py_HashDouble(temp);
#endif
}

static Py_hash_t
Pympfr_hash(PympfrObject *self)
{
    if (self->hash_cache == -1)
        self->hash_cache = _mpfr_hash(self->f);
    return self->hash_cache;
}

/* This function is used in gmpy_mpany. */

static PyObject *
Pympfr_pow(PyObject *base, PyObject *exp, PyObject *m)
{
    PympfrObject *tempb, *tempe, *result;
#ifdef WITHMPC
    PympcObject *mpc_result;
#endif

    if (m != Py_None) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
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
    result->rc = mpfr_pow(result->f, tempb->f,
                          tempe->f, context->now.mpfr_round);
#ifdef WITHMPC
    if (result && mpfr_nanflag_p() && context->now.allow_complex) {
        /* If we don't get a valid result, or the result is a nan, then just
         * return the original mpfr value. */
        if (!(mpc_result = (PympcObject*)Pympc_pow(base, exp, m)) ||
            MPC_IS_NAN_P(mpc_result)) {

            Py_XDECREF((PyObject*)mpc_result);
            context->now.invalid = 1;
            GMPY_INVALID("invalid operation in 'mpfr' pow()");
            goto done;
        }
        /* return a valid complex result */
        Py_DECREF(result);
        result = (PympfrObject*)mpc_result;
        goto done;
    }
#endif

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

static PyObject *
Pympfr_sqrt(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("sqrt() requires 'mpfr' argument");

#ifdef WITHMPC
    if (mpfr_sgn(Pympfr_AS_MPFR(self)) < 0 && context->now.allow_complex) {
        Py_DECREF(self);
        return Pympc_sqrt(self, other);
    }
#endif

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpfr_sqrt(result->f, Pympfr_AS_MPFR(self),
                           context->now.mpfr_round);

    MPFR_CLEANUP_SELF("sqrt()");
}

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

    mpfr_clear_flags();
    result->rc = mpfr_rec_sqrt(result->f, Pympfr_AS_MPFR(self),
                               context->now.mpfr_round);

    MPFR_CLEANUP_SELF("rec_sqrt()");
}

PyDoc_STRVAR(doc_mpfr_root,
"root(x,n) -> number\n\n"
"Return n-th root of x. The result always an 'mpfr'.");

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

    PARSE_TWO_MPFR_ARGS(other, "reldiff() requires 'mpfr','mpfr' arguments");

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
    if (!(r = Pympfr_new(0))) \
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
x.__ceil__(): returns an 'mpfr' that is the smallest integer >= x\n\
";
static char doc_g_mpfr_ceil[]="\
ceil(x): returns an 'mpfr' that is the smallest integer >= x\n\
x must be an 'mpfr', or else gets coerced to one.\n\
";

MPFR_UNIOP_NOROUND(ceil)

static char doc_mpfr_floor[]="\
x.__floor__(): returns an 'mpfr' that is the smallest integer <= x\n\
";
static char doc_g_mpfr_floor[]="\
floor(x): returns an 'mpfr' that is the smallest integer <= x\n\
x must be an 'mpfr', or else gets coerced to one.\n\
";

MPFR_UNIOP_NOROUND(floor);

static char doc_mpfr_trunc[]="\
x.__trunc__(): returns an 'mpf'r that is truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0).\n\
";
static char doc_g_mpfr_trunc[]="\
trunc(x): returns an 'mpfr' that is x truncated towards 0\n\
(same as x.floor() if x>=0, x.ceil() if x<0). x must be an\n\
'mpfr', or else gets coerced to one.\n\
";

MPFR_UNIOP_NOROUND(trunc)

static char doc_g_mpfr_round2[]="\
round2(x): returns an 'mpfr' that is x rounded to the nearest\n\
integer, with ties rounded away from 0.";

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

PyDoc_STRVAR(doc_g_mpfr_rint,
"rint(x) -> mpfr\n\n"
"Return x rounded to the nearest integer using the current rounding\n"
"mode.");

MPFR_UNIOP(rint)

PyDoc_STRVAR(doc_g_mpfr_rint_ceil,
"rint_ceil(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next higher or equal integer and then, if needed, using the current\n"
"rounding mode.");

MPFR_UNIOP(rint_ceil)

PyDoc_STRVAR(doc_g_mpfr_rint_floor,
"rint_floor(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"next lower or equal integer and then, if needed, using the current\n"
"rounding mode.");

MPFR_UNIOP(rint_floor)

PyDoc_STRVAR(doc_g_mpfr_rint_round,
"rint_round(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding to the\n"
"nearest integer (ties away from 0) and then, if needed, using\n"
"the current rounding mode.");

MPFR_UNIOP(rint_round)

PyDoc_STRVAR(doc_g_mpfr_rint_trunc,
"rint_trunc(x) -> mpfr\n\n"
"Return x rounded to the nearest integer by first rounding towards\n"
"zero and then, if needed, using the current rounding mode.");

MPFR_UNIOP(rint_trunc)

PyDoc_STRVAR(doc_g_mpfr_frac,
"frac(x) -> mpfr\n\n"
"Return fractional part of x.");

MPFR_UNIOP(frac)

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

/* Needed for square() in mpz_mpany.c */

MPFR_UNIOP(sqr)

PyDoc_STRVAR(doc_g_mpfr_cbrt,
"cbrt(x) -> mpfr\n\n"
"Return the cube root of x.");

MPFR_UNIOP(cbrt)

/* Called via gmpy_mpany so doc-string is there. */

MPFR_UNIOP(log)

static char doc_g_mpfr_log2[]="\
log2(x): returns base-2 logarithm of x.\n\
";

MPFR_UNIOP(log2)

static char doc_g_mpfr_log10[]="\
log10(x): returns base-10 logarithm of x.\n\
";

MPFR_UNIOP(log10)

/* Called via gmpy_mpany so doc-string is there. */

MPFR_UNIOP(exp)

static char doc_g_mpfr_exp2[]="\
exp2(x): returns 2**x.\n\
";

MPFR_UNIOP(exp2)

static char doc_g_mpfr_exp10[]="\
exp10(x): returns 10**x.\n\
";

MPFR_UNIOP(exp10)

MPFR_UNIOP(sin)

MPFR_UNIOP(cos)

MPFR_UNIOP(tan)

static char doc_g_mpfr_sec[]="\
sec(x): returns secant of x; x in radians.\n\
";

MPFR_UNIOP(sec)

static char doc_g_mpfr_csc[]="\
csc(x): returns cosecant of x; x in radians.\n\
";

MPFR_UNIOP(csc)

static char doc_g_mpfr_cot[]="\
cot(x): returns cotangent of x; x in radians.\n\
";

MPFR_UNIOP(cot)

static PyObject *
Pympfr_acos(PyObject* self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("acos() requires 'mpfr' argument");

#ifdef WITHMPC
    if (!mpfr_nan_p(Pympfr_AS_MPFR(self)) &&
            (mpfr_cmp_si(Pympfr_AS_MPFR(self), 1) > 0 ||
            mpfr_cmp_si(Pympfr_AS_MPFR(self), -1) < 0) &&
            context->now.allow_complex) {
        Py_DECREF(self);
        return Pympc_acos(self, other);
    }
#endif

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_clear_flags();
    result->rc = mpfr_acos(result->f, Pympfr_AS_MPFR(self),
                           context->now.mpfr_round);
    MPFR_CLEANUP_SELF("acos()");
}

static PyObject *
Pympfr_asin(PyObject* self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("asin() requires 'mpfr' argument");

#ifdef WITHMPC
    if (!mpfr_nan_p(Pympfr_AS_MPFR(self)) &&
            (mpfr_cmp_si(Pympfr_AS_MPFR(self), 1) > 0 ||
            mpfr_cmp_si(Pympfr_AS_MPFR(self), -1) < 0) &&
            context->now.allow_complex) {
        Py_DECREF(self);
        return Pympc_asin(self, other);
    }
#endif

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_clear_flags();
    result->rc = mpfr_asin(result->f, Pympfr_AS_MPFR(self),
                           context->now.mpfr_round);
    MPFR_CLEANUP_SELF("asin()");
}

MPFR_UNIOP(atan)

MPFR_UNIOP(cosh)

MPFR_UNIOP(sinh)

MPFR_UNIOP(tanh)

static char doc_g_mpfr_sech[]="\
sech(x): returns hyperbolic secant of x.\n\
";

MPFR_UNIOP(sech)

static char doc_g_mpfr_csch[]="\
csch(x): returns hyperbolic cosecant of x.\n\
";

MPFR_UNIOP(csch)

static char doc_g_mpfr_coth[]="\
coth(x): returns hyperbolic cotangent of x.\n\
";

MPFR_UNIOP(coth)

MPFR_UNIOP(acosh)

MPFR_UNIOP(asinh)

static PyObject *
Pympfr_atanh(PyObject* self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPFR_OTHER("atanh() requires 'mpfr' argument");

#ifdef WITHMPC
    if (!mpfr_nan_p(Pympfr_AS_MPFR(self)) &&
            (mpfr_cmp_si(Pympfr_AS_MPFR(self), 1) > 0 ||
            mpfr_cmp_si(Pympfr_AS_MPFR(self), -1) < 0) &&
            context->now.allow_complex) {
        Py_DECREF(self);
        return Pympc_atanh(self, other);
    }
#endif

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        return NULL;
    }
    mpfr_clear_flags();
    result->rc = mpfr_asin(result->f, Pympfr_AS_MPFR(self),
                           context->now.mpfr_round);
    MPFR_CLEANUP_SELF("atanh()");
}


static char doc_g_mpfr_log1p[]="\
log1p(x): returns logarithm of (1+x).\n\
";

MPFR_UNIOP(log1p)

static char doc_g_mpfr_expm1[]="\
expm1(x): returns exponential(x) - 1.\n\
";

MPFR_UNIOP(expm1)

static char doc_g_mpfr_eint[]="\
eint(x): returns exponential integral of x.\n\
";

MPFR_UNIOP(eint)

static char doc_g_mpfr_li2[]="\
li2(x): returns real part of dilogarithm of x.\n\
";

MPFR_UNIOP(li2)

static char doc_g_mpfr_gamma[]="\
gamma(x): returns gamma of x.\n\
";

MPFR_UNIOP(gamma)

static char doc_g_mpfr_lngamma[]="\
lngamma(x): returns logarithm of gamma(x).\n\
";

MPFR_UNIOP(lngamma)

PyDoc_STRVAR(doc_g_mpfr_lgamma,
"lgamma(x) -> (mpfr, int)\n\n"
"Return a tuple containing the logarithm of the absolute value of\n"
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

static char doc_g_mpfr_digamma[]="\
digamma(x): returns digamma of x.\n\
";

MPFR_UNIOP(digamma)

static char doc_g_mpfr_zeta[]="\
zeta(x): returns Riemann zeta of x.\n\
";

MPFR_UNIOP(zeta)

static char doc_g_mpfr_erf[]="\
erf(x): returns error function of x.\n\
";

MPFR_UNIOP(erf)

static char doc_g_mpfr_erfc[]="\
erfc(x): returns complementary error function of x.\n\
";

MPFR_UNIOP(erfc)

static char doc_g_mpfr_j0[]="\
j0(x): returns first kind Bessel function of order 0 of x.\n\
";

MPFR_UNIOP(j0)

static char doc_g_mpfr_j1[]="\
j1(x): returns first kind Bessel function of order 1 of x.\n\
";

MPFR_UNIOP(j1)

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

static char doc_g_mpfr_y0[]="\
y0(x): returns second kind Bessel function of order 0 of x.\n\
";

MPFR_UNIOP(y0)

PyDoc_STRVAR(doc_g_mpfr_y1,
"y1(x) -> mpfr\n\n"
"Return second kind Bessel function of order 1 of x.");

MPFR_UNIOP(y1)

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

PyDoc_STRVAR(doc_g_mpfr_ai,
"ai(x) -> mpfr\n\n"
"Return Airy function of x.");

MPFR_UNIOP(ai)

static PyObject *
Pympfr_add(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "add() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_add(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("add()");
}

static PyObject *
Pympfr_sub(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "sub() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_sub(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("sub()");
}

static PyObject *
Pympfr_mul(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "mul() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_mul(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("mul()");
}

static PyObject *
Pympfr_div(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "div() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_div(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("div()");
}

PyDoc_STRVAR(doc_g_mpfr_fmod,
"fmod(x, y) -> mpfr\n\n"
"Return x - n * y where n is the integer quotient of x/y, rounded to 0.");

static PyObject *
Pympfr_fmod(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "fmod() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_fmod(result->f, Pympfr_AS_MPFR(self),
                           Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("fmod()");
}

PyDoc_STRVAR(doc_g_mpfr_remainder,
"remainder(x, y) -> mpfr\n\n"
"Return x - n * y where n is the integer quotient of x/y, rounded to\n"
"the nearest integer and ties rounded to even.");

static PyObject *
Pympfr_remainder(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "remainder() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_remainder(result->f, Pympfr_AS_MPFR(self),
                                Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("remainder()");
}

PyDoc_STRVAR(doc_g_mpfr_remquo,
"remquo(x, y) -> (mpfr, int)\n\n"
"Return a tuple containing the remainder(x,y) and the low bits of the\n"
"quotient.");

static PyObject *
Pympfr_remquo(PyObject* self, PyObject *args)
{
    PyObject *result, *other;
    PympfrObject *value;
    long quobits = 0;

    PARSE_TWO_MPFR_ARGS(other, "remquo() requires 'mpfr', 'mpfr' argument");

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

PyDoc_STRVAR(doc_g_mpfr_frexp,
"frexp(x) -> (int, mpfr)\n\n"
"Return a tuple containing the exponent and mantissa of x quotient.");

static PyObject *
Pympfr_frexp(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympfrObject *value;
    mpfr_exp_t exp = 0;

    PARSE_ONE_MPFR_OTHER("frexp() requires 'mpfr' argument");

    value = Pympfr_new(0);
    result = PyTuple_New(2);
    if (!value || !result)
        goto done;

    mpfr_clear_flags();
    value->rc = mpfr_frexp(&exp, value->f, Pympfr_AS_MPFR(self),
                           context->now.mpfr_round);
    MERGE_FLAGS;
    CHECK_FLAGS("frexp()");

  done:
    Py_DECREF(self);
    Py_DECREF(other);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)value);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, PyIntOrLong_FromSsize_t((Py_ssize_t)exp));
        PyTuple_SET_ITEM(result, 1, (PyObject*)value);
    }
    return result;
}

PyDoc_STRVAR(doc_g_mpfr_atan2,
"atan2(y, x) -> mpfr\n\n"
"Return arc-tangent of (y/x).");

static PyObject *
Pympfr_atan2(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "atan2() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_atan2(result->f, Pympfr_AS_MPFR(self),
                            Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("atan2()");
}

PyDoc_STRVAR(doc_g_mpfr_agm,
"agm(x, y) -> mpfr\n\n"
"Return arithmetic-geometric mean of x and y.");

static PyObject *
Pympfr_agm(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "agm() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_agm(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("agm()");
}

PyDoc_STRVAR(doc_g_mpfr_hypot,
"hypot(y, x) -> mpfr\n\n"
"Return square root of (x**2 + y**2).");

static PyObject *
Pympfr_hypot(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "hypot() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_hypot(result->f, Pympfr_AS_MPFR(self),
                            Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("hypot()");
}

PyDoc_STRVAR(doc_g_mpfr_max,
"max(y, x) -> mpfr\n\n"
"Return maximum of x and y.");

static PyObject *
Pympfr_max(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "max() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_max(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("max()");
}

PyDoc_STRVAR(doc_g_mpfr_min,
"min(y, x) -> mpfr\n\n"
"Return minimum of x and y.");

static PyObject *
Pympfr_min(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "min() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(0)))
        goto done;

    mpfr_clear_flags();
    result->rc = mpfr_min(result->f, Pympfr_AS_MPFR(self),
                          Pympfr_AS_MPFR(other), context->now.mpfr_round);
    MPFR_CLEANUP_SELF_OTHER("min()");
}

PyDoc_STRVAR(doc_g_mpfr_nexttoward,
"next_toward(y, x) -> mpfr\n\n"
"Return the next 'mpfr' from x in the direction of y.");

static PyObject *
Pympfr_nexttoward(PyObject *self, PyObject *args)
{
    PympfrObject *result;
    PyObject *other;

    PARSE_TWO_MPFR_ARGS(other, "next_toward() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(self)))))
        goto done;

    mpfr_clear_flags();
    mpfr_set(result->f, Pympfr_AS_MPFR(self), context->now.mpfr_round);
    mpfr_nexttoward(result->f, Pympfr_AS_MPFR(other));
    result->rc = 0;
    MPFR_CLEANUP_SELF_OTHER("next_toward()");
}

PyDoc_STRVAR(doc_g_mpfr_nextabove,
"next_above(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward +Infinity.");

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

PyDoc_STRVAR(doc_g_mpfr_nextbelow,
"next_below(x) -> mpfr\n\n"
"Return the next 'mpfr' from x toward -Infinity.");

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

    PARSE_TWO_MPFR_ARGS(other, "is_lessgreater() requires 'mpfr','mpfr' arguments");

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

    PARSE_TWO_MPFR_ARGS(other, "unordered() requires 'mpfr','mpfr' arguments");

    temp = mpfr_unordered_p(Pympfr_AS_MPFR(self), Pympfr_AS_MPFR(other));
    Py_DECREF(self);
    Py_DECREF(other);
    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

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

PyDoc_STRVAR(doc_g_mpfr_fsum,
"fsum(iterable) -> mpfr\n\n"
"Return an accurate sum of the values in the iterable.");

static PyObject *
Pympfr_fsum(PyObject *self, PyObject *other)
{
    PympfrObject *temp, *result;
    mpfr_ptr *tab;
    int errcode;
    Py_ssize_t i, seq_length = 0;
    
    if (!(result = Pympfr_new(0)))
        return NULL;
        
    if (!(other = PySequence_List(other))) {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("argument must be an iterable");
        return NULL;
    }
        
    /* other contains a new list containing all the values from the
     * iterable. Now make sure each item in the list is an mpfr.
     */
     
    seq_length = PyList_GET_SIZE(other);
    for (i=0; i < seq_length; i++) {
        if (!(temp = Pympfr_From_Real(PyList_GET_ITEM(other, i), 0))) {
            Py_DECREF(other);
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("all items in iterable must be real numbers");
            return NULL;
        }
      
        errcode = PyList_SetItem(other, i,(PyObject*)temp);
        if (errcode < 0) {
            Py_DECREF(other);
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("all items in iterable must be real numbers");
            return NULL;
        }
    }
    
    /* create an array of pointers to the mpfr_t field of a Pympfr object */
    
    if (!(tab = (mpfr_ptr *)GMPY_MALLOC((sizeof(mpfr_srcptr) * seq_length)))) {
        Py_DECREF(other);
        Py_DECREF((PyObject*)result);
        return PyErr_NoMemory();
    }
    for (i=0; i < seq_length; i++) {
        temp = (PympfrObject*)PyList_GET_ITEM(other, i);
        tab[i] = temp->f;
    }
    result->rc = mpfr_sum(result->f, tab, seq_length, context->now.mpfr_round);
    Py_DECREF(other);
    GMPY_FREE(tab);

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpfr_degrees,
"degrees(x) -> mpfr\n\n"
"Convert angle x from radians to degrees.");

static PyObject *
Pympfr_degrees(PyObject *self, PyObject *other)
{
    PympfrObject *result, *temp;

    PARSE_ONE_MPFR_OTHER("degrees() requires 'mpfr' argument");

    result = Pympfr_new(0);
    temp = Pympfr_new(context->now.mpfr_prec + 20);
    if (!result || !temp) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_clear_flags();
    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_ui_div(temp->f, 180, temp->f, MPFR_RNDN);
    mpfr_mul(result->f, temp->f, Pympfr_AS_MPFR(self), MPFR_RNDN);
    Py_DECREF((PyObject*)temp);
    MPFR_CLEANUP_SELF("degrees()");
}

PyDoc_STRVAR(doc_g_mpfr_radians,
"radians(x) -> mpfr\n\n"
"Convert angle x from degrees to radians.");

static PyObject *
Pympfr_radians(PyObject *self, PyObject *other)
{
    PympfrObject *result, *temp;

    PARSE_ONE_MPFR_OTHER("radians() requires 'mpfr' argument");

    result = Pympfr_new(0);
    temp = Pympfr_new(context->now.mpfr_prec + 20);
    if (!result || !temp) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_clear_flags();
    mpfr_const_pi(temp->f, MPFR_RNDN);
    mpfr_div_ui(temp->f, temp->f, 180, MPFR_RNDN);
    mpfr_mul(result->f, Pympfr_AS_MPFR(self), temp->f, MPFR_RNDN);
    Py_DECREF((PyObject*)temp);
    MPFR_CLEANUP_SELF("radians()");
}

PyDoc_STRVAR(doc_mpfr_format,
"x.__format__(fmt) -> string\n\n"
"Return a Python string by formatting 'x' using the format string\n"
"'fmt'. A valid format string consists of:\n"
"     optional alignment code:\n"
"        '<' -> left shifted in field\n"
"        '>' -> right shifted in field\n"
"        '^' -> centered in field\n"
"     optional leading sign code\n"
"        '+' -> always display leading sign\n"
"        '-' -> only display minus for negative values\n"
"        ' ' -> minus for negative values, space for positive values\n"
"     optional width.precision\n"
"     optional rounding mode:\n"
"        'U' -> round toward plus Infinity\n"
"        'D' -> round toward minus Infinity\n"
"        'Y' -> round away from zero\n"
"        'Z' -> round toward zero\n"
"        'N' -> round to nearest\n"
"     optional conversion code:\n"
"        'a','A' -> hex format\n"
"        'b'     -> binary format\n"
"        'e','E' -> scientific format\n"
"        'f','F' -> fixed point format\n"
"        'g','G' -> fixed or float format\n\n"
"The default format is '.6f'.");

static PyObject *
Pympfr_format(PyObject *self, PyObject *args)
{
    PyObject *result = 0, *mpfrstr = 0;
    char *buffer = 0, *newbuf = 0, *fmtcode = 0, *p1, *p2, *p3;
    char mpfrfmt[100], fmt[30];
    int buflen;
    int seensign = 0, seenalign = 0, seendecimal = 0, seendigits = 0;
    int seenround = 0, seenconv = 0;

    if (!Pympfr_Check(self)) {
        TYPE_ERROR("requires mpfr type");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "s", &fmtcode))
        return NULL;

    p2 = mpfrfmt;
    p3 = fmt;
    *(p2++) = '%';

    for (p1 = fmtcode; *p1 != '\00'; p1++) {
        if (*p1 == '<' || *p1 == '>' || *p1 == '^') {
            if (seenalign || seensign || seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p3++) = *p1;
                seenalign = 1;
                continue;
            }
        }
        if (*p1 == '+' || *p1 == ' ') {
            if (seensign || seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '-') {
            if (seensign || seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '.') {
            if (seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seendecimal = 1;
                continue;
            }
        }
        if (isdigit(*p1)) {
            if (seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else if (seendecimal) {
                *(p2++) = *p1;
                continue;
            }
            else {
                if (p3 == fmt) {
                    *(p3++) = '>';
                    seenalign = 1;
                }
                *(p3++) = *p1;
                continue;
            }
        }
        if (!seendigits) {
            seendigits = 1;
            *(p2++) = 'R';
        }
        if (*p1 == 'U' || *p1 == 'D' || *p1 == 'Y' || *p1 == 'Z' ||
            *p1 == 'N' ) {
            if (seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seenround = 1;
                continue;
            }
        }
        if (*p1 == 'a' || *p1 == 'A' || *p1 == 'b' || *p1 == 'e' ||
            *p1 == 'E' || *p1 == 'f' || *p1 == 'F' || *p1 == 'g' ||
            *p1 == 'G' ) {
            *(p2++) = *p1;
            seenconv = 1;
            break;
        }
        VALUE_ERROR("Invalid conversion specification");
        return NULL;
    }

    if (!seendigits)
        *(p2++) = 'R';
    if (!seenconv)
        *(p2++) = 'f';

    *(p2) = '\00';
    *(p3) = '\00';

    buflen = mpfr_asprintf(&buffer, mpfrfmt, Pympfr_AS_MPFR(self));

    /* If there isn't a decimal point in the output and the output
     * only consists of digits, then append .0 */
    if (strlen(buffer) == strspn(buffer, "+- 0123456789")) {
        newbuf = GMPY_MALLOC(buflen + 3);
        if (!newbuf) {
            mpfr_free_str(buffer);
            return PyErr_NoMemory();
        }
        *newbuf = '\0';
        strcat(newbuf, buffer);
        strcat(newbuf, ".0");
        mpfr_free_str(buffer);
        mpfrstr = Py_BuildValue("s", newbuf);
        GMPY_FREE(newbuf);
    }
    else {
        mpfrstr = Py_BuildValue("s", buffer);
        mpfr_free_str(buffer);
    }
    if (!mpfrstr) {
        return NULL;
    }

    result = PyObject_CallMethod(mpfrstr, "__format__", "(s)", fmt);
    Py_DECREF(mpfrstr);
    return result;
}

#ifdef PY3
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,              /* nb_negative             */
    (unaryfunc) Pympfr_pos,              /* nb_positive             */
    (unaryfunc) Pympfr_abs,              /* nb_absolute             */
    (inquiry) Pympfr_nonzero,            /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympfr2PyLong,           /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympfr2PyFloat,          /* nb_float                */
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
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympfr_neg,              /* nb_negative             */
    (unaryfunc) Pympfr_pos,              /* nb_positive             */
    (unaryfunc) Pympfr_abs,              /* nb_absolute             */
    (inquiry) Pympfr_nonzero,            /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympfr2PyInt,            /* nb_int                  */
    (unaryfunc) Pympfr2PyLong,           /* nb_long                 */
    (unaryfunc) Pympfr2PyFloat,          /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add          */
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

static PyGetSetDef Pympfr_getseters[] =
{
    {"precision", (getter)Pympfr_getprec_attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)Pympfr_getrc_attrib, NULL, "return code", NULL},
    {"imag", (getter)Pympfr_getimag_attrib, NULL, "imaginary component", NULL},
    {"real", (getter)Pympfr_getreal_attrib, NULL, "real component", NULL},
    {NULL}
};

static PyMethodDef Pympfr_methods [] =
{
    { "__ceil__", Pympfr_ceil, METH_NOARGS, doc_mpfr_ceil },
    { "__floor__", Pympfr_floor, METH_NOARGS, doc_mpfr_floor },
    { "__format__", Pympfr_format, METH_VARARGS, doc_mpfr_format },
    { "__trunc__", Pympfr_trunc, METH_NOARGS, doc_mpfr_trunc },
    { "as_integer_ratio", Pympfr_integer_ratio, METH_NOARGS, doc_mpfr_integer_ratio },
    { "as_mantissa_exp", Pympfr_mantissa_exp, METH_NOARGS, doc_mpfr_mantissa_exp },
    { "as_simple_fraction", (PyCFunction)Pympfr_simple_fraction, METH_VARARGS | METH_KEYWORDS, doc_mpfr_simple_fraction },
    { "binary", Pympany_binary, METH_NOARGS, doc_binarym },
    { "conjugate", Pympfr_conjugate, METH_NOARGS, doc_mpfr_conjugate },
    { "digits", Pympfr_digits, METH_VARARGS, doc_mpfr_digits },
    { "is_inf", Pympfr_is_inf, METH_NOARGS, doc_mpfr_is_inf },
    { "is_integer", Pympfr_is_integer, METH_NOARGS, doc_mpfr_is_integer },
    { "is_lessgreater", Pympfr_is_lessgreater, METH_VARARGS, doc_mpfr_is_lessgreater },
    { "is_nan", Pympfr_is_nan, METH_NOARGS, doc_mpfr_is_nan },
    { "is_number", Pympfr_is_number, METH_NOARGS, doc_mpfr_is_number },
    { "is_regular", Pympfr_is_regular, METH_NOARGS, doc_mpfr_is_regular },
    { "is_signed", Pympfr_is_signed, METH_NOARGS, doc_mpfr_is_signed },
    { "is_unordered", Pympfr_is_unordered, METH_VARARGS, doc_mpfr_is_unordered },
    { "is_zero", Pympfr_is_zero, METH_NOARGS, doc_mpfr_is_zero },
    { NULL, NULL, 1 }
};

static PyTypeObject Pympfr_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpfr",                                 /* tp_name          */
    sizeof(PympfrObject),                   /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympfr_dealloc,            /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympfr2repr,                 /* tp_repr          */
    &mpfr_number_methods,                   /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympfr_hash,                 /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympfr2str,                  /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "GNU Multi Precision floating point",   /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympfr_methods,                         /* tp_methods       */
        0,                                  /* tp_members       */
    Pympfr_getseters,                       /* tp_getset        */
};

