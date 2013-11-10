/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013 Case Van Horsen            *
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

/* This file contains all the conversion functions for gmpy2.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. gmpy2 uses a LBYL (Look Before You Leap) approach and
 * identifies the numeric type before conversion before conversion to a gmpy2
 * type. The basic operations (+, -, *, /) are optimized to directly work with
 * some basic types such as C longs or doubles.
 *
 * Support for the Decimal type is a challenge. For the basic operations, it
 * is most accurate to convert a Decimal instance into an mpq and then use
 * MPFR's functions to accurately operate on an mpfr and mpq. This approach is
 * challenging because (1) a large exponent can create a very large mpq and
 * (2) the changes made to C-coded version of Decimal in Python 3.3.
 *
 */

/* Functions that operate strictly on mpfr. */

/* Make a copy of an mpfr object. If bits is 0, the new object will have
 * the same precision as the original object. If the requested precision
 * is less than the precision of the original object, the new object
 * will be rounded to requested precision using the current rounding mode.
 */

static MPFR_Object *
Pympfr_From_Pympfr(PyObject *self, mpfr_prec_t bits)
{
    MPFR_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (bits == 0)
        bits = mpfr_get_prec(MPFR(self));

    if ((result = (MPFR_Object*)Pympfr_new(bits))) {
        result->rc = mpfr_set(result->f,
                              MPFR(self),
                              context->ctx.mpfr_round);
    }

    return result;
}

/* Return a copy of an mpfr, using the precision of the context argument. */

static MPFR_Object *
Pympfr_From_Pympfr_context(PyObject *self, GMPyContextObject *context)
{
    MPFR_Object *result;

    if ((result = (MPFR_Object*)Pympfr_new_context(context))) {
        result->rc = mpfr_set(result->f,
                              MPFR(self),
                              context->ctx.mpfr_round);
    }

    return result;
}

static MPFR_Object *
Pympfr_From_PyFloat(PyObject *self, mpfr_prec_t bits)
{
    MPFR_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPFR_Object*)Pympfr_new(bits))) {
        result->rc = mpfr_set_d(result->f,
                                PyFloat_AS_DOUBLE(self),
                                context->ctx.mpfr_round);
    }

    return result;
}

static MPFR_Object *
Pympfr_From_PyFloat_bits_context(PyObject *self,
                                 mpfr_prec_t bits,
                                 GMPyContextObject *context)
{
    MPFR_Object *result;

    if ((result = (MPFR_Object*)Pympfr_new_bits_context(bits, context))) {
        result->rc = mpfr_set_d(result->f,
                                PyFloat_AS_DOUBLE(self),
                                GET_MPFR_ROUND(context));
    }

    return result;
}

static MPFR_Object *
Pympfr_From_Pympz(PyObject *self, mpfr_prec_t bits)
{
    MPFR_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPFR_Object*)Pympfr_new(bits))) {
        result->rc = mpfr_set_z(result->f,
                                MPZ(self),
                                context->ctx.mpfr_round);
    }

    return result;
}

static MPFR_Object *
Pympfr_From_Pympz_context(PyObject *self,
                          mpfr_prec_t bits,
                          GMPyContextObject *context)
{
    MPFR_Object *result;

    if ((result = (MPFR_Object*)Pympfr_new_bits_context(bits, context))) {
        result->rc = mpfr_set_z(result->f,
                                MPZ(self),
                                GET_MPFR_ROUND(context));
    }

    return result;
}

#define Pympfr_From_Pyxmpz Pympfr_From_Pympz
#define Pympfr_From_Pyxmpz_context Pympfr_From_Pympz_context

static MPZ_Object *
Pympfr_To_Pympz(PyObject *self)
{
    MPZ_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = GMPy_MPZ_New())) {
        if (mpfr_nan_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, MPFR(self), context->ctx.mpfr_round);
    }

    return result;
}

static XMPZ_Object *
Pympfr_To_Pyxmpz(PyObject *self)
{
    XMPZ_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = GMPy_XMPZ_New())) {
        if (mpfr_nan_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, MPFR(self), context->ctx.mpfr_round);
    }

    return result;
}

/* Return the simpliest rational number that approximates 'self' to the
 * requested precision 'err'. If 'err' is negative, then the requested
 * precision is -2**abs(int(err)). If 'err' is NULL, then the requested
 * precision is -2**prec. If 'prec' is 0, then the requested precision is
 * the precision of 'self'.
 */

static MPQ_Object *
stern_brocot(MPFR_Object* self, MPFR_Object *err, mpfr_prec_t prec, int mayz)
{
    MPQ_Object *result = 0;
    int i, negative, errsign;
    mpfr_t f, al, a, r1[3], r2[3], minerr, curerr, newerr, temp;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

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
        prec = (mpfr_prec_t)(-mpfr_get_si(err->f, context->ctx.mpfr_round));

    if (errsign <= 0 && (prec < 2 || prec > mpfr_get_prec(self->f))) {
        VALUE_ERROR("Requested precision out-of-bounds.");
        return NULL;
    }

    if (!(result = GMPy_MPQ_New()))
        return NULL;

    mpfr_init2(minerr, F2Q_PREC);
    if (errsign <= 0) {
        mpfr_set_ui(minerr, 1, MPFR_RNDN);
        mpfr_div_2si(minerr, minerr, prec, context->ctx.mpfr_round);
    }
    else {
        mpfr_set(minerr, err->f, context->ctx.mpfr_round);
    }

    mpfr_init2(f, prec);
    if (mpfr_sgn(self->f) < 0) {
        negative = 1;
        mpfr_abs(f, self->f, context->ctx.mpfr_round);
    }
    else {
        negative = 0;
        mpfr_set(f, self->f, context->ctx.mpfr_round);
    }

    mpfr_init2(al, prec);
    mpfr_set(al, f, context->ctx.mpfr_round);
    mpfr_init2(a, prec);
    mpfr_floor(a, al);
    mpfr_init2(temp, prec);
    for (i=0; i<3; ++i) {
        mpfr_init2(r1[i], prec);
        mpfr_init2(r2[i], prec);
    }
    mpfr_set_si(r1[0], 0, MPFR_RNDN);
    mpfr_set_si(r1[1], 0, MPFR_RNDN);
    mpfr_set_si(r1[2], 1, MPFR_RNDN);
    mpfr_set_si(r2[0], 0, MPFR_RNDN);
    mpfr_set_si(r2[1], 1, MPFR_RNDN);
    mpfr_set(r2[2], a, context->ctx.mpfr_round);
    mpfr_init2(curerr, F2Q_PREC);
    mpfr_init2(newerr, F2Q_PREC);
    mpfr_reldiff(curerr, f, a, context->ctx.mpfr_round);
    while (mpfr_cmp(curerr, minerr) > 0) {
        mpfr_sub(temp, al, a, context->ctx.mpfr_round);
        mpfr_ui_div(al, 1, temp, context->ctx.mpfr_round);
        mpfr_floor(a, al);
        mpfr_swap(r1[0], r1[1]);
        mpfr_swap(r1[1], r1[2]);
        mpfr_mul(r1[2], r1[1], a, context->ctx.mpfr_round);
        mpfr_add(r1[2], r1[2], r1[0], context->ctx.mpfr_round);
        mpfr_swap(r2[0], r2[1]);
        mpfr_swap(r2[1], r2[2]);
        mpfr_mul(r2[2], r2[1], a, context->ctx.mpfr_round);
        mpfr_add(r2[2], r2[2], r2[0], context->ctx.mpfr_round);
        mpfr_div(temp, r2[2], r1[2], context->ctx.mpfr_round);
        mpfr_reldiff(newerr, f, temp, context->ctx.mpfr_round);
        if (mpfr_cmp(curerr, newerr) <= 0) {
            mpfr_swap(r1[1],r1[2]);
            mpfr_swap(r2[1],r2[2]);
            break;
        }
        mpfr_swap(curerr, newerr);
    }

    if (mayz && (mpfr_cmp_ui(r1[2],1) == 0)) {
        Py_DECREF((PyObject*)result);
        result = (MPQ_Object*)GMPy_MPZ_New();
        mpfr_get_z(MPZ(result), r2[2], context->ctx.mpfr_round);
        if (negative)
            mpz_neg(MPZ(result), MPZ(result));
    }
    else {
        mpfr_get_z(mpq_numref(result->q), r2[2], context->ctx.mpfr_round);
        mpfr_get_z(mpq_denref(result->q), r1[2], context->ctx.mpfr_round);
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

static MPQ_Object *
Pympfr_To_Pympq(PyObject *self)
{
    return stern_brocot((MPFR_Object*)self, 0, 0, 0);
}

static MPFR_Object *
Pympfr_From_Pympq(PyObject *self, mpfr_prec_t bits)
{
    MPFR_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPFR_Object*)Pympfr_new(bits)))
        result->rc = mpfr_set_q(result->f, MPQ(self),
                                context->ctx.mpfr_round);
    return result;
}

static MPFR_Object *
Pympfr_From_Pympq_bits_context(PyObject *self, mpfr_prec_t bits,
                               GMPyContextObject *context)
{
    MPFR_Object *result;

    if ((result = (MPFR_Object*)Pympfr_new_bits_context(bits, context)))
        result->rc = mpfr_set_q(result->f, MPQ(self),
                                context->ctx.mpfr_round);
    return result;
}

static MPFR_Object *
Pympfr_From_PyLong(PyObject *self, mpfr_prec_t bits)
{
    MPFR_Object *result;
    PyObject *temp = (PyObject*)GMPy_MPZ_From_PyIntOrLong(self);

    if (!temp)
        return NULL;
    result = Pympfr_From_Pympz(temp, bits);
    Py_DECREF(temp);
    return result;
}

static MPFR_Object *
Pympfr_From_PyLong_context(PyObject *self, mpfr_prec_t bits,
                          GMPyContextObject *context)
{
    MPFR_Object *result;
    PyObject *temp = (PyObject*)GMPy_MPZ_From_PyIntOrLong(self);

    if (!temp)
        return NULL;
    result = Pympfr_From_Pympz_context(temp, bits, context);
    Py_DECREF(temp);
    return result;
}

#ifdef PY2
static MPFR_Object *
Pympfr_From_PyInt(PyObject *self, mpfr_prec_t bits)
{
    MPFR_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPFR_Object*)Pympfr_new(bits)))
        result->rc = mpfr_set_si(result->f, PyInt_AsLong(self),
                                 context->ctx.mpfr_round);
    return result;
}

static MPFR_Object *
Pympfr_From_PyInt_bits_context(PyObject *self, mpfr_prec_t bits,
                          GMPyContextObject *context)
{
    MPFR_Object *result;

    if ((result = (MPFR_Object*)Pympfr_new_bits_context(bits, context)))
        result->rc = mpfr_set_si(result->f, PyInt_AsLong(self),
                                 context->ctx.mpfr_round);
    return result;
}

static PyObject *
Pympfr_To_PyInt(MPFR_Object *self)
{
    PyObject *result;
    MPZ_Object *temp = Pympfr_To_Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = GMPy_PyIntOrLong_From_MPZ(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}
#endif

static MPFR_Object *
Pympfr_From_PyStr(PyObject *s, int base, mpfr_prec_t bits)
{
    MPFR_Object *result;
    char *cp, *endptr;
    mpfr_prec_t prec;
    Py_ssize_t len;
    PyObject *ascii_str = NULL;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = PyBytes_AsString(s);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = PyBytes_AsString(ascii_str);
    }

    if (bits > 0)
        prec = bits;
    else
        prec = context->ctx.mpfr_prec;

    if (!(result = (MPFR_Object*)Pympfr_new(prec))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* delegate the rest to MPFR */
    result->rc = mpfr_strtofr(result->f, cp, &endptr, base,
                              context->ctx.mpfr_round);

    if (len != (Py_ssize_t)(endptr - cp)) {
        VALUE_ERROR("invalid digits");
        Py_DECREF((PyObject*)result);
        Py_XDECREF(ascii_str);
        return NULL;
    }
    Py_XDECREF(ascii_str);

    return result;
}

static MPFR_Object *
Pympfr_From_PyStr_context(PyObject *s, int base, mpfr_prec_t bits,
                          GMPyContextObject *context)
{
    MPFR_Object *result;
    char *cp, *endptr;
    mpfr_prec_t prec;
    Py_ssize_t len;
    PyObject *ascii_str = NULL;

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = PyBytes_AsString(s);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = PyBytes_AsString(ascii_str);
    }

    if (bits > 0)
        prec = bits;
    else
        prec = context->ctx.mpfr_prec;

    if (!(result = (MPFR_Object*)Pympfr_new(prec))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* delegate the rest to MPFR */
    result->rc = mpfr_strtofr(result->f, cp, &endptr, base,
                              context->ctx.mpfr_round);

    if (len != (Py_ssize_t)(endptr - cp)) {
        VALUE_ERROR("invalid digits");
        Py_DECREF((PyObject*)result);
        Py_XDECREF(ascii_str);
        return NULL;
    }
    Py_XDECREF(ascii_str);

    return result;
}

static PyObject *
Pympfr_To_PyLong(MPFR_Object *self)
{
    PyObject *result;
    MPZ_Object *temp = Pympfr_To_Pympz((PyObject*)self);

    if (!temp) return NULL;

    result = GMPy_PyLong_From_MPZ(temp);
    Py_DECREF((PyObject*)temp);

    return result;
}

static PyObject *
Pympfr_To_PyFloat(MPFR_Object *self)
{
    double res;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    res = mpfr_get_d(self->f, context->ctx.mpfr_round);

    return PyFloat_FromDouble(res);
}

static PyObject*
Pympfr_To_PyStr(MPFR_Object *self, int base, int digits)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    /* check arguments are valid */
    assert(MPFR_Check((PyObject*)self));
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
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, context->ctx.mpfr_round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in Pympfr_To_PyStr");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
    mpfr_free_str(buffer);
    return result;
}

static MPFR_Object *
Pympfr_From_Decimal(PyObject* obj, mpfr_prec_t bits)
{
    MPFR_Object *result;
    MPQ_Object *temp;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    result = (MPFR_Object*)Pympfr_new_bits_context(bits, context);
    temp = Pympq_From_DecimalRaw(obj);

    if (!temp || !result) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    if (!mpz_cmp_si(mpq_numref(temp->q), 0)) {
        if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
            mpfr_set_nan(result->f);
        }
        else {
            mpfr_set_zero(result->f, mpz_sgn(mpq_denref(temp->q)));
        }
    }
    else if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
        if (mpz_cmp_si(mpq_numref(temp->q), 0) < 0) {
            mpfr_set_inf(result->f, -1);
        }
        else {
            mpfr_set_inf(result->f, 1);
        }
    }
    else {
        Py_DECREF((PyObject*)result);
        result = Pympfr_From_Pympq_bits_context((PyObject*)temp, bits, context);
    }
    Py_DECREF((PyObject*)temp);
    return result;
}

static MPFR_Object *
Pympfr_From_Decimal_context(PyObject* obj,
                            mpfr_prec_t bits,
                            GMPyContextObject *context)
{
    MPFR_Object *result;
    MPQ_Object *temp;

    result = (MPFR_Object*)Pympfr_new_bits_context(bits, context);
    temp = Pympq_From_DecimalRaw(obj);

    if (!temp || !result) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    if (!mpz_cmp_si(mpq_numref(temp->q), 0)) {
        if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
            mpfr_set_nan(result->f);
        }
        else {
            mpfr_set_zero(result->f, mpz_sgn(mpq_denref(temp->q)));
        }
    }
    else if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
        if (mpz_cmp_si(mpq_numref(temp->q), 0) < 0) {
            mpfr_set_inf(result->f, -1);
        }
        else {
            mpfr_set_inf(result->f, 1);
        }
    }
    else {
        Py_DECREF((PyObject*)result);
        result = Pympfr_From_Pympq_bits_context((PyObject*)temp, bits, context);
    }
    Py_DECREF((PyObject*)temp);
    return result;
}

/*
 * If obj is a Pympfr and bits is 0 or bits is the same as the precision of
 * obj, then a new reference is created.
 *
 * For all other numerical types with bits = 0, the conversion is rounded to
 * context->ctx.mpfr_prec.
 */

static MPFR_Object *
Pympfr_From_Real(PyObject* obj, mpfr_prec_t bits)
{
    MPFR_Object* newob = 0;
    MPQ_Object* temp = 0;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (MPFR_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpfr is still
         * valid in the current context. */
        if (!bits || mpfr_get_prec(MPFR(obj)) == bits) {
            newob = (MPFR_Object*) obj;
            Py_INCREF(obj);
        }
        else {
            newob = Pympfr_From_Pympfr((PyObject*)obj, bits);
        }
    }
    else if (MPFR_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer valid
         * and mpfr_check_range needs to be called. */
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpfr' incompatible with current context");
            return NULL;
        }
        if ((newob = (MPFR_Object*)Pympfr_new(mpfr_get_prec(MPFR(obj))))) {
            mpfr_set(newob->f, MPFR(obj), context->ctx.mpfr_round);
            newob->round_mode = ((MPFR_Object*)obj)->round_mode;
            newob->rc = ((MPFR_Object*)obj)->rc;
            newob->rc = mpfr_check_range(newob->f, newob->rc, newob->round_mode);
        }
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympfr_From_PyFloat(obj, bits);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympfr_From_PyInt(obj, bits);
#endif
    }
    else if (MPQ_Check(obj)) {
        newob = Pympfr_From_Pympq(obj, bits);
    }
    else if (MPZ_Check(obj)) {
        newob = Pympfr_From_Pympz(obj, bits);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympfr_From_PyLong(obj, bits);
    }
    else if (XMPZ_Check(obj)) {
        newob = Pympfr_From_Pyxmpz(obj, bits);
    }
    else if (IS_DECIMAL(obj)) {
        newob = Pympfr_From_Decimal(obj, bits);
    }
    else if (IS_FRACTION(obj)) {
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympfr_From_Pympq((PyObject*)temp, bits);
            Py_DECREF((PyObject*)temp);
        }
    }
    if (!newob)
        TYPE_ERROR("object could not be converted to 'mpfr'");
    return newob;
}

/* GMPy_MPFR_From_Real_Temp() converts a real number to an mpfr. When
 * converting values that can be converted exactly (i.e. floating-point using
 * radix-2 represetnation), the conversion is done with the maximum possible
 * precision. Regardless of the context's precision, the precision of the
 * returned value will not be decreased. This is done to minimize rounding
 * error. This value returned by this function is primarily intended for
 * internal use. See GMPy_MPFR_From_Real_Prec() to convert a real number to an
 * mpfr with precision and rounding controlled by the context.
 *
 * Note: Even though the precision of the value returned by ..._Temp() is
 *       not be constrained by the context, the exponent of the returned
 *       value is guaranteed to be valid as per the context.
 */

static MPFR_Object *
GMPy_MPFR_From_Real_Temp(PyObject *obj, GMPyContextObject *context)
{
    MPFR_Object *result = NULL;

    /* Check if obj is an mpfr and exponent is valid. */

    if (MPFR_CheckAndExp(obj)) {
        /* Return a new reference with the precision of the input. */
        result = (MPFR_Object*) obj;
        Py_INCREF(obj);
        return result;
    }

    /* The exponent is not valid. */

    if (MPFR_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer valid
         * and mpfr_check_range needs to be called. */
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpfr' incompatible with current context");
            return NULL;
        }
        if ((result = (MPFR_Object*)Pympfr_new_bits_context(mpfr_get_prec(MPFR(obj)), context))) {
            mpfr_set(result->f, MPFR(obj), GET_MPFR_ROUND(context));
            result->round_mode = ((MPFR_Object*)obj)->round_mode;
            result->rc = mpfr_check_range(result->f, ((MPFR_Object*)obj)->rc, result->round_mode);
        }
        return result;
    }

    /* To prevent losing precision when converting a standard Python float
     * to an temporary mpfr, we specify 53 bits of precision.
     */
    if (PyFloat_Check(obj))
        return Pympfr_From_PyFloat_bits_context(obj, 53, context);

#ifdef PY2
    if (PyInt_Check(obj))
        return Pympfr_From_PyInt_bits_context(obj, 0, context);
#endif

    if (MPQ_Check(obj))
        return Pympfr_From_Pympq_bits_context(obj, 0, context);

    if (MPZ_Check(obj))
        return Pympfr_From_Pympz_context(obj, 0, context);

    if (PyLong_Check(obj))
        return Pympfr_From_PyLong_context(obj, 0, context);

    if (XMPZ_Check(obj))
        return Pympfr_From_Pyxmpz_context(obj, 0, context);

    if (IS_DECIMAL(obj))
        return Pympfr_From_Decimal_context(obj, 0, context);

    if (IS_FRACTION(obj)) {
        MPQ_Object *tempq = NULL;

        if ((tempq = Pympq_From_Fraction(obj))) {
            result = Pympfr_From_Pympq_bits_context((PyObject*)tempq, 0, context);
            Py_DECREF((PyObject*)tempq);
        }
        return result;
    }

    TYPE_ERROR("object could not be converted to 'mpfr'");
    return NULL;
}

/* This function should eventually go away. */

static MPFR_Object *
Pympfr_From_Real_bits_context(PyObject* obj, mpfr_prec_t bits, GMPyContextObject *context)
{
    MPFR_Object* newob = 0;
    MPQ_Object* temp = 0;

    if (MPFR_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpfr is still
         * valid in the current context. */
        if (!bits || mpfr_get_prec(MPFR(obj)) == bits) {
            newob = (MPFR_Object*) obj;
            Py_INCREF(obj);
        }
        else {
            newob = Pympfr_From_Pympfr((PyObject*)obj, bits);
        }
    }
    else if (MPFR_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer valid
         * and mpfr_check_range needs to be called. */
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpfr' incompatible with current context");
            return NULL;
        }
        if ((newob = (MPFR_Object*)Pympfr_new_bits_context(mpfr_get_prec(MPFR(obj)),
                                                            context))) {
            mpfr_set(newob->f, MPFR(obj), GET_MPFR_ROUND(context));
            newob->round_mode = ((MPFR_Object*)obj)->round_mode;
            newob->rc = ((MPFR_Object*)obj)->rc;
            newob->rc = mpfr_check_range(newob->f, newob->rc, newob->round_mode);
        }
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympfr_From_PyFloat_bits_context(obj, bits, context);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympfr_From_PyInt_bits_context(obj, bits, context);
#endif
    }
    else if (MPQ_Check(obj)) {
        newob = Pympfr_From_Pympq_bits_context(obj, bits, context);
    }
    else if (MPZ_Check(obj)) {
        newob = Pympfr_From_Pympz_context(obj, bits, context);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympfr_From_PyLong_context(obj, bits, context);
    }
    else if (XMPZ_Check(obj)) {
        newob = Pympfr_From_Pyxmpz_context(obj, bits, context);
    }
    else if (IS_DECIMAL(obj)) {
        newob = Pympfr_From_Decimal_context(obj, bits, context);
    }
    else if (IS_FRACTION(obj)) {
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympfr_From_Pympq_bits_context((PyObject*)temp, bits, context);
            Py_DECREF((PyObject*)temp);
        }
    }
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
    MPFR_Object* newob = Pympfr_From_Real(arg, 0);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpfr'");
        return 0;
    }
}

/* str and repr implementations for mpfr */
static PyObject *
Pympfr_To_Str(MPFR_Object *self)
{
    PyObject *result, *temp;
    long precision;
    char fmtstr[30];

    precision = (long)(log10(2) * (double)mpfr_get_prec(MPFR(self))) + 2;

    sprintf(fmtstr, "{0:.%ldg}", precision);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympfr_To_Repr(MPFR_Object *self)
{
    PyObject *result, *temp;
    long precision, bits;
    char fmtstr[30];

    bits = mpfr_get_prec(MPFR(self));
    precision = (long)(log10(2) * (double)bits) + 2;

    if (mpfr_number_p(MPFR(self)) && bits != DBL_MANT_DIG)
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

static PyObject *
raw_mpfr_ascii(mpfr_t self, int base, int digits, int round)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;

    /* Process special cases first */
    if (!(mpfr_regular_p(self))) {
        if (mpfr_nan_p(self)) {
            return Py_BuildValue("(sii)", "nan", 0, 0);
        }
        else if (mpfr_inf_p(self) && !mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "inf", 0, 0);
        }
        else if (mpfr_inf_p(self) && mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "-inf", 0, 0);
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "-0", 0, mpfr_get_prec(self));
        }
        else {
            return Py_BuildValue("(sii)", "0", 0, mpfr_get_prec(self));
        }
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self, round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in raw_mpfr_ascii");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self));
    mpfr_free_str(buffer);
    return result;
}
