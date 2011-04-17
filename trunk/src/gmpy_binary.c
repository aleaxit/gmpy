/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_binary.c                                                           *
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

/* Conversion routines between GMPY2 objects and a compact, portable
 * binary representation. The binary format of GMPY2 is not compatible
 * with GMPY 1.x. Methods to read the old format are provided.
 */

/* Provide functions to access the old binary formats. */

PyDoc_STRVAR(doc_g_mpz_from_old_binary,
"mpz_from_old_binary(string) -> mpz\n\n"
"Return an mpz from a GMPY 1.x binary format.");
static PyObject *
Pympz_From_Old_Binary(PyObject *self, PyObject *other)
{
    unsigned char *cp;
    Py_ssize_t len;
    int negative = 0;
    PympzObject *result;

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("mpz_from_old_binary() requires bytes argument");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;

    len = PyBytes_Size(other);
    cp = (unsigned char*)PyBytes_AsString(other);

    if (cp[len-1] == 0xFF) {
        negative = 1;
        --len;
    }
    mpz_import(result->z, len, -1, sizeof(char), 0, 0, cp);
    if (negative)
        mpz_neg(result->z, result->z);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_g_mpq_from_old_binary,
"mpq_from_old_binary(string) -> mpq\n\n"
"Return an mpq from a GMPY 1.x binary format.");
static PyObject *
Pympq_From_Old_Binary(PyObject *self, PyObject *other)
{
    unsigned char *cp;
    Py_ssize_t len;
    int topper, negative, numlen;
    mpz_t numerator, denominator;
    PympqObject *result;

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("mpq_from_old_binary() requires bytes argument");
        return NULL;
    }

    if (!(result = Pympq_new()))
        return NULL;

    len = PyBytes_Size(other);
    cp = (unsigned char*)PyBytes_AsString(other);

    if (len < 6) {
        VALUE_ERROR("invalid mpq binary (too short)");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    topper = cp[3] & 0x7f;
    negative = cp[3] & 0x80;
    numlen = cp[0] + 256 * (cp[1] + 256 * (cp[2] + 256 * topper));
    if (len < (4 + numlen + 1)) {
        VALUE_ERROR("invalid mpq binary (num len)");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_inoc(numerator);
    mpz_inoc(denominator);
    mpz_import(numerator, numlen, -1, sizeof(char), 0, 0, cp+4);
    mpz_import(denominator, len-4-numlen, -1, sizeof(char), 0, 0, cp+4+numlen);
    if (negative)
        mpz_neg(numerator, numerator);

    mpq_set_num(result->q, numerator);
    mpq_set_den(result->q, denominator);
    mpq_canonicalize(result->q);
    mpz_cloc(numerator);
    mpz_cloc(denominator);
    return (PyObject*)result;
}

#ifdef WITHMPFR
PyDoc_STRVAR(doc_g_mpfr_from_old_binary,
"mpfr_from_old_binary(string) -> mpfr\n\n"
"Return an mpfr from a GMPY 1.x binary mpf format.");
static PyObject *
Pympfr_From_Old_Binary(PyObject *self, PyObject *other)
{
    unsigned char *cp;
    Py_ssize_t len;
    PympfrObject *result;
    mpfr_t digit;
    mpfr_prec_t prec;
    int i, codebyte, resusign, exposign, resuzero, precilen;
    unsigned int expomag = 0;

    if (!(PyBytes_Check(other))) {
        TYPE_ERROR("mpfr_from_old_binary() requires bytes argument");
        return NULL;
    }

    len = PyBytes_Size(other);
    cp = (unsigned char*)PyBytes_AsString(other);

    if (len == 1) {
        prec = 0;
    }
    else {
        prec = (mpfr_prec_t)(8 * (len - 5));
        if ((len>=5) && (cp[0]&8)) {
            prec = 0;
            for (i=4; i>0; --i) {
                prec = (prec << 8) | cp[i];
            }
        }
    }

    /*
     * binary format for MP floats: first, a code-byte, then, a LSB
     * 4-byte unsigned int (exponent magnitude), then the "mantissa"
     * (actually, "significand", but "mantissa" is the usual term...)
     * in MSB form.
     *
     * The codebyte encodes both the signs, exponent and result, or
     * also the zeroness of the result (in which case, nothing more).
     */
    codebyte = cp[0];
    resusign = codebyte & 1;
    exposign = codebyte & 2;
    resuzero = codebyte & 4;
    precilen = (codebyte & 8)?4:0;

    /* mpfr zero has a very compact (1-byte) binary encoding!-) */
    if (resuzero) {
        if (!(result = Pympfr_new(prec)))
            return NULL;
        result->rc = mpfr_set_ui(result->f, 0, context->now.mpfr_round);
        return (PyObject*)result;
    }

    /* all other numbers are 6+ bytes: codebyte, 4-byte exp, 1+
     * bytes for the mantissa; check this string is 6+ bytes
     */
    if (len < 6 + precilen) {
        VALUE_ERROR("invalid mpf binary encoding (too short)");
        return NULL;
    }

    if (!(result = Pympfr_new(prec)))
        return NULL;

    /* reconstruct exponent */
    for (i = 4 + precilen; i > precilen; --i) {
        expomag = (expomag<<8) | cp[i];
    }

    /* reconstruct 'mantissa' (significand) */
    mpfr_set_si(result->f, 0, context->now.mpfr_round);
    mpfr_init2(digit, prec);
    for (i = 5 + precilen; i<len; i++) {
        mpfr_set_ui(digit, cp[i], context->now.mpfr_round);
        mpfr_div_2ui(digit, digit, (unsigned long)((i-4-precilen) * 8),
                     context->now.mpfr_round);
        mpfr_add(result->f, result->f, digit, context->now.mpfr_round);
    }
    mpfr_clear(digit);
    /* apply exponent, with its appropriate sign */
    if (exposign)
        mpfr_div_2ui(result->f, result->f, 8*expomag, context->now.mpfr_round);
    else
        mpfr_mul_2ui(result->f, result->f, 8*expomag, context->now.mpfr_round);
    /* apply significand-sign (sign of the overall number) */
    if (resusign)
        mpfr_neg(result->f, result->f, context->now.mpfr_round);

    return (PyObject*)result;
}
#endif


#if 0
/* Format of binary representation of an mpfr.
 *
 * bytes[0..3]         "mpfr"
 * byte[4.0]           if set, mantissa is negative
 * byte[4.1]           if set, exponent is negative
 * byte[4.2]           if set, value is 0, see [20.0] for sign
 * byte[4.3]           if set, value is Infinity, see [20.0] for sign
 * byte[4.4]           if set, value is NaN
 * byte[4.5...7]       # of bytes - 1 in precision and exponent
 * bytes[5...]         precision
 * bytes[5+plen...]    absolute value of exponent
 * bytes[5+2*plen...]  mantissa
 */

static PympfrObject *
Pympfr_From_Binary(PyObject *s)
{
    PympfrObject *newob;
    unsigned char *cp;
    mpfr_prec_t prec;
    mpfr_exp_t expt;
    Py_ssize_t i, len;
    PyObject *ascii_str = NULL;
    int codebyte, msgn, esgn, plen;

    /* Accept either a sequence of bytes or Unicode string. */
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

    /* Verify this is really an mpfr binary string. */
    if (!(strcmp("mpfr", cp)) {
        VALUE_ERROR("invalid byte sequence");
        Py_XDECREF(ascii_str);
        return NULL;
    }

    if (!(newob = Pympfr_new(0))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    codebyte = cp[4];
    msgn = codebyte & 1;
    esgn = codebyte & 2;
    plen = ((codebyte & 224) >> 5) + 1;

    if (codebyte & 4) {
        if (msgn)
            mpfr_set_zero(Pympfr_AS_MPFR(newob), -1);
        else
            mpfr_set_zero(Pympfr_AS_MPFR(newob), 1);
        Py_XDECREF(ascii_str);
        return newob;
    }

    if (codebyte & 8) {
        if (msgn)
            mpfr_set_inf(Pympfr_AS_MPFR(newob), -1);
        else
            mpfr_set_inf(Pympfr_AS_MPFR(newob), 1);
        Py_XDECREF(ascii_str);
        return newob;
    }

    if (codebyte & 16) {
        mpfr_set_nan(Pympfr_AS_MPFR(newob));
        Py_XDECREF(ascii_str);
        return newob;
    }

    /* Check for correct string length. */
    if (len < 22) {
        VALUE_ERROR("string too short to be a gmpy2.mpfr binary encoding");
        Py_DECREF((PyObject*)newob);
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Read the next 8 bytes for the precision. */
    prec = 0;
    for (i = 5; i < 13; ++i) {
        prec = (prec << 8) | cp[i];
    }

    /* Read the next 8 bytes for the exponent. */
    expt = 0;
    for (i = 13; i < 21; ++i) {
        expt = (expt << 8) | cp[i];
    }

    /* Set the precision. */
    if (prec > MPFR_PREC_MIN && prec < MPFR_PREC_MAX)
        mpfr_set_prec(Pympfr_AS_MPFR(newob), prec);
    else {
        VALUE_ERROR("value for precision is too large");
        Py_DECREF((PyObject*)newob);
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* reconstruct 'mantissa' (significand) */
    mpfr_set_si(newob->f, 0, context.mpfr_round);
    mpfr_init2(digit, prec);
    for (i = 5 + precilen; i<len; i++) {
        mpfr_set_ui(digit, cp[i], context.mpfr_round);
        mpfr_div_2ui(digit, digit, (unsigned long)((i-4-precilen) * 8),
                     context.mpfr_round);
        mpfr_add(newob->f, newob->f, digit, context.mpfr_round);
    }
    mpfr_clear(digit);
    /* apply exponent, with its appropriate sign */
    if (exposign)
        mpfr_div_2ui(newob->f, newob->f, 8*expomag, context.mpfr_round);
    else
        mpfr_mul_2ui(newob->f, newob->f, 8*expomag, context.mpfr_round);
    /* apply significand-sign (sign of the overall number) */
    if (resusign)
        mpfr_neg(newob->f, newob->f, context.mpfr_round);

    Py_XDECREF(ascii_str);
    return newob;
}
#endif
