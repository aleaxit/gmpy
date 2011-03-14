/*
 * The contents of this file are currently not used. This file will be used
 * to support gmpy2.to_binary() and gmpy2.from_binary(). When those two
 * functions are implemented, binary support will be removed from the object
 * constructors.
 */


/* Conversion routines between GMPY2 objects and a compact, portable
 * binary representation. The binary format of GMPY2 is not compatible
 * with GMPY 1.x.
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
