/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq_conv.c                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen             *
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

#ifdef PY2
static PympqObject *
PyInt2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympq_new()))
        mpq_set_si(newob->q, PyInt_AsLong(self), 1);
    return newob;
}
#endif

static PympqObject *
Pympz2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pympq_new()))
        mpq_set_z(newob->q, Pympz_AS_MPZ(self));
    return newob;
}

static PympqObject *
Pyxmpz2Pympq(PyObject * obj)
{
    PympqObject *newob;

    assert(Pyxmpz_Check(obj));
    if ((newob = Pympq_new()))
        mpq_set_z(newob->q, Pyxmpz_AS_MPZ(obj));
    return newob;
}

static PympzObject *
Pympq2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pympq_Check(self));
    if ((newob = Pympz_new()))
        mpz_set_q(newob->z, Pympq_AS_MPQ(self));
    return newob;
}

static PyxmpzObject *
Pympq2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pympq_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set_q(newob->z, Pympq_AS_MPQ(self));
    return newob;
}

static PympqObject *
PyLong2Pympq(PyObject *self)
{
    PympqObject *newob;
    PyObject *temp = (PyObject*)PyLong2Pympz(self);

    if (!temp)
        return NULL;
    newob = Pympz2Pympq(temp);
    Py_DECREF(temp);
    return newob;
}

static PympqObject *
PyFloat2Pympq(PyObject *self)
{
    PympqObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pympq_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpq' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpq' does not support Infinity");
            return NULL;
        }
        mpq_set_d(newob->q, d);
    }
    return newob;
}

/*
 * mpq conversion from string includes from-binary (base-256 LSB string
 * of bytes) and 'true' from-string (bases 2 to 62; bases 8 and 16 are
 * special -- decorations of leading 0/0x are allowed (but not required).
 * For 'true-bases' 2..62 there is a '/' separator between numerator and
 * denominator (if none, just numerator!); decimal point NOT allowed.
 *
 * Added in gmpy 1.02: also support a string of the form '12.34', i.e.,
 * WITH a decimal point and WITHOUT a slash
 *
 * Binary-form: 4-byte numerator length (upper bit set if <0), then
 * numerator (as above for mpz), then denominator (ditto).
 */
static PympqObject *
PyStr2Pympq(PyObject *stringarg, long base)
{
    PympqObject *newob;
    unsigned char *cp;
    Py_ssize_t len;
    int i;
    PyObject *ascii_str = NULL;

    assert(PyStrOrUnicode_Check(stringarg));
    if (!(newob = Pympq_new()))
        return NULL;
    if (PyBytes_Check(stringarg)) {
        len = PyBytes_Size(stringarg);
        cp = (unsigned char*)PyBytes_AsString(stringarg);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(stringarg);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            Py_DECREF((PyObject*)newob);
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }

    /* Don't allow NULL characters */
    for (i=0; i<len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }
    }
    /* trickily delegate the rest to GMP avoiding allocations/copies */
    {
        char *whereslash = strchr((char*)cp, '/');
        char *wheredot = strchr((char*)cp, '.');
        if (whereslash && wheredot) {
            VALUE_ERROR("illegal string: both . and / found");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }

        if (wheredot) {
            char *counter;
            unsigned long digits = 0;
            if (base != 10) {
                VALUE_ERROR("illegal string: embedded . requires base=10");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }

            counter = wheredot;
            digits = 0;
            *wheredot = ' ';
            while (*++counter != '\0') {
                if (isdigit(*counter))
                    digits++;
            }
            if (-1 == mpz_set_str(mpq_numref(newob->q), (char*)cp, base)) {
                if (wheredot)
                    *wheredot = '.';
                VALUE_ERROR("invalid digits");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }
            mpz_ui_pow_ui(mpq_denref(newob->q), 10, digits);
            mpq_canonicalize(newob->q);
            *wheredot = '.';
            return (PympqObject*)newob;
        }

        if (whereslash)
            *whereslash = 0;
        if (-1 == mpz_set_str(mpq_numref(newob->q), (char*)cp, base)) {
            if (whereslash)
                *whereslash = '/';
            VALUE_ERROR("invalid digits");
            Py_DECREF((PyObject*)newob);
            Py_XDECREF(ascii_str);
            return NULL;
        }
        if (whereslash) {
            *whereslash = '/';
            if (-1==mpz_set_str(mpq_denref(newob->q), whereslash+1, base)) {
                VALUE_ERROR("invalid digits");
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                return NULL;
            }
            if (0==mpz_sgn(mpq_denref(newob->q))) {
                Py_DECREF((PyObject*)newob);
                Py_XDECREF(ascii_str);
                ZERO_ERROR("zero denominator in 'mpq'");
                return NULL;
            }
            mpq_canonicalize(newob->q);
        }
        else {
            mpz_set_ui(mpq_denref (newob->q), 1);
        }
    }
    Py_XDECREF(ascii_str);
    return newob;
}

static PyObject *
Pympq2PyLong(PympqObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympq2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz2PyLong(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}

#ifdef PY2
static PyObject *
Pympq2PyInt(PympqObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympq2Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz_To_Integer(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}
#endif

static PyObject *
Pympq2PyFloat(PympqObject *self)
{
    double res = mpq_get_d(self->q);

    return PyFloat_FromDouble(res);
}

/*
 *  build binary representation of mpq (base-256 little-endian
 *  for num, then den; before those, 4 bytes for _length_ of
 *  numerator, which also encode sign as the single top bit).
 */

static PyObject *
Pympq2binary(PympqObject *self)
{
    size_t sizenum, sizeden, size, sizetemp;
    int negative = 0;
    char *buffer;
    int i;
    PyObject *result;

    if (mpq_sgn(self->q) < 0) {
        negative = 1;
        mpz_abs(mpq_numref(self->q), mpq_numref(self->q));
    }
    else {
        negative = 0;
    }
    assert(mpz_sgn(mpq_denref(self->q)) > 0);

    sizenum = (mpz_sizeinbase(mpq_numref(self->q), 2) + 7) / 8;
    sizeden = (mpz_sizeinbase(mpq_denref(self->q), 2) + 7) / 8;
    size = sizenum + sizeden + 4;

    TEMP_ALLOC(buffer, size);

    sizetemp = sizenum;
    for (i=0; i<4; i++) {
        buffer[i] = (char)(sizetemp & 0xff);
        sizetemp >>= 8;
    }
    if (negative)
        buffer[3] |= 0x80;
    buffer[4] = 0x00;

    mpz_export(buffer+4, NULL, -1, sizeof(char), 0, 0, mpq_numref(self->q));
    mpz_export(buffer+sizenum+4, NULL, -1, sizeof(char), 0, 0, mpq_denref(self->q));
    if (negative) {
        mpz_neg( mpq_numref(self->q), mpq_numref(self->q));
    }
    result = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return result;
}

static int qden_1(mpq_t q)
{
    return 0 == mpz_cmp_ui(mpq_denref(q),1);
}

static PyObject *
Pympq_ascii(PympqObject *self, int base, int option)
{
    PyObject *result = 0, *numstr = 0, *denstr = 0;
    char buffer[50], *p;

    numstr = mpz_ascii(mpq_numref(self->q), base, 0);
    if (!numstr)
        return NULL;

    /* Check if denominator is 1 and no tag is requested. If so, just
     * return the numerator.
     */
    if (!(option & 1) && qden_1(self->q))
        return numstr;

    denstr = mpz_ascii(mpq_denref(self->q), base, 0);
    if (!denstr) {
        Py_DECREF(numstr);
        return NULL;
    }

    /* Build the format string. */
    p = buffer;
    if (option & 1) {
        *(p++) = 'm';
        *(p++) = 'p';
        *(p++) = 'q';
        *(p++) = '(';
    }
#ifdef PY2
    *(p++) = '%';
    *(p++) = 's';
    if (!mpz_fits_slong_p(mpq_numref(self->q)))
        *(p++) = 'L';
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 's';
    if (!mpz_fits_slong_p(mpq_denref(self->q)))
        *(p++) = 'L';
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';
    result = PyString_FromFormat(buffer, PyString_AS_STRING(numstr),
                                 PyString_AS_STRING(denstr));
#else
    *(p++) = '%';
    *(p++) = 'U';
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 'U';
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';
    result = PyUnicode_FromFormat(buffer, numstr, denstr);
#endif
    Py_DECREF(numstr);
    Py_DECREF(denstr);
    return result;
}

static int isRational(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

/* NOTE: Pympq_From_Decimal returns an invalid mpq object when attempting to
 *       convert a NaN or inifinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 */

static PympqObject*
Pympq_From_Decimal(PyObject* obj)
{
    PympqObject *result;
    PyObject *d_exp, *d_int, *d_sign, *d_is_special;
    gmp_si exp;
    mpz_t temp;
    const char *string;

    if (!(result = Pympq_new()))
        return NULL;
    mpq_set_si(result->q, 0, 1);

    d_exp = PyObject_GetAttrString(obj, "_exp");
    d_int = PyObject_GetAttrString(obj, "_int");
    d_sign = PyObject_GetAttrString(obj, "_sign");
    d_is_special = PyObject_GetAttrString(obj, "_is_special");
    if (!d_exp || !d_int || !d_sign || !d_is_special) {
        SYSTEM_ERROR("Object does not appear to be Decimal");
        Py_XDECREF(d_exp);
        Py_XDECREF(d_int);
        Py_XDECREF(d_sign);
        Py_XDECREF(d_is_special);
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    if (PyObject_IsTrue(d_is_special)) {
        string = Py2or3String_AsString(d_exp);
        if (string[0] == 'N' || string[0] == 'n') {
            mpz_set_si(mpq_denref(result->q), 0);
            Py_DECREF(d_exp);
            Py_DECREF(d_int);
            Py_DECREF(d_sign);
            Py_DECREF(d_is_special);
            return result;
        }
        if (string[0] == 'F') {
            if (PyObject_IsTrue(d_sign))
                mpq_set_si(result->q, -1, 0);
            else
                mpq_set_si(result->q, 1, 0);
            Py_DECREF(d_exp);
            Py_DECREF(d_int);
            Py_DECREF(d_sign);
            Py_DECREF(d_is_special);
            return result;
        }
        SYSTEM_ERROR("Cannot convert Decimal to mpq");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_set_PyStr(mpq_numref(result->q), d_int, 10);

    if (PyObject_IsTrue(d_sign))
        mpz_mul_si(mpq_numref(result->q), mpq_numref(result->q), -1);

    exp = PyIntOrLong_AsGmp_si(d_exp);
    if (exp == -1 && PyErr_Occurred()) {
        SYSTEM_ERROR("Decimal _exp is not valid or overflow occurred");
        Py_DECREF((PyObject*)result);
        Py_DECREF(d_exp);
        Py_DECREF(d_int);
        Py_DECREF(d_sign);
        Py_DECREF(d_is_special);
        return NULL;
    }

    mpz_inoc(temp);
    if (exp <= 0)
        mpz_ui_pow_ui(mpq_denref(result->q), 10, (gmp_ui)(-exp));
    else {
        mpz_inoc(temp);
        mpz_ui_pow_ui(temp, 10, (gmp_si)(exp));
        mpz_mul(mpq_numref(result->q), mpq_numref(result->q), temp);
        mpz_cloc(temp);
    }

    mpq_canonicalize(result->q);
    Py_DECREF(d_exp);
    Py_DECREF(d_int);
    Py_DECREF(d_sign);
    Py_DECREF(d_is_special);

    return result;
}

static PympqObject*
Pympq_From_Real(PyObject* obj)
{
    PympqObject* newob = 0;

    if (Pympq_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympqObject *) obj;
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympq(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympq(obj);
#endif
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pympq(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympq(obj);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympq(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympq(obj);
    }
    else if (isDecimal(obj)) {
        if ((newob = Pympq_From_Decimal(obj))) {
            if (!mpz_cmp_si(mpq_denref(newob->q), 0)) {
                if (mpz_get_si(mpq_numref(newob->q)) == 0)
                    VALUE_ERROR("'mpq' does not support NaN");
                else
                    VALUE_ERROR("'mpq' does not support Infinity");
                Py_DECREF((PyObject*)newob);
                newob = NULL;
            }
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }

    return newob;
}

/* Convert an integer or mpz to mpq. */

static PympqObject*
Pympq_From_Rational(PyObject* obj)
{
    PympqObject* newob = 0;

    if (Pympq_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympqObject *) obj;
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympq(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympq(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympq(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympq(obj);
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympq(s, 10);
            Py_DECREF(s);
        }
    }

    return newob;
}

/*
 * coerce any number to a mpq
 */

int
Pympq_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympqObject* newob = Pympq_From_Rational(arg);
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        if (!PyErr_Occurred()) {
            TYPE_ERROR("argument can not be converted to 'mpq'");
        }
        return 0;
    }
}

/* str and repr implementations for mpq */
static PyObject *
Pympq2str(PympqObject *self)
{
    /* base-10, no tag */
    return Pympq_ascii(self, 10, 0);
}

static PyObject *
Pympq2repr(PympqObject *self)
{
    /* base-10, with tag */
    return Pympq_ascii(self, 10, 1);
}





