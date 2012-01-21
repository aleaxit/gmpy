/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_conv.c                                                         *
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

static int isFraction(PyObject* obj)
{
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isDecimal(PyObject* obj)
{
    if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) return 1;

    return 0;
}

static PyxmpzObject *
Pyxmpz2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pyxmpz_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PympzObject *
Pyxmpz2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(Pyxmpz_Check(self));
    if ((newob = Pympz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PyxmpzObject *
Pympz2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(Pympz_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

#ifdef PY2
static PympzObject *
PyInt2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pympz_new()))
        mpz_set_si(newob->z, PyInt_AS_LONG(self));
    return newob;
}

static PyxmpzObject *
PyInt2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(PyInt_Check(self));
    if ((newob = Pyxmpz_new()))
        mpz_set_si(newob->z, PyInt_AsLong(self));
    return newob;
}
#endif

static PympzObject *
PyFloat2Pympz(PyObject *self)
{
    PympzObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pympz_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(newob->z, d);
    }
    return newob;
}

static PyxmpzObject *
PyFloat2Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    assert(PyFloat_Check(self));
    if ((newob = Pyxmpz_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(newob->z, d);
    }
    return newob;
}

/* For fast conversion between PyLong and mpz, we use code located in
 * mpz_pylong.c.
 */
static PympzObject *
PyLong2Pympz(PyObject * obj)
{
    PympzObject *newob;
    if (!(newob = Pympz_new()))
        return NULL;
    mpz_set_PyLong(Pympz_AS_MPZ(newob), obj);
    return newob;
}

static PyxmpzObject *
PyLong2Pyxmpz(PyObject * obj)
{
    PyxmpzObject *newob;
    if (!(newob = Pyxmpz_new()))
        return NULL;
    mpz_set_PyLong(Pyxmpz_AS_MPZ(newob), obj);
    return newob;
}

/* mpz conversion from string includes from-binary (base-256 LSB string
 * of bytes) and 'true' from-string (bases 2 to 62; bases 8 and 16 are
 * special -- decorations of leading 0/0x are allowed (not required).
 *
 * Binary form was previously (0.6) limited to >=0 values; now (0.7)
 * extended to <0 values, by adding one last sign-byte, 0xFF for <0,
 * 0x00 for >0 (the latter only if the #bits is exact multiple of 8).
 */

/* mpz_set_PyStr returns -1 on error, 1 if successful. */

static int
mpz_set_PyStr(mpz_ptr z, PyObject *s, long base)
{
    unsigned char *cp;
    Py_ssize_t len;
    int i;
    PyObject *ascii_str = NULL;

    assert(PyStrOrUnicode_Check(s));

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = (unsigned char*)PyBytes_AsString(s);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return -1;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }

    /* Don't allow NULL characters */
    for (i=0; i<len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            Py_XDECREF(ascii_str);
            return -1;
        }
    }
    /* delegate rest to GMP's _set_str function */
    if (base==0) {
        if (cp[0]=='0') {
            if (cp[1]=='b') {
                base = 2;
                cp+=2;
            }
            else if (cp[1]=='o') {
                base = 8;
                cp+=2;
            }
            else if (cp[1]=='x') {
                base = 16;
                cp+=2;
            }
            else {
                base = 10;
            }
        }
        else {
            base = 10;
        }
    }
    if (-1 == mpz_set_str(z, (char*)cp, base)) {
        VALUE_ERROR("invalid digits");
        Py_XDECREF(ascii_str);
        return -1;
    }
    Py_XDECREF(ascii_str);
    return 1;
}

static PympzObject *
PyStr2Pympz(PyObject *s, long base)
{
    PympzObject *newob;

    assert(PyStrOrUnicode_Check(s));
    if (!(newob = Pympz_new()))
        return NULL;
    if (mpz_set_PyStr(newob->z, s, base) == -1) {
        Py_DECREF((PyObject*)newob);
        return NULL;
    }
    return newob;
}

static PyxmpzObject *
PyStr2Pyxmpz(PyObject *s, long base)
{
    PyxmpzObject *newob;

    assert(PyStrOrUnicode_Check(s));
    if (!(newob = Pyxmpz_new()))
        return NULL;
    if (mpz_set_PyStr(newob->z, s, base) == -1) {
        Py_DECREF((PyObject*)newob);
        return NULL;
    }
    return newob;
}

/* For fast mpz to PyLong conversion, we use code located in mpz_pylong.
 */
static PyObject *
Pympz2PyLong(PympzObject *self)
{
    return mpz_get_PyLong(Pympz_AS_MPZ(self));
}

static PyObject *
Pyxmpz2PyLong(PyxmpzObject *self)
{
    return mpz_get_PyLong(Pyxmpz_AS_MPZ(self));
}

static PyObject *
Pympz_To_Integer(PympzObject *self)
{
#ifdef PY3
    return Pympz2PyLong(self);
#else
    if (mpz_fits_slong_p(self->z))
        return PyInt_FromLong(mpz_get_si(self->z));
    else
        return Pympz2PyLong(self);
#endif
}

static PyObject *
Pyxmpz_To_Integer(PyxmpzObject *self)
{
#ifdef PY3
    return Pyxmpz2PyLong(self);
#else
    if (mpz_fits_slong_p(self->z))
        return PyInt_FromLong(mpz_get_si(self->z));
    else
        return Pyxmpz2PyLong(self);
#endif
}

static PyObject *
Pympz2PyFloat(PympzObject *self)
{
    double res = mpz_get_d(self->z);

    return PyFloat_FromDouble(res);
}

/*
 *  build binary representation of mpz (base-256 little-endian)
 *  Note: design limitation used to forbid binary repr of <0 mpz;
 *  this has now been remedied, but at the price of full compatibility
 *  with files saved in gmpy releases 0.6 and earlier.
 */
static PyObject *
mpz2binary(mpz_ptr z)
{
    size_t size, usize;
    int negative, needtrail;
    char *buffer;
    PyObject *s;

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z); /* Change the sign temporarily! */
    }
    else {
        negative = 0;
    }

    size = mpz_sizeinbase(z, 2);
    needtrail = (size%8)==0;
    usize = size = (size + 7) / 8;
    if (negative || needtrail)
        ++size;

    TEMP_ALLOC(buffer, size);
    buffer[0] = 0x00;
    mpz_export(buffer, NULL, -1, sizeof(char), 0, 0, z);
    if (usize < size) {
        buffer[usize] = negative?0xff:0x00;
    }
    if (negative) {
        mpz_neg(z, z);
    }
    s = PyBytes_FromStringAndSize(buffer, size);
    TEMP_FREE(buffer, size);
    return s;
}

static PyObject *
Pympz2binary(PympzObject *self)
{
    return mpz2binary(self->z);
}

static PyObject *
Pyxmpz2binary(PyxmpzObject *self)
{
    return mpz2binary(self->z);
}

/* Format an mpz into any base (2 to 62). Bits in the option parameter
 * control various behaviors:
 *   bit 0: if set, output is wrapped with 'mpz(...)'
 *   bit 1: if set, a '+' is included for positive numbers
 *   bit 2: if set, a ' ' is included for positive nubmers
 *   bit 3: if set, a '0b', '0o', or '0x' is included for binary, octal, hex
 *   bit 4: if set, no prefix is included for binary, octal, hex
 *
 * Note: if neither bit 3 or 4 is set, prefixes that match the platform default
 * are included.
 *
 * If base < 0, capital letters are used
 */
static char* ztag = "mpz(";
static PyObject *
mpz_ascii(mpz_t z, int base, int option)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if (!((base == 0) || ((base >= -36) && (base <= -2)) ||
        ((base >= 2) && (base <= 62)) )) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 62");
        return NULL;
    }

    /* Allocate extra space for:
     *
     * minus sign and trailing NULL byte (2)
     * 'mpz()' tag                       (5)
     * '0x' prefix                       (2)
     * 'L' suffix                        (1)
     *                                  -----
     *                                   10
     */
    size = mpz_sizeinbase(z, base) + 11;
    TEMP_ALLOC(buffer, size);

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z);
    }

    p = buffer;
    if (option & 1) {
       strcpy(p, ztag);
       p += strlen(p);
    }

    if (negative)
        *(p++) = '-';
    else if (option & 2)
        *(p++) = '+';
    else if (option & 4)
        *(p++) = ' ';

    if (option & 8) {
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }
    else if (!(option & 24)) {
    #ifdef PY2
        if (base == 8) {
            *(p++) = '0';
        }
    #else
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
    #endif
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }

    mpz_get_str(p, base, z);     /* Doesn't return number of characters */
    p = buffer + strlen(buffer); /* Address of NULL byte */
#ifdef PY2
    if ((option & 1) && !mpz_fits_slong_p(z))
        *(p++) = 'L';
#endif
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';

    result = Py_BuildValue("s", buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}

/* Format an xmpz into any base (2 to 62). If with_tag != 0, the the output
 * is wrapped with 'xmpz(...)'. If with_sign > 0, a plus or minus leading
 * sign is always included. If with_sign < 0, a space is included instead of
 * a plus sign.
 */
static char* xztag = "xmpz(";
static PyObject *
xmpz_ascii(mpz_t z, int base, int option)
{
    PyObject *result;
    char *buffer, *p;
    int negative = 0;
    size_t size;

    if (!((base == 0) || ((base >= -36) && (base <= -2)) ||
        ((base >= 2) && (base <= 62)) )) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 62");
        return NULL;
    }

    /* Allocate extra space for:
     *
     * minus sign and trailing NULL byte (2)
     * 'gmpy2.xmpz()' tag                (6)
     * '0x' prefix                       (2)
     * 'L' suffix                        (1)
     *                                  -----
     *                                   11
     */
    size = mpz_sizeinbase(z, base) + 12;
    TEMP_ALLOC(buffer, size);

    if (mpz_sgn(z) < 0) {
        negative = 1;
        mpz_neg(z, z);
    }

    p = buffer;
    if (option & 1) {
       strcpy(p, xztag);
       p += strlen(p);
    }

    if (negative)
        *(p++) = '-';
    else if (option & 2)
        *(p++) = '+';
    else if (option & 4)
        *(p++) = ' ';

    if (option & 8) {
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }
    else if (!(option & 24)) {
    #ifdef PY2
        if (base == 8) {
            *(p++) = '0';
        }
    #else
        if (base == 2) {
            *(p++) = '0';
            *(p++) = 'b';
        }
        else if (base == 8) {
            *(p++) = '0';
            *(p++) = 'o';
        }
    #endif
        else if (base == 16) {
            *(p++) = '0';
            *(p++) = 'x';
        }
        else if (base == -16) {
            *(p++) = '0';
            *(p++) = 'X';
        }
    }

    mpz_get_str(p, base, z);     /* Doesn't return number of characters */
    p = buffer + strlen(buffer); /* Address of NULL byte */
#ifdef PY2
    if ((option & 1) && !mpz_fits_slong_p(z))
        *(p++) = 'L';
#endif
    if (option & 1)
        *(p++) = ')';
    *(p++) = '\00';

    result = Py_BuildValue("s", buffer);
    if (negative == 1) {
        mpz_neg(z, z);
    }
    TEMP_FREE(buffer, size);
    return result;
}

static PyObject *
Pympz_ascii(PympzObject *self, int base, int option)
{
    return mpz_ascii(self->z, base, option);
}

static PyObject *
Pyxmpz_ascii(PyxmpzObject *self, int base, int option)
{
    return xmpz_ascii(self->z, base, option);
}

static int isInteger(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pyxmpz_Check(obj))      return 1;

    return 0;
}

/* Number conversion routines
 *
 * The routines anynum2mpX will attempt to convert any number-like object into
 * into a gmpy object. These routines are intended for construction of mpXs.
 * The accepted number-like objects are:
 *      1) int (Python 2.x)
 *      2) long (Python 2.x and 3.x)
 *      3) float
 *      4) Decimal
 *      5) Fraction
 *      6) other gmpy objects
 *
 * The routine Pympz_From_Integer will only convert integer-like objects into to a
 * gmpy mpz. The accepted integer-like objects are:
 *      1) int
 *      2) long
 *      3) mpz
 *      4) xmpz
 *
 * The routine Pympq_From_Rational will convert an integer- and rational-like
 * object into a gmpy mpq. The accepted objects are:
 *      1) int
 *      2) long
 *      3) Fraction
 *      4) mpz
 *      5) mpq
 *      6) xmpz
 */

/* NOTE: Pympq_From_Decimal returns an invalid mpq object when attempting to
 *       convert a NaN or inifinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 */

static PympzObject*
anynum2Pympz(PyObject* obj)
{
    PympzObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympz_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympzObject *) obj;
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympz(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympz(obj);
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pympz(obj);
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pympz(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympz(obj);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = PyLong2Pympz(s);
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pympz((PyObject *)temp);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }

    return newob;
}

static PyxmpzObject*
anynum2Pyxmpz(PyObject* obj)
{
    PyxmpzObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympz_Check(obj)) {
        newob = Pympz2Pyxmpz(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pyxmpz(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pyxmpz(obj);
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pyxmpz(obj);
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr2Pyxmpz(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pyxmpz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pyxmpz(obj);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = PyLong2Pyxmpz(s);
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pyxmpz((PyObject *)temp);
            Py_DECREF(s);
            Py_DECREF((PyObject*)temp);
        }
    }

    return newob;
}

/* Convert an Integer-like object (as determined by isInteger) to
 * a Pympz. Returns NULL and raises a TypeError if obj is not an
 * Integer-like object.
 */

static PympzObject*
Pympz_From_Integer(PyObject* obj)
{
    PympzObject* newob = 0;

    if (Pympz_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympzObject*) obj;
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympz(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympz(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympz(obj);
    }
    if (!newob)
        TYPE_ERROR("conversion error in Pympz_From_Integer");
    return newob;
}

/* Convert an Integer-like object (as determined by isInteger) to a 
 * C long or C unsigned long.
 */

static long
clong_From_Integer(PyObject *obj)
{
    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in clong_From_Integer");
            return -1;
        }
    }
    TYPE_ERROR("conversion error in clong_From_Integer");
    return -1;
}

/*
 * Convert an Integer-like object (as determined by isInteger) to
 * a gmp_si. On all platforms except 64-bit Windows, gmp_si is the same
 * as a C long. Returns -1 and raises OverflowError if the the number is
 * too large. Returns -1 and raises TypeError if obj was not an
 * Integer-like object.
 */

#ifndef _WIN64

/* Working with C long. */

static gmp_si
gmp_si_From_Integer(PyObject *obj)
{
    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in gmp_si_From_Integer");
            return -1;
        }
    }
    TYPE_ERROR("conversion error in gmp_si_From_Integer");
    return -1;
}

/* Working with C unsigned long. */

static gmp_ui
gmp_ui_From_Integer(PyObject *obj)
{
    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsUnsignedLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_ulong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_ui(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in gmp_ui_From_Integer");
            return (gmp_ui)-1;
        }
    }
    TYPE_ERROR("conversion error in gmp_ui_From_Integer");
    return (gmp_ui)-1;
}

#define MP_BITCNT_FROM_INTEGER(obj) gmp_ui_From_Integer(obj)
#define PyLong_AsGmp_siAndOverflow(a,b) PyLong_AsLongAndOverflow(a,b)

#else

/* Working with C long long. Can also assume using > MPIR 2.5. */

static gmp_si
gmp_si_From_Integer(PyObject *obj)
{
    PyObject* temp;

    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsLongLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_si(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in gmp_si_From_Integer");
            return -1;
        }
    }
    TYPE_ERROR("conversion error in gmp_si_From_Integer");
    return -1;
}

static gmp_ui
gmp_ui_From_Integer(PyObject *obj)
{
    PyObject* temp;

    if (PyIntOrLong_Check(obj)) {
        return PyLong_AsUnsignedLongLong(obj);
    }
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_ui(Pympz_AS_MPZ(obj))) {
            return mpz_get_ui(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in gmp_ui_From_Integer");
            return (gmp_ui)-1;
        }
    }
    TYPE_ERROR("conversion error in gmp_ui_From_Integer");
    return -1;
}

#define MP_BITCNT_FROM_INTEGER(obj) gmp_ui_From_Integer(obj)
#define PyLong_AsGmp_siAndOverflow(a,b) PyLong_AsLongLongAndOverflow(a,b)

#endif

/*
 * Convert an Integer-like object (as determined by isInteger) to
 * a Py_ssize_t. Returns -1 and raises OverflowError if the the number is
 * too large. Returns -1 and raises TypeError if obj was not an
 * Integer-like object.
 */

static Py_ssize_t
ssize_t_From_Integer(PyObject *obj)
{
    Py_ssize_t val;
    PyObject* temp;

    if (PyLong_Check(obj)) {
        return PyLong_AsSsize_t(obj);
    }
#ifdef PY2
    else if (PyInt_Check(obj)) {
        return PyInt_AsSsize_t(obj);
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return (Py_ssize_t)mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            /* This section should only be called on Win64. */
            temp = mpz_get_PyLong(Pympz_AS_MPZ(obj));
            if (!temp) {
                TYPE_ERROR("conversion error in ssize_t_From_Integer");
                return -1;
            }
            else {
                val = PyLong_AsSsize_t(temp);
                Py_DECREF(temp);
                return val;
            }
        }
    }
    TYPE_ERROR("conversion error in ssize_t_From_Integer");
    return -1;
}

/*
 * coerce any number to a mpz
 */
int
Pympz_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympzObject* newob = Pympz_From_Integer(arg);
    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpz'");
        return 0;
    }
}

/* str and repr implementations for mpz */
static PyObject *
Pympz2str(PympzObject *self)
{
    /* base-10, no tag */
    return Pympz_ascii(self, 10, 0);
}

static PyObject *
Pympz2repr(PympzObject *self)
{
    /* base-10, with tag */
    return Pympz_ascii(self, 10, 1);
}

/* str and repr implementations for xmpz */
static PyObject *
Pyxmpz2str(PyxmpzObject *self)
{
    /* base-10, no tag */
    return Pyxmpz_ascii(self, 10, 0);
}

static PyObject *
Pyxmpz2repr(PyxmpzObject *self)
{
    /* base-10, with tag */
    return Pyxmpz_ascii(self, 10, 1);
}




