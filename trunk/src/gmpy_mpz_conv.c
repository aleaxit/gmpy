/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_conv.c                                                         *
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

static int isInteger(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pyxmpz_Check(obj))      return 1;

    return 0;
}

static int isFraction(PyObject* obj)
{
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isDecimal(PyObject* obj)
{
#if PY_VERSION_HEX < 0x03030000
    if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) return 1;
#else
    if (!strcmp(Py_TYPE(obj)->tp_name, "decimal.Decimal")) return 1;
#endif

    return 0;
}

static PyxmpzObject *
Pyxmpz_From_Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    if ((newob = (PyxmpzObject*)Pyxmpz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PympzObject *
Pympz_From_Pyxmpz(PyObject *self)
{
    PympzObject *newob;

    if ((newob = (PympzObject*)Pympz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

static PyxmpzObject *
Pyxmpz_From_Pympz(PyObject *self)
{
    PyxmpzObject *newob;

    if ((newob = (PyxmpzObject*)Pyxmpz_new()))
        mpz_set(newob->z, Pyxmpz_AS_MPZ(self));
    return newob;
}

#ifdef PY2
static PympzObject *
Pympz_From_PyInt(PyObject *self)
{
    PympzObject *newob;

    if ((newob = (PympzObject*)Pympz_new()))
        mpz_set_si(newob->z, PyInt_AS_LONG(self));
    return newob;
}

static PyxmpzObject *
Pyxmpz_From_PyInt(PyObject *self)
{
    PyxmpzObject *newob;

    if ((newob = (PyxmpzObject*)Pyxmpz_new()))
        mpz_set_si(newob->z, PyInt_AsLong(self));
    return newob;
}
#endif

static PympzObject *
Pympz_From_PyFloat(PyObject *self)
{
    PympzObject *newob;

    if ((newob = (PympzObject*)Pympz_new())) {
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
Pyxmpz_From_PyFloat(PyObject *self)
{
    PyxmpzObject *newob;

    if ((newob = (PyxmpzObject*)Pyxmpz_new())) {
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
Pympz_From_PyLong(PyObject * obj)
{
    PympzObject *newob;
    if (!(newob = (PympzObject*)Pympz_new()))
        return NULL;
    mpz_set_PyLong(newob->z, obj);
    return newob;
}

static PyxmpzObject *
Pyxmpz_From_PyLong(PyObject * obj)
{
    PyxmpzObject *newob;
    if (!(newob = (PyxmpzObject*)Pyxmpz_new()))
        return NULL;
    mpz_set_PyLong(newob->z, obj);
    return newob;
}

/* mpz_set_PyStr returns -1 on error, 1 if successful. */

static int
mpz_set_PyStr(mpz_ptr z, PyObject *s, long base)
{
    unsigned char *cp;
    Py_ssize_t len;
    size_t i;
    PyObject *ascii_str = NULL;

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
Pympz_From_PyStr(PyObject *s, long base)
{
    PympzObject *newob;

    if (!(newob = (PympzObject*)Pympz_new()))
        return NULL;

    if (mpz_set_PyStr(newob->z, s, base) == -1) {
        Py_DECREF((PyObject*)newob);
        return NULL;
    }
    return newob;
}

static PyxmpzObject *
Pyxmpz_From_PyStr(PyObject *s, long base)
{
    PyxmpzObject *newob;

    if (!(newob = (PyxmpzObject*)Pyxmpz_new()))
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
Pympz_To_PyLong(PympzObject *self)
{
    return mpz_get_PyLong(Pympz_AS_MPZ(self));
}

static PyObject *
Pyxmpz_To_PyLong(PyxmpzObject *self)
{
    return mpz_get_PyLong(Pyxmpz_AS_MPZ(self));
}

/* The _To_PyIntOrLong functions should be used when converting a number back
 * to a Python value since is automatically returns an "int" or "long" when
 * using Python 2.x
 */

static PyObject *
Pympz_To_PyIntOrLong(PympzObject *self)
{
#ifdef PY3
    return Pympz_To_PyLong(self);
#else
    if (mpz_fits_slong_p(self->z))
        /* cast is safe since we know if first in a signed long */
        return PyInt_FromLong((long)mpz_get_si(self->z));
    else
        return Pympz_To_PyLong(self);
#endif
}

static PyObject *
Pyxmpz_To_PyIntOrLong(PyxmpzObject *self)
{
#ifdef PY3
    return Pyxmpz_To_PyLong(self);
#else
    if (mpz_fits_slong_p(self->z))
        /* cast is safe since we know if first in a signed long */
        return PyInt_FromLong((long)mpz_get_si(self->z));
    else
        return Pyxmpz_To_PyLong(self);
#endif
}

static PyObject *
Pympz_To_PyFloat(PympzObject *self)
{
    double res = mpz_get_d(self->z);

    return PyFloat_FromDouble(res);
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
        VALUE_ERROR("base must be in the interval 2 ... 62");
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
        VALUE_ERROR("base must be in the interval 2 ... 62");
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
Pympz_To_PyStr(PympzObject *self, int base, int option)
{
    return mpz_ascii(self->z, base, option);
}

static PyObject *
Pyxmpz_To_PyStr(PyxmpzObject *self, int base, int option)
{
    return xmpz_ascii(self->z, base, option);
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
Pympz_From_Number(PyObject* obj)
{
    PympzObject* newob = 0;
    PympqObject* temp = 0;

    if (Pympz_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympzObject *) obj;
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympz_From_PyInt(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = Pympz_From_PyLong(obj);
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq_To_Pympz(obj);
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr_To_Pympz(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = Pympz_From_PyFloat(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pympz_From_Pyxmpz(obj);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = Pympz_From_PyLong(s);
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = Pympq_From_PyStr(s, 10);
            if (temp) {
                newob = Pympq_To_Pympz((PyObject *)temp);
                Py_DECREF((PyObject*)temp);
            }
            Py_DECREF(s);
        }
    }

    return newob;
}

static PyxmpzObject*
Pyxmpz_From_Number(PyObject* obj)
{
    PyxmpzObject* newob = NULL;
    PympqObject* temp = NULL;

    if (Pympz_Check(obj)) {
        newob = Pyxmpz_From_Pympz(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pyxmpz_From_PyInt(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = Pyxmpz_From_PyLong(obj);
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq_To_Pyxmpz(obj);
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr_To_Pyxmpz(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = Pyxmpz_From_PyFloat(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz_From_Pyxmpz(obj);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = Pyxmpz_From_PyLong(s);
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = Pympq_From_PyStr(s, 10);
            if (temp) {
                newob = Pympq_To_Pyxmpz((PyObject *)temp);
                Py_DECREF((PyObject*)temp);
            }
            Py_DECREF(s);
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
        newob = Pympz_From_PyInt(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = Pympz_From_PyLong(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pympz_From_Pyxmpz(obj);
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
            return (long)mpz_get_si(Pympz_AS_MPZ(obj));
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
 * a mpir_si. On all platforms except 64-bit Windows, mpir_si is the same
 * as a C long. Returns -1 and raises OverflowError if the the number is
 * too large. Returns -1 and raises TypeError if obj was not an
 * Integer-like object.
 */

#ifndef _WIN64

/* Working with C long. */

static mpir_si
SI_From_Integer(PyObject *obj)
{
    if (PyLong_Check(obj)) {
        return PyLong_AsLong(obj);
    }
#ifdef PY2
    else if (PyInt_Check(obj)) {
        return PyInt_AsLong(obj);
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in SI_From_Integer");
            return -1;
        }
    }
    TYPE_ERROR("conversion error in SI_From_Integer");
    return -1;
}

/* Working with C unsigned long. */

static mpir_ui
UI_From_Integer(PyObject *obj)
{
    if (PyLong_Check(obj)) {
        return PyLong_AsUnsignedLong(obj);
    }
#ifdef PY2
    if (PyInt_Check(obj)) {
        long temp = PyInt_AsLong(obj);
        /* Create an OverflowError for negative values. */
        if (temp < 0) {
            OVERFLOW_ERROR("can't convert negative value to unsigned int");
            return (mpir_ui)-1;
        }
        return temp;
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_ulong_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_ui(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in UI_From_Integer");
            return (mpir_ui)-1;
        }
    }
    TYPE_ERROR("conversion error in UI_From_Integer");
    return (mpir_ui)-1;
}

#define MP_BITCNT_FROM_INTEGER(obj) UI_From_Integer(obj)
#define PyLong_AsSIAndOverflow(a,b) PyLong_AsLongAndOverflow(a,b)

#else

/* Working with C long long. Can also assume using > MPIR 2.5. */

static mpir_si
SI_From_Integer(PyObject *obj)
{
    if (PyLong_Check(obj)) {
        return PyLong_AsLongLong(obj);
    }
#ifdef PY2
    else if (PyInt_Check(obj)) {
        return (mpir_si)PyInt_AsLong(obj);
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_si_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_si(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in SI_From_Integer");
            return -1;
        }
    }
    TYPE_ERROR("conversion error in SI_From_Integer");
    return -1;
}

static mpir_ui
UI_From_Integer(PyObject *obj)
{
    if (PyLong_Check(obj)) {
        /* Returns an OverflowError for negative values. */
        return PyLong_AsUnsignedLongLong(obj);
    }
#ifdef PY2
    else if (PyInt_Check(obj)) {
        long temp = PyInt_AsLong(obj);
        /* Create an OverflowError for negative values. */
        if (temp < 0) {
            OVERFLOW_ERROR("can't convert negative value to unsigned int");
            return (mpir_ui)-1;
        }
        return (mpir_ui)temp;
    }
#endif
    else if (CHECK_MPZANY(obj)) {
        if (mpz_fits_ui_p(Pympz_AS_MPZ(obj))) {
            return mpz_get_ui(Pympz_AS_MPZ(obj));
        }
        else {
            OVERFLOW_ERROR("overflow in UI_From_Integer");
            return (mpir_ui)-1;
        }
    }
    TYPE_ERROR("conversion error in UI_From_Integer");
    return (mpir_ui)-1;
}

#define MP_BITCNT_FROM_INTEGER(obj) UI_From_Integer(obj)
#define PyLong_AsSIAndOverflow(a,b) PyLong_AsLongLongAndOverflow(a,b)

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

/* currently not in use */
#if 0
static int
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
#endif

/* str and repr implementations for mpz */
static PyObject *
Pympz_To_Str(PympzObject *self)
{
    /* base-10, no tag */
    return Pympz_To_PyStr(self, 10, 0);
}

static PyObject *
Pympz_To_Repr(PympzObject *self)
{
    /* base-10, with tag */
    return Pympz_To_PyStr(self, 10, 1);
}

/* str and repr implementations for xmpz */
static PyObject *
Pyxmpz_To_Str(PyxmpzObject *self)
{
    /* base-10, no tag */
    return Pyxmpz_To_PyStr(self, 10, 0);
}

static PyObject *
Pyxmpz_To_Repr(PyxmpzObject *self)
{
    /* base-10, with tag */
    return Pyxmpz_To_PyStr(self, 10, 1);
}




