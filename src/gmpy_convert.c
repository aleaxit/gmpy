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

static int isRational(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

#ifdef WITHMPFR
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
#endif

#ifdef WITHMPC
static int isComplex(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (Pympc_Check(obj))       return 1;
    if (PyFloat_Check(obj))     return 1;
    if (PyComplex_Check(obj))   return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}
#endif

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
            OVERFLOW_ERROR("'mpz' does not support Infinity");
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
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
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
    mpz_set_PyIntOrLong(newob->z, obj);
    return newob;
}

static PyxmpzObject *
Pyxmpz_From_PyLong(PyObject * obj)
{
    PyxmpzObject *newob;
    if (!(newob = (PyxmpzObject*)Pyxmpz_new()))
        return NULL;
    mpz_set_PyIntOrLong(newob->z, obj);
    return newob;
}

/* mpz_set_PyStr returns -1 on error, 1 if successful. */

static int
mpz_set_PyStr(mpz_ptr z, PyObject *s, int base)
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
Pympz_From_PyStr(PyObject *s, int base)
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
Pyxmpz_From_PyStr(PyObject *s, int base)
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
 * using Python 2.x. The _To_PyLong functions (above) should only be used
 * when a PyLong is specifically needed for Python 2.x.
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

    if (Py_IS_INFINITY(res)) {
        OVERFLOW_ERROR("'mpz' too large to convert to float");
        return NULL;
    }

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
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympq_To_Pympz((PyObject *)temp);
            Py_DECREF((PyObject*)temp);
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
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympq_To_Pyxmpz((PyObject *)temp);
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

#ifdef PY2
static PympqObject *
Pympq_From_PyInt(PyObject *self)
{
    PympqObject *newob;

    if ((newob = (PympqObject*)Pympq_new()))
        mpq_set_si(newob->q, PyInt_AsLong(self), 1);

    return newob;
}
#endif

static PympqObject *
Pympq_From_Pympz(PyObject *self)
{
    PympqObject *newob;

    if ((newob = (PympqObject*)Pympq_new()))
        mpq_set_z(newob->q, Pympz_AS_MPZ(self));

    return newob;
}

static PympqObject *
Pympq_From_Pyxmpz(PyObject * obj)
{
    PympqObject *newob;

    if ((newob = (PympqObject*)Pympq_new()))
        mpq_set_z(newob->q, Pyxmpz_AS_MPZ(obj));

    return newob;
}

static PympzObject *
Pympq_To_Pympz(PyObject *self)
{
    PympzObject *newob;

    if ((newob = (PympzObject*)Pympz_new()))
        mpz_set_q(newob->z, Pympq_AS_MPQ(self));

    return newob;
}

static PyxmpzObject *
Pympq_To_Pyxmpz(PyObject *self)
{
    PyxmpzObject *newob;

    if ((newob = (PyxmpzObject*)Pyxmpz_new()))
        mpz_set_q(newob->z, Pympq_AS_MPQ(self));

    return newob;
}

static PympqObject *
Pympq_From_PyLong(PyObject *self)
{
    PympqObject *newob;
    PyObject *temp = (PyObject*)Pympz_From_PyLong(self);

    if (!temp)
        return NULL;

    newob = Pympq_From_Pympz(temp);
    Py_DECREF(temp);
    return newob;
}

static PympqObject *
Pympq_From_PyFloat(PyObject *self)
{
    PympqObject *newob;

    if ((newob = (PympqObject*)Pympq_new())) {
        double d = PyFloat_AsDouble(self);
        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)newob);
            VALUE_ERROR("'mpq' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)newob);
            OVERFLOW_ERROR("'mpq' does not support Infinity");
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
Pympq_From_PyStr(PyObject *stringarg, int base)
{
    PympqObject *newob;
    unsigned char *cp;
    Py_ssize_t len;
    int i;
    PyObject *ascii_str = NULL;
    mpz_t temp;
    long expt = 0;

    if (!(newob = (PympqObject*)Pympq_new()))
        return NULL;

    if (PyBytes_Check(stringarg)) {
        len = PyBytes_Size(stringarg);
        cp = (unsigned char*)PyBytes_AsString(stringarg);
    }
    else {
        ascii_str = PyUnicode_AsASCIIString(stringarg);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            goto error;
        }
        len = PyBytes_Size(ascii_str);
        cp = (unsigned char*)PyBytes_AsString(ascii_str);
    }

    /* Don't allow NULL characters */
    for (i=0; i<len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            goto error;
        }
    }
    /* trickily delegate the rest to GMP avoiding allocations/copies */
    {
        char *whereslash = strchr((char*)cp, '/');
        char *wheredot = strchr((char*)cp, '.');
        char *whereexp = strchr((char*)cp, 'E');

        if (whereslash && wheredot) {
            VALUE_ERROR("illegal string: both . and / found");
            goto error;
        }

        if (wheredot && (base != 10)) {
            VALUE_ERROR("illegal string: embedded . requires base=10");
            goto error;
        }

        /* If base=10, no slash is found, and an exponent symbol is found, then
         * assume we have decimal number in scientific format.
         */
        if (whereexp && !whereslash && (base == 10)) {
            /* Temporarily shorten the string and continue processing as
             * normal. We'll deal with the possible exponent later.
             */
            *whereexp = '\0';
            expt = atol(whereexp+1);
        }

        if (wheredot) {
            char *counter;
            long digits = 0;

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
                /* Restore the exponent! */
                if (whereexp && (base == 10))
                    *whereexp = '\0';
                VALUE_ERROR("invalid digits");
                goto error;
            }
            /* Process the exponent. */
            digits = expt - digits;
            mpz_inoc(temp);
            if (digits < 0) {
                mpz_ui_pow_ui(mpq_denref(newob->q), 10, (unsigned long)(-digits));
            }
            else {
                mpz_ui_pow_ui(temp, 10, (unsigned long)digits);
                mpz_mul(mpq_numref(newob->q), mpq_numref(newob->q), temp);
                mpz_set_ui(mpq_denref(newob->q), 1);
            }
            mpz_cloc(temp);
            mpq_canonicalize(newob->q);

            /* Restore the decimal point. */
            *wheredot = '.';

            /* Restore the exponent! */
            if (whereexp && !whereslash && (base == 10))
                *whereexp = '\0';

            goto finish;
        }

        if (whereslash)
            *whereslash = 0;
        if (-1 == mpz_set_str(mpq_numref(newob->q), (char*)cp, base)) {
            if (whereslash)
                *whereslash = '/';
            VALUE_ERROR("invalid digits");
            goto error;
        }
        if (whereslash) {
            *whereslash = '/';
            if (-1 == mpz_set_str(mpq_denref(newob->q), whereslash+1, base)) {
                VALUE_ERROR("invalid digits");
                goto error;
            }
            if (0==mpz_sgn(mpq_denref(newob->q))) {
                ZERO_ERROR("zero denominator in 'mpq'");
                goto error;
            }
            mpq_canonicalize(newob->q);
        }
        else {
            mpz_inoc(temp);
            if (expt < 0) {
                mpz_ui_pow_ui(mpq_denref(newob->q), 10, (unsigned long)(-expt));
            }
            else {
                mpz_ui_pow_ui(temp, 10, (unsigned long)expt);
                mpz_mul(mpq_numref(newob->q), mpq_numref(newob->q), temp);
                mpz_set_ui(mpq_denref(newob->q), 1);
            }
            mpz_cloc(temp);
            mpq_canonicalize(newob->q);
            if (whereexp && (base == 10))
                *whereexp = 'E';
        }
    }

  finish:
    Py_XDECREF(ascii_str);
    return newob;

  error:
    Py_DECREF((PyObject*)newob);
    Py_XDECREF(ascii_str);
    return NULL;
}

static PyObject *
Pympq_To_PyLong(PympqObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympq_To_Pympz((PyObject*)self);

    if (!temp)
        return NULL;

    result = Pympz_To_PyLong(temp);
    Py_DECREF((PyObject*)temp);

    return result;
}

#ifdef PY2
static PyObject *
Pympq_To_PyInt(PympqObject *self)
{
    PyObject* result;
    PympzObject *temp = Pympq_To_Pympz((PyObject*)self);

    if (!temp)
        return NULL;

    result = Pympz_To_PyIntOrLong(temp);
    Py_DECREF((PyObject*)temp);

    return result;
}
#endif

static PyObject *
Pympq_To_PyFloat(PympqObject *self)
{
    double res = mpq_get_d(self->q);

    return PyFloat_FromDouble(res);
}

static int qden_1(mpq_t q)
{
    return 0 == mpz_cmp_ui(mpq_denref(q),1);
}

static PyObject *
Pympq_To_PyStr(PympqObject *self, int base, int option)
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

/* NOTE: Pympq_From_DecimalRaw returns an invalid mpq object when attempting to
 *       convert a NaN or inifinity. If the denominator is 0, then interpret
 *       the numerator as:
 *         -1: -Infinity
 *          0: Nan
 *          1: Infinity
 *
 *       If the numerator is 0 and the denominator is negative, then the value
 *       is -0.
 *
 *       These conventions are not supported by GMP/MPIR, but are used by
 *       MPFR.
 */

#if PY_VERSION_HEX < 0x03030000
static PympqObject*
Pympq_From_DecimalRaw(PyObject* obj)
{
    PympqObject *result;
    PyObject *d_exp, *d_int, *d_sign, *d_is_special;
    mpir_si exp;
    mpz_t temp;
    const char *string;

    if (!(result = (PympqObject*)Pympq_new()))
        return NULL;
    mpq_set_si(result->q, 0, 1);

    d_exp = PyObject_GetAttrString(obj, "_exp");
    d_int = PyObject_GetAttrString(obj, "_int");
    d_sign = PyObject_GetAttrString(obj, "_sign");
    d_is_special = PyObject_GetAttrString(obj, "_is_special");
    if (!d_exp || !d_int || !d_sign || !d_is_special) {
        SYSTEM_ERROR("Object does not appear to be Decimal");
        goto error;
    }

    if (PyObject_IsTrue(d_is_special)) {
        string = Py2or3String_AsString(d_exp);
        if (string[0] == 'N' || string[0] == 'n') {
            mpz_set_si(mpq_denref(result->q), 0);
            goto okay;
        }
        if (string[0] == 'F') {
            if (PyObject_IsTrue(d_sign))
                mpq_set_si(result->q, -1, 0);
            else
                mpq_set_si(result->q, 1, 0);
            goto okay;
        }
        SYSTEM_ERROR("Cannot convert Decimal to mpq");
        goto error;
    }

    if (mpz_set_PyStr(mpq_numref(result->q), d_int, 10) == -1) {
        SYSTEM_ERROR("Cannot convert Decimal to mpq");
        goto error;
    }

    exp = PyIntOrLong_AsSI(d_exp);
    if (exp == -1 && PyErr_Occurred()) {
        SYSTEM_ERROR("Decimal _exp is not valid or overflow occurred");
        goto error;
    }

    mpz_inoc(temp);
    if (exp <= 0)
        mpz_ui_pow_ui(mpq_denref(result->q), 10, (mpir_ui)(-exp));
    else {
        mpz_inoc(temp);
        mpz_ui_pow_ui(temp, 10, (mpir_ui)(exp));
        mpz_mul(mpq_numref(result->q), mpq_numref(result->q), temp);
        mpz_cloc(temp);
    }

    mpq_canonicalize(result->q);

    /* For -0, we need a negative denominator. */
    if (PyObject_IsTrue(d_sign)) {
        if (!mpz_cmp_si(mpq_numref(result->q), 0))
            mpz_set_si(mpq_denref(result->q), -1);
        else
            mpz_mul_si(mpq_numref(result->q), mpq_numref(result->q), -1);
    }

  okay:
    Py_DECREF(d_exp);
    Py_DECREF(d_int);
    Py_DECREF(d_sign);
    Py_DECREF(d_is_special);
    return result;

  error:
    Py_XDECREF(d_exp);
    Py_XDECREF(d_int);
    Py_XDECREF(d_sign);
    Py_XDECREF(d_is_special);
    Py_DECREF((PyObject*)result);
    return NULL;

}
#else
static PympqObject*
Pympq_From_DecimalRaw(PyObject* obj)
{
    PympqObject *result;
    PyObject *temp = NULL, *d_is_inf = NULL, *d_is_nan = NULL;
    PyObject *d_is_zero = NULL, *d_is_signed = NULL, *s = NULL;

    if (!(result = (PympqObject*)Pympq_new()))
        return NULL;

    d_is_inf = PyObject_GetAttrString(obj, "is_infinite");
    d_is_nan = PyObject_GetAttrString(obj, "is_nan");
    d_is_zero = PyObject_GetAttrString(obj, "is_zero");
    d_is_signed = PyObject_GetAttrString(obj, "is_signed");
    if (!d_is_inf || !d_is_nan || !d_is_zero || !d_is_signed) {
        SYSTEM_ERROR("Object does not appear to be Decimal");
        goto error;
    }

    if (!(temp = PyObject_CallFunctionObjArgs(d_is_nan, NULL)))
        goto error;
    if (PyObject_IsTrue(temp)) {
        mpz_set_si(mpq_numref(result->q), 0);
        mpz_set_si(mpq_denref(result->q), 0);
        goto okay;
    }
    Py_DECREF(temp);

    if (!(temp = PyObject_CallFunctionObjArgs(d_is_inf, NULL)))
        goto error;
    if (PyObject_IsTrue(temp)) {
        Py_DECREF(temp);
        if (!(temp = PyObject_CallFunctionObjArgs(d_is_signed, NULL)))
            goto error;
        if (PyObject_IsTrue(temp)) {
            mpz_set_si(mpq_numref(result->q), -1);
            mpz_set_si(mpq_denref(result->q), 0);
        }
        else {
            mpz_set_si(mpq_numref(result->q), 1);
            mpz_set_si(mpq_denref(result->q), 0);
        }
        goto okay;
    }
    Py_DECREF(temp);

    if (!(temp = PyObject_CallFunctionObjArgs(d_is_zero, NULL)))
        goto error;
    if (PyObject_IsTrue(temp)) {
        Py_DECREF(temp);
        if (!(temp = PyObject_CallFunctionObjArgs(d_is_signed, NULL)))
            goto error;
        if (PyObject_IsTrue(temp)) {
            mpz_set_si(mpq_numref(result->q), 0);
            mpz_set_si(mpq_denref(result->q), -1);
        }
        else {
            mpz_set_si(mpq_numref(result->q), 0);
            mpz_set_si(mpq_denref(result->q), 1);
        }
        goto okay;
    }

    Py_DECREF(result);

    s = PyObject_Str(obj);
    if (s) {
        result = Pympq_From_PyStr(s, 10);
        Py_DECREF(s);
    }

  okay:
    Py_DECREF(temp);
    Py_DECREF(d_is_inf);
    Py_DECREF(d_is_nan);
    Py_DECREF(d_is_zero);
    Py_DECREF(d_is_signed);
    return result;

  error:
    Py_XDECREF(temp);
    Py_XDECREF(d_is_inf);
    Py_XDECREF(d_is_nan);
    Py_XDECREF(d_is_zero);
    Py_XDECREF(d_is_signed);
    Py_DECREF((PyObject*)result);
    return NULL;
}
#endif

static PympqObject*
Pympq_From_Decimal(PyObject* obj)
{
    PympqObject *result;

    if ((result = Pympq_From_DecimalRaw(obj))) {
        if (!mpz_cmp_si(mpq_numref(result->q), 0)) {
            if (mpz_cmp_si(mpq_denref(result->q), 0) < 0) {
                VALUE_ERROR("'mpq' does not support -0");
                goto error;
            }
            else if (mpz_cmp_si(mpq_denref(result->q), 0) == 0) {
                VALUE_ERROR("'mpq' does not support NaN");
                goto error;
            }
        }
        else if (!mpz_cmp_si(mpq_denref(result->q), 0)) {
            OVERFLOW_ERROR("'mpq' does not support Infinity");
            goto error;
        }
    }
    return result;

  error:
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PympqObject*
Pympq_From_Fraction(PyObject* obj)
{
    PympqObject *result;
    PyObject *num, *den;

    if (!(result = (PympqObject*)Pympq_new()))
        return NULL;
    mpq_set_si(result->q, 0, 1);

    num = PyObject_GetAttrString(obj, "numerator");
    den = PyObject_GetAttrString(obj, "denominator");
    if (!num || !PyIntOrLong_Check(num) || !den || !PyIntOrLong_Check(den)) {
        SYSTEM_ERROR("Object does not appear to be Fraction");
        Py_XDECREF(num);
        Py_XDECREF(den);
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    mpz_set_PyIntOrLong(mpq_numref(result->q), num);
    mpz_set_PyIntOrLong(mpq_denref(result->q), den);
    Py_DECREF(num);
    Py_DECREF(den);
    return result;
}

static PympqObject*
Pympq_From_Number(PyObject* obj)
{
    PympqObject* newob = 0;

    if (Pympq_Check(obj)) {
        Py_INCREF(obj);
        newob = (PympqObject *) obj;
    }
    else if (Pympz_Check(obj)) {
        newob = Pympq_From_Pympz(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympq_From_PyInt(obj);
#endif
    }
#ifdef WITHMPFR
    else if (Pympfr_Check(obj)) {
        newob = Pympfr_To_Pympq(obj);
    }
#endif
    else if (PyFloat_Check(obj)) {
        newob = Pympq_From_PyFloat(obj);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympq_From_PyLong(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pympq_From_Pyxmpz(obj);
    }
    else if (isDecimal(obj)) {
        newob = Pympq_From_Decimal(obj);
    }
    else if (isFraction(obj)) {
        newob = Pympq_From_Fraction(obj);
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
        newob = Pympq_From_Pympz(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympq_From_PyInt(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = Pympq_From_PyLong(obj);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pympq_From_Pyxmpz(obj);
    }
    else if (isFraction(obj)) {
        newob = Pympq_From_Fraction(obj);
    }

    return newob;
}

/*
 * coerce any number to a mpq
 */

int
Pympq_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympqObject* newob = Pympq_From_Number(arg);
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
Pympq_To_Str(PympqObject *self)
{
    /* base-10, no tag */
    return Pympq_To_PyStr(self, 10, 0);
}

static PyObject *
Pympq_To_Repr(PympqObject *self)
{
    /* base-10, with tag */
    return Pympq_To_PyStr(self, 10, 1);
}

#ifdef WITHMPFR
/* Functions that operate strictly on mpfr. */

/* Make a copy of an mpfr object. If bits is 0, the new object will have
 * the same precision as the original object. If the requested precision
 * is less than the precision of the original object, the new object
 * will be rounded to requested precision using the current rounding mode.
 */

static PympfrObject *
Pympfr_From_Pympfr(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if (bits == 0)
        bits = mpfr_get_prec(Pympfr_AS_MPFR(self));

    if ((result = (PympfrObject*)Pympfr_new(bits))) {
        result->rc = mpfr_set(result->f,
                              Pympfr_AS_MPFR(self),
                              context->ctx.mpfr_round);
    }

    return result;
}

static PympfrObject *
Pympfr_From_PyFloat(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = (PympfrObject*)Pympfr_new(bits))) {
        result->rc = mpfr_set_d(result->f,
                                PyFloat_AS_DOUBLE(self),
                                context->ctx.mpfr_round);
    }

    return result;
}

static PympfrObject *
Pympfr_From_Pympz(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = (PympfrObject*)Pympfr_new(bits))) {
        result->rc = mpfr_set_z(result->f,
                                Pympz_AS_MPZ(self),
                                context->ctx.mpfr_round);
    }

    return result;
}

#define Pympfr_From_Pyxmpz Pympfr_From_Pympz

static PympzObject *
Pympfr_To_Pympz(PyObject *self)
{
    PympzObject *result;

    if ((result = (PympzObject*)Pympz_new())) {
        if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, Pympfr_AS_MPFR(self), context->ctx.mpfr_round);
    }

    return result;
}

static PyxmpzObject *
Pympfr_To_Pyxmpz(PyObject *self)
{
    PyxmpzObject *result;

    if ((result = (PyxmpzObject*)Pyxmpz_new())) {
        if (mpfr_nan_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(Pympfr_AS_MPFR(self))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, Pympfr_AS_MPFR(self), context->ctx.mpfr_round);
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
        prec = (mpfr_prec_t)(-mpfr_get_si(err->f, context->ctx.mpfr_round));

    if (errsign <= 0 && (prec < 2 || prec > mpfr_get_prec(self->f))) {
        VALUE_ERROR("Requested precision out-of-bounds.");
        return NULL;
    }

    if (!(result = (PympqObject*)Pympq_new()))
        return NULL;

    mpfr_init2(minerr, F2Q_PREC);
    if (errsign <= 0) {
        mpfr_set_ui(minerr, 1, context->ctx.mpfr_round);
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
    mpfr_set_si(r1[0], 0, context->ctx.mpfr_round);
    mpfr_set_si(r1[1], 0, context->ctx.mpfr_round);
    mpfr_set_si(r1[2], 1, context->ctx.mpfr_round);
    mpfr_set_si(r2[0], 0, context->ctx.mpfr_round);
    mpfr_set_si(r2[1], 1, context->ctx.mpfr_round);
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
        result = (PympqObject*)Pympz_new();
        mpfr_get_z(Pympz_AS_MPZ(result), r2[2], context->ctx.mpfr_round);
        if (negative)
            mpz_neg(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result));
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

static PympqObject *
Pympfr_To_Pympq(PyObject *self)
{
    return stern_brocot((PympfrObject*)self, 0, 0, 0);
}

static PympfrObject *
Pympfr_From_Pympq(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = (PympfrObject*)Pympfr_new(bits)))
        result->rc = mpfr_set_q(result->f, Pympq_AS_MPQ(self),
                                context->ctx.mpfr_round);
    return result;
}

static PympfrObject *
Pympfr_From_PyLong(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;
    PyObject *temp = (PyObject*)Pympz_From_PyLong(self);

    if (!temp)
        return NULL;
    result = Pympfr_From_Pympz(temp, bits);
    Py_DECREF(temp);
    return result;
}

#ifdef PY2
static PympfrObject *
Pympfr_From_PyInt(PyObject *self, mpfr_prec_t bits)
{
    PympfrObject *result;

    if ((result = (PympfrObject*)Pympfr_new(bits)))
        result->rc = mpfr_set_si(result->f, PyInt_AsLong(self),
                                 context->ctx.mpfr_round);
    return result;
}

static PyObject *
Pympfr_To_PyInt(PympfrObject *self)
{
    PyObject *result;
    PympzObject *temp = Pympfr_To_Pympz((PyObject*)self);

    if (!temp)
        return NULL;
    result = Pympz_To_PyIntOrLong(temp);
    Py_DECREF((PyObject*)temp);
    return result;
}
#endif

static PympfrObject *
Pympfr_From_PyStr(PyObject *s, int base, mpfr_prec_t bits)
{
    PympfrObject *result;
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

    if (!(result = (PympfrObject*)Pympfr_new(prec))) {
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
Pympfr_To_PyLong(PympfrObject *self)
{
    PyObject *result;
    PympzObject *temp = Pympfr_To_Pympz((PyObject*)self);

    if (!temp) return NULL;

    result = Pympz_To_PyLong(temp);
    Py_DECREF((PyObject*)temp);

    return result;
}

static PyObject *
Pympfr_To_PyFloat(PympfrObject *self)
{
    double res = mpfr_get_d(self->f, context->ctx.mpfr_round);

    return PyFloat_FromDouble(res);
}

static PyObject*
Pympfr_To_PyStr(PympfrObject *self, int base, int digits)
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
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, context->ctx.mpfr_round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in Pympfr_To_PyStr");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
    mpfr_free_str(buffer);
    return result;
}

static PympfrObject *
Pympfr_From_Decimal(PyObject* obj, mpfr_prec_t bits)
{
    PympfrObject *result = (PympfrObject*)Pympfr_new(0);
    PympqObject *temp = Pympq_From_DecimalRaw(obj);

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
        result = Pympfr_From_Pympq((PyObject*)temp, bits);
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
            newob = Pympfr_From_Pympfr((PyObject*)obj, bits);
        }
    }
    else if (Pympfr_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer valid
         * and mpfr_check_range needs to be called. */
        if (context->ctx.trap_expbound) {
            GMPY_EXPBOUND("exponent of existing 'mpfr' incompatible with current context");
            return NULL;
        }
        if ((newob = (PympfrObject*)Pympfr_new(mpfr_get_prec(Pympfr_AS_MPFR(obj))))) {
            mpfr_set(newob->f, Pympfr_AS_MPFR(obj), context->ctx.mpfr_round);
            newob->round_mode = ((PympfrObject*)obj)->round_mode;
            newob->rc = ((PympfrObject*)obj)->rc;
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
    else if (Pympq_Check(obj)) {
        newob = Pympfr_From_Pympq(obj, bits);
    }
    else if (Pympz_Check(obj)) {
        newob = Pympfr_From_Pympz(obj, bits);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympfr_From_PyLong(obj, bits);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pympfr_From_Pyxmpz(obj, bits);
    }
    else if (isDecimal(obj)) {
        newob = Pympfr_From_Decimal(obj, bits);
    }
    else if (isFraction(obj)) {
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

/* str and repr implementations for mpfr */
static PyObject *
Pympfr_To_Str(PympfrObject *self)
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
Pympfr_To_Repr(PympfrObject *self)
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
#endif

#ifdef WITHMPC
static PympcObject *
Pympc_From_Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, Pympc_AS_MPC(self));
    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        mpc_set(result->c, Pympc_AS_MPC(self), GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
Pympc_From_PyComplex(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(self),
                    PyComplex_ImagAsDouble(self), GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
Pympc_From_Pympfr(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (!rprec)
        rprec = mpfr_get_prec(Pympfr_AS_MPFR(self));
    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_fr(result->c, Pympfr_AS_MPFR(self),
                                GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
Pympc_From_PyFloat(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (!rprec)
        rprec = DBL_MANT_DIG;
    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(self),
                               GET_MPC_ROUND(context));
    return result;
}

static PyObject *
Pympc_To_PyFloat(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'float'");
    return NULL;
}

PyDoc_STRVAR(doc_mpc_complex, "Convert 'mpc' to 'complex'.");

static PyObject *
Pympc_To_PyComplex(PyObject *self, PyObject *other)
{
    double real = mpfr_get_d(mpc_realref(Pympc_AS_MPC(self)),
                             GET_REAL_ROUND(context));
    double imag = mpfr_get_d(mpc_imagref(Pympc_AS_MPC(self)),
                             GET_IMAG_ROUND(context));

    return PyComplex_FromDoubles(real, imag);
}

static PympcObject *
Pympc_From_Pympz(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_z(result->c, Pympz_AS_MPZ(self),
                                GET_MPC_ROUND(context));
    return result;
}

#define Pympc_From_Pyxmpz Pympc_From_Pympz

static PympcObject *
Pympc_From_Pympq(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_q(result->c, Pympq_AS_MPQ(self),
                               GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
Pympc_From_PyLong(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;
    PyObject *temp = (PyObject*)Pympz_From_PyLong(self);

    if (!temp)
        return NULL;
    result = Pympc_From_Pympz(temp, rprec, iprec);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympc_To_PyLong(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'long'");
    return NULL;
}

#ifdef PY2
static PympcObject *
Pympc_From_PyInt(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = (PympcObject*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_si(result->c, PyInt_AsLong(self),
                                GET_MPC_ROUND(context));
    return result;
}
static PyObject *
Pympc_To_PyIntOrLong(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'int'");
    return NULL;
}
#endif

/* Conversion to/from MPC
 * Python's string representation of a complex number differs from the format
 * used by MPC. Both MPC and Python surround the complex number with '(' and
 * ')' but Python adds a 'j' after the imaginary component and MPC requires a
 * space between the real and imaginery components. PyStr2Pympc tries to work
 * around the differences as follows reading two MPFR-compatible numbers from
 * the string and storing into the real and imaginary components respectively.
 */

static PympcObject *
Pympc_From_PyStr(PyObject *s, int base, mpfr_prec_t rbits, mpfr_prec_t ibits)
{
    PympcObject *newob;
    PyObject *ascii_str = NULL;
    Py_ssize_t len;
    char *cp, *unwind, *tempchar, *lastchar;
    int firstp = 0, lastp = 0, real_rc = 0, imag_rc = 0;

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = (char*)PyBytes_AsString(s);
    }
    else if (PyUnicode_Check(s)) {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = (char*)PyBytes_AsString(ascii_str);
    }
    else {
        TYPE_ERROR("string required for PyStr2Pympc");
        return NULL;
    }

    if (!(newob = (PympcObject*)Pympc_new(rbits, ibits))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Don't allow NULL characters */
    if (strlen(cp) != len) {
        VALUE_ERROR("string without NULL characters expected");
        Py_DECREF((PyObject*)newob);
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Get a pointer to the last valid character (ignoring trailing
     * whitespace.) */
    lastchar = cp + len - 1;
    while (isspace(*lastchar))
        lastchar--;

    /* Skip trailing ). */
    if (*lastchar == ')') {
        lastp = 1;
        lastchar--;
    }

    /* Skip trailing j. */
    if (*lastchar == 'j')
        lastchar--;

    /* Skip leading whitespace. */
    while (isspace(*cp))
        cp++;

    /* Skip a leading (. */
    if (*cp == '(') {
        firstp = 1;
        cp++;
    }

    if (firstp != lastp) goto invalid_string;

    /* Read the real component first. */
    unwind = cp;
    real_rc = mpfr_strtofr(mpc_realref(newob->c), cp, &tempchar, base,
                           GET_REAL_ROUND(context));
    /* Verify that at least one valid character was read. */
    if (cp == tempchar) goto invalid_string;
    /* If the next character is a j, then the real component is 0 and
     * we just read the imaginary componenet.
     */
    if (*tempchar == 'j') {
        mpfr_set_zero(mpc_realref(newob->c), +1);
        cp = unwind;
    }
    else {
        /* Read the imaginary component next. */
        cp = tempchar;
    }
    imag_rc = mpfr_strtofr(mpc_imagref(newob->c), cp, &tempchar, base,
                           GET_IMAG_ROUND(context));

    if (cp == tempchar && tempchar > lastchar)
        goto valid_string;

    if (*tempchar != 'j' && *cp != ' ')
        goto invalid_string;

    if (tempchar <= lastchar)
        goto invalid_string;

  valid_string:
    Py_XDECREF(ascii_str);
    newob->rc = MPC_INEX(real_rc, imag_rc);
    return newob;

  invalid_string:
    VALUE_ERROR("invalid string in mpc()");
    Py_DECREF((PyObject*)newob);
    Py_XDECREF(ascii_str);
    return NULL;
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

static PyObject *
Pympc_To_PyStr(PympcObject *self, int base, int digits)
{
    PyObject *tempreal = 0, *tempimag = 0, *result;

    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    tempreal = raw_mpfr_ascii(mpc_realref(self->c), base, digits,
                            MPC_RND_RE(GET_MPC_ROUND(context)));
    tempimag = raw_mpfr_ascii(mpc_imagref(self->c), base, digits,
                            MPC_RND_IM(GET_MPC_ROUND(context)));

    if (!tempreal || !tempimag) {
        Py_XDECREF(tempreal);
        Py_XDECREF(tempimag);
        return NULL;
    }

    result = Py_BuildValue("(NN)", tempreal, tempimag);
    if (!result) {
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
    }
    return result;
}

/*
 * If obj is a Pympc and rprec/iprec are 0/0 or the same as the precision of
 * obj, then a new reference is created.
 *
 * For all other numerical types with bits = 0, the conversion is rounded
 * according to the context.
 */

static PympcObject *
Pympc_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject* newob = 0;
    PympqObject* temp = 0;
    mpfr_prec_t pr = 0, pi = 0;
    int rr, ri, dr, di;

    if (Pympc_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpc is still
         * valid in the current context. */
        if (!rprec && !iprec) {
            Py_INCREF(obj);
            newob = (PympcObject*)obj;
        }
        else {
            mpc_get_prec2(&pr, &pi, Pympc_AS_MPC(obj));
            if (rprec == pr && iprec == pi) {
                Py_INCREF(obj);
                newob = (PympcObject*)obj;
            }
            else {
                newob = Pympc_From_Pympc((PyObject*)obj, rprec, iprec);
            }
        }
    }
    else if (Pympc_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer
         * valid and mpfr_check_range needs to be called. */
        if (context->ctx.trap_expbound) {
            GMPY_EXPBOUND("exponent of existing 'mpc' incompatible with current context");
            return NULL;
        }
        /* Get the real and imaginary precisions. */
        mpc_get_prec2(&pr, &pi, Pympc_AS_MPC(obj));

        /* Get the real and imaginary inexact codes. */
        rr = MPC_INEX_RE( ((PympcObject*)obj)->rc );
        ri = MPC_INEX_IM( ((PympcObject*)obj)->rc );

        /* Get the real and imaginary rounding modes. */
        dr = MPC_RND_RE( ((PympcObject*)obj)->round_mode );
        di = MPC_RND_IM( ((PympcObject*)obj)->round_mode );

        if ((newob = (PympcObject*)Pympc_new(pr, pi))) {
            mpc_set(newob->c, Pympc_AS_MPC(obj), GET_MPC_ROUND(context));
            newob->round_mode = ((PympcObject*)obj)->round_mode;
            rr = mpfr_check_range(mpc_realref(newob->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(newob->c), ri, di);
            newob->rc = MPC_INEX(rr, ri);
        }
    }
    else if (Pympfr_Check(obj)) {
            newob = Pympc_From_Pympfr((PyObject*)obj, rprec, iprec);
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympc_From_PyFloat(obj, rprec, iprec);
    }
    else if (PyComplex_Check(obj)) {
            newob = Pympc_From_PyComplex(obj, rprec, iprec);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympc_From_PyInt(obj, rprec, iprec);
#endif
    }
    else if (Pympq_Check(obj)) {
        newob = Pympc_From_Pympq(obj, rprec, iprec);
    }
    else if (Pympz_Check(obj)) {
        newob = Pympc_From_Pympz(obj, rprec, iprec);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympc_From_PyLong(obj, rprec, iprec);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pympc_From_Pyxmpz(obj, rprec, iprec);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = Pympc_From_PyStr(s, 10, rprec, iprec);
            if (!newob) {
                Py_DECREF(s);
                return NULL;
            }
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympc_From_Pympq((PyObject *)temp, rprec, iprec);
            Py_DECREF((PyObject*)temp);
        }
    }
    return newob;
}

/*
 * coerce any number to a mpc
 */

int
Pympc_convert_arg(PyObject *arg, PyObject **ptr)
{
    PympcObject* newob = Pympc_From_Complex(arg, 0, 0);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("can't convert argument 'mpc'");
        return 0;
    }
}

/* str and repr implementations for mpc */
static PyObject *
Pympc_To_Str(PympcObject *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[30];

    mpc_get_prec2(&rbits, &ibits, Pympc_AS_MPC(self));
    rprec = (long)(log10(2) * (double)rbits) + 2;
    iprec = (long)(log10(2) * (double)ibits) + 2;

    sprintf(fmtstr, "{0:.%ld.%ldg}", rprec, iprec);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympc_To_Repr(PympcObject *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[30];

    mpc_get_prec2(&rbits, &ibits, Pympc_AS_MPC(self));
    rprec = (long)(log10(2) * (double)rbits) + 2;
    iprec = (long)(log10(2) * (double)ibits) + 2;

    if (rbits != DBL_MANT_DIG || ibits !=DBL_MANT_DIG)
        sprintf(fmtstr, "mpc('{0:.%ld.%ldg}',(%ld,%ld))",
                rprec, iprec, rbits, ibits);
    else
        sprintf(fmtstr, "mpc('{0:.%ld.%ldg}')", rprec, iprec);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}
#endif



