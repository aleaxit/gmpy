/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_gmp.c                                                     *
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

/* This file contains all the conversion functions for GMP data types.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. gmpy2 uses a LBYL (Look Before You Leap) approach and
 * identifies the numeric type before conversion before conversion to a gmpy2
 * type. The basic operations (+, -, *, /) are optimized to directly work with
 * some basic types such as C longs or doubles.
 */

static XMPZ_Object *
GMPy_XMPZ_From_XMPZ(PyObject *self)
{
    XMPZ_Object *result;

    if (!XMPZ_Check(self)) {
        TYPE_ERROR("object with type xmpz expected");
        return NULL;
    }

    if ((result = GMPy_XMPZ_New()))
        mpz_set(result->z, MPZ(self));

    return result;
}

static MPZ_Object *
GMPy_MPZ_From_XMPZ(PyObject *self)
{
    MPZ_Object *result;

    if (!XMPZ_Check(self)) {
        TYPE_ERROR("object with type xmpz expected");
        return NULL;
    }

    if ((result = GMPy_MPZ_New()))
        mpz_set(result->z, MPZ(self));

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPZ(PyObject *self)
{
    XMPZ_Object *result;

    if (!MPZ_Check(self)) {
        TYPE_ERROR("object with type mpz expected");
        return NULL;
    }

    if ((result = GMPy_XMPZ_New()))
        mpz_set(result->z, MPZ(self));

    return result;
}

#ifdef PY2
static MPZ_Object *
GMPy_MPZ_From_PyInt(PyObject *self)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New()))
        mpz_set_si(result->z, PyInt_AS_LONG(self));
    return result;
}

static XMPZ_Object *
Pyxmpz_From_PyInt(PyObject *self)
{
    XMPZ_Object *newob;

    if ((newob = GMPy_XMPZ_New()))
        mpz_set_si(newob->z, PyInt_AsLong(self));
    return newob;
}
#endif

static MPZ_Object *
Pympz_From_PyFloat(PyObject *self)
{
    MPZ_Object *newob;

    if ((newob = GMPy_MPZ_New())) {
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

static XMPZ_Object *
Pyxmpz_From_PyFloat(PyObject *self)
{
    XMPZ_Object *newob;

    if ((newob = GMPy_XMPZ_New())) {
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

static MPZ_Object *
GMPy_MPZ_From_PyLong(PyObject *obj)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    if (-1 == mpz_set_PyIntOrLong(result->z, obj)) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    else {
        return result;
    }
}

static XMPZ_Object *
GMPy_XMPZ_From_PyLong(PyObject *obj)
{
    XMPZ_Object *result;

    if (!(result = GMPy_XMPZ_New()))
        return NULL;

    if (-1 == mpz_set_PyIntOrLong(result->z, obj)) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    else {
        return result;
    }
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

static MPZ_Object *
Pympz_From_PyStr(PyObject *s, int base)
{
    MPZ_Object *newob;

    if (!(newob = GMPy_MPZ_New()))
        return NULL;

    if (mpz_set_PyStr(newob->z, s, base) == -1) {
        Py_DECREF((PyObject*)newob);
        return NULL;
    }
    return newob;
}

static XMPZ_Object *
Pyxmpz_From_PyStr(PyObject *s, int base)
{
    XMPZ_Object *newob;

    if (!(newob = GMPy_XMPZ_New()))
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
Pympz_To_PyLong(MPZ_Object *self)
{
    return mpz_get_PyLong(MPZ(self));
}

static PyObject *
Pyxmpz_To_PyLong(XMPZ_Object *self)
{
    return mpz_get_PyLong(MPZ(self));
}

/* The _To_PyIntOrLong functions should be used when converting a number back
 * to a Python value since is automatically returns an "int" or "long" when
 * using Python 2.x. The _To_PyLong functions (above) should only be used
 * when a PyLong is specifically needed for Python 2.x.
 */

static PyObject *
Pympz_To_PyIntOrLong(MPZ_Object *self)
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
Pyxmpz_To_PyIntOrLong(XMPZ_Object *self)
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
Pympz_To_PyFloat(MPZ_Object *self)
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
Pympz_To_PyStr(MPZ_Object *self, int base, int option)
{
    return mpz_ascii(self->z, base, option);
}

static PyObject *
Pyxmpz_To_PyStr(XMPZ_Object *self, int base, int option)
{
    return xmpz_ascii(self->z, base, option);
}

static MPZ_Object *
GMPy_MPZ_From_Number_Temp(PyObject *obj)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        Py_INCREF(obj);
        return (MPZ_Object*)obj;
    }

#ifdef PY2
    if (PyInt_Check(obj))
        return Pympz_From_PyInt(obj);
#endif

    if (PyLong_Check(obj))
        return GMPy_MPZ_From_PyLong(obj);

    if (MPQ_Check(obj))
        return Pympq_To_Pympz(obj);

    if (MPFR_Check(obj))
        return Pympfr_To_Pympz(obj);

    if (PyFloat_Check(obj))
        return Pympz_From_PyFloat(obj);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ(obj);

    if (IS_DECIMAL(obj)) {
        PyObject *temp = PyNumber_Long(obj);

        if (temp) {
            result = GMPy_MPZ_From_PyLong(temp);
            Py_DECREF(temp);
        }
        return result;
    }

    if (IS_FRACTION(obj)) {
        MPQ_Object *temp = Pympq_From_Fraction(obj);

        if (temp) {
            result = Pympq_To_Pympz((PyObject*)temp);
            Py_DECREF((PyObject*)temp);
        }
        return result;
    }

    TYPE_ERROR("cannot convert object to mpz");
    return result;
}

static MPZ_Object *
GMPy_MPZ_From_Number_New(PyObject *obj)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        if ((result = GMPy_MPZ_New()))
            mpz_set(result->z, MPZ(obj));

        return result;
    }

#ifdef PY2
    if (PyInt_Check(obj))
        return Pympz_From_PyInt(obj);
#endif

    if (PyLong_Check(obj))
        return GMPy_MPZ_From_PyLong(obj);

    if (MPQ_Check(obj))
        return Pympq_To_Pympz(obj);

    if (MPFR_Check(obj))
        return Pympfr_To_Pympz(obj);

    if (PyFloat_Check(obj))
        return Pympz_From_PyFloat(obj);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ(obj);

    if (IS_DECIMAL(obj)) {
        PyObject *temp = PyNumber_Long(obj);

        if (temp) {
            result = GMPy_MPZ_From_PyLong(temp);
            Py_DECREF(temp);
        }
        return result;
    }

    if (IS_FRACTION(obj)) {
        MPQ_Object *temp = Pympq_From_Fraction(obj);

        if (temp) {
            result = Pympq_To_Pympz((PyObject*)temp);
            Py_DECREF((PyObject*)temp);
        }
        return result;
    }

    TYPE_ERROR("cannot convert object to mpz");
    return result;
}


static XMPZ_Object*
Pyxmpz_From_Number(PyObject* obj)
{
    XMPZ_Object* newob = NULL;
    MPQ_Object* temp = NULL;

    if (MPZ_Check(obj)) {
        newob = GMPy_XMPZ_From_MPZ(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pyxmpz_From_PyInt(obj);
#endif
    }
    else if (PyLong_Check(obj)) {
        newob = GMPy_XMPZ_From_PyLong(obj);
    }
    else if (MPQ_Check(obj)) {
        newob = Pympq_To_Pyxmpz(obj);
    }
    else if (MPFR_Check(obj)) {
        newob = Pympfr_To_Pyxmpz(obj);
    }
    else if (PyFloat_Check(obj)) {
        newob = Pyxmpz_From_PyFloat(obj);
    }
    else if (XMPZ_Check(obj)) {
        newob = GMPy_XMPZ_From_XMPZ(obj);
    }
    else if (IS_DECIMAL(obj)) {
        PyObject *s = PyNumber_Long(obj);
        if (s) {
            newob = GMPy_XMPZ_From_PyLong(s);
            Py_DECREF(s);
        }
    }
    else if (IS_FRACTION(obj)) {
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

static MPZ_Object *
GMPy_MPZ_From_Integer_Temp(PyObject *obj)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        Py_INCREF(obj);
        return (MPZ_Object*)obj;
    }

#ifdef PY2
    if (PyInt_Check(obj))
        return Pympz_From_PyInt(obj);
#endif

    if (PyLong_Check(obj))
        return GMPy_MPZ_From_PyLong(obj);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ(obj);

    TYPE_ERROR("cannot convert object to mpz");
    return result;
}

static MPZ_Object *
GMPy_MPZ_From_Integer_New(PyObject *obj)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        if ((result = GMPy_MPZ_New()))
            mpz_set(result->z, MPZ(obj));
        return result;
    }

#ifdef PY2
    if (PyInt_Check(obj))
        return Pympz_From_PyInt(obj);
#endif

    if (PyLong_Check(obj))
        return GMPy_MPZ_From_PyLong(obj);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ(obj);

    TYPE_ERROR("cannot convert object to mpz");
    return result;
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
        if (mpz_fits_slong_p(MPZ(obj))) {
            return (long)mpz_get_si(MPZ(obj));
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
        if (mpz_fits_slong_p(MPZ(obj))) {
            return mpz_get_si(MPZ(obj));
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
        if (mpz_fits_ulong_p(MPZ(obj))) {
            return mpz_get_ui(MPZ(obj));
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
        if (mpz_fits_si_p(MPZ(obj))) {
            return mpz_get_si(MPZ(obj));
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
        if (mpz_fits_ui_p(MPZ(obj))) {
            return mpz_get_ui(MPZ(obj));
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
        if (mpz_fits_slong_p(MPZ(obj))) {
            return (Py_ssize_t)mpz_get_si(MPZ(obj));
        }
        else {
            /* This section should only be called on Win64. */
            temp = mpz_get_PyLong(MPZ(obj));
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
    MPZ_Object* newob = Pympz_From_Integer(arg);

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
Pympz_To_Str(MPZ_Object *self)
{
    /* base-10, no tag */
    return Pympz_To_PyStr(self, 10, 0);
}

static PyObject *
Pympz_To_Repr(MPZ_Object *self)
{
    /* base-10, with tag */
    return Pympz_To_PyStr(self, 10, 1);
}

/* str and repr implementations for xmpz */
static PyObject *
Pyxmpz_To_Str(XMPZ_Object *self)
{
    /* base-10, no tag */
    return Pyxmpz_To_PyStr(self, 10, 0);
}

static PyObject *
Pyxmpz_To_Repr(XMPZ_Object *self)
{
    /* base-10, with tag */
    return Pyxmpz_To_PyStr(self, 10, 1);
}

#ifdef PY2
static MPQ_Object *
Pympq_From_PyInt(PyObject *self)
{
    MPQ_Object *newob;

    if ((newob = (MPQ_Object*)Pympq_new()))
        mpq_set_si(newob->q, PyInt_AsLong(self), 1);

    return newob;
}
#endif

static MPQ_Object *
Pympq_From_Pympz(PyObject *self)
{
    MPQ_Object *newob;

    if ((newob = (MPQ_Object*)Pympq_new()))
        mpq_set_z(newob->q, MPZ(self));

    return newob;
}

static MPQ_Object *
Pympq_From_Pyxmpz(PyObject * obj)
{
    MPQ_Object *newob;

    if ((newob = (MPQ_Object*)Pympq_new()))
        mpq_set_z(newob->q, MPZ(obj));

    return newob;
}

static MPZ_Object *
Pympq_To_Pympz(PyObject *self)
{
    MPZ_Object *newob;

    if ((newob = GMPy_MPZ_New()))
        mpz_set_q(newob->z, MPQ(self));

    return newob;
}

static XMPZ_Object *
Pympq_To_Pyxmpz(PyObject *self)
{
    XMPZ_Object *newob;

    if ((newob = GMPy_XMPZ_New()))
        mpz_set_q(newob->z, MPQ(self));

    return newob;
}

static MPQ_Object *
Pympq_From_PyLong(PyObject *self)
{
    MPQ_Object *newob;
    PyObject *temp = (PyObject*)GMPy_MPZ_From_PyLong(self);

    if (!temp)
        return NULL;

    newob = Pympq_From_Pympz(temp);
    Py_DECREF(temp);
    return newob;
}

static MPQ_Object *
Pympq_From_PyFloat(PyObject *self)
{
    MPQ_Object *newob;

    if ((newob = (MPQ_Object*)Pympq_new())) {
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
static MPQ_Object *
Pympq_From_PyStr(PyObject *stringarg, int base)
{
    MPQ_Object *newob;
    unsigned char *cp;
    Py_ssize_t len;
    int i;
    PyObject *ascii_str = NULL;
    mpz_t temp;
    long expt = 0;

    if (!(newob = (MPQ_Object*)Pympq_new()))
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
Pympq_To_PyLong(MPQ_Object *self)
{
    PyObject* result;
    MPZ_Object *temp = Pympq_To_Pympz((PyObject*)self);

    if (!temp)
        return NULL;

    result = Pympz_To_PyLong(temp);
    Py_DECREF((PyObject*)temp);

    return result;
}

#ifdef PY2
static PyObject *
Pympq_To_PyInt(MPQ_Object *self)
{
    PyObject* result;
    MPZ_Object *temp = Pympq_To_Pympz((PyObject*)self);

    if (!temp)
        return NULL;

    result = Pympz_To_PyIntOrLong(temp);
    Py_DECREF((PyObject*)temp);

    return result;
}
#endif

static PyObject *
Pympq_To_PyFloat(MPQ_Object *self)
{
    double res = mpq_get_d(self->q);

    return PyFloat_FromDouble(res);
}

static int qden_1(mpq_t q)
{
    return 0 == mpz_cmp_ui(mpq_denref(q),1);
}

static PyObject *
Pympq_To_PyStr(MPQ_Object *self, int base, int option)
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
static MPQ_Object*
Pympq_From_DecimalRaw(PyObject* obj)
{
    MPQ_Object *result;
    PyObject *d_exp, *d_int, *d_sign, *d_is_special;
    mpir_si exp;
    mpz_t temp;
    const char *string;

    if (!(result = (MPQ_Object*)Pympq_new()))
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
static MPQ_Object*
Pympq_From_DecimalRaw(PyObject* obj)
{
    MPQ_Object *result;
    PyObject *temp = NULL, *d_is_inf = NULL, *d_is_nan = NULL;
    PyObject *d_is_zero = NULL, *d_is_signed = NULL, *s = NULL;

    if (!(result = (MPQ_Object*)Pympq_new()))
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

static MPQ_Object*
Pympq_From_Decimal(PyObject* obj)
{
    MPQ_Object *result;

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

static MPQ_Object*
Pympq_From_Fraction(PyObject* obj)
{
    MPQ_Object *result;
    PyObject *num, *den;

    if (!(result = (MPQ_Object*)Pympq_new()))
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

static MPQ_Object*
Pympq_From_Number(PyObject* obj)
{
    MPQ_Object* newob = 0;

    if (MPQ_Check(obj)) {
        Py_INCREF(obj);
        newob = (MPQ_Object *) obj;
    }
    else if (MPZ_Check(obj)) {
        newob = Pympq_From_Pympz(obj);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympq_From_PyInt(obj);
#endif
    }
    else if (MPFR_Check(obj)) {
        newob = Pympfr_To_Pympq(obj);
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympq_From_PyFloat(obj);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympq_From_PyLong(obj);
    }
    else if (XMPZ_Check(obj)) {
        newob = Pympq_From_Pyxmpz(obj);
    }
    else if (IS_DECIMAL(obj)) {
        newob = Pympq_From_Decimal(obj);
    }
    else if (IS_FRACTION(obj)) {
        newob = Pympq_From_Fraction(obj);
    }

    return newob;
}

/* Convert an integer or mpz to mpq. */

static MPQ_Object*
Pympq_From_Rational(PyObject* obj)
{
    MPQ_Object* newob = 0;

    if (MPQ_Check(obj)) {
        Py_INCREF(obj);
        newob = (MPQ_Object *) obj;
    }
    else if (MPZ_Check(obj)) {
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
    else if (XMPZ_Check(obj)) {
        newob = Pympq_From_Pyxmpz(obj);
    }
    else if (IS_FRACTION(obj)) {
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
    MPQ_Object* newob = Pympq_From_Number(arg);
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
Pympq_To_Str(MPQ_Object *self)
{
    /* base-10, no tag */
    return Pympq_To_PyStr(self, 10, 0);
}

static PyObject *
Pympq_To_Repr(MPQ_Object *self)
{
    /* base-10, with tag */
    return Pympq_To_PyStr(self, 10, 1);
}

