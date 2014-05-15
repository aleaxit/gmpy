/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_gmp.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Case Van Horsen      *
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

/* ======================================================================== *
 * Conversion between native Python objects and MPZ.                        *
 * ======================================================================== */

static MPZ_Object *
GMPy_MPZ_From_PyIntOrLong(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    assert(PyIntOrLong_Check(obj));

    if(!(result = GMPy_MPZ_New(context)))
        return NULL;

    mpz_set_PyIntOrLong(result->z, obj);
    return result;
}

static MPZ_Object *
GMPy_MPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    if (mpz_set_PyStr(result->z, s, base) == -1) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    return result;
}

static MPZ_Object *
GMPy_MPZ_From_PyFloat(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    assert(PyFloat_Check(obj));

    if ((result = GMPy_MPZ_New(context))) {
        double d = PyFloat_AsDouble(obj);

        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(result->z, d);
    }
    return result;
}

#ifdef PY2
static PyObject *
GMPy_PyLong_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    assert(MPZ_Check(obj));

    return mpz_get_PyLong(obj->z);
}

static PyObject *
GMPy_MPZ_Long_Slot(MPZ_Object *self)
{
    return GMPy_PyLong_From_MPZ(self, NULL);
}
#endif

/* The PyIntOrLong functions should be used when converting a number back
 * to a Python value since is automatically returns an "int" or "long" when
 * using Python 2.x. The PyLong_From functions (above) should only be used
 * when a PyLong is specifically needed for Python 2.x.
 */

static PyObject *
GMPy_PyIntOrLong_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    assert(MPZ_Check(obj));

#ifdef PY3
    return mpz_get_PyLong(obj->z);
#else
    if (mpz_fits_slong_p(obj->z))
        /* cast is safe since we know it fits in a signed long */
        return PyInt_FromLong((long)mpz_get_si(obj->z));
    else
        return mpz_get_PyLong(obj->z);
#endif
}

static PyObject *
GMPy_MPZ_Int_Slot(MPZ_Object *self)
{
    return GMPy_PyIntOrLong_From_MPZ(self, NULL);
}

static PyObject *
GMPy_PyFloat_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    double res;

    assert(MPZ_Check(obj));

    res = mpz_get_d(obj->z);

    if (Py_IS_INFINITY(res)) {
        OVERFLOW_ERROR("'mpz' too large to convert to float");
        return NULL;
    }

    return PyFloat_FromDouble(res);
}

static PyObject *
GMPy_MPZ_Float_Slot(MPZ_Object *self)
{
    return GMPy_PyFloat_From_MPZ(self, NULL);
}

static PyObject *
GMPy_PyStr_From_MPZ(MPZ_Object *obj, int base, int option, CTXT_Object *context)
{
    assert(MPZ_Check(obj));

    return mpz_ascii(obj->z, base, option, 0);
}

static MPZ_Object *
GMPy_MPZ_From_Number(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        Py_INCREF(obj);
        return (MPZ_Object*)obj;
    }

    if (PyIntOrLong_Check(obj))
        return GMPy_MPZ_From_PyIntOrLong(obj, context);

    if (MPQ_Check(obj))
        return GMPy_MPZ_From_MPQ((MPQ_Object*)obj, context);

    if (MPFR_Check(obj))
        return GMPy_MPZ_From_MPFR((MPFR_Object*)obj, context);

    if (PyFloat_Check(obj))
        return GMPy_MPZ_From_PyFloat(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_DECIMAL(obj)) {
        PyObject *temp = PyNumber_Long(obj);

        if (temp) {
            result = GMPy_MPZ_From_PyIntOrLong(temp, context);
            Py_DECREF(temp);
        }
        return result;
    }

    if (IS_FRACTION(obj)) {
        MPQ_Object *temp = GMPy_MPQ_From_Fraction(obj, context);

        if (temp) {
            result = GMPy_MPZ_From_MPQ(temp, context);
            Py_DECREF((PyObject*)temp);
        }
        return result;
    }

    TYPE_ERROR("cannot convert object to mpz");
    return result;
}

static MPZ_Object *
GMPy_MPZ_From_Integer(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (MPZ_Check(obj)) {
        Py_INCREF(obj);
        return (MPZ_Object*)obj;
    }

    if (PyIntOrLong_Check(obj))
        return GMPy_MPZ_From_PyIntOrLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPZ_From_XMPZ((XMPZ_Object*)obj, context);

    TYPE_ERROR("cannot convert object to mpz");
    return result;
}

/* str and repr implementations for mpz */
static PyObject *
GMPy_MPZ_Str_Slot(MPZ_Object *self)
{
    /* base-10, no tag */
    return GMPy_PyStr_From_MPZ(self, 10, 0, NULL);
}

static PyObject *
GMPy_MPZ_Repr_Slot(MPZ_Object *self)
{
    /* base-10, with tag */
    return GMPy_PyStr_From_MPZ(self, 10, 1, NULL);
}

/* Helper function for argument parsing. Not currently used. */

#if 0
static int
GMPy_MPZ_convert_arg(PyObject *arg, PyObject **ptr)
{
    MPZ_Object *result = GMPy_MPZ_From_Integer(arg, NULL);

    if (result) {
        *ptr = (PyObject*)result;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpz'");
        return 0;
    }
}
#endif
/* ======================================================================== *
 * Conversion between native Python objects/MPZ and XMPZ.                   *
 * ======================================================================== */

static XMPZ_Object *
GMPy_XMPZ_From_PyIntOrLong(PyObject *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    assert(PyIntOrLong_Check(obj));

    if(!(result = GMPy_XMPZ_New(context)))
        return NULL;

    mpz_set_PyIntOrLong(result->z, obj);
    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    XMPZ_Object *result;

    if (!(result = GMPy_XMPZ_New(context)))
        return NULL;

    if (mpz_set_PyStr(result->z, s, base) == -1) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_PyFloat(PyObject *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    assert(PyFloat_Check(obj));

    if ((result = GMPy_XMPZ_New(context))) {
        double d = PyFloat_AsDouble(obj);

        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        mpz_set_d(result->z, d);
    }
    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    assert(MPZ_Check(obj));

    if ((result = GMPy_XMPZ_New(context)))
        mpz_set(result->z, obj->z);

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    assert(XMPZ_Check(obj));

    if ((result = GMPy_XMPZ_New(context)))
        mpz_set(result->z, obj->z);

    return result;
}

static PyObject *
GMPy_PyStr_From_XMPZ(XMPZ_Object *obj, int base, int option, CTXT_Object *context)
{
    assert(XMPZ_Check(obj));

    return mpz_ascii(obj->z, base, option, 1);
}

static MPZ_Object *
GMPy_MPZ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    assert(XMPZ_Check(obj));

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, obj->z);

    return result;
}

static XMPZ_Object*
GMPy_XMPZ_From_Number_New(PyObject *obj, CTXT_Object *context)
{
    XMPZ_Object *result = NULL;

    if (MPZ_Check(obj))
        return GMPy_XMPZ_From_MPZ((MPZ_Object*)obj, context);

    if (PyIntOrLong_Check(obj))
        return GMPy_XMPZ_From_PyIntOrLong(obj, context);

    if (MPQ_Check(obj))
        return GMPy_XMPZ_From_MPQ((MPQ_Object*)obj, context);

    if (MPFR_Check(obj))
        return GMPy_XMPZ_From_MPFR((MPFR_Object*)obj, context);

    if (PyFloat_Check(obj))
        return GMPy_XMPZ_From_PyFloat(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_XMPZ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_DECIMAL(obj)) {
        PyObject *temp = PyNumber_Long(obj);

        if (temp) {
            result = GMPy_XMPZ_From_PyIntOrLong(temp, context);
            Py_DECREF(temp);
        }
        return result;
    }

    if (IS_FRACTION(obj)) {
        MPQ_Object *temp = GMPy_MPQ_From_Fraction(obj, context);

        if (temp) {
            result = GMPy_XMPZ_From_MPQ(temp, context);
            Py_DECREF((PyObject*)temp);
        }
        return result;
    }

    TYPE_ERROR("cannot convert object to xmpz");
    return result;
}

/* str and repr implementations for xmpz */
static PyObject *
GMPy_XMPZ_Str_Slot(XMPZ_Object *self)
{
    /* base-10, no tag */
    return GMPy_PyStr_From_XMPZ(self, 10, 0, NULL);
}

static PyObject *
GMPy_XMPZ_Repr_Slot(XMPZ_Object *self)
{
    /* base-10, with tag */
    return GMPy_PyStr_From_XMPZ(self, 10, 1, NULL);
}

/* ======================================================================== *
 * Conversion between Integer objects and C types.                          *
 * ======================================================================== */

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

/* ======================================================================== *
 * Conversion between native Python objects/MPZ/XMPZ and MPQ.               *
 * ======================================================================== */

static MPQ_Object *
GMPy_MPQ_From_PyIntOrLong(PyObject *obj, CTXT_Object *context)
{
    MPQ_Object *result;
    MPZ_Object *temp;

    assert(PyIntOrLong_Check(obj));

    temp = GMPy_MPZ_From_PyIntOrLong(obj, context);

    if (!temp)
        return NULL;

    if ((result = GMPy_MPQ_New(context))) {
        mpq_set_z(result->q, temp->z);
        Py_DECREF((PyObject*)temp);
    }
    return result;
}

static MPQ_Object *
GMPy_MPQ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    MPQ_Object *result;
    char *cp;
    Py_ssize_t len, i;
    PyObject *ascii_str = NULL;
    long expt = 0;

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = PyBytes_AsString(s);
    }
    else if (PyUnicode_Check(s)) {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            goto error;
        }
        len = PyBytes_Size(ascii_str);
        cp = PyBytes_AsString(ascii_str);
    }
    else {
        TYPE_ERROR("object is not string or Unicode");
        goto error;
    }

    /* Don't allow NULL characters */
    for (i = 0; i < len; i++) {
        if (cp[i] == '\0') {
            VALUE_ERROR("string contains NULL characters");
            goto error;
        }
    }

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
             * normal. We'll deal with the exponent later.
             */
            *whereexp = '\0';
            expt = atol(whereexp+1);
        }

        /* Handle the case where an embedded decimal point exists. An optional
         * exponent is also allowed. We take advantage of the fact that
         * mpz_set_str() ignores embedded space characters. We count the
         * number of digits after the decimal point and then replace the '.'
         * with a ' '. We've already inserted a NUL byte to terminate the
         * string in the case when an exponent was entered. The value for
         * the exponent has already been read.
         */

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
            if (-1 == mpz_set_str(mpq_numref(result->q), (char*)cp, base)) {
                if (wheredot)
                    *wheredot = '.';
                /* Restore the exponent! */
                if (whereexp && !whereslash && (base == 10))
                    *whereexp = 'E';
                VALUE_ERROR("invalid digits");
                goto error;
            }
            /* Process the exponent. */
            digits = expt - digits;
            if (digits < 0) {
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(-digits));
            }
            else {
                /* Use the denominator instead of a temporary variable. */
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(digits));
                mpz_mul(mpq_numref(result->q), mpq_numref(result->q), mpq_denref(result->q));
                mpz_set_ui(mpq_denref(result->q), 1);
            }
            mpq_canonicalize(result->q);

            /* Restore the decimal point. */
            *wheredot = '.';

            /* Restore the exponent! */
            if (whereexp && (base == 10))
                *whereexp = 'E';

            goto finish;
        }

        if (whereslash)
            *whereslash = '\0';

        /* Read the numerator. */
        if (-1 == mpz_set_str(mpq_numref(result->q), (char*)cp, base)) {
            if (whereslash)
                *whereslash = '/';
            VALUE_ERROR("invalid digits");
            goto error;
        }

        /* If a slash was present, read the denominator. */
        if (whereslash) {
            *whereslash = '/';
            if (-1 == mpz_set_str(mpq_denref(result->q), whereslash+1, base)) {
                VALUE_ERROR("invalid digits");
                goto error;
            }
            if (0 == mpz_sgn(mpq_denref(result->q))) {
                ZERO_ERROR("zero denominator in mpq()");
                goto error;
            }
            mpq_canonicalize(result->q);
        }
        else {
            /* Since no slash was present, either the denominator is 1 or a
             * power of 10.
             */

            if (expt <= 0) {
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(-expt));
            }
            else {
                mpz_ui_pow_ui(mpq_denref(result->q), 10, (unsigned long)(expt));
                mpz_mul(mpq_numref(result->q), mpq_numref(result->q), mpq_denref(result->q));
                mpz_set_ui(mpq_denref(result->q), 1);
            }
            mpq_canonicalize(result->q);
            if (whereexp && (base == 10))
                *whereexp = 'E';
        }
    }

  finish:
    Py_XDECREF(ascii_str);
    return result;

  error:
    Py_DECREF((PyObject*)result);
    Py_XDECREF(ascii_str);
    return NULL;
}

static MPQ_Object *
GMPy_MPQ_From_PyFloat(PyObject *obj, CTXT_Object *context)
{
    MPQ_Object *result;

    assert(PyFloat_Check(obj));

    if ((result = GMPy_MPQ_New(context))) {
        double d = PyFloat_AsDouble(obj);

        if (Py_IS_NAN(d)) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpq' does not support NaN");
            return NULL;
        }
        if (Py_IS_INFINITY(d)) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpq' does not support Infinity");
            return NULL;
        }
        mpq_set_d(result->q, d);
    }

    return result;
}

static MPQ_Object *
GMPy_MPQ_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    MPQ_Object *result;

    assert(MPZ_Check(obj));

    if ((result = GMPy_MPQ_New(context)))
        mpq_set_z(result->q, obj->z);

    return result;
}

static MPQ_Object *
GMPy_MPQ_From_XMPZ(XMPZ_Object *obj, CTXT_Object *context)
{
    MPQ_Object *result;

    assert(XMPZ_Check(obj));

    if ((result = GMPy_MPQ_New(context)))
        mpq_set_z(result->q, obj->z);

    return result;
}

static MPZ_Object *
GMPy_MPZ_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    assert(MPQ_Check(obj));

    if ((result = GMPy_MPZ_New(context)))
        mpz_set_q(result->z, obj->q);

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    XMPZ_Object *result;

    assert(MPQ_Check(obj));

    if ((result = GMPy_XMPZ_New(context)))

        mpz_set_q(result->z, obj->q);

    return result;
}
#ifdef PY2
static PyObject *
GMPy_PyLong_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    PyObject *result;
    MPZ_Object *temp;

    assert(MPQ_Check(obj));

    temp = GMPy_MPZ_From_MPQ(obj, context);

    if (!temp)
        return NULL;

    result = GMPy_PyLong_From_MPZ(temp, context);
    Py_DECREF((PyObject*)temp);

    return result;
}

static PyObject *
GMPy_MPQ_Long_Slot(MPQ_Object *self)
{
    return GMPy_PyLong_From_MPQ(self, NULL);
}
#endif

static PyObject *
GMPy_PyIntOrLong_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    PyObject *result;
    MPZ_Object *temp;

    assert(MPQ_Check(obj));

    temp = GMPy_MPZ_From_MPQ(obj, context);

    if (!temp)
        return NULL;

    result = GMPy_PyIntOrLong_From_MPZ(temp, context);
    Py_DECREF((PyObject*)temp);

    return result;
}

static PyObject *
GMPy_MPQ_Int_Slot(MPQ_Object *self)
{
    return GMPy_PyIntOrLong_From_MPQ(self, NULL);
}

static char* _qtag = "mpq(";

static PyObject *
GMPy_PyStr_From_MPQ(MPQ_Object *obj, int base, int option, CTXT_Object *context)
{
    PyObject *result = NULL, *numstr = NULL, *denstr = NULL;
    char buffer[50], *p;

    numstr = mpz_ascii(mpq_numref(obj->q), base, 0, 0);
    if (!numstr)
        return NULL;

    /* Check if denominator is 1 and no tag is requested. If so, just
     * return the numerator.
     */
    if (!(option & 1) && (0 == mpz_cmp_ui(mpq_denref(obj->q),1)))
        return numstr;

    denstr = mpz_ascii(mpq_denref(obj->q), base, 0, 0);
    if (!denstr) {
        Py_DECREF(numstr);
        return NULL;
    }

    /* Build the format string. */
    p = buffer;
    if (option & 1) {
        strcpy(p, _qtag);
        p += strlen(p);
    }

#ifdef PY2
    *(p++) = '%';
    *(p++) = 's';
    if (!mpz_fits_slong_p(mpq_numref(obj->q)))
        *(p++) = 'L';
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 's';
    if (!mpz_fits_slong_p(mpq_denref(obj->q)))
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

static PyObject *
GMPy_PyFloat_From_MPQ(MPQ_Object *obj, CTXT_Object *context)
{
    double res;

    assert(MPQ_Check(obj));

    res = mpq_get_d(obj->q);

    if (Py_IS_INFINITY(res)) {
        OVERFLOW_ERROR("'mpq' too large to convert to float");
        return NULL;
    }

    return PyFloat_FromDouble(res);
}

static PyObject *
GMPy_MPQ_Float_Slot(MPQ_Object *self)
{
    return GMPy_PyFloat_From_MPQ(self, NULL);
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
GMPy_MPQ_From_DecimalRaw(PyObject* obj, CTXT_Object *context)
{
    MPQ_Object *result;
    PyObject *d_exp, *d_int, *d_sign, *d_is_special;
    mpir_si exp;
    const char *string;

    if (!(result = GMPy_MPQ_New(context)))
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

    if (exp <= 0)
        mpz_ui_pow_ui(mpq_denref(result->q), 10, (mpir_ui)(-exp));
    else {
        mpz_ui_pow_ui(mpq_denref(result->q), 10, (mpir_ui)(exp));
        mpz_mul(mpq_numref(result->q), mpq_numref(result->q), mpq_denref(result->q));
        mpz_set_si(mpq_denref(result->q), 1);
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
GMPy_MPQ_From_DecimalRaw(PyObject* obj, CTXT_Object *context)
{
    MPQ_Object *result;
    PyObject *temp = NULL, *d_is_inf = NULL, *d_is_nan = NULL;
    PyObject *d_is_zero = NULL, *d_is_signed = NULL, *s = NULL;

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;
    mpq_set_si(result->q, 0, 1);

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
            mpq_set_si(result->q, -1, 0);
        }
        else {
            mpq_set_si(result->q, 1, 0);
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
            mpz_set_si(mpq_denref(result->q), -1);
        }
        goto okay;
    }

    Py_DECREF((PyObject*)result);

    s = PyObject_Str(obj);
    if (s) {
        result = GMPy_MPQ_From_PyStr(s, 10, context);
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
GMPy_MPQ_From_Decimal(PyObject* obj, CTXT_Object *context)
{
    MPQ_Object *result;

    if ((result = GMPy_MPQ_From_DecimalRaw(obj, context))) {
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
GMPy_MPQ_From_Fraction(PyObject* obj, CTXT_Object *context)
{
    MPQ_Object *result;
    PyObject *num, *den;

    if (!(result = GMPy_MPQ_New(context)))
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
GMPy_MPQ_From_Number(PyObject *obj, CTXT_Object *context)
{
    if (MPQ_Check(obj)) {
        Py_INCREF(obj);
        return (MPQ_Object*)obj;
    }

    if (MPZ_Check(obj))
        return GMPy_MPQ_From_MPZ((MPZ_Object*)obj, context);

    if (MPFR_Check(obj))
        return GMPy_MPQ_From_MPFR((MPFR_Object*)obj, context);

    if (PyFloat_Check(obj))
        return GMPy_MPQ_From_PyFloat(obj, context);

    if (PyIntOrLong_Check(obj))
        return GMPy_MPQ_From_PyIntOrLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPQ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_DECIMAL(obj))
        return GMPy_MPQ_From_Decimal(obj, context);

    if (IS_FRACTION(obj))
        return GMPy_MPQ_From_Fraction(obj, context);

    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

static MPQ_Object*
GMPy_MPQ_From_Rational(PyObject *obj, CTXT_Object *context)
{
    if (MPQ_Check(obj)) {
        Py_INCREF(obj);
        return (MPQ_Object*)obj;
    }

    if (MPZ_Check(obj))
        return GMPy_MPQ_From_MPZ((MPZ_Object*)obj, context);

    if (PyIntOrLong_Check(obj))
        return GMPy_MPQ_From_PyIntOrLong(obj, context);

    if (XMPZ_Check(obj))
        return GMPy_MPQ_From_XMPZ((XMPZ_Object*)obj, context);

    if (IS_FRACTION(obj))
        return GMPy_MPQ_From_Fraction(obj, context);

    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

/*
 * coerce any number to a mpq
 */

int
GMPy_MPQ_convert_arg(PyObject *arg, PyObject **ptr)
{
    MPQ_Object* result = GMPy_MPQ_From_Number(arg, NULL);

    if (result) {
        *ptr = (PyObject*)result;
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
GMPy_MPQ_Str_Slot(MPQ_Object *self)
{
    /* base-10, no tag */
    return GMPy_PyStr_From_MPQ(self, 10, 0, NULL);
}

static PyObject *
GMPy_MPQ_Repr_Slot(MPQ_Object *self)
{
    /* base-10, with tag */
    return GMPy_PyStr_From_MPQ(self, 10, 1, NULL);
}

