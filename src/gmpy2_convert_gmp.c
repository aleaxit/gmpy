/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_gmp.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018, 2019 Case Van Horsen                  *
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
    int negative;
    Py_ssize_t len;
    PyLongObject *templong = (PyLongObject*)obj;

    assert(PyIntOrLong_Check(obj));

    if(!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

#ifdef PY2
    if (PyInt_Check(obj)) {
        mpz_set_si(result->z, PyInt_AS_LONG(obj));
        return result;
    }
#endif

    switch (Py_SIZE(templong)) {
    case -1:
        mpz_set_si(result->z, -(sdigit)templong->ob_digit[0]);
        break;
    case 0:
        mpz_set_si(result->z, 0);
        break;
    case 1:
        mpz_set_si(result->z, templong->ob_digit[0]);
        break;
    default:
        mpz_set_si(result->z, 0);

        if (Py_SIZE(templong) < 0) {
            len = - Py_SIZE(templong);
            negative = 1;
        } else {
            len = Py_SIZE(templong);
            negative = 0;
        }

        mpz_import(result->z, len, -1, sizeof(templong->ob_digit[0]), 0,
                   sizeof(templong->ob_digit[0])*8 - PyLong_SHIFT, templong->ob_digit);

        if (negative) {
            mpz_neg(result->z, result->z);
        }
    }
    return result;
}

/* To support creation of temporary mpz objects. */
static void
mpz_set_PyIntOrLong(mpz_t z, PyObject *obj)
{
    int negative;
    Py_ssize_t len;
    PyLongObject *templong = (PyLongObject*)obj;

#ifdef PY2
    if (PyInt_Check(obj)) {
        mpz_set_si(z, PyInt_AS_LONG(obj));
        return;
    }
#endif

    switch (Py_SIZE(templong)) {
    case -1:
        mpz_set_si(z, -(sdigit)templong->ob_digit[0]);
        break;
    case 0:
        mpz_set_si(z, 0);
        break;
    case 1:
        mpz_set_si(z, templong->ob_digit[0]);
        break;
    default:
        mpz_set_si(z, 0);

        if (Py_SIZE(templong) < 0) {
            len = - Py_SIZE(templong);
            negative = 1;
        } else {
            len = Py_SIZE(templong);
            negative = 0;
        }

        mpz_import(z, len, -1, sizeof(templong->ob_digit[0]), 0,
                   sizeof(templong->ob_digit[0])*8 - PyLong_SHIFT, templong->ob_digit);

        if (negative) {
            mpz_neg(z, z);
        }
    }
    return;
}

static MPZ_Object *
GMPy_MPZ_From_PyStr(PyObject *s, int base, CTXT_Object *context)
{
    MPZ_Object *result;

    if (!(result = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

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

static PyObject *
GMPy_PyLong_From_MPZ(MPZ_Object *obj, CTXT_Object *context)
{
    int negative;
    size_t count, size;
    PyLongObject *result;

    assert(CHECK_MPZANY(obj));

    /* Assume gmp uses limbs as least as large as the builtin longs do */
    assert(mp_bits_per_limb >= PyLong_SHIFT);

    if (mpz_sgn(obj->z) < 0) {
        negative = 1;
    } else {
        negative = 0;
    }

    size = (mpz_sizeinbase(obj->z, 2) + PyLong_SHIFT - 1) / PyLong_SHIFT;

    if (!(result = _PyLong_New(size))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_export(result->ob_digit, &count, -1, sizeof(result->ob_digit[0]), 0,
               sizeof(result->ob_digit[0])*8 - PyLong_SHIFT, obj->z);

    if (count == 0) {
        result->ob_digit[0] = 0;
    }

    /* long_normalize() is file-static so we must reimplement it */
    /* longobjp = long_normalize(longobjp); */
    while ((size>0) && (result->ob_digit[size-1] == 0)) {
        size--;
    }
    Py_SIZE(result) = size;

    if (negative) {
        Py_SIZE(result) = - Py_SIZE(result);
    }
    return (PyObject*)result;
}

#ifdef PY2
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
    assert(CHECK_MPZANY(obj));

#ifdef PY2
    if (mpz_fits_slong_p(obj->z)) {
        /* cast is safe since we know it fits in a signed long */
        return PyInt_FromLong((long)mpz_get_si(obj->z));
    }
#endif

    return GMPy_PyLong_From_MPZ(obj, context);
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

    assert(CHECK_MPZANY(obj));

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
    assert(CHECK_MPZANY(obj));

    return mpz_ascii(obj->z, base, option, 0);
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

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        result = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (result != NULL && MPZ_Check(result)) {
            return result;
        }
        else {
            Py_XDECREF((PyObject*)result);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpz");
    return NULL;
}



static MPZ_Object *
GMPy_MPZ_From_IntegerAndCopy(PyObject *obj, CTXT_Object *context)
{
    MPZ_Object *result = NULL, *temp = NULL;

    result = GMPy_MPZ_From_Integer(obj, context);

    if (result == NULL)
        return result;

    if (Py_REFCNT(result) == 1)
        return result;

    if (!(temp = GMPy_MPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_set(temp->z, result->z);
    Py_DECREF((PyObject*)result);
    return temp;
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

#ifdef SHARED
/* Helper function for argument parsing. Not used in static build. */

static int
GMPy_MPZ_ConvertArg(PyObject *arg, PyObject **ptr)
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
    int negative;
    Py_ssize_t len;
    PyLongObject *templong = (PyLongObject*)obj;

    assert(PyIntOrLong_Check(obj));

    if(!(result = GMPy_XMPZ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

#ifdef PY2
    if (PyInt_Check(obj)) {
        mpz_set_si(result->z, PyInt_AS_LONG(obj));
        return result;
    }
#endif

    switch (Py_SIZE(templong)) {
    case -1:
        mpz_set_si(result->z, -(sdigit)templong->ob_digit[0]);
        break;
    case 0:
        mpz_set_si(result->z, 0);
        break;
    case 1:
        mpz_set_si(result->z, templong->ob_digit[0]);
        break;
    default:
        mpz_set_si(result->z, 0);

        if (Py_SIZE(templong) < 0) {
            len = - Py_SIZE(templong);
            negative = 1;
        } else {
            len = Py_SIZE(templong);
            negative = 0;
        }

        mpz_import(result->z, len, -1, sizeof(templong->ob_digit[0]), 0,
                   sizeof(templong->ob_digit[0])*8 - PyLong_SHIFT, templong->ob_digit);

        if (negative) {
            mpz_neg(result->z, result->z);
        }
    }
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

    assert(CHECK_MPZANY(obj));

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
    char exp_char = 'E';
    Py_ssize_t len, i;
    PyObject *ascii_str = NULL;
    long expt = 0;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

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

        if (!whereexp) {
            whereexp = strchr((char*)cp, 'e');
            exp_char = 'e';
        }

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
                    *whereexp = exp_char;
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
                *whereexp = exp_char;

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
                *whereexp = exp_char;
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

    assert(CHECK_MPZANY(obj));

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

    if (!temp) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

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

    if (!temp) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

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
    if (!numstr) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    /* Check if denominator is 1 and no tag is requested. If so, just
     * return the numerator.
     */
    if (!(option & 1) && (0 == mpz_cmp_ui(mpq_denref(obj->q),1)))
        return numstr;

    denstr = mpz_ascii(mpq_denref(obj->q), base, 0, 0);
    if (!denstr) {
        /* LCOV_EXCL_START */
        Py_DECREF(numstr);
        return NULL;
        /* LCOV_EXCL_STOP */
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
    if (option & 1)
        *(p++) = ',';
    else
        *(p++) = '/';
    *(p++) = '%';
    *(p++) = 's';
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

static MPQ_Object*
GMPy_MPQ_From_Fraction(PyObject* obj, CTXT_Object *context)
{
    MPQ_Object *result;
    PyObject *num, *den;

    if (!(result = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }
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

    if (IS_FRACTION(obj))
        return GMPy_MPQ_From_Fraction(obj, context);

    if (HAS_MPQ_CONVERSION(obj)) {
        MPQ_Object * res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (HAS_MPZ_CONVERSION(obj)) {
        MPZ_Object * res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPQ_Object * temp = GMPy_MPQ_From_MPZ(res, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

static MPQ_Object*
GMPy_MPQ_From_RationalAndCopy(PyObject *obj, CTXT_Object *context)
{
    MPQ_Object *result = NULL, *temp = NULL;

    result = GMPy_MPQ_From_Rational(obj, context);

    if (result == NULL)
        return result;

    if (Py_REFCNT(result) == 1)
        return result;

    if (!(temp = GMPy_MPQ_New(context))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpq_set(temp->q, result->q);
    Py_DECREF((PyObject*)result);
    return temp;
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

    if (HAS_MPQ_CONVERSION(obj)) {
        MPQ_Object * res = (MPQ_Object *) PyObject_CallMethod(obj, "__mpq__", NULL);

        if (res != NULL && MPQ_Check(res)) {
            return res;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

    if (HAS_MPZ_CONVERSION(obj)) {
        MPZ_Object * res = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (res != NULL && MPZ_Check(res)) {
            MPQ_Object * temp = GMPy_MPQ_From_MPZ(res, context);
            Py_DECREF(res);
            return temp;
        }
        else {
            Py_XDECREF((PyObject*)res);
            goto error;
        }
    }

  error:
    TYPE_ERROR("cannot convert object to mpq");
    return NULL;
}

/*
 * coerce any number to a mpq
 */

#ifdef SHARED
/* Helper function for argument parsing. Not used in static build. */

int
GMPy_MPQ_ConvertArg(PyObject *arg, PyObject **ptr)
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

#endif

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

