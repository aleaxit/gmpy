/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_utils.c                                                   *
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

/* ======================================================================== *
 * Conversion between Integer objects and C types.                          *
 * ======================================================================== *
 *
 * Optimized routines for converting an Integer object (Python's integer
 * type(s), mpz, plus types defining __mpz__) to various C types.
 *
 * Note: These functions will not set any exceptions!
 *
 */

#include "longintrepr.h"

#define PY_ABS_LONG_MIN     (0-(unsigned long)LONG_MIN)

/* The various ...AndError functions indicate exceptions as follows:
 *   1) If error is -1, then the argument is negative and the result could
 *      not fit into the desired return type.
 *   2) If error is 1, then the argument is positive and the result could
 *      not fit into the desired return type.
 *   3) If error is 2, then the argument is not a recognized type.
 *   4) If error is 0, then a valid result has been returned.
 *
 * The functions special-case the most common paths and then call the
 * corresponding Python functions to support the general case.
 *
 */

/* The following functions test the most common cases as quickly as
 * possible and then defer to PyLong_AsLongAndOverflow().
 */

static long
GMPy_Integer_AsLongAndError(PyObject *obj, int *error)
{
    register PyLongObject *v;
    MPZ_Object *temp_mpz = NULL;
    unsigned long x, prev;
    long res = 0;
    Py_ssize_t i;
    int sign;

    *error = 0;

#ifdef PY2
    if (PyInt_Check(obj)) {
        return PyInt_AS_LONG(obj);
    }
#endif

    if (PyLong_Check(obj)) {
        v = (PyLongObject *)obj;
        i = Py_SIZE(v);

        switch (i) {
        case -1:
            res = -(sdigit)v->ob_digit[0];
            break;
        case 0:
            break;
        case 1:
            res = v->ob_digit[0];
            break;
        default:
            sign = 1;
            x = 0;
            if (i < 0) {
                sign = -1;
                i = -(i);
            }
            while (--i >= 0) {
                prev = x;
                x = (x << PyLong_SHIFT) + v->ob_digit[i];
                if ((x >> PyLong_SHIFT) != prev) {
                    *error = sign;
                    return res;
                }
            }
            /* Haven't lost any bits, but casting to long requires extra care.
             */
            if (x <= (unsigned long)LONG_MAX) {
                res = (long)x * sign;
            }
            else if (sign < 0 && x == PY_ABS_LONG_MIN) {
                res = LONG_MIN;
            }
            else {
                *error = sign;
            }
        }
        return res;
    }

    if (CHECK_MPZANY(obj)) {
        if (mpz_fits_slong_p(MPZ(obj))) {
            res = (long) mpz_get_si(MPZ(obj));
        }
        else {
            *error = mpz_sgn(MPZ(obj));
            res = 0;
        }
        return res;
    }

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        temp_mpz = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (temp_mpz != NULL && MPZ_Check(temp_mpz)) {
            if (mpz_fits_slong_p(MPZ(temp_mpz))) {
                res = (long) mpz_get_si(MPZ(temp_mpz));
            }
            else {
                *error = mpz_sgn(MPZ(temp_mpz));
                res = 0;
            }
        }
        Py_XDECREF((PyObject*)temp_mpz);
        return res;
    }

    *error = 2;
    return 0;
    /* Fall back to using PyLong_AsLongAndOverFlow. This allows gmpy2 to
     * follow the same semantics as Python for conversion.
     */

/*    res = PyLong_AsLongAndOverflow(obj, error);

    if ((res == -1) && PyErr_Occurred()) {
        PyErr_Clear();
        *error = 2;
        return 0;
    }
    else {
        if (*error) {
            return 0;
        }
        else {
            return res;
        }
    }
*/
}

static unsigned long
GMPy_Integer_AsUnsignedLongAndError(PyObject *obj, int *error)
{
    register PyLongObject *v;
    MPZ_Object *temp_mpz = NULL;
    unsigned long x, prev;
    unsigned long res = 0;
    Py_ssize_t i;

    *error = 0;

#ifdef PY2
    if (PyInt_Check(obj)) {
        long temp = PyInt_AS_LONG(obj);
        if (temp < 0) {
            *error = -1;
            return 0;
        }
        else {
            return (unsigned long)temp;
        }
    }
#endif

    if (PyLong_Check(obj)) {
        v = (PyLongObject *)obj;
        i = Py_SIZE(v);

        if (i < 0) {
            *error = -1;
            return res;
        }

        switch (i) {
        case 0:
            break;
        case 1:
            res = v->ob_digit[0];
            break;
        default:
            x = 0;
            while (--i >= 0) {
                prev = x;
                x = (x << PyLong_SHIFT) + v->ob_digit[i];
                if ((x >> PyLong_SHIFT) != prev) {
                    *error = 1;
                    return res;
                }
            }
            res = x;
        }
        return res;
    }

    if (CHECK_MPZANY(obj)) {
        if (mpz_fits_ulong_p(MPZ(obj))) {
            res = (unsigned long) mpz_get_ui(MPZ(obj));
        }
        else {
            *error = mpz_sgn(MPZ(obj));
            res = 0;
        }
        return res;
    }

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        temp_mpz = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (temp_mpz != NULL && MPZ_Check(temp_mpz)) {
            if (mpz_fits_ulong_p(MPZ(temp_mpz))) {
                res = (unsigned long) mpz_get_ui(MPZ(temp_mpz));
            }
            else {
                *error = mpz_sgn(MPZ(temp_mpz));
                res = 0;
            }
        }
        Py_XDECREF((PyObject*)temp_mpz);
        return res;
    }

    *error = 2;
    return 0;
}

static long
c_long_From_Integer(PyObject *obj)
{
    long result;
    int error;

    result = GMPy_Integer_AsLongAndError(obj, &error);
    if (!error) {
        return result;
    }
    else {
        if (error == 2) {
            TYPE_ERROR("could not convert object to integer");
        }
        else {
            OVERFLOW_ERROR("value too large to convert to C long");
        }
        return -1;
    }
}

static unsigned long
c_ulong_From_Integer(PyObject *obj)
{
    unsigned long result;
    int error;

    result = GMPy_Integer_AsUnsignedLongAndError(obj, &error);
    if (!error) {
        return result;
    }
    else {
        if (error == 2) {
            TYPE_ERROR("could not convert object to integer");
        }
        else if (error == 1) {
            OVERFLOW_ERROR("value too large to convert to C unsigned long");
        }
        else if (error < 0) {
            VALUE_ERROR("a non-negative value is required");
        }
        return (unsigned long)(-1);
    }
}

/* The follow code is only used on Windows x64 platform. We check for _WIN64
 * but we then assume the PY_LONG_LONG is defined.
 */

#ifdef _WIN64

#define PY_ABS_LLONG_MIN (0-(unsigned PY_LONG_LONG)PY_LLONG_MIN)

static PY_LONG_LONG
GMPy_Integer_AsLongLongAndError(PyObject *obj, int *error)
{
    register PyLongObject *v;
    MPZ_Object *temp_mpz = NULL;
    unsigned PY_LONG_LONG x, prev;
    PY_LONG_LONG res = 0;
    Py_ssize_t i;
    int sign;

    *error = 0;

#ifdef PY2
    if (PyInt_Check(obj)) {
        return (PY_LONG_LONG)PyInt_AS_LONG(obj);
    }
#endif

    if (PyLong_Check(obj)) {
        v = (PyLongObject *)obj;
        i = Py_SIZE(v);

        switch (i) {
        case -1:
            res = -(sdigit)v->ob_digit[0];
            break;
        case 0:
            break;
        case 1:
            res = v->ob_digit[0];
            break;
        default:
            sign = 1;
            x = 0;
            if (i < 0) {
                sign = -1;
                i = -(i);
            }
            while (--i >= 0) {
                prev = x;
                x = (x << PyLong_SHIFT) + v->ob_digit[i];
                if ((x >> PyLong_SHIFT) != prev) {
                    *error = sign;
                    return res;
                }
            }
            /* Haven't lost any bits, but casting to long requires extra care.
             */
            if (x <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
                res = (PY_LONG_LONG)x * sign;
            }
            else if (sign < 0 && x == PY_ABS_LLONG_MIN) {
                res = PY_LLONG_MIN;
            }
            else {
                *error = sign;
            }
        }
        return res;
    }

    if (CHECK_MPZANY(obj)) {
        sign = mpz_sgn(MPZ(obj));
        if (sign) {
            if (mpz_sizeinbase(MPZ(obj), 256) <= sizeof(x)) {
                x = 0;
                mpz_export(&x, NULL, 1, sizeof(x), 0, 0, MPZ(obj));
            }
            if (x <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
                res = (PY_LONG_LONG)x * sign;
            }
            else if (sign < 0 && x == PY_ABS_LLONG_MIN) {
                res = PY_LLONG_MIN;
            }
            else {
                *error = sign;
            }
        return res;
        }
    }

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        temp_mpz = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (temp_mpz != NULL && MPZ_Check(temp_mpz)) {
            sign = mpz_sgn(MPZ(obj));
            if (sign) {
                if (mpz_sizeinbase(MPZ(obj), 256) <= sizeof(x)) {
                    x = 0;
                    mpz_export(&x, NULL, 1, sizeof(x), 0, 0, MPZ(obj));
                }
                if (x <= (unsigned PY_LONG_LONG)PY_LLONG_MAX) {
                    res = (PY_LONG_LONG)x * sign;
                }
                else if (sign < 0 && x == PY_ABS_LLONG_MIN) {
                    res = PY_LLONG_MIN;
                }
                else {
                    *error = sign;
                }
            }
            Py_XDECREF((PyObject*)temp_mpz);
            return res;
        }
    }

    *error = 2;
    return 0;
}

static unsigned PY_LONG_LONG
GMPy_Integer_AsUnsignedLongLongAndError(PyObject *obj, int *error)
{
    register PyLongObject *v;
    MPZ_Object *temp_mpz = NULL;
    unsigned PY_LONG_LONG x, prev, res = 0;
    Py_ssize_t i;
    int sign;

    *error = 0;

#ifdef PY2
    if (PyInt_Check(obj)) {
        long temp = PyInt_AS_LONG(obj);
        if (temp < 0) {
            *error = -1;
            return res;
        }
        else {
            return (unsigned PY_LONG_LONG)temp;
        }
    }
#endif

    if (PyLong_Check(obj)) {
        v = (PyLongObject *)obj;
        i = Py_SIZE(v);

        if (i < 0) {
            *error = -1;
            return res;
        }

        switch (i) {
        case 0:
            break;
        case 1:
            res = v->ob_digit[0];
            break;
        default:
            x = 0;
            while (--i >= 0) {
                prev = x;
                x = (x << PyLong_SHIFT) + v->ob_digit[i];
                if ((x >> PyLong_SHIFT) != prev) {
                    *error = 1;
                    return res;
                }
            }
            res = x;
        }
        return res;
    }

    if (CHECK_MPZANY(obj)) {
        sign = mpz_sgn(MPZ(obj));
        if (sign < 0) {
            *error = -1;
            return res;
        }
        else if (sign) {
            if (mpz_sizeinbase(MPZ(obj), 256) <= sizeof(res)) {
                mpz_export(&res, NULL, 1, sizeof(x), 0, 0, MPZ(obj));
            }
            else {
                *error = 1;
            }
        return res;
        }
    }

    if (HAS_STRICT_MPZ_CONVERSION(obj)) {
        temp_mpz = (MPZ_Object *) PyObject_CallMethod(obj, "__mpz__", NULL);

        if (temp_mpz != NULL && MPZ_Check(temp_mpz)) {
            sign = mpz_sgn(MPZ(obj));
            if (sign < 0) {
                *error = -1;
            }
            else if (sign) {
                if (mpz_sizeinbase(MPZ(obj), 256) <= sizeof(res)) {
                    mpz_export(&res, NULL, 1, sizeof(x), 0, 0, MPZ(obj));
                }
            }
            else {
                *error = 1;
            }
            Py_XDECREF((PyObject*)temp_mpz);
            return res;
        }
    }

    *error = 2;
    return 0;
}

static  PY_LONG_LONG
c_longlong_From_Integer(PyObject *obj)
{
    PY_LONG_LONG result;
    int error;

    result = GMPy_Integer_AsLongLongAndError(obj, &error);
    if (!error) {
        return result;
    }
    else {
        if (error == 2) {
            TYPE_ERROR("could not convert object to integer");
        }
        else {
            OVERFLOW_ERROR("value too large to convert to C long long");
        }
        return -1;
    }
}

static unsigned PY_LONG_LONG
c_ulonglong_From_Integer(PyObject *obj)
{
    unsigned PY_LONG_LONG result;
    int error;

    result = GMPy_Integer_AsUnsignedLongLongAndError(obj, &error);
    if (!error) {
        return result;
    }
    else {
        if (error == 2) {
            TYPE_ERROR("could not convert object to integer");
        }
        else if (error == 1) {
            OVERFLOW_ERROR("value too large to convert to C unsigned long long");
        }
        else if (error < 0) {
            VALUE_ERROR("a non-negative value is required");
        }
        return (unsigned PY_LONG_LONG)(-1);
    }
}

#endif

