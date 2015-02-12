/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpmath.c                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015 Case Van Horsen                                          *
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


/* Internal helper function for mpmath. */

static PyObject *
mpmath_build_mpf(long sign, PympzObject *man, PyObject *exp, mpir_si bc)
{
    PyObject *tup, *tsign, *tbc;

    if (!(tup = PyTuple_New(4))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        return NULL;
    }

    if (!(tsign = PyIntOrLong_FromLong(sign))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        Py_DECREF(tup);
        return NULL;
    }

    if (!(tbc = PyIntOrLong_FromSI(bc))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        Py_DECREF(tup);
        Py_DECREF(tsign);
        return NULL;
    }

    PyTuple_SET_ITEM(tup, 0, tsign);
    PyTuple_SET_ITEM(tup, 1, (PyObject*)man);
    PyTuple_SET_ITEM(tup, 2, (exp)?exp:PyIntOrLong_FromLong(0));
    PyTuple_SET_ITEM(tup, 3, tbc);
    return tup;
}

PyDoc_STRVAR(doc_mpmath_normalizeg,
"_mpmath_normalize(...): helper function for mpmath.");

static PyObject *
Pympz_mpmath_normalize(PyObject *self, PyObject *args)
{
    long sign = 0;
    mpir_si bc = 0, prec = 0, shift, zbits, carry = 0;
    PyObject *exp = 0, *newexp = 0, *newexp2 = 0, *tmp = 0, *rndstr = 0;
    PympzObject *man = 0, *upper = 0, *lower = 0;
    char rnd = 0;

    if (PyTuple_GET_SIZE(args) == 6) {
        /* Need better error-checking here. Under Python 3.0, overflow into
           C-long is possible. */
        sign = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        man = (PympzObject *)PyTuple_GET_ITEM(args, 1);
        exp = PyTuple_GET_ITEM(args, 2);
        bc = SI_From_Integer(PyTuple_GET_ITEM(args, 3));
        prec = SI_From_Integer(PyTuple_GET_ITEM(args, 4));
        rndstr = PyTuple_GET_ITEM(args, 5);
        if (PyErr_Occurred()) {
            TYPE_ERROR("arguments long, PympzObject*, PyObject*, long, long, char needed");
            return NULL;
        }
    }
    else {
        TYPE_ERROR("6 arguments required");
        return NULL;
    }

    if (!Pympz_Check(man)) {
        TYPE_ERROR("argument is not an mpz");
        return NULL;
    }

    /* If rndstr really is a string, extract the first character. */
    if (Py2or3String_Check(rndstr)) {
        rnd = Py2or3String_AsString(rndstr)[0];
    }
    else {
        VALUE_ERROR("invalid rounding mode specified");
        return NULL;
    }

    /* If the mantissa is 0, return the normalized representation. */
    if (!mpz_sgn(man->z)) {
        Py_INCREF((PyObject*)man);
        return mpmath_build_mpf(0, man, 0, 0);
    }

    /* if bc <= prec and the number is odd return it */
    if ((bc <= prec) && mpz_odd_p(man->z)) {
        Py_INCREF((PyObject*)man);
        Py_INCREF((PyObject*)exp);
        return mpmath_build_mpf(sign, man, exp, bc);
    }

    if (!(upper = (PympzObject*)Pympz_new()) ||
        !(lower = (PympzObject*)Pympz_new())) {
        Py_XDECREF((PyObject*)upper);
        Py_XDECREF((PyObject*)lower);
    }

    shift = bc - prec;
    if (shift>0) {
        switch (rnd) {
            case 'f':
                if(sign) {
                    mpz_cdiv_q_2exp(upper->z, man->z, shift);
                }
                else {
                    mpz_fdiv_q_2exp(upper->z, man->z, shift);
                }
                break;
            case 'c':
                if(sign) {
                    mpz_fdiv_q_2exp(upper->z, man->z, shift);
                }
                else {
                    mpz_cdiv_q_2exp(upper->z, man->z, shift);
                }
                break;
            case 'd':
                mpz_fdiv_q_2exp(upper->z, man->z, shift);
                break;
            case 'u':
                mpz_cdiv_q_2exp(upper->z, man->z, shift);
                break;
            case 'n':
            default:
                mpz_tdiv_r_2exp(lower->z, man->z, shift);
                mpz_tdiv_q_2exp(upper->z, man->z, shift);
                if (mpz_sgn(lower->z)) {
                    /* lower is not 0 so it must have at least 1 bit set */
                    if (mpz_sizeinbase(lower->z, 2) == shift) {
                        /* lower is >= 1/2 */
                        if (mpz_scan1(lower->z, 0) == shift-1) {
                            /* lower is exactly 1/2 */
                            if (mpz_odd_p(upper->z))
                                carry = 1;
                        }
                        else {
                            carry = 1;
                        }
                    }
                }
                if (carry)
                    mpz_add_ui(upper->z, upper->z, 1);
        }

        if (!(tmp = PyIntOrLong_FromSI(shift))) {
            Py_DECREF((PyObject*)upper);
            Py_DECREF((PyObject*)lower);
            return NULL;
        }

        if (!(newexp = PyNumber_Add(exp, tmp))) {
            Py_DECREF((PyObject*)upper);
            Py_DECREF((PyObject*)lower);
            Py_DECREF(tmp);
            return NULL;
        }
        Py_DECREF(tmp);
        bc = prec;
    }
    else {
        mpz_set(upper->z, man->z);
        newexp = exp;
        Py_INCREF(newexp);
    }

    /* Strip trailing 0 bits. */
    if ((zbits = mpz_scan1(upper->z, 0)))
        mpz_tdiv_q_2exp(upper->z, upper->z, zbits);

    if (!(tmp = PyIntOrLong_FromSI(zbits))) {
        Py_DECREF((PyObject*)upper);
        Py_DECREF((PyObject*)lower);
        Py_DECREF(newexp);
        return NULL;
    }
    if (!(newexp2 = PyNumber_Add(newexp, tmp))) {
        Py_DECREF((PyObject*)upper);
        Py_DECREF((PyObject*)lower);
        Py_DECREF(tmp);
        Py_DECREF(newexp);
        return NULL;
    }
    Py_DECREF(newexp);
    Py_DECREF(tmp);

    bc -= zbits;
    /* Check if one less than a power of 2 was rounded up. */
    if (!mpz_cmp_ui(upper->z, 1))
        bc = 1;

    Py_DECREF((PyObject*)lower);
    return mpmath_build_mpf(sign, upper, newexp2, bc);
}

PyDoc_STRVAR(doc_mpmath_createg,
"_mpmath_create(...): helper function for mpmath.");

static PyObject *
Pympz_mpmath_create(PyObject *self, PyObject *args)
{
    long sign;
    mpir_si bc, shift, zbits, carry = 0, prec = 0;
    PyObject *exp = 0, *newexp = 0, *newexp2 = 0, *tmp = 0;
    PympzObject *man = 0, *upper = 0, *lower = 0;

    const char *rnd = "f";

    if (PyTuple_GET_SIZE(args) < 2) {
        TYPE_ERROR("mpmath_create() expects 'mpz','int'[,'int','str'] arguments");
        return NULL;
    }

    switch (PyTuple_GET_SIZE(args)) {
        case 4:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 3));
        case 3:
            prec = SI_From_Integer(PyTuple_GET_ITEM(args, 2));
            if (prec == -1 && PyErr_Occurred())
                return NULL;
            prec = ABS(prec);
        case 2:
            exp = PyTuple_GET_ITEM(args, 1);
        case 1:
            man = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
            if (!man) {
                TYPE_ERROR("mpmath_create() expects 'mpz','int'[,'int','str'] arguments");
                return NULL;
            }
    }

    /* If the mantissa is 0, return the normalized representation. */
    if (!mpz_sgn(man->z)) {
        return mpmath_build_mpf(0, man, 0, 0);
    }

    upper = (PympzObject*)Pympz_new();
    lower = (PympzObject*)Pympz_new();
    if (!upper || !lower) {
        Py_DECREF((PyObject*)man);
        Py_XDECREF((PyObject*)upper);
        Py_XDECREF((PyObject*)lower);
        return NULL;
    }

    /* Extract sign, make man positive, and set bit count */
    sign = (mpz_sgn(man->z) == -1);
    mpz_abs(upper->z, man->z);
    bc = mpz_sizeinbase(upper->z, 2);

    if (!prec) prec = bc;

    shift = bc - prec;
    if (shift > 0) {
        switch (rnd[0]) {
            case 'f':
                if(sign) {
                    mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                }
                else {
                    mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                }
                break;
            case 'c':
                if(sign) {
                    mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                }
                else {
                    mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                }
                break;
            case 'd':
                mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                break;
            case 'u':
                mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                break;
            case 'n':
            default:
                mpz_tdiv_r_2exp(lower->z, upper->z, shift);
                mpz_tdiv_q_2exp(upper->z, upper->z, shift);
                if (mpz_sgn(lower->z)) {
                    /* lower is not 0 so it must have at least 1 bit set */
                    if (mpz_sizeinbase(lower->z, 2)==shift) {
                        /* lower is >= 1/2 */
                        if (mpz_scan1(lower->z, 0)==shift-1) {
                            /* lower is exactly 1/2 */
                            if (mpz_odd_p(upper->z))
                                carry = 1;
                        }
                        else {
                            carry = 1;
                        }
                    }
                }
                if (carry)
                    mpz_add_ui(upper->z, upper->z, 1);
        }
        if (!(tmp = PyIntOrLong_FromSI(shift))) {
            Py_DECREF((PyObject*)upper);
            Py_DECREF((PyObject*)lower);
            return NULL;
        }
        if (!(newexp = PyNumber_Add(exp, tmp))) {
            Py_DECREF((PyObject*)man);
            Py_DECREF((PyObject*)upper);
            Py_DECREF((PyObject*)lower);
            Py_DECREF(tmp);
            return NULL;
        }
        Py_DECREF(tmp);
        bc = prec;
    }
    else {
        newexp = exp;
        Py_INCREF(newexp);
    }

    /* Strip trailing 0 bits. */
    if ((zbits = mpz_scan1(upper->z, 0)))
        mpz_tdiv_q_2exp(upper->z, upper->z, zbits);

    if (!(tmp = PyIntOrLong_FromSI(zbits))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF((PyObject*)upper);
        Py_DECREF((PyObject*)lower);
        Py_DECREF(newexp);
        return NULL;
    }
    if (!(newexp2 = PyNumber_Add(newexp, tmp))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF((PyObject*)upper);
        Py_DECREF((PyObject*)lower);
        Py_DECREF(tmp);
        Py_DECREF(newexp);
        return NULL;
    }
    Py_DECREF(newexp);
    Py_DECREF(tmp);

    bc -= zbits;
    /* Check if one less than a power of 2 was rounded up. */
    if (!mpz_cmp_ui(upper->z, 1))
        bc = 1;

    Py_DECREF((PyObject*)lower);
    Py_DECREF((PyObject*)man);
    return mpmath_build_mpf(sign, upper, newexp2, bc);
}
