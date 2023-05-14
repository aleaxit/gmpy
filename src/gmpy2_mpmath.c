/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpmath.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2023 Case Van Horsen                                   *
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
mpmath_build_mpf(long sign, MPZ_Object *man, PyObject *exp, mp_bitcnt_t bc)
{
    PyObject *tup, *tsign, *tbc;

    if (!(tup = PyTuple_New(4))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        return NULL;
    }

    if (!(tsign = PyLong_FromLong(sign))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        Py_DECREF(tup);
        return NULL;
    }

    if (!(tbc = PyLong_FromMpBitCnt(bc))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        Py_DECREF(tup);
        Py_DECREF(tsign);
        return NULL;
    }

    PyTuple_SET_ITEM(tup, 0, tsign);
    PyTuple_SET_ITEM(tup, 1, (PyObject*)man);
    PyTuple_SET_ITEM(tup, 2, (exp)?exp:PyLong_FromLong(0));
    PyTuple_SET_ITEM(tup, 3, tbc);
    return tup;
}

static long
mpmath_get_sign(PyObject *x)
{
    /* Assume properly normalized arguments. A value of 0 implies positive, and
     * value of 1 implies negative.
     */

    if (PyLong_Check(x)) {
        return PyLong_AsLong(x);
    }

    if (MPZ_Check(x)) {
        if (mpz_sgn(MPZ(x)) < 0) {
            return 1;
        }
        else {
            return 0;
        }
    }

    TYPE_ERROR("could not convert object to integer");
    return (long)-1;
}

PyDoc_STRVAR(doc_mpmath_normalizeg,
"_mpmath_normalize(...): helper function for mpmath.");

static PyObject *
Pympz_mpmath_normalize_fast(PyObject *self, PyObject *const *args, Py_ssize_t nargs, PyObject *kwnames)
{
    Py_ssize_t sign = 0, zbits, bc = 0, prec = 0, carry = 0;
    mp_bitcnt_t shift = 0;
    PyObject *exp = NULL, *newexp = NULL, *newexp2 = NULL, *tmp = NULL, *rndstr = NULL;
    MPZ_Object *man = NULL, *upper = NULL, *lower = NULL;
    Py_UCS4 rnd = 0;

    if (PyVectorcall_NARGS(nargs) == 6) {
        /* Need better error-checking here. Under Python 3.0, overflow into
           C-long is possible. */
        sign = mpmath_get_sign(args[0]);
        man = (MPZ_Object*)args[1];
        exp = args[2];
        bc = PyLong_AsSsize_t(args[3]);
        prec = PyLong_AsSsize_t(args[4]);
        rndstr = args[5];

        if ((sign == -1) || (bc == -1) || (prec == -1)) {
            TYPE_ERROR("arguments long, MPZ_Object*, PyObject*, long, long, char needed");
            return NULL;
        }
    }
    else {
        TYPE_ERROR("6 arguments required");
        return NULL;
    }

    if (!MPZ_Check(man)) {
        /* Try to convert to an mpz... */
        if (!(man = GMPy_MPZ_From_Integer((PyObject*)man, NULL))) {
            TYPE_ERROR("argument is not an mpz");
            return NULL;
        }
    }

    /* If rndstr really is a string, extract the first character. */
    if (PyUnicode_Check(rndstr)) {
        rnd = PyString_1Char(rndstr);
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

    if (!(upper = GMPy_MPZ_New(NULL)) || !(lower = GMPy_MPZ_New(NULL))) {
        Py_XDECREF((PyObject*)upper);
        Py_XDECREF((PyObject*)lower);
        return NULL;
    }

    if (bc > prec) {
        shift = bc - prec;
        switch (rnd) {
            case (Py_UCS4)'f':
                if(sign) {
                    mpz_cdiv_q_2exp(upper->z, man->z, shift);
                }
                else {
                    mpz_fdiv_q_2exp(upper->z, man->z, shift);
                }
                break;
            case (Py_UCS4)'c':
                if(sign) {
                    mpz_fdiv_q_2exp(upper->z, man->z, shift);
                }
                else {
                    mpz_cdiv_q_2exp(upper->z, man->z, shift);
                }
                break;
            case (Py_UCS4)'d':
                mpz_fdiv_q_2exp(upper->z, man->z, shift);
                break;
            case (Py_UCS4)'u':
                mpz_cdiv_q_2exp(upper->z, man->z, shift);
                break;
            case (Py_UCS4)'n':
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

        if (!(tmp = PyLong_FromUnsignedLong((unsigned long)shift))) {
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

    if (!(tmp = PyLong_FromLong(zbits))) {
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

PyDoc_STRVAR(doc_mpmath_create,
"_mpmath_create(...): helper function for mpmath.");

static PyObject *
Pympz_mpmath_create_fast(PyObject *self, PyObject *const *args, Py_ssize_t nargs, PyObject *kwnames)
{
    long sign, zbits, bc = 0, prec = 0, carry = 0;
    mp_bitcnt_t shift = 0;
    Py_ssize_t n = 0;
    PyObject *exp = NULL, *newexp = NULL, *newexp2 = NULL, *tmp = NULL;
    MPZ_Object *man = NULL, *upper = NULL, *lower = NULL;

    Py_UCS4 rnd = (Py_UCS4)'f';

    n = PyVectorcall_NARGS(nargs);

    if (n < 2) {
        TYPE_ERROR("mpmath_create() expects 'mpz','int'[,'int','str'] arguments");
        return NULL;
    }

    switch (n) {
        case 4:
            rnd = PyString_1Char(args[3]);
        case 3:
            prec = GMPy_Integer_AsLong(args[2]);
            if (prec == -1) {
                VALUE_ERROR("could not convert prec to positive int");
                return NULL;
            }
        case 2:
            exp = args[1];
        case 1:
            man = GMPy_MPZ_From_Integer(args[0], NULL);
            if (!man) {
                TYPE_ERROR("mpmath_create() expects 'mpz','int'[,'int','str'] arguments");
                return NULL;
            }
    }

    /* If the mantissa is 0, return the normalized representation. */
    if (!mpz_sgn(man->z)) {
        return mpmath_build_mpf(0, man, 0, 0);
    }

    upper = GMPy_MPZ_New(NULL);
    lower = GMPy_MPZ_New(NULL);
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

    if (!prec) {
        prec = bc;
    }

    if (bc > prec) {
        shift = bc - prec;
        switch (rnd) {
            case (Py_UCS4)'f':
                if (sign) {
                    mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                }
                else {
                    mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                }
                break;
            case (Py_UCS4)'c':
                if (sign) {
                    mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                }
                else {
                    mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                }
                break;
            case (Py_UCS4)'d':
                mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                break;
            case (Py_UCS4)'u':
                mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                break;
            case (Py_UCS4)'n':
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
                if (carry) {
                    mpz_add_ui(upper->z, upper->z, 1);
                }
        }
        if (!(tmp = PyLong_FromUnsignedLong((unsigned long)shift))) {
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

    if (!(tmp = PyLong_FromLong(zbits))) {
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
