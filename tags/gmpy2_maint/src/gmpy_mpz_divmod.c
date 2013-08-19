/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_divmod.c                                                       *
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


/* This file contains functions related to division and remainder. */

/*
 **************************************************************************
 * Ceiling division and remainder.
 **************************************************************************
 */

PyDoc_STRVAR(doc_gmpy_c_divmod,
"c_divmod(x, y) -> (quotient, remainder)\n\n"
"Return the quotient and remainder of x divided by y. The quotient\n"
"is rounded towards +Inf (ceiling rounding) and the remainder will\n"
"have the opposite sign of y. x and y must be integers.");

static PyObject *
Pygmpy_c_divmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *q, *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_divmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE_TWO_MPZ_TUPLE(q, r, result);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("c_divmod() division by 0");
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_qr(q->z, r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("c_divmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("c_divmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_qr(q->z, r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

PyDoc_STRVAR(doc_gmpy_c_div,
"c_div(x, y) -> quotient\n\n"
"Return the quotient of x divided by y. The quotient is rounded\n"
"towards +Inf (ceiling rounding). x and y must be integers.");

static PyObject *
Pygmpy_c_div(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *q, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_div() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(q = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("c_div() division by 0");
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_cdiv_q(q->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("c_div() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("c_div() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_cdiv_q(q->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)q;
}

PyDoc_STRVAR(doc_gmpy_c_mod,
"c_mod(x, y) -> remainder\n\n"
"Return the remainder of x divided by y. The remainder will have\n"
"the opposite sign of y. x and y must be integers.");

static PyObject *
Pygmpy_c_mod(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_mod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(r = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("c_mod() division by 0");
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_cdiv_r(r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("c_mod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("c_mod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_cdiv_r(r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)r;
}

/*
 **************************************************************************
 * Floor division and remainder.
 **************************************************************************
 */

PyDoc_STRVAR(doc_gmpy_f_divmod,
"f_divmod(x, y) -> (quotient, remainder)\n\n"
"Return the quotient and remainder of x divided by y. The quotient\n"
"is rounded towards -Inf (floor rounding) and the remainder will\n"
"have the same sign as y. x and y must be integers.");

static PyObject *
Pygmpy_f_divmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *q, *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_divmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE_TWO_MPZ_TUPLE(q, r, result);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("f_divmod() division by 0");
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(q->z, r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("f_divmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("f_divmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(q->z, r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

PyDoc_STRVAR(doc_gmpy_f_div,
"f_div(x, y) -> quotient\n\n"
"Return the quotient of x divided by y. The quotient is rounded\n"
"towards -Inf (floor rounding). x and y must be integers.");

static PyObject *
Pygmpy_f_div(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *q, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_div() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(q = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("f_div() division by 0");
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_fdiv_q(q->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("f_div() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("f_div() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_fdiv_q(q->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)q;
}

PyDoc_STRVAR(doc_gmpy_f_mod,
"f_mod(x, y) -> remainder\n\n"
"Return the remainder of x divided by y. The remainder will have\n"
"the same sign as y. x and y must be integers.");

static PyObject *
Pygmpy_f_mod(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_mod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(r = (PympzObject*)Pympz_new()))
        return NULL;

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("f_mod() division by 0");
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_fdiv_r(r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("f_mod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("f_mod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_fdiv_r(r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)r;
}

/*
 **************************************************************************
 * Truncating division and remainder.
 **************************************************************************
 */

PyDoc_STRVAR(doc_gmpy_t_divmod,
"t_divmod(x, y) -> (quotient, remainder)\n\n"
"Return the quotient and remainder of x divided by y. The quotient\n"
"is rounded towards zero (truncation) and the remainder will have\n"
"the same sign as x. x and y must be integers.");

static PyObject *
Pygmpy_t_divmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *q, *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_divmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE_TWO_MPZ_TUPLE(q, r, result);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("t_divmod() division by 0");
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_qr(q->z, r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("t_divmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("t_divmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_qr(q->z, r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

PyDoc_STRVAR(doc_gmpy_t_div,
"t_div(x, y) -> quotient\n\n"
"Return the quotient of x divided by y. The quotient is rounded\n"
"towards 0. x and y must be integers.");

static PyObject *
Pygmpy_t_div(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *q, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_div() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(q = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("t_div() division by 0");
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_tdiv_q(q->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("t_div() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("t_div() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_tdiv_q(q->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)q;
}

PyDoc_STRVAR(doc_gmpy_t_mod,
"t_mod(x, y) -> remainder\n\n"
"Return the remainder of x divided by y. The remainder will have\n"
"the same sign as x. x and y must be integers.");

static PyObject *
Pygmpy_t_mod(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_mod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(r = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("t_mod() division by 0");
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_tdiv_r(r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("t_mod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("t_mod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_tdiv_r(r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)r;
}
