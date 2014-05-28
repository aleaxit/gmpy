/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_divmod.c                                                      *
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


/* This file contains functions related to division and remainder. */

/*
 **************************************************************************
 * Ceiling division and remainder.
 **************************************************************************
 */

PyDoc_STRVAR(doc_c_divmod,
"c_divmod(x, y) -> (quotient, remainder)\n\n"
"Return the quotient and remainder of x divided by y. The quotient\n"
"is rounded towards +Inf (ceiling rounding) and the remainder will\n"
"have the opposite sign of y. x and y must be integers.");

static PyObject *
GMPy_MPZ_c_divmod(PyObject *self, PyObject *args)
{
    PyObject *result;
    MPZ_Object *q, *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_divmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_New(NULL);
    r = GMPy_MPZ_New(NULL);
    result = PyTuple_New(2);
    if (!tempx || !tempy || !q || !r || !result)
        goto err;

    if(mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("c_divmod() division by 0");
        goto err;
    }

    mpz_cdiv_qr(q->z, r->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;

  err:
    Py_XDECREF(result);
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)r);
    return NULL;
}

PyDoc_STRVAR(doc_c_div,
"c_div(x, y) -> quotient\n\n"
"Return the quotient of x divided by y. The quotient is rounded\n"
"towards +Inf (ceiling rounding). x and y must be integers.");

static PyObject *
GMPy_MPZ_c_div(PyObject *self, PyObject *args)
{
    MPZ_Object *q, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_div() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_New(NULL);
    if (!tempx || !tempy || !q)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("c_div() division by 0");
        goto err;
    }

    mpz_cdiv_q(q->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)q;

  err:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)q);
    return NULL;
}

PyDoc_STRVAR(doc_c_mod,
"c_mod(x, y) -> remainder\n\n"
"Return the remainder of x divided by y. The remainder will have\n"
"the opposite sign of y. x and y must be integers.");

static PyObject *
GMPy_MPZ_c_mod(PyObject *self, PyObject *args)
{
    MPZ_Object *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("c_mod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    r = GMPy_MPZ_New(NULL);
    if (!tempx || !tempy || !r)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("c_mod() division by 0");
        goto err;
    }

    mpz_cdiv_r(r->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)r;

  err:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)r);
    return NULL;
}

/*
 **************************************************************************
 * Floor division and remainder.
 **************************************************************************
 */

PyDoc_STRVAR(doc_f_divmod,
"f_divmod(x, y) -> (quotient, remainder)\n\n"
"Return the quotient and remainder of x divided by y. The quotient\n"
"is rounded towards -Inf (floor rounding) and the remainder will\n"
"have the same sign as y. x and y must be integers.");

static PyObject *
GMPy_MPZ_f_divmod(PyObject *self, PyObject *args)
{
    PyObject *result;
    MPZ_Object *q, *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_divmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_New(NULL);
    r = GMPy_MPZ_New(NULL);
    result = PyTuple_New(2);
    if (!tempx || !tempy || !q || !r || !result)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("f_divmod() division by 0");
        goto err;
    }

    mpz_fdiv_qr(q->z, r->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;

  err:
    Py_XDECREF(result);
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)r);
    return NULL;
}

PyDoc_STRVAR(doc_f_div,
"f_div(x, y) -> quotient\n\n"
"Return the quotient of x divided by y. The quotient is rounded\n"
"towards -Inf (floor rounding). x and y must be integers.");

static PyObject *
GMPy_MPZ_f_div(PyObject *self, PyObject *args)
{
    MPZ_Object *q, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_div() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_New(NULL);
    if (!tempx || !tempy || !q)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("f_div() division by 0");
        goto err;
    }

    mpz_fdiv_q(q->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)q;

  err:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)q);
    return NULL;
}

PyDoc_STRVAR(doc_f_mod,
"f_mod(x, y) -> remainder\n\n"
"Return the remainder of x divided by y. The remainder will have\n"
"the same sign as y. x and y must be integers.");

static PyObject *
GMPy_MPZ_f_mod(PyObject *self, PyObject *args)
{
    MPZ_Object *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("f_mod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    r = GMPy_MPZ_New(NULL);
    if (!tempx || !tempy || !r)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("f_mod() division by 0");
        goto err;
    }

    mpz_fdiv_r(r->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)r;

  err:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)r);
    return NULL;
}

/*
 **************************************************************************
 * Truncating division and remainder.
 **************************************************************************
 */

PyDoc_STRVAR(doc_t_divmod,
"t_divmod(x, y) -> (quotient, remainder)\n\n"
"Return the quotient and remainder of x divided by y. The quotient\n"
"is rounded towards zero (truncation) and the remainder will have\n"
"the same sign as x. x and y must be integers.");

static PyObject *
GMPy_MPZ_t_divmod(PyObject *self, PyObject *args)
{
    PyObject *result;
    MPZ_Object *q, *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_divmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_New(NULL);
    r = GMPy_MPZ_New(NULL);
    result = PyTuple_New(2);
    if (!tempx || !tempy || !q || !r || !result)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("t_divmod() division by 0");
        goto err;
    }

    mpz_tdiv_qr(q->z, r->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;

  err:
    Py_XDECREF(result);
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)q);
    Py_XDECREF((PyObject*)r);
    return NULL;
}

PyDoc_STRVAR(doc_t_div,
"t_div(x, y) -> quotient\n\n"
"Return the quotient of x divided by y. The quotient is rounded\n"
"towards 0. x and y must be integers.");

static PyObject *
GMPy_MPZ_t_div(PyObject *self, PyObject *args)
{
    MPZ_Object *q, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_div() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    q = GMPy_MPZ_New(NULL);
    if (!tempx || !tempy || !q)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("t_div() division by 0");
        goto err;
    }

    mpz_tdiv_q(q->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)q;

  err:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)q);
    return NULL;
}

PyDoc_STRVAR(doc_t_mod,
"t_mod(x, y) -> remainder\n\n"
"Return the remainder of x divided by y. The remainder will have\n"
"the same sign as x. x and y must be integers.");

static PyObject *
GMPy_MPZ_t_mod(PyObject *self, PyObject *args)
{
    MPZ_Object *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("t_mod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL);
    tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL);
    r = GMPy_MPZ_New(NULL);
    if (!tempx || !tempy || !r)
        goto err;

    if (mpz_sgn(tempy->z) == 0) {
        ZERO_ERROR("t_mod() division by 0");
        goto err;
    }

    mpz_tdiv_r(r->z, tempx->z, tempy->z);

    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return (PyObject*)r;

  err:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)r);
    return NULL;
}
