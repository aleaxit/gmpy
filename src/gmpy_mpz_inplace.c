/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_inplace.c                                                      *
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


/* Provides inplace operations for mpz. */

#include <math.h>

static PyObject *
Pympz_inplace_add(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(other)) {
        mpz_add(rz->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, other);
            mpz_add(rz->z, Pympz_AS_MPZ(self), tempz);
            mpz_cloc(tempz);
        }
        else if(temp_si >= 0) {
            mpz_add_ui(rz->z, Pympz_AS_MPZ(self), temp_si);
        }
        else {
            mpz_sub_ui(rz->z, Pympz_AS_MPZ(self), -temp_si);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_inplace_sub(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(other)) {
        mpz_sub(rz->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, other);
            mpz_sub(rz->z, Pympz_AS_MPZ(self), tempz);
            mpz_cloc(tempz);
        }
        else if(temp_si >= 0) {
            mpz_sub_ui(rz->z, Pympz_AS_MPZ(self), temp_si);
        }
        else {
            mpz_add_ui(rz->z, Pympz_AS_MPZ(self), -temp_si);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_inplace_mul(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(other)) {
        mpz_mul(rz->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, other);
            mpz_mul(rz->z, Pympz_AS_MPZ(self), tempz);
            mpz_cloc(tempz);
        }
        else {
            mpz_mul_si(rz->z, Pympz_AS_MPZ(self), temp_si);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
Pympz_inplace_floordiv(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        mpz_fdiv_q(rz->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, other);
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(self), tempz);
            mpz_cloc(tempz);
        }
        else if(temp_si == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        else if(temp_si > 0) {
            mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(self), temp_si);
        }
        else {
            mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(self), -temp_si);
            mpz_neg(rz->z, rz->z);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_inplace_rem(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

     if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("mpz modulo by zero");
            return NULL;
        }
        mpz_fdiv_r(rz->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
        return (PyObject*)rz;
    }

   if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, other);
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(self), tempz);
            mpz_cloc(tempz);
        }
        else if(temp_si > 0) {
            mpz_fdiv_r_ui(rz->z, Pympz_AS_MPZ(self), temp_si);
        }
        else if(temp_si == 0) {
            ZERO_ERROR("mpz modulo by zero");
            return NULL;
        }
        else {
            mpz_cdiv_r_ui(rz->z, Pympz_AS_MPZ(self), -temp_si);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_inplace_rshift(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) < 0) {
            VALUE_ERROR("negative shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
        if (!mpz_fits_si_p(Pympz_AS_MPZ(other))) {
            OVERFLOW_ERROR("outrageous shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
        temp_si = mpz_get_si(Pympz_AS_MPZ(other));
        mpz_fdiv_q_2exp(rz->z, Pympz_AS_MPZ(self), temp_si);
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            VALUE_ERROR("outrageous shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
        else if(temp_si >= 0) {
            mpz_fdiv_q_2exp(rz->z, Pympz_AS_MPZ(self), temp_si);
            return (PyObject *)rz;
        }
        else {
            VALUE_ERROR("negative shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_inplace_lshift(PyObject *self, PyObject *other)
{
    PympzObject *rz;
    mpir_si temp_si;
    int overflow;

    if (!(rz = (PympzObject*)Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) < 0) {
            VALUE_ERROR("negative shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
        if (!mpz_fits_si_p(Pympz_AS_MPZ(other))) {
            OVERFLOW_ERROR("outrageous shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
        temp_si = mpz_get_si(Pympz_AS_MPZ(other));
        mpz_mul_2exp(rz->z, Pympz_AS_MPZ(self), temp_si);
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        temp_si = PyLong_AsSIAndOverflow(other, &overflow);
        if (overflow) {
            VALUE_ERROR("outrageous shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
        else if(temp_si >= 0) {
            mpz_mul_2exp(rz->z, Pympz_AS_MPZ(self), temp_si);
            return (PyObject *)rz;
        }
        else {
            VALUE_ERROR("negative shift count");
            Py_DECREF((PyObject*)rz);
            return NULL;
        }
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_inplace_pow(PyObject *self, PyObject *other, PyObject *mod)
{
    PympzObject *r, *e = 0;
    mpir_ui el;

    if (mod != Py_None) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    if (!(e = Pympz_From_Integer(other))) {
        PyErr_Clear();
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (mpz_sgn(e->z) < 0) {
        PyErr_Clear();
        Py_DECREF((PyObject*)e);
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (!mpz_fits_ui_p(e->z)) {
        PyErr_Clear();
        Py_DECREF((PyObject*)e);
        Py_RETURN_NOTIMPLEMENTED;
    }
    if (!(r = (PympzObject*)Pympz_new())) {
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    el = mpz_get_ui(e->z);
    mpz_pow_ui(r->z, Pympz_AS_MPZ(self), el);
    Py_DECREF((PyObject*)e);
    return (PyObject*)r;
}

