/* * * * * * * * * * * * * * * * * * * * ** * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_inplace.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018 Case Van Horsen                        *
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


/* Provides inplace mutating operations for xmpz. */

#include <math.h>

/* Inplace xmpz addition. */

static PyObject *
GMPy_XMPZ_IAdd_Slot(PyObject *self, PyObject *other)
{
    /* Try to make mpz + small_int faster */
    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp >= 0) {
                mpz_add_ui(MPZ(self), MPZ(self), temp);
            }
            else {
                mpz_sub_ui(MPZ(self), MPZ(self), -temp);
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_add(MPZ(self), MPZ(self), global.tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (CHECK_MPZANY(other)) {
        mpz_add(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz subtraction.
 */

static PyObject *
GMPy_XMPZ_ISub_Slot(PyObject *self, PyObject *other)
{
    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp >= 0) {
                mpz_sub_ui(MPZ(self), MPZ(self), temp);
            }
            else {
                mpz_add_ui(MPZ(self), MPZ(self), -temp);
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_sub(MPZ(self), MPZ(self), global.tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (CHECK_MPZANY(other)) {
        mpz_sub(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz multiplication.
 */

static PyObject *
GMPy_XMPZ_IMul_Slot(PyObject *self, PyObject *other)
{
    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            mpz_mul_si(MPZ(self), MPZ(self), temp);
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_mul(MPZ(self), MPZ(self), global.tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (CHECK_MPZANY(other)) {
        mpz_mul(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
GMPy_XMPZ_IFloorDiv_Slot(PyObject *self, PyObject *other)
{
    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp == 0) {
                ZERO_ERROR("xmpz division by zero");
                return NULL;
            }
            else if(temp > 0) {
                mpz_fdiv_q_ui(MPZ(self), MPZ(self), temp);
            }
            else {
                mpz_cdiv_q_ui(MPZ(self), MPZ(self), -temp);
                mpz_neg(MPZ(self), MPZ(self));
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_fdiv_q(MPZ(self), MPZ(self), global.tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) == 0) {
            ZERO_ERROR("xmpz division by zero");
            return NULL;
        }
        mpz_fdiv_q(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz remainder.
 */

static PyObject *
GMPy_XMPZ_IRem_Slot(PyObject *self, PyObject *other)
{
    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp > 0) {
                mpz_fdiv_r_ui(MPZ(self), MPZ(self), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("xmpz modulo by zero");
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(MPZ(self), MPZ(self), -temp);
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_fdiv_r(MPZ(self), MPZ(self), global.tempz);
        }
        Py_INCREF(self);
        return self;
    }

    if (CHECK_MPZANY(other)) {
        if(mpz_sgn(MPZ(other)) == 0) {
            ZERO_ERROR("xmpz modulo by zero");
            return NULL;
        }
        mpz_fdiv_r(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz rshift.
 */

static PyObject *
GMPy_XMPZ_IRshift_Slot(PyObject *self, PyObject *other)
{
    mp_bitcnt_t shift;

    if (IS_INTEGER(other)) {
        shift = mp_bitcnt_t_From_Integer(other);
        if (shift == (mp_bitcnt_t)(-1) && PyErr_Occurred())
            return NULL;

        mpz_fdiv_q_2exp(MPZ(self), MPZ(self), shift);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz lshift.
 */

static PyObject *
GMPy_XMPZ_ILshift_Slot(PyObject *self, PyObject *other)
{
    mp_bitcnt_t shift;

    if (IS_INTEGER(other)) {
        shift = mp_bitcnt_t_From_Integer(other);
        if (shift == (mp_bitcnt_t)(-1) && PyErr_Occurred())
            return NULL;

        mpz_mul_2exp(MPZ(self), MPZ(self), shift);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz_pow.
 */

static PyObject *
GMPy_XMPZ_IPow_Slot(PyObject *self, PyObject *other, PyObject *mod)
{
    mp_bitcnt_t exp;

    exp = mp_bitcnt_t_From_Integer(other);
    if (exp == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        PyErr_Clear();
        Py_RETURN_NOTIMPLEMENTED;
    }

    mpz_pow_ui(MPZ(self), MPZ(self), exp);
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

/* Inplace xmpz and.
 */

static PyObject *
GMPy_XMPZ_IAnd_Slot(PyObject *self, PyObject *other)
{
    if (CHECK_MPZANY(other)) {
        mpz_and(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    if (PyIntOrLong_Check(other)) {
        mpz_set_PyIntOrLong(global.tempz, other);
        mpz_and(MPZ(self), MPZ(self), global.tempz);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz xor.
 */

static PyObject *
GMPy_XMPZ_IXor_Slot(PyObject *self, PyObject *other)
{
    if(CHECK_MPZANY(other)) {
        mpz_xor(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    if(PyIntOrLong_Check(other)) {
        mpz_set_PyIntOrLong(global.tempz, other);
        mpz_xor(MPZ(self), MPZ(self), global.tempz);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz or.
 */

static PyObject *
GMPy_XMPZ_IIor_Slot(PyObject *self, PyObject *other)
{
    if(CHECK_MPZANY(other)) {
        mpz_ior(MPZ(self), MPZ(self), MPZ(other));
        Py_INCREF(self);
        return self;
    }

    if(PyIntOrLong_Check(other)) {
        mpz_set_PyIntOrLong(global.tempz, other);
        mpz_ior(MPZ(self), MPZ(self), global.tempz);
        Py_INCREF(self);
        return self;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

