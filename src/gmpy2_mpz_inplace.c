/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_inplace.c                                                     *
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


/* Provides inplace operations for mpz. */

#include <math.h>

static PyObject *
GMPy_MPZ_IAdd_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;

    if (CHECK_MPZANY(other)) {
        if ((result = GMPy_MPZ_New(NULL))) {
            mpz_add(result->z, MPZ(self), MPZ(other));
        }
        return (PyObject*)result;
    }

    if (PyLong_CheckExact(other)) {
        if ((result = GMPy_MPZ_New(NULL))) {
            switch (Py_SIZE((PyLongObject*)other)) {
            case -1:
                mpz_sub_ui(result->z, MPZ(self), ((PyLongObject*)other)->ob_digit[0]);
                return (PyObject*)result;
            case 0:
            case 1:
                mpz_add_ui(result->z, MPZ(self), ((PyLongObject*)other)->ob_digit[0]);
                return (PyObject*)result;
            default:
                break;
            }
        }
        else {
            return NULL;
        }
    }

    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if ((result = GMPy_MPZ_New(NULL))) {
            if (!error) {
                if (temp >= 0) {
                    mpz_add_ui(result->z, MPZ(self), temp);
                }
                else {
                    mpz_sub_ui(result->z, MPZ(self), -temp);
                }
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, other);
                mpz_add(result->z, MPZ(self), global.tempz);
            }
        }
        return (PyObject*)result;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPZ_ISub_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *rz;

    if (!(rz =  GMPy_MPZ_New(NULL)))
        return NULL;

    if (CHECK_MPZANY(other)) {
        mpz_sub(rz->z, MPZ(self), MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp >= 0) {
                mpz_sub_ui(rz->z, MPZ(self), temp);
            }
            else {
                mpz_add_ui(rz->z, MPZ(self), -temp);
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_sub(rz->z, MPZ(self), global.tempz);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPZ_IMul_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *rz;

    if (!(rz =  GMPy_MPZ_New(NULL)))
        return NULL;

    if (CHECK_MPZANY(other)) {
        mpz_mul(rz->z, MPZ(self), MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            mpz_mul_si(rz->z, MPZ(self), temp);
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_mul(rz->z, MPZ(self), global.tempz);
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
GMPy_MPZ_IFloorDiv_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *rz;

    if (!(rz =  GMPy_MPZ_New(NULL)))
        return NULL;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        mpz_fdiv_q(rz->z, MPZ(self), MPZ(other));
        return (PyObject*)rz;
    }

    if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp == 0) {
                ZERO_ERROR("mpz division by zero");
                return NULL;
            }
            else if(temp > 0) {
                mpz_fdiv_q_ui(rz->z, MPZ(self), temp);
            }
            else {
                mpz_cdiv_q_ui(rz->z, MPZ(self), -temp);
                mpz_neg(rz->z, rz->z);
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_fdiv_q(rz->z, MPZ(self), global.tempz);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPZ_IRem_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *rz;

    if (!(rz =  GMPy_MPZ_New(NULL)))
        return NULL;

     if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) == 0) {
            ZERO_ERROR("mpz modulo by zero");
            return NULL;
        }
        mpz_fdiv_r(rz->z, MPZ(self), MPZ(other));
        return (PyObject*)rz;
    }

   if (PyIntOrLong_Check(other)) {
        int error;
        native_si temp = GMPy_Integer_AsNative_siAndError(other, &error);

        if (!error) {
            if (temp > 0) {
                mpz_fdiv_r_ui(rz->z, MPZ(self), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("mpz modulo by zero");
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(rz->z, MPZ(self), -temp);
            }
        }
        else {
            mpz_set_PyIntOrLong(global.tempz, other);
            mpz_fdiv_r(rz->z, MPZ(self), global.tempz);
        }
        return (PyObject*)rz;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPZ_IRshift_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *rz;
    mp_bitcnt_t shift;

    if (IS_INTEGER(other)) {
        shift = mp_bitcnt_t_From_Integer(other);
        if (shift == (mp_bitcnt_t)(-1) && PyErr_Occurred())
            return NULL;

        if (!(rz =  GMPy_MPZ_New(NULL)))
            return NULL;

        mpz_fdiv_q_2exp(rz->z, MPZ(self), shift);
        return (PyObject *)rz;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPZ_ILshift_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *rz;
    mp_bitcnt_t shift;

    if (IS_INTEGER(other)) {
        shift = mp_bitcnt_t_From_Integer(other);
        if (shift == (mp_bitcnt_t)(-1) && PyErr_Occurred())
            return NULL;

        if (!(rz =  GMPy_MPZ_New(NULL)))
            return NULL;

        mpz_mul_2exp(rz->z, MPZ(self), shift);
        return (PyObject *)rz;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
GMPy_MPZ_IPow_Slot(PyObject *self, PyObject *other, PyObject *mod)
{
    MPZ_Object *r;
    mp_bitcnt_t exp;

    exp = mp_bitcnt_t_From_Integer(other);
    if (exp == (mp_bitcnt_t)(-1) && PyErr_Occurred()) {
        PyErr_Clear();
        Py_RETURN_NOTIMPLEMENTED;
    }

    if (!(r =  GMPy_MPZ_New(NULL)))
        return NULL;

    mpz_pow_ui(r->z, MPZ(self), exp);
    return (PyObject*)r;
}

