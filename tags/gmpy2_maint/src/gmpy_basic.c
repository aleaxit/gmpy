/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_basic.c                                                            *
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

/* Generic arithmetic methods for gmpy types.
 * The routines are highly optimized for performance.
 */

#include <math.h>

/* Generic addition
 *
 * Support addition for gmpy types with automatic conversion of Python types.
 *
 * The following conversion logic is used:
 *  1) 'mpz' combined with an integer type returns an 'mpz'
 *  2) 'mpz' combined with a rational type returns an 'mpq'
 *  3) 'mpz' combined with a floating-point type returns an 'mpf'
 *  4) 'mpq' combined with an integer or rational type returns an 'mpq'
 *  5) 'mpq' combines with a floating-point type returns an 'mpf'
 *  6) Any type combines with 'mpc' returns an 'mpc.
 *
 * The most common inputs are processed as efficiently as possible.
 */

/* To optimize the performance of mpfr op mpfr and mpc op mpc, dedicated
 * functions are provided for the mpfr and mpc types. The optimal path assumes
 * both arguments are of the same type. If not, the general-purpose Pybasic
 * function is called.
 *
 * These functions are located in either gmpy_mpfr.c or gmpy_mpc.c. Ther are
 * call Pympfr_op_fast or Pympc_op_fast.
 */

static PyObject *
Pybasic_add(PyObject *a, PyObject *b)
{
    mpz_t tempz;
#ifdef WITHMPFR
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
#else
    PympzObject *rz = 0;
#endif
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *rf = 0, *paf = 0, *pbf = 0;
#else
    double tempdouble;
#endif
#ifdef WITHMPC
    PympcObject *rc = 0, *pac = 0, *pbc = 0;
#endif
    mpir_si temp_si;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_add(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si >= 0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), -temp_si);
            }
            return (PyObject*)rz;
        }
        if (CHECK_MPZANY(b)) {
            mpz_add(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (CHECK_MPZANY(b)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(a)) {
            temp_si = PyLong_AsSIAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, a);
                mpz_add(rz->z, Pympz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si >0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(b), temp_si);
            }
            else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(b), -temp_si);
            }
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

#ifdef WITHMPFR

    /* Pympfr_add_fast already handles the case when both a and b are valid
     * mpfr instances.
     */

    if (Pympfr_CheckAndExp(a)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_add_z(rf->f, Pympfr_AS_MPFR(a), pbz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbz);
            MPFR_CLEANUP_RF(addition);
        }
        if (isRational(b) || isDecimal(b)) {
            if (!(pbq = Pympq_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Rational or Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_add_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(addition);
        }
        if (PyFloat_Check(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_add_d(rf->f, Pympfr_AS_MPFR(a), PyFloat_AS_DOUBLE(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(addition);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (Pympfr_CheckAndExp(b)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (isInteger(a)) {
            if (!(paz = Pympz_From_Number(a))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_add_z(rf->f, Pympfr_AS_MPFR(b), paz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paz);
            MPFR_CLEANUP_RF(addition);
        }
        if (isRational(a) || isDecimal(a)) {
            if (!(paq = Pympq_From_Number(a))) {
                SYSTEM_ERROR("Can not convert Rational or Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_add_q(rf->f, Pympfr_AS_MPFR(b), paq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paq);
            MPFR_CLEANUP_RF(addition);
        }
        if (PyFloat_Check(a)) {
            mpfr_clear_flags();
            rf->rc = mpfr_add_d(rf->f, Pympfr_AS_MPFR(b), PyFloat_AS_DOUBLE(a),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(addition);
        }
        Py_DECREF((PyObject*)rf);
    }
#endif

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_add(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*)rq;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        rf->rc = mpfr_add(rf->f, paf->f, pbf->f, context->ctx.mpfr_round);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(addition);
    }
#else
    /* Support mpz+float and float+mpz. */
    if (CHECK_MPZANY(a) && PyFloat_Check(b)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(a));
        tempdouble = tempdouble + PyFloat_AsDouble(b);
        return PyFloat_FromDouble(tempdouble);
    }
    if (CHECK_MPZANY(b) && PyFloat_Check(a)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(b));
        tempdouble = PyFloat_AsDouble(a) + tempdouble;
        return PyFloat_FromDouble(tempdouble);
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        pac = Pympc_From_Complex(a, 0, 0);
        pbc = Pympc_From_Complex(b, 0, 0);
        if (!pac || !pbc) {
            SYSTEM_ERROR("Can not convert Complex to 'mpc'");
            Py_XDECREF((PyObject*)pac);
            Py_XDECREF((PyObject*)pbc);
            return NULL;
        }
        if (!(rc = (PympcObject*)Pympc_new(0, 0))) {
            Py_DECREF((PyObject*)pac);
            Py_DECREF((PyObject*)pbc);
            return NULL;
        }
        rc->rc = mpc_add(rc->c, pac->c, pbc->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)pac);
        Py_DECREF((PyObject*)pbc);
        MPC_CLEANUP(rc, "addition");
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

/* Generic Subtraction
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pybasic_sub(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *rf = 0, *paf = 0, *pbf = 0;
#else
    double tempdouble;
#endif
#ifdef WITHMPC
    PympcObject *rc = 0, *pac = 0, *pbc = 0;
#endif
    mpir_si temp_si;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_sub(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si >= 0) {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), -temp_si);
            }
            return (PyObject*)rz;
        }
        if (Pympz_Check(b)) {
            mpz_sub(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (CHECK_MPZANY(b)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(a)) {
            temp_si = PyLong_AsSIAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, a);
                mpz_sub(rz->z, tempz, Pympz_AS_MPZ(b));
                mpz_cloc(tempz);
            }
            else if (temp_si >= 0) {
                mpz_ui_sub(rz->z, temp_si, Pympz_AS_MPZ(b));
            }
            else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(b), -temp_si);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

#ifdef WITHMPFR

    /* Pympfr_sub_fast already handles the case when both a and b are valid
     * mpfr instances.
     */

    if (Pympfr_CheckAndExp(a)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_sub_z(rf->f, Pympfr_AS_MPFR(a), pbz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbz);
            MPFR_CLEANUP_RF(subtraction);
        }
        if (isRational(b)) {
            if (!(pbq = Pympq_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_sub_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(subtraction);
        }
        if (isDecimal(b)) {
            if (!(pbq = Pympq_From_Decimal(b))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_sub_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(subtraction);
        }
        if (PyFloat_Check(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_sub_d(rf->f, Pympfr_AS_MPFR(a), PyFloat_AS_DOUBLE(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(subtraction);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (Pympfr_CheckAndExp(b)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (isInteger(a)) {
            if (!(paz = Pympz_From_Number(a))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_sub_z(rf->f, Pympfr_AS_MPFR(b), paz->z, context->ctx.mpfr_round);
            mpfr_neg(rf->f, rf->f, context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paz);
            MPFR_CLEANUP_RF(subtraction);
        }
        if (isRational(a)) {
            if (!(paq = Pympq_From_Number(a))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_sub_q(rf->f, Pympfr_AS_MPFR(b), paq->q, context->ctx.mpfr_round);
            mpfr_neg(rf->f, rf->f, context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paq);
            MPFR_CLEANUP_RF(subtraction);
        }
        if (isDecimal(a)) {
            if (!(paq = Pympq_From_Decimal(a))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_sub_q(rf->f, Pympfr_AS_MPFR(b), paq->q, context->ctx.mpfr_round);
            mpfr_neg(rf->f, rf->f, context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paq);
            MPFR_CLEANUP_RF(subtraction);
        }
        if (PyFloat_Check(a)) {
            mpfr_clear_flags();
            rf->rc = mpfr_sub_d(rf->f, Pympfr_AS_MPFR(b), PyFloat_AS_DOUBLE(a), context->ctx.mpfr_round);
            mpfr_neg(rf->f, rf->f, context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(subtraction);
        }
        Py_DECREF((PyObject*)rf);
    }
#endif

    if (isInteger(a) && isInteger(b)) {
        paz = Pympz_From_Number(a);
        pbz = Pympz_From_Number(b);
        if (!paz || !pbz) {
            SYSTEM_ERROR("Can not convert Integer to 'mpz'");
            Py_XDECREF((PyObject*)paz);
            Py_XDECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = (PympzObject*)Pympz_new())) {
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_sub(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz);
        Py_DECREF((PyObject*)pbz);
        return (PyObject*)rz;
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_sub(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*)rq;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        rf->rc = mpfr_sub(rf->f, paf->f, pbf->f, context->ctx.mpfr_round);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(subtraction);
    }
#else
    /* Support mpz-float and float-mpz. */
    if (CHECK_MPZANY(a) && PyFloat_Check(b)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(a));
        tempdouble = tempdouble-PyFloat_AsDouble(b);
        return PyFloat_FromDouble(tempdouble);
    }
    if (CHECK_MPZANY(b) && PyFloat_Check(a)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(b));
        tempdouble = PyFloat_AsDouble(a)-tempdouble;
        return PyFloat_FromDouble(tempdouble);
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        pac = Pympc_From_Complex(a, 0, 0);
        pbc = Pympc_From_Complex(b, 0, 0);
        if (!pac || !pbc) {
            SYSTEM_ERROR("Can not convert Complex to 'mpc'");
            Py_XDECREF((PyObject*)pac);
            Py_XDECREF((PyObject*)pbc);
            return NULL;
        }
        if (!(rc = (PympcObject*)Pympc_new(0, 0))) {
            Py_DECREF((PyObject*)pac);
            Py_DECREF((PyObject*)pbc);
            return NULL;
        }
        rc->rc = mpc_sub(rc->c, pac->c, pbc->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)pac);
        Py_DECREF((PyObject*)pbc);
        MPC_CLEANUP(rc, "subtraction");
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

/* Generic Multiplication
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pybasic_mul(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *rf = 0, *paf = 0, *pbf = 0;
#else
    double tempdouble;
#endif
#ifdef WITHMPC
    PympcObject *rc = 0, *pac = 0, *pbc = 0;
#endif
    mpir_si temp_si;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_mul(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            return (PyObject*)rz;
        }
        if (Pympz_Check(b)) {
            mpz_mul(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (CHECK_MPZANY(b)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(a)) {
            temp_si = PyLong_AsSIAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, a);
                mpz_mul(rz->z, Pympz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(b), temp_si);
            }
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

#ifdef WITHMPFR
    if (Pympfr_CheckAndExp(a)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (Pympfr_CheckAndExp(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_mul(rf->f, Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(b),
                              context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_mul_z(rf->f, Pympfr_AS_MPFR(a), pbz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbz);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (isRational(b)) {
            if (!(pbq = Pympq_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_mul_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (isDecimal(b)) {
            if (!(pbq = Pympq_From_Decimal(b))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_mul_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (PyFloat_Check(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_mul_d(rf->f, Pympfr_AS_MPFR(a), PyFloat_AS_DOUBLE(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(multiplication);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (Pympfr_CheckAndExp(b)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (isInteger(a)) {
            if (!(paz = Pympz_From_Number(a))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_mul_z(rf->f, Pympfr_AS_MPFR(b), paz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paz);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (isRational(a)) {
            if (!(paq = Pympq_From_Number(a))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_mul_q(rf->f, Pympfr_AS_MPFR(b), paq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paq);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (isDecimal(a)) {
            if (!(paq = Pympq_From_Decimal(a))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_mul_q(rf->f, Pympfr_AS_MPFR(b), paq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)paq);
            MPFR_CLEANUP_RF(multiplication);
        }
        if (PyFloat_Check(a)) {
            mpfr_clear_flags();
            rf->rc = mpfr_mul_d(rf->f, Pympfr_AS_MPFR(b), PyFloat_AS_DOUBLE(a),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(multiplication);
        }
        Py_DECREF((PyObject*)rf);
    }
#endif

    if (isInteger(a) && isInteger(b)) {
        paz = Pympz_From_Number(a);
        pbz = Pympz_From_Number(b);
        if (!paz || !pbz) {
            SYSTEM_ERROR("Can not convert Integer to 'mpz'");
            Py_XDECREF((PyObject*)paz);
            Py_XDECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = (PympzObject*)Pympz_new())) {
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_mul(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz);
        Py_DECREF((PyObject*)pbz);
        return (PyObject*)rz;
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_mul(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*)rq;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        rf->rc = mpfr_mul(rf->f, paf->f, pbf->f, context->ctx.mpfr_round);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(multiplication);
    }
#else
    /* Support mpz*float and float*mpz. */
    if (CHECK_MPZANY(a) && PyFloat_Check(b)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(a));
        if (Py_IS_INFINITY(tempdouble)) {
            OVERFLOW_ERROR("'mpz' too large to convert to float");
            return NULL;
        }
        tempdouble = tempdouble*PyFloat_AsDouble(b);
        return PyFloat_FromDouble(tempdouble);
    }
    if (CHECK_MPZANY(b) && PyFloat_Check(a)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(b));
        if (Py_IS_INFINITY(tempdouble)) {
            OVERFLOW_ERROR("'mpz' too large to convert to float");
            return NULL;
        }
        tempdouble = PyFloat_AsDouble(a)*tempdouble;
        return PyFloat_FromDouble(tempdouble);
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        pac = Pympc_From_Complex(a, 0, 0);
        pbc = Pympc_From_Complex(b, 0, 0);
        if (!pac || !pbc) {
            SYSTEM_ERROR("Can not convert Complex to 'mpc'");
            Py_XDECREF((PyObject*)pac);
            Py_XDECREF((PyObject*)pbc);
            return NULL;
        }
        if (!(rc = (PympcObject*)Pympc_new(0, 0))) {
            Py_DECREF((PyObject*)pac);
            Py_DECREF((PyObject*)pbc);
            return NULL;
        }
        rc->rc = mpc_mul(rc->c, pac->c, pbc->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)pac);
        Py_DECREF((PyObject*)pbc);
        MPC_CLEANUP(rc, "multiplication");
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpfr when
 * the arguments are mpfr.
 */

static PyObject *
Pybasic_floordiv(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *rf = 0, *paf = 0, *pbf = 0;
#endif
    mpir_si temp_si;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            else if (temp_si > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), -temp_si);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }
        if (CHECK_MPZANY(b)) {
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (CHECK_MPZANY(b)) {
        if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            return NULL;
        }
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(a)) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, a);
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

#ifdef WITHMPFR
    if (Pympfr_CheckAndExp(a)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (Pympfr_CheckAndExp(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_div(rf->f, Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(b),
                              MPFR_RNDD);
            rf->rc = mpfr_floor(rf->f, rf->f);
            MPFR_CLEANUP_RF(division);
        }
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_z(rf->f, Pympfr_AS_MPFR(a), pbz->z,
                                MPFR_RNDD);
            rf->rc = mpfr_floor(rf->f, rf->f);
            Py_DECREF((PyObject*)pbz);
            MPFR_CLEANUP_RF(division);
        }
        if (isRational(b)) {
            if (!(pbq = Pympq_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                MPFR_RNDD);
            rf->rc = mpfr_floor(rf->f, rf->f);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(division);
        }
        if (isDecimal(b)) {
            if (!(pbq = Pympq_From_Decimal(b))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                MPFR_RNDD);
            rf->rc = mpfr_floor(rf->f, rf->f);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(division);
        }
        if (PyFloat_Check(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_div_d(rf->f, Pympfr_AS_MPFR(a), PyFloat_AS_DOUBLE(b),
                                MPFR_RNDD);
            rf->rc = mpfr_floor(rf->f, rf->f);
            MPFR_CLEANUP_RF(division);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (Pympfr_CheckAndExp(b)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;

        /* Need an mpfr_z_div() to provide optimal support for isInteger(), and
         * need an mpfr_q_div() to provide optimal support for isRational() and
         * isDecimal().
         */

        if (PyFloat_Check(a)) {
            mpfr_clear_flags();
            rf->rc = mpfr_d_div(rf->f, PyFloat_AS_DOUBLE(a), Pympfr_AS_MPFR(b),
                                MPFR_RNDD);
            rf->rc = mpfr_floor(rf->f, rf->f);
            MPFR_CLEANUP_RF(division);
        }
        Py_DECREF((PyObject*)rf);
    }
#endif

    if (isInteger(a) && isInteger(b)) {
        paz = Pympz_From_Number(a);
        pbz = Pympz_From_Number(b);
        if (!paz || !pbz) {
            SYSTEM_ERROR("Can not convert Integer to 'mpz'");
            Py_XDECREF((PyObject*)paz);
            Py_XDECREF((PyObject*)pbz);
            return NULL;
        }
        if (mpz_sgn(pbz->z)==0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = (PympzObject*)Pympz_new())) {
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_fdiv_q(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz);
        Py_DECREF((PyObject*)pbz);
        return (PyObject*)rz;
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new()) || !(rz = (PympzObject*)Pympz_new())) {
            Py_XDECREF((PyObject*)rq);
            Py_XDECREF((PyObject*)rz);
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(rz->z, mpq_numref(rq->q), mpq_denref(rq->q));
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        Py_DECREF((PyObject*)rq);
        return (PyObject*)rz;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_XDECREF((PyObject*)rf);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        rf->rc = mpfr_div(rf->f, paf->f, pbf->f, MPFR_RNDD);
        rf->rc = mpfr_floor(rf->f, rf->f);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(division);
    }
#else
    /* Handle mixed mpz and PyFloat floor division. */
    if (CHECK_MPZANY(a) && PyFloat_Check(b)) {
        PyObject *temp = NULL, *result = NULL;

        if ((temp = Pympz_To_PyLong((PympzObject*)a))) {
            result = PyNumber_FloorDivide(temp, b);
            Py_DECREF(temp);
        }
        return result;
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        TYPE_ERROR("can't take floor of complex number.");
        return NULL;
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_truediv follows the / semantics from Python 3.x. The result types
 * are:
 *   mpz  / mpz  -> mpfr
 *   mpq  / mpq  -> mpq
 *   mpfr / mpfr -> mpfr
 *
 * The behavior of mpq now mimics the behavior of fractions.Fraction.
 */

static PyObject *
Pybasic_truediv(PyObject *a, PyObject *b)
{
    PympzObject *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *rf = 0, *paf = 0, *pbf = 0;
#else
    double tempdouble;
#endif
#ifdef WITHMPC
    PympcObject *rc = 0, *pac = 0, *pbc = 0;
#endif
    mpq_t tempq;

#ifdef WITHMPFR
    if (Pympfr_CheckAndExp(a)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (Pympfr_CheckAndExp(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_div(rf->f, Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(b),
                              context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(division);
        }
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_z(rf->f, Pympfr_AS_MPFR(a), pbz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbz);
            MPFR_CLEANUP_RF(division);
        }
        if (isRational(b)) {
            if (!(pbq = Pympq_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(division);
        }
        if (isDecimal(b)) {
            if (!(pbq = Pympq_From_Decimal(b))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(division);
        }
        if (PyFloat_Check(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_div_d(rf->f, Pympfr_AS_MPFR(a), PyFloat_AS_DOUBLE(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(division);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (Pympfr_CheckAndExp(b)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;

        /* Need an mpfr_z_div() to provide optimal support for isInteger(), and
         * need an mpfr_q_div() to provide optimal support for isRational() and
         * isDecimal().
         */

        if (PyFloat_Check(a)) {
            mpfr_clear_flags();
            rf->rc = mpfr_d_div(rf->f, PyFloat_AS_DOUBLE(a), Pympfr_AS_MPFR(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(division);
        }
        Py_DECREF((PyObject*)rf);
    }
#endif

   if (isInteger(a) && isInteger(b)) {
        paz = Pympz_From_Number(a);
        pbz = Pympz_From_Number(b);
        if (!paz || !pbz) {
            SYSTEM_ERROR("Can not convert Integer to 'mpz'");
            Py_XDECREF((PyObject*)paz);
            Py_XDECREF((PyObject*)pbz);
            return NULL;
        }
        if (mpz_sgn(pbz->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
#ifdef WITHMPFR
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
#endif
        mpq_init(tempq);
        mpq_set_num(tempq, paz->z);
        mpq_set_den(tempq, pbz->z);
        mpq_canonicalize(tempq);
#ifdef WITHMPFR
        mpfr_clear_flags();
        rf->rc = mpfr_set_q(rf->f, tempq, context->ctx.mpfr_round);
#else
        tempdouble = mpq_get_d(tempq);
#endif
        mpq_clear(tempq);
        Py_DECREF((PyObject*)paz);
        Py_DECREF((PyObject*)pbz);
#ifdef WITHMPFR
        MPFR_CLEANUP_RF(division);
#else
        return PyFloat_FromDouble(tempdouble);
#endif
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*) rq;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        rf->rc = mpfr_div(rf->f, paf->f, pbf->f, context->ctx.mpfr_round);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(division);
    }
#else
    /* Support mpz/float and float/mpz. */
    if (CHECK_MPZANY(a) && PyFloat_Check(b)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(a));
        if (Py_IS_INFINITY(tempdouble)) {
            OVERFLOW_ERROR("'mpz' too large to convert to float");
            return NULL;
        }
        tempdouble = tempdouble/PyFloat_AsDouble(b);
        return PyFloat_FromDouble(tempdouble);
    }
    if (CHECK_MPZANY(b) && PyFloat_Check(a)) {
        tempdouble = mpz_get_d(Pympz_AS_MPZ(b));
        if (Py_IS_INFINITY(tempdouble)) {
            OVERFLOW_ERROR("'mpz' too large to convert to float");
            return NULL;
        }
        tempdouble = PyFloat_AsDouble(a)/tempdouble;
        return PyFloat_FromDouble(tempdouble);
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        pac = Pympc_From_Complex(a, 0, 0);
        pbc = Pympc_From_Complex(b, 0, 0);
        if (!pac || !pbc) {
            SYSTEM_ERROR("Can not convert Complex to 'mpc'");
            Py_XDECREF((PyObject*)pac);
            Py_XDECREF((PyObject*)pbc);
            return NULL;
        }
        if (MPC_IS_ZERO_P(pbc)) {
            context->ctx.divzero = 1;
            if (context->ctx.trap_divzero) {
                GMPY_DIVZERO("'mpc' division by zero");
                Py_DECREF((PyObject*)pac);
                Py_DECREF((PyObject*)pbc);
                return NULL;
            }
        }
        if (!(rc = (PympcObject*)Pympc_new(0, 0))) {
            Py_DECREF((PyObject*)pac);
            Py_DECREF((PyObject*)pbc);
            return NULL;
        }
        rc->rc = mpc_div(rc->c, pac->c, pbc->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)pac);
        Py_DECREF((PyObject*)pbc);
        MPC_CLEANUP(rc, "division");
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

#ifdef PY2
/* Pympany_div2 follows the conversions rules for Python 2.x. The behavior is
 * a mix of floordiv and truediv. The type conversion behavior is:
 *   mpz  / mpz  -> mpz
 *   mpq  / mpq  -> mpq
 *   mpfr / mpfr -> mpfr
 *
 * A division operator with these properties is not available with Python 3.x.
 */

static PyObject *
Pybasic_div2(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *rf = 0, *paf = 0, *pbf = 0;
#endif
#ifdef WITHMPC
    PympcObject *rc = 0, *pac = 0, *pbc = 0;
#endif
    mpir_si temp_si;
    int overflow;

    /* Use floordiv for integer types. */

    if (CHECK_MPZANY(a)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            else if (temp_si > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), -temp_si);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }
        if (CHECK_MPZANY(b)) {
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (isInteger(a) && isInteger(b)) {
        paz = Pympz_From_Number(a);
        pbz = Pympz_From_Number(b);
        if (!paz || !pbz) {
            SYSTEM_ERROR("Can not convert Integer to 'mpz'");
            Py_XDECREF((PyObject*)paz);
            Py_XDECREF((PyObject*)pbz);
            return NULL;
        }
        if (mpz_sgn(pbz->z)==0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = (PympzObject*)Pympz_new())) {
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_fdiv_q(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz);
        Py_DECREF((PyObject*)pbz);
        return (PyObject*)rz;
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*) rq;
    }

    /* Use truediv for floating-point types. */

#ifdef WITHMPFR
    if (Pympfr_CheckAndExp(a)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;
        if (Pympfr_CheckAndExp(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_div(rf->f, Pympfr_AS_MPFR(a), Pympfr_AS_MPFR(b),
                              context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(division);
        }
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Integer to 'mpz'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_z(rf->f, Pympfr_AS_MPFR(a), pbz->z,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbz);
            MPFR_CLEANUP_RF(division);
        }
        if (isRational(b)) {
            if (!(pbq = Pympq_From_Number(b))) {
                SYSTEM_ERROR("Can not convert Rational to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(division);
        }
        if (isDecimal(b)) {
            if (!(pbq = Pympq_From_Decimal(b))) {
                SYSTEM_ERROR("Can not convert Decimal to 'mpq'");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_clear_flags();
            rf->rc = mpfr_div_q(rf->f, Pympfr_AS_MPFR(a), pbq->q,
                                context->ctx.mpfr_round);
            Py_DECREF((PyObject*)pbq);
            MPFR_CLEANUP_RF(division);
        }
        if (PyFloat_Check(b)) {
            mpfr_clear_flags();
            rf->rc = mpfr_div_d(rf->f, Pympfr_AS_MPFR(a), PyFloat_AS_DOUBLE(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(division);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (Pympfr_CheckAndExp(b)) {
        if (!(rf = (PympfrObject*)Pympfr_new(0)))
            return NULL;

        /* Need an mpfr_z_div() to provide optimal support for isInteger(), and
         * need an mpfr_q_div() to provide optimal support for isRational() and
         * isDecimal().
         */

        if (PyFloat_Check(a)) {
            mpfr_clear_flags();
            rf->rc = mpfr_d_div(rf->f, PyFloat_AS_DOUBLE(a), Pympfr_AS_MPFR(b),
                                context->ctx.mpfr_round);
            MPFR_CLEANUP_RF(division);
        }
        Py_DECREF((PyObject*)rf);
    }

    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        rf->rc = mpfr_div(rf->f, paf->f, pbf->f, context->ctx.mpfr_round);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(division);
    }
#else
    /* Handle mixed mpz and PyFloat true division. */
    if (CHECK_MPZANY(a) && PyFloat_Check(b)) {
        PyObject *temp = NULL, *result = NULL;

        if ((temp = Pympz_To_PyLong((PympzObject*)a))) {
            result = PyNumber_TrueDivide(temp, b);
            Py_DECREF(temp);
        }
        return result;
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        pac = Pympc_From_Complex(a, 0, 0);
        pbc = Pympc_From_Complex(b, 0, 0);
        if (!pac || !pbc) {
            SYSTEM_ERROR("Can not convert Complex to 'mpc'");
            Py_XDECREF((PyObject*)pac);
            Py_XDECREF((PyObject*)pbc);
            return NULL;
        }
        if (MPC_IS_ZERO_P(pbc)) {
            context->ctx.divzero = 1;
            if (context->ctx.trap_divzero) {
                GMPY_DIVZERO("'mpc' division by zero");
                Py_DECREF((PyObject*)pac);
                Py_DECREF((PyObject*)pbc);
                return NULL;
            }
        }
        if (!(rc = (PympcObject*)Pympc_new(0, 0))) {
            Py_DECREF((PyObject*)pac);
            Py_DECREF((PyObject*)pbc);
            return NULL;
        }
        rc->rc = mpc_div(rc->c, pac->c, pbc->c, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)pac);
        Py_DECREF((PyObject*)pbc);
        MPC_CLEANUP(rc, "division");
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}
#endif

/* Pympany_rem follows the % semantics from Python 3.x. The result types
 * are:
 *   mpz  % mpz  -> mpz
 *   mpq  % mpq  -> mpq
 *   mpfr % mpfr -> mpfr
 *
 * The behavior of mpq now mimics the behavior of fractions.Fraction.
 */

static PyObject *
Pybasic_rem(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *qf = 0, *rf = 0, *paf = 0, *pbf = 0;
#endif
    mpir_si temp_si;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!(rz = (PympzObject*)Pympz_new()))
            return NULL;
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_r_ui(rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(rz->z, Pympz_AS_MPZ(a), -temp_si);
            }
            return (PyObject*)rz;
        }
        if (CHECK_MPZANY(b)) {
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (CHECK_MPZANY(b)) {
        if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            return NULL;
        }
        if (!(rz = (PympzObject*)Pympz_new())) return NULL;
        if (PyIntOrLong_Check(a)) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, a);
            mpz_fdiv_r(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }
        Py_DECREF((PyObject*)rz);
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to 'mpq'");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = (PympqObject*)Pympq_new())) {
            Py_XDECREF((PyObject*)rq);
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpz_inoc(tempz);
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(tempz, mpq_numref(rq->q), mpq_denref(rq->q));
        /* Need to calculate paq - rz * pbq */
        mpq_set_z(rq->q, tempz);
        mpq_mul(rq->q, rq->q, pbq->q);
        mpq_sub(rq->q, paq->q, rq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        mpz_cloc(tempz);
        return (PyObject*)rq;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (mpfr_zero_p(pbf->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.trap_divzero) {
                GMPY_DIVZERO("'mpfr' division by zero in modulo");
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)pbf);
                return NULL;
            }
        }
        if (!(rf = (PympfrObject*)Pympfr_new(0)) ||
            !(qf = (PympfrObject*)Pympfr_new(0))) {
            Py_XDECREF((PyObject*)rf);
            Py_XDECREF((PyObject*)qf);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        if (mpfr_nan_p(paf->f) || mpfr_nan_p(pbf->f) || mpfr_inf_p(paf->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.trap_invalid) {
                GMPY_INVALID("'mpfr' invalid operation in modulo");
                Py_DECREF((PyObject*)rf);
                Py_DECREF((PyObject*)qf);
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)pbf);
                return NULL;
            }
            else {
                mpfr_set_nan(rf->f);
            }
        }
        else if (mpfr_inf_p(pbf->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.trap_invalid) {
                GMPY_INVALID("'mpfr' invalid operation in modulo");
                Py_DECREF((PyObject*)rf);
                Py_DECREF((PyObject*)qf);
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)pbf);
                return NULL;
            }
            if (mpfr_signbit(pbf->f)) {
                mpfr_set_inf(rf->f, -1);
            }
            else {
                rf->rc = mpfr_set(rf->f, paf->f, context->ctx.mpfr_round);
            }
        }
        else {
            mpfr_div(qf->f, paf->f, pbf->f, MPFR_RNDD);
            mpfr_floor(qf->f, qf->f);
            rf->rc = mpfr_fms(rf->f, qf->f, pbf->f, paf->f, context->ctx.mpfr_round);
            mpfr_neg(rf->f, rf->f, context->ctx.mpfr_round);
        }
        Py_XDECREF((PyObject*)qf);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        MPFR_CLEANUP_RF(rem);
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        TYPE_ERROR("can't mod complex numbers");
        return NULL;
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_divmod follows the semantics from Python 3.x. The result types
 * are:
 *   divmod(mpz, mpz)   -> (mpz, mpz)
 *   divmod(mpq, mpq)   -> (mpz, mpq)
 *   divmod(mpfr, mpfr) -> (mpfr, mpfr)
 *
 * The behavior of mpq now mimics the behavior of fractions.Fraction.
 */

static PyObject *
Pybasic_divmod(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    mpz_t tempz;
    PympzObject *qz = 0, *rz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
#ifdef WITHMPFR
    PympfrObject *qf = 0, *rf = 0, *paf = 0, *pbf = 0;
#endif
    mpir_si temp_si;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!( r= PyTuple_New(2)) ||
            !(rz = (PympzObject*)Pympz_new()) ||
            !(qz = (PympzObject*)Pympz_new())) {
            Py_XDECREF((PyObject*)rz);
            Py_XDECREF((PyObject*)qz);
            Py_XDECREF(r);
            return NULL;
        }
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                mpz_fdiv_qr(qz->z, rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_qr_ui(qz->z, rz->z, Pympz_AS_MPZ(a), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                Py_DECREF((PyObject*)qz);
                Py_DECREF(r);
                return NULL;
            }
            else {
                mpz_cdiv_qr_ui(qz->z, rz->z, Pympz_AS_MPZ(a), -temp_si);
                mpz_neg(qz->z, qz->z);
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
        if (CHECK_MPZANY(b)) {
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rz);
                Py_DECREF((PyObject*)qz);
                Py_DECREF(r);
                return NULL;
            }
            mpz_fdiv_qr(qz->z, rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
        Py_DECREF((PyObject*)rz);
        Py_DECREF((PyObject*)qz);
        Py_DECREF(r);
    }

    if (CHECK_MPZANY(b) && PyIntOrLong_Check(a)) {
        if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            return NULL;
        }
        if (!(r = PyTuple_New(2)) ||
            !(rz = (PympzObject*)Pympz_new()) ||
            !(qz = (PympzObject*)Pympz_new())) {
            Py_XDECREF((PyObject*)rz);
            Py_XDECREF((PyObject*)qz);
            Py_XDECREF(r);
            return NULL;
        }

        mpz_inoc(tempz);
        mpz_set_PyIntOrLong(tempz, a);
        mpz_fdiv_qr(qz->z, rz->z, tempz, Pympz_AS_MPZ(b));
        mpz_cloc(tempz);
        PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
        PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
        return r;
    }

    if (isRational(a) && isRational(b)) {
        paq = Pympq_From_Number(a);
        pbq = Pympq_From_Number(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert Rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("'mpq' division or modulo by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(r = PyTuple_New(2)) ||
            !(rq = (PympqObject*)Pympq_new()) ||
            !(qz = (PympzObject*)Pympz_new())) {
            Py_XDECREF((PyObject*)r);
            Py_XDECREF((PyObject*)rq);
            Py_XDECREF((PyObject*)qz);
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(qz->z, mpq_numref(rq->q), mpq_denref(rq->q));
        /* Need to calculate paq - rz * pbq */
        mpq_set_z(rq->q, qz->z);
        mpq_mul(rq->q, rq->q, pbq->q);
        mpq_sub(rq->q, paq->q, rq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
        PyTuple_SET_ITEM(r, 1, (PyObject*)rq);
        return r;
    }

#ifdef WITHMPFR
    if (isReal(a) && isReal(b)) {
        paf = Pympfr_From_Real(a, 0);
        pbf = Pympfr_From_Real(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert Real to 'mpfr'");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (mpfr_zero_p(pbf->f)) {
            context->ctx.divzero = 1;
            if (context->ctx.trap_divzero) {
                GMPY_DIVZERO("'mpfr' division by zero in divmod");
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)pbf);
                return NULL;
            }
        }
        if (!(r = PyTuple_New(2)) ||
            !(qf = (PympfrObject*)Pympfr_new(0)) ||
            !(rf = (PympfrObject*)Pympfr_new(0))) {
            Py_XDECREF(r);
            Py_XDECREF((PyObject*)qf);
            Py_XDECREF((PyObject*)rf);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_clear_flags();
        if (mpfr_nan_p(paf->f) || mpfr_nan_p(pbf->f) || mpfr_inf_p(paf->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.trap_invalid) {
                GMPY_INVALID("'mpfr' invalid operation in divmod");
                Py_DECREF(r);
                Py_DECREF((PyObject*)qf);
                Py_DECREF((PyObject*)rf);
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)pbf);
                return NULL;
            }
            else {
                mpfr_set_nan(qf->f);
                mpfr_set_nan(rf->f);
            }
        }
        else if (mpfr_inf_p(pbf->f)) {
            context->ctx.invalid = 1;
            if (context->ctx.trap_invalid) {
                GMPY_INVALID("'mpfr' invalid operation in divmod");
                Py_DECREF(r);
                Py_DECREF((PyObject*)qf);
                Py_DECREF((PyObject*)rf);
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)pbf);
                return NULL;
            }
            else {
                if (mpfr_zero_p(paf->f)) {
                    mpfr_set_zero(qf->f, mpfr_sgn(pbf->f));
                    mpfr_set_zero(rf->f, mpfr_sgn(pbf->f));
                }
                else if ((mpfr_signbit(paf->f)) != (mpfr_signbit(pbf->f))) {
                    mpfr_set_si(qf->f, -1, context->ctx.mpfr_round);
                    mpfr_set_inf(rf->f, mpfr_sgn(pbf->f));
                }
                else {
                    mpfr_set_si(qf->f, 0, context->ctx.mpfr_round);
                    rf->rc = mpfr_set(rf->f, paf->f, context->ctx.mpfr_round);
                }
            }
        }
        else {
            mpfr_div(qf->f, paf->f, pbf->f, MPFR_RNDD);
            mpfr_floor(qf->f, qf->f);
            rf->rc = mpfr_fms(rf->f, qf->f, pbf->f, paf->f, context->ctx.mpfr_round);
            mpfr_neg(rf->f, rf->f, context->ctx.mpfr_round);
        }
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        SUBNORMALIZE(rf);
        SUBNORMALIZE(qf);
        MERGE_FLAGS;
        if (mpfr_underflow_p() && context->ctx.trap_underflow) {
            GMPY_UNDERFLOW("'mpfr' underflow in divmod");
            Py_DECREF((PyObject*)rf);
            Py_DECREF((PyObject*)qf);
            Py_DECREF(r);
            return NULL;
        }
        if (mpfr_overflow_p() && context->ctx.trap_overflow) {
            GMPY_OVERFLOW("'mpfr' overflow in divmod");
            Py_DECREF((PyObject*)rf);
            Py_DECREF((PyObject*)qf);
            Py_DECREF(r);
            return NULL;
        }
        if (mpfr_inexflag_p() && context->ctx.trap_inexact) {
            GMPY_INEXACT("'mpfr' inexact result in divmod");
            Py_DECREF((PyObject*)rf);
            Py_DECREF((PyObject*)qf);
            Py_DECREF(r);
            return NULL;
        }
        PyTuple_SET_ITEM(r, 0, (PyObject*)qf);
        PyTuple_SET_ITEM(r, 1, (PyObject*)rf);
        return r;
    }
#endif

#ifdef WITHMPC
    if (isComplex(a) && isComplex(b)) {
        TYPE_ERROR("can't take floor or mod of complex number.");
        return NULL;
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}
