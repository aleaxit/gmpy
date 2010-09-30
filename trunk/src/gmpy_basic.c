/* gmpy_basic.c
 *
 * Generic methods for gmpy types. The routines are highly optimized.
 *
 * This file should be considered part of gmpy.c
 */

#include <math.h>

/* Generic addition
 *
 * Support addition for gmpy types with automatic conversion of Python types.
 *
 * When adding an 'mpz' and 'xmpz', the result type is determined by checking
 * options.prefer_mutable.
 *
 * The following conversion logic is used:
 *  1) 'mpz' combined with an integer type returns an 'mpz'
 *  2) 'mpz' combined with a rational type returns an 'mpq'
 *  3) 'mpz' combined with a floating-point type returns an 'mpf'
 *  4) 'mpq' combined with an integer or rational type returns an 'mpq'
 *  5) 'mpq' combines with a floating-point type returns an 'mpf'
 *
 * The most common inputs are processed as efficiently as possible.
 */

static PyObject *
Pympany_add(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;

    if (CHECK_MPZANY(a)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* mpz + integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Adding (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_add(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp >= 0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), temp);
            }
            else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rz;
        }

        /* mpz + mpz */
        if (CHECK_MPZANY(b)) {
            TRACE("Adding (mpz,mpz)\n");
            mpz_add(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (CHECK_MPZANY(b)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* integer + mpz */
        if (PyIntOrLong_Check(a)) {
            TRACE("Adding (long,mpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_add(rz->z, Pympz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            }
            else if (temp >=0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(b), temp);
            }
            else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(b), -temp);
            }
            return (PyObject*)rz;
        }
    }

    if (Pympf_Check(a)) {
        TRACE("Adding (mpf,number)\n");
        if (!(rf = Pympf_new(0))) {
            return NULL;
        }
        if (Pympf_Check(b)) {
            mpfr_add(rf->f, Pympf_AS_MPF(a), Pympf_AS_MPF(b), options.rounding);
        }
        else if (isInteger(b)) {
            if (!(pbz = Pympz_From_Integer(b))) {
                SYSTEM_ERROR("Can not convert number to mpz");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_add_z(rf->f, Pympf_AS_MPF(a), pbz->z, options.rounding);
            Py_DECREF((PyObject*)pbz);
            return (PyObject*)rf;
        }
        else if (isRational(b)) {
            if (!(pbq = Pympq_From_Rational(b))) {
                SYSTEM_ERROR("Can not convert number to mpq");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_add_q(rf->f, Pympf_AS_MPF(a), pbq->q, options.rounding);
            Py_DECREF((PyObject*)pbq);
            return (PyObject*)rf;
        }
        else if (PyFloat_Check(b)) {
            mpfr_add_d(rf->f, Pympf_AS_MPF(a), PyFloat_AS_DOUBLE(b), options.rounding);
            return (PyObject*)rf;
        }
    }

    if (Pympf_Check(b)) {
        TRACE("Adding (number,mpf)\n");
        if (!(rf = Pympf_new(0))) {
            return NULL;
        }
        if (isInteger(a)) {
            if (!(paz = Pympz_From_Integer(a))) {
                SYSTEM_ERROR("Can not convert number to mpz");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_add_z(rf->f, Pympf_AS_MPF(b), paz->z, options.rounding);
            Py_DECREF((PyObject*)paz);
            return (PyObject*)rf;
        }
        else if (isRational(a)) {
            if (!(paq = Pympq_From_Rational(a))) {
                SYSTEM_ERROR("Can not convert number to mpq");
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_add_q(rf->f, Pympf_AS_MPF(b), paq->q, options.rounding);
            Py_DECREF((PyObject*)paq);
            return (PyObject*)rf;
        }
        else if (PyFloat_Check(a)) {
            mpfr_add_d(rf->f, Pympf_AS_MPF(b), PyFloat_AS_DOUBLE(a), options.rounding);
            return (PyObject*)rf;
        }
    }

    if (isInteger(a) && isInteger(b)) {
        TRACE("Adding (integer,integer)\n");
        paz = Pympz_From_Integer(a);
        pbz = Pympz_From_Integer(b);
        if (!paz || !pbz) {
            SYSTEM_ERROR("Can not convert integer to mpz");
            Py_XDECREF((PyObject*)paz);
            Py_XDECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz);
            Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_add(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz);
        Py_DECREF((PyObject*)pbz);
        return (PyObject*)rz;
    }

    if (isRational(a) && isRational(b)) {
        TRACE("Adding (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_add(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*)rq;
    }

    if (isFloat(a) && isFloat(b)) {
        TRACE("Adding (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(0))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_add(rf->f, paf->f, pbf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*)rf;
    }

    Py_RETURN_NOTIMPLEMENTED;
}

/* Generic Subtraction
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pympany_sub(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PyxmpzObject *rxz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;

    if (Pympz_Check(a)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* mpz - integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Subtracting (mpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_sub(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp >= 0) {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), temp);
            }
            else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rz;
        }

        /* mpz - mpz */
        if (Pympz_Check(b)) {
            TRACE("Subtracting (mpz,mpz)\n");
            mpz_sub(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }

        /* mpz - xmpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            TRACE("Subtracting (mpz,xmpz)\n");
            mpz_sub(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pympz_Check(b)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* integer - mpz */
        if (PyIntOrLong_Check(a)) {
            TRACE("Subtracting (long,mpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_sub(rz->z, tempz, Pympz_AS_MPZ(b));
                mpz_cloc(tempz);
            }
            else if (temp >= 0) {
                mpz_ui_sub(rz->z, temp, Pympz_AS_MPZ(b));
            }
            else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(b), -temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }

        /* xmpz - mpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            TRACE("Subtracting (xmpz,mpz)\n");
            mpz_sub(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pyxmpz_Check(a)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz - integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Subtracting (xmpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_sub(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp >= 0) {
                mpz_sub_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            else {
                mpz_add_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rxz;
        }

        /* xmpz - xmpz */
        if (Pyxmpz_Check(b)) {
            TRACE("Subtracting (xmpz,xmpz)\n");
            mpz_sub(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }

        /* xmpz - mpz */
        if ((options.prefer_mutable) && Pympz_Check(b)) {
            TRACE("Subtracting (xmpz,mpz)\n");
            mpz_sub(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (Pyxmpz_Check(b)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer - xmpz */
        if (PyIntOrLong_Check(a)) {
            TRACE("Subtracting (long,xmpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_sub(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
                mpz_cloc(tempz);
            }
            else if (temp >= 0) {
                mpz_ui_sub(rxz->z, temp, Pyxmpz_AS_MPZ(b));
            }
            else {
                mpz_add_ui(rxz->z, Pyxmpz_AS_MPZ(b), -temp);
                mpz_neg(rxz->z, rxz->z);
            }
            return (PyObject*)rxz;
        }

        /* mpz - xmpz */
        if ((options.prefer_mutable) && Pympz_Check(a)) {
            TRACE("Subtracting (mpz,xmpz)\n");
            mpz_sub(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (isRational(a) && isRational(b)) {
        TRACE("Subtracting (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_sub(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*)rq;
    }

    if (Pympf_Check(a) || PyFloat_Check(a)) {
        TRACE("Subtracting (number,number)\n");
        if (!(paf = Pympf_From_Float(a, 0))) {
            SYSTEM_ERROR("Can not convert number to mpf");
            return NULL;
        }
        if (!(rf = Pympf_new(0))) {
            Py_DECREF((PyObject*)paf);
            return NULL;
        }
        if (isInteger(b)) {
            if (!(pbz = Pympz_From_Integer(b))) {
                SYSTEM_ERROR("Can not convert number to mpz");
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_sub_z(rf->f, paf->f, pbz->z, options.rounding);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbz);
            return (PyObject*)rf;
        }
        else if (isRational(b)) {
            if (!(pbq = Pympq_From_Rational(b))) {
                SYSTEM_ERROR("Can not convert number to mpq");
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_sub_q(rf->f, paf->f, pbq->q, options.rounding);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbq);
            return (PyObject*)rf;
        }
        else if (PyFloat_Check(b)) {
            mpfr_sub_d(rf->f, paf->f, PyFloat_AS_DOUBLE(b), options.rounding);
            Py_DECREF((PyObject*)paf);
            return (PyObject*)rf;
        }
        else if (isFloat(b)) {
            if (!(pbf = Pympf_From_Float(b, options.precision))) {
                SYSTEM_ERROR("Can not convert number to mpf");
                Py_DECREF((PyObject*)paf);
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_sub(rf->f, paf->f, pbf->f, options.rounding);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return (PyObject*)rf;
        }
    }

    if (Pympf_Check(b) || PyFloat_Check(b)) {
        TRACE("Subtracting (number,number)\n");
        if (!(pbf = Pympf_From_Float(b, 0))) {
            SYSTEM_ERROR("Can not convert number to mpf");
            return NULL;
        }
        if (!(rf = Pympf_new(0))) {
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (isInteger(a)) {
            if (!(paz = Pympz_From_Integer(a))) {
                SYSTEM_ERROR("Can not convert number to mpz");
                Py_DECREF((PyObject*)pbf);
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_sub_z(rf->f, pbf->f, paz->z, options.rounding);
            mpfr_neg(rf->f, rf->f, options.rounding);
            Py_DECREF((PyObject*)pbf);
            Py_DECREF((PyObject*)paz);
            return (PyObject*)rf;
        }
        else if (isRational(a)) {
            if (!(paq = Pympq_From_Rational(a))) {
                SYSTEM_ERROR("Can not convert number to mpq");
                Py_DECREF((PyObject*)pbf);
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_add_q(rf->f, pbf->f, paq->q, options.rounding);
            mpfr_neg(rf->f, rf->f, options.rounding);
            Py_DECREF((PyObject*)pbf);
            Py_DECREF((PyObject*)paq);
            return (PyObject*)rf;
        }
        else if (PyFloat_Check(a)) {
            mpfr_sub_d(rf->f, pbf->f, PyFloat_AS_DOUBLE(a), options.rounding);
            mpfr_neg(rf->f, rf->f, options.rounding);
            Py_DECREF((PyObject*)pbf);
            return (PyObject*)rf;
        }
        else if (isFloat(a)) {
            if (!(paf = Pympf_From_Float(a, options.precision))) {
                SYSTEM_ERROR("Can not convert number to mpf");
                Py_DECREF((PyObject*)pbf);
                Py_DECREF((PyObject*)rf);
                return NULL;
            }
            mpfr_add(rf->f, pbf->f, paf->f, options.rounding);
            mpfr_neg(rf->f, rf->f, options.rounding);
            Py_DECREF((PyObject*)pbf);
            Py_DECREF((PyObject*)paf);
            return (PyObject*)rf;
        }
    }
    Py_RETURN_NOTIMPLEMENTED;
}

/* Generic Multiplication
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pympany_mul(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0;
    PyxmpzObject *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;

    if (Pympz_Check(a)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* mpz * integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Multiplying (mpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_mul(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(a), temp);
            }
            return (PyObject*)rz;
        }

        /* mpz * mpz */
        if (Pympz_Check(b)) {
            TRACE("Multiplying (mpz,mpz)\n");
            mpz_mul(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }

        /* mpz * xmpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            TRACE("Multiplying (mpz,xmpz)\n");
            mpz_mul(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pympz_Check(b)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* integer * mpz */
        if (PyIntOrLong_Check(a)) {
            TRACE("Multiplying (long,mpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_mul(rz->z, Pympz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(b), temp);
            }
            return (PyObject*)rz;
        }

        /* xmpz * mpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            TRACE("Multiplying (mpz,xmpz)\n");
            mpz_mul(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pyxmpz_Check(a)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz * integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Multiplying (xmpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_mul(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            return (PyObject*)rxz;
        }

        /* xmpz * xmpz */
        if (Pyxmpz_Check(b)) {
            TRACE("Multiplying (xmpz,xmpz)\n");
            mpz_mul(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }

        /* xmpz * mpz */
        if ((options.prefer_mutable) && Pympz_Check(b)) {
            TRACE("Multiplying (xmpz,mpz)\n");
            mpz_mul(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (Pyxmpz_Check(b)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer * xmpz */
        if (PyIntOrLong_Check(a)) {
            TRACE("Multiplying (long,xmpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_mul(rxz->z, Pyxmpz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            }
            else {
                mpz_mul_si(rxz->z, Pyxmpz_AS_MPZ(b), temp);
            }
            return (PyObject*)rxz;
        }

        /* mpz * xmpz */
        if ((options.prefer_mutable) && Pympz_Check(a)) {
            TRACE("Multiplying (mpz,xmpz)\n");
            mpz_mul(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (isRational(a) && isRational(b)) {
        TRACE("Multiplying (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_mul(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*)rq;
    }

    if (isFloat(a) && isFloat(b)) {
        TRACE("Multiplying (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(options.precision))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_mul(rf->f, paf->f, pbf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*)rf;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
Pympany_floordiv(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0;
    PyxmpzObject *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;

    if (Pympz_Check(a)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* mpz // integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Floor divide (mpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp == 0) {
                ZERO_ERROR("mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            else if (temp > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
            }
            else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), -temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }

        /* mpz // mpz */
        if (Pympz_Check(b)) {
            TRACE("Floor divide (mpz,mpz)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }

        /* mpz // xmpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            TRACE("Floor divide (mpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pympz_Check(b)) {
        if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        if (!(rz = Pympz_new()))
            return NULL;

        /* integer // mpz */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("Floor divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }
#endif

        if (PyLong_Check(a)) {
            TRACE("Floor divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }

        /* xmpz // mpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            TRACE("Floor divide (xmpz,mpz)\n");
            mpz_fdiv_q(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pyxmpz_Check(a)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz // integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Floor divide (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp == 0) {
                ZERO_ERROR("xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            else if (temp > 0) {
                mpz_fdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            else {
                mpz_cdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
                mpz_neg(rxz->z, rxz->z);
            }
            return (PyObject*)rxz;
        }

        /* xmpz // xmpz */
        if (Pyxmpz_Check(b)) {
            TRACE("Floor divide (xmpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }

        /* xmpz // mpz */
        if ((options.prefer_mutable) && Pympz_Check(b)) {
            TRACE("Floor divide (xmpz,mpz)\n");
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (Pyxmpz_Check(b)) {
        if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("xmpz division by zero");
            return NULL;
        }
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer // xmpz */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("Floor divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }
#endif

        if (PyLong_Check(a)) {
            TRACE("Floor divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }

        /* mpz // xmpz */
        if ((options.prefer_mutable) && Pympz_Check(a)) {
            TRACE("Floor divide (mpz,xmpz)\n");
            mpz_fdiv_q(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (isRational(a) && isRational(b)) {
        TRACE("Floor divide (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("mpq division by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new()) || !(rz = Pympz_new())) {
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

    if (isFloat(a) && isFloat(b)) {
        TRACE("Floor divide (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (mpfr_zero_p(pbf->f)) {
            ZERO_ERROR("mpf division by zero");
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(options.precision))) {
            Py_XDECREF((PyObject*)rf);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_div(rf->f, paf->f, pbf->f, options.rounding);
        mpfr_floor(rf->f, rf->f);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*) rf;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_truediv follows the / semantics from Python 3.x. The result types
 * are:
 *   mpz / mpz -> mpf
 *   mpq / mpq -> mpq
 *   mpf / mpf -> mpf
 *
 * The behavior of mpq now mimics the behavior of fractions.Fraction.
 */

static PyObject *
Pympany_truediv(PyObject *a, PyObject *b)
{
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;

    if (CHECK_MPZANY(b) && (mpz_sgn(Pympz_AS_MPZ(b)) == 0)) {
        ZERO_ERROR("mpz division by zero");
        return NULL;
    }

    if (Pympq_Check(b) && (mpq_sgn(Pympq_AS_MPQ(b)) == 0)) {
        ZERO_ERROR("mpq division by zero");
        return NULL;
    }

    if (isInteger(a) && isInteger(b)) {
        TRACE("True divide (integer,integer)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (mpfr_sgn(pbf->f)==0) {
            ZERO_ERROR("mpz division by zero");
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(options.precision))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_div(rf->f, paf->f, pbf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*) rf;
    }

    if (isRational(a) && isRational(b)) {
        TRACE("True divide (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("mpq division by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq);
        Py_DECREF((PyObject*)pbq);
        return (PyObject*) rq;
    }

    if (isFloat(a) && isFloat(b)) {
        TRACE("True divide (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            if (PyErr_Occurred()) {
                PyErr_Clear();
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                Py_RETURN_NOTIMPLEMENTED;
            }
            else {
                SYSTEM_ERROR("Internal error status is confused.");
                return NULL;
            }
        }
        if (mpfr_zero_p(pbf->f)) {
            ZERO_ERROR("mpf division by zero");
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(options.precision))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_div(rf->f, paf->f, pbf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*) rf;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

#ifdef PY2
/* Pympany_div2 follows the conversions rules for Python 2.x. The behavior is
 * a mix of floordiv and truediv. The type conversion behavior is:
 *   mpz / mpz -> mpz
 *   mpq / mpq -> mpq
 *   mpf / mpf -> mpf
 *
 * A division operator with these properties is not available with Python 3.x.
 */

static PyObject *
Pympany_div2(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    mpz_t tempz;
    PympzObject *rz = 0;
    PyxmpzObject *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;
    int overflow;

    /* Use floordiv for integer types. */

    if (Pympz_Check(a)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* mpz / integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("True divide (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }

        /* mpz / mpz */
        if (Pympz_Check(b)) {
            TRACE("True divide (mpz,mpz)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }

        /* mpz / xmpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            TRACE("True divide (mpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pympz_Check(b)) {
        if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        if (!(rz = Pympz_new()))
            return NULL;

        /* integer / mpz */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("True divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }
#endif
        if (PyLong_Check(a)) {
            TRACE("True divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }

        /* xmpz / mpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            TRACE("True divide (xmpz,mpz)\n");
            mpz_fdiv_q(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pyxmpz_Check(a)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz / integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("True divide (xmpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            else {
                mpz_cdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
                mpz_neg(rxz->z, rxz->z);
            }
            return (PyObject*)rxz;
        }

        /* xmpz / xmpz */
        if (Pyxmpz_Check(b)) {
            TRACE("True divide (xmpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }

        /* xmpz / mpz */
        if ((options.prefer_mutable) && Pympz_Check(b)) {
            TRACE("True divide (xmpz,mpz)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (Pyxmpz_Check(b)) {
        if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("xmpz division by zero");
            return NULL;
        }
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer / xmpz */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("True divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }
#endif
        if (PyLong_Check(a)) {
            TRACE("True divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }

        /* mpz / xmpz */
        if ((options.prefer_mutable) && Pympz_Check(a)) {
            TRACE("True divide (mpz,xmpz)\n");
            mpz_fdiv_q(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    /* Use truediv for rational types. */

    if (isRational(a) && isRational(b)) {
        TRACE("True divide (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("mpq division by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
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

    if (isFloat(a) && isFloat(b)) {
        TRACE("True divide (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((pyObject*)pbf);
            return NULL;
        }
        if (mpfr_zero_p(pbf->f)) {
            ZERO_ERROR("mpf division by zero");
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(options.precision))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_div(rf->f, paf->f, pbf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*) rf;
    }
    Py_RETURN_NOTIMPLEMENTED;
}
#endif

/* Pympany_rem follows the % semantics from Python 3.x. The result types
 * are:
 *   mpz % mpz -> mpz
 *   mpq % mpq -> mpq
 *   mpf % mpf -> mpf
 *
 * The behavior of mpq now mimics the behavior of fractions.Fraction.
 */

static PyObject *
Pympany_rem(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    PympzObject *rz = 0;
    PyxmpzObject *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;

    if (Pympz_Check(a)) {
        if (!(rz = Pympz_new()))
            return NULL;

        /* mpz % integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Modulo (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_r_ui(rz->z, Pympz_AS_MPZ(a), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("mpz modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rz;
        }

        /* mpz % mpz */
        if (Pympz_Check(b)) {
            TRACE("Modulo (integer,integer)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }

        /* mpz % xmpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            TRACE("Modulo (mpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pympz_Check(b)) {
        if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        if (!(rz = Pympz_new()))
            return NULL;

        /* integer % mpz */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("Modulo (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_r(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }
#endif
        if (PyLong_Check(a)) {
            TRACE("Modulo (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_r(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }

        /* xmpz % mpz */
        if ((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            TRACE("Modulo (xmpz,mpz)\n");
            mpz_fdiv_r(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if (Pyxmpz_Check(a)) {
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz % integer */
        if (PyIntOrLong_Check(b)) {
            TRACE("Modulo (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_r(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_r_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("xmpz modulo by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rxz;
        }

        /* xmpz % xmpz */
        if (Pyxmpz_Check(b)) {
            TRACE("Modulo (xmpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz modulo by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_r(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }

        /* xmpz % mpz */
        if ((options.prefer_mutable) && Pympz_Check(b)) {
            TRACE("Modulo (xmpz,mpz)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz modulo by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_r(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (Pyxmpz_Check(b)) {
        if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            ZERO_ERROR("mpz division by zero");
            return NULL;
        }
        if (!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer % xmpz */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("Modulo (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_r(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }
#endif
        if (PyLong_Check(a)) {
            TRACE("Modulo (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_r(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }

        /* mpz % xmpz */
        if ((options.prefer_mutable) && Pympz_Check(a)) {
            TRACE("Modulo (mpz,xmpz)\n");
            mpz_fdiv_r(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if (isRational(a) && isRational(b)) {
        TRACE("Modulo (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("mpq modulo by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
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

    if (isFloat(a) && isFloat(b)) {
        TRACE("Modulo (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (mpfr_zero_p(pbf->f)) {
            ZERO_ERROR("mpf modulo by zero");
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(options.rounding))) {
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_div(rf->f, paf->f, pbf->f, options.rounding);
        mpfr_floor(rf->f, rf->f);
        mpfr_mul(rf->f, pbf->f, rf->f, options.rounding);
        mpfr_sub(rf->f, paf->f, rf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return (PyObject*)rf;
    }
    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_divmod follows the semantics from Python 3.x. The result types
 * are:
 *   divmod(mpz, mpz) -> (mpz, mpz)
 *   divmod(mpq, mpq) -> (mpz, mpq)
 *   divmod(mpf, mpf) -> (mpf, mpf)
 *
 * The behavior of mpq now mimics the behavior of fractions.Fraction.
 */

static PyObject *
Pympany_divmod(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    mpz_t tempz;
    PympzObject *qz = 0, *rz = 0;
    PyxmpzObject *qxz = 0, *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *qf = 0, *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;

    if (Pympz_Check(a)) {
        if (!(r=PyTuple_New(2)) || !(rz=Pympz_new()) || !(qz=Pympz_new())) {
            Py_XDECREF((PyObject*)rz);
            Py_XDECREF((PyObject*)qz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(mpz,integer) */
        if (PyIntOrLong_Check(b)) {
            TRACE("divmod (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_qr(qz->z, rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_qr_ui(qz->z, rz->z, Pympz_AS_MPZ(a), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("mpz divmod by zero");
                Py_DECREF((PyObject*)rz);
                Py_DECREF((PyObject*)qz);
                return NULL;
            }
            else {
                mpz_cdiv_qr_ui(qz->z, rz->z, Pympz_AS_MPZ(a), -temp);
                mpz_neg(qz->z, qz->z);
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }

        /* divmod(mpz,mpz) */
        if (Pympz_Check(b)) {
            TRACE("divmod (mpz,mpz)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz divmod by zero");
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

        /* divmod(mpz,xmpz) */
        if ((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            TRACE("divmod (mpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("mpz divmod by zero");
                Py_DECREF((PyObject*)rz);
                Py_DECREF((PyObject*)qz);
                Py_DECREF(r);
                return NULL;
            }
            mpz_fdiv_qr(qz->z, rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
   }

    if (Pympz_Check(b)) {
        if (mpz_sgn(Pympz_AS_MPZ(b))==0) {
            ZERO_ERROR("mpz modulo by zero");
            return NULL;
        }
        if (!(r=PyTuple_New(2)) || !(rz=Pympz_new()) || !(qz=Pympz_new())) {
            Py_XDECREF((PyObject*)rz);
            Py_XDECREF((PyObject*)qz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(integer,mpz) */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("divmod (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_qr(qz->z, rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
#endif
        if (PyLong_Check(a)) {
            TRACE("divmod (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_qr(qz->z, rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }

        /* divmod(xmpz,mpz) */
        if ((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            TRACE("divmod (xmpz,mpz)\n");
            mpz_fdiv_qr(qz->z, rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
    }

    if (Pyxmpz_Check(a)) {
        if (!(r=PyTuple_New(2)) || !(rxz=Pyxmpz_new()) || !(qxz=Pyxmpz_new())) {
            Py_XDECREF((PyObject*)rxz);
            Py_XDECREF((PyObject*)qxz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(xmpz,integer) */
        if (PyIntOrLong_Check(b)) {
            TRACE("divmod (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_qr(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else if (temp > 0) {
                mpz_fdiv_qr_ui(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            else if (temp == 0) {
                ZERO_ERROR("xmpz divmod by zero");
                Py_DECREF((PyObject*)rxz);
                Py_DECREF((PyObject*)qxz);
                return NULL;
            }
            else {
                mpz_cdiv_qr_ui(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), -temp);
                mpz_neg(qxz->z, qxz->z);
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }

        /* divmod(xmpz,xmpz) */
        if (Pyxmpz_Check(b)) {
            TRACE("divmod (xmpz,xmpz)\n");
            if (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz divmod by zero");
                Py_DECREF((PyObject*)rxz);
                Py_DECREF((PyObject*)qxz);
                Py_DECREF(r);
                return NULL;
            }
            mpz_fdiv_qr(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }

        /* divmod(xmpz,mpz) */
        if ((options.prefer_mutable) && Pympz_Check(b)) {
            TRACE("divmod (xmpz,mpz)\n");
            if (mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                ZERO_ERROR("xmpz divmod by zero");
                Py_DECREF((PyObject*)rxz);
                Py_DECREF((PyObject*)qxz);
                Py_DECREF(r);
                return NULL;
            }
            mpz_fdiv_qr(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }
    }

    if (Pyxmpz_Check(b)) {
        if (mpz_sgn(Pyxmpz_AS_MPZ(b))==0) {
            ZERO_ERROR("xmpz divmod by zero");
            return NULL;
        }
        if (!(r=PyTuple_New(2)) || !(rxz=Pyxmpz_new()) || !(qxz=Pyxmpz_new())) {
            Py_XDECREF((PyObject*)rxz);
            Py_XDECREF((PyObject*)qxz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(integer,xmpz) */
#ifdef PY2
        if (PyInt_Check(a)) {
            TRACE("divmod (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_qr(qxz->z, rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }
#endif
        if (PyLong_Check(a)) {
            TRACE("divmod (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_qr(qxz->z, rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }

        /* divmod(mpz,xmpz) */
        if ((options.prefer_mutable) && Pympz_Check(a)) {
            TRACE("divmod (xmpz,mpz)\n");
            mpz_fdiv_qr(qxz->z, rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }
    }

    if (isRational(a) && isRational(b)) {
        TRACE("Divmod (rational,rational)\n");
        paq = Pympq_From_Rational(a);
        pbq = Pympq_From_Rational(b);
        if (!paq || !pbq) {
            SYSTEM_ERROR("Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq);
            Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (mpq_sgn(pbq->q)==0) {
            ZERO_ERROR("mpq divmod by zero");
            Py_DECREF((PyObject*)paq);
            Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new()) || !(qz = Pympz_new())) {
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
        return Py_BuildValue("(NN)", qz, rq);
    }

    if (isFloat(a) && isFloat(b)) {
        TRACE("Divmod (number,number)\n");
        paf = Pympf_From_Float(a, 0);
        pbf = Pympf_From_Float(b, 0);
        if (!paf || !pbf) {
            SYSTEM_ERROR("Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf);
            Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if (mpfr_zero_p(pbf->f)) {
            ZERO_ERROR("mpf divmod by zero");
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(qf = Pympf_new(options.precision)) ||
            !(rf = Pympf_new(options.precision))) {
            Py_XDECREF((PyObject*)qf);
            Py_XDECREF((PyObject*)rf);
            Py_DECREF((PyObject*)paf);
            Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpfr_div(qf->f, paf->f, pbf->f, options.rounding);
        mpfr_floor(qf->f, qf->f);
        mpfr_mul(rf->f, pbf->f, qf->f, options.rounding);
        mpfr_sub(rf->f, paf->f, rf->f, options.rounding);
        Py_DECREF((PyObject*)paf);
        Py_DECREF((PyObject*)pbf);
        return Py_BuildValue("(NN)", qf, rf);
    }
    Py_RETURN_NOTIMPLEMENTED;
}
