/* gmpy_basic.c
 *
 * Generic methods for gmpy types.
 *
 * This file should be considered part of gmpy.c
 */

#include <math.h>

/* Generic addition
 *
 * Support addition for gmpy types with automatic conversion of Python types.
*
 * The following conversion logic is used:
 *  1) 'mpz' combined with an integer type returns an 'mpz'
 *  2) 'mpz' combined with an integer or rational type returns an 'mpq'
 *  3) 'mpz' combined with a floating-point type returns an 'mpf'
 *  4) 'mpq' combined with an integer or rational type returns an 'mpq'
 *  5) 'mpq' combines with a floating-point type returns an 'mpf'
 *
 * The most common inputs are processed as efficiently as possible.
 */

static PyObject *
Pympany_add(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    /* Try to make mpz + small_int faster */

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Adding (mpz,small_int)\n");
        if((temp=Py2or3Int_AsLong(b))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_sub_ui(rz->z, ((PympzObject*)a)->z, -temp);
                return (PyObject *) rz;
            }
        } else {
            if (!(rz = Pympz_new())) return NULL;
            mpz_add_ui(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(Pympz_Check(b) && Py2or3Int_Check(a)) {
        if (options.debug) fprintf(stderr, "Adding (small_int,mpz)\n");
        if((temp=Py2or3Int_AsLong(a))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_sub_ui(rz->z, ((PympzObject*)b)->z, -temp);
                return (PyObject *) rz;
            }
        } else {
            if (!(rz = Pympz_new())) return NULL;
            mpz_add_ui(rz->z, ((PympzObject*)b)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Adding (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_add(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return (PyObject *) rz;
    }

    if(isRational(a) && isRational(b)) {
        if (options.debug) fprintf(stderr, "Adding (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_add(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        return (PyObject *) rq;
    }

    if(isNumber(a) && isNumber(b)) {
        if (options.debug) fprintf(stderr, "Adding (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        if (!(rf = Pympf_new(bits))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_add(rf->f, paf->f, pbf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        Pympf_normalize(rf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Generic Subtraction
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pympany_sub(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Subtracting (mpz,small_int)\n");
        if((temp=Py2or3Int_AsLong(b))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_add_ui(rz->z, ((PympzObject*)a)->z, -temp);
                return (PyObject *) rz;
            }
        } else {
            if (!(rz = Pympz_new())) return NULL;
            mpz_sub_ui(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(Pympz_Check(b) && Py2or3Int_Check(a)) {
        if (options.debug) fprintf(stderr, "Subtracting (small_int,mpz)\n");
        if((temp=Py2or3Int_AsLong(a))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_add_ui(rz->z, ((PympzObject*)b)->z, -temp);
                mpz_neg(rz->z, rz->z);
                return (PyObject *) rz;
            }
        } else {
            if (!(rz = Pympz_new())) return NULL;
            mpz_sub_ui(rz->z, ((PympzObject*)b)->z, temp);
            mpz_neg(rz->z, rz->z);
            return (PyObject *) rz;
        }
    }

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Subtractinging (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_sub(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return (PyObject *) rz;
    }

    if(isRational(a) && isRational(b)) {
        if (options.debug) fprintf(stderr, "Subtracting (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_sub(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        return (PyObject *) rq;
    }

    if(isNumber(a) && isNumber(b)) {
        if (options.debug) fprintf(stderr, "Subtracting (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        if (!(rf = Pympf_new(bits))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_sub(rf->f, paf->f, pbf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        Pympf_normalize(rf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Generic Multiplication
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pympany_mul(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Multiplying (mpz,small_int)\n");
        temp=Py2or3Int_AsLong(b);
        if(PyErr_Occurred()) {
            PyErr_Clear();
        } else {
            if (!(rz = Pympz_new())) return NULL;
            mpz_mul_si(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(Pympz_Check(b) && Py2or3Int_Check(a)) {
        if (options.debug) fprintf(stderr, "Multiplying (mpz,small_int)\n");
        temp=Py2or3Int_AsLong(a);
        if(PyErr_Occurred()) {
            PyErr_Clear();
        } else {
            if (!(rz = Pympz_new())) return NULL;
            mpz_mul_si(rz->z, ((PympzObject*)b)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Multiplying (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_mul(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return (PyObject *) rz;
    }

    if(isRational(a) && isRational(b)) {
        if (options.debug) fprintf(stderr, "Multiplying (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_mul(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        return (PyObject *) rq;
    }

    if(isNumber(a) && isNumber(b)) {
        if (options.debug) fprintf(stderr, "Multiplying (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        if (!(rf = Pympf_new(bits))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_mul(rf->f, paf->f, pbf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        Pympf_normalize(rf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
Pympany_floordiv(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Dividing (mpz,small_int)\n");
        if((temp=Py2or3Int_AsLong(b))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_cdiv_q_ui(rz->z, ((PympzObject*)a)->z, -temp);
                mpz_neg(rz->z, rz->z);
                return (PyObject *) rz;
            }
        } else {
            if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                return NULL;
            }
            if (!(rz = Pympz_new())) return NULL;
            mpz_fdiv_q_ui(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Floor divide (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if(mpz_sgn(pbz->z)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_fdiv_q(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return (PyObject *) rz;
    }

    if(isRational(a) && isRational(b)) {
        if(options.debug) fprintf(stderr, "Floor divide (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if(mpq_sgn(pbq->q)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpq division by zero");
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new()) || !(rz = Pympz_new())) {
            Py_XDECREF((PyObject*)rq); Py_XDECREF((PyObject*)rz);
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(rz->z, mpq_numref(rq->q), mpq_denref(rq->q));
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        Py_DECREF((PyObject*)rq);
        return (PyObject *) rz;
    }

    if(isNumber(a) && isNumber(b)) {
        if(options.debug) fprintf(stderr, "Floor divide (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if(mpf_sgn(pbf->f)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpf division by zero");
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        if (!(rf = Pympf_new(bits))) {
            Py_XDECREF((PyObject*)rf);
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(rf->f, paf->f, pbf->f);
        mpf_floor(rf->f, rf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
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
    PyObject *r = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    unsigned int bits;

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "True divide (integer,integer)\n");
        paf = anynum2Pympf(a, 0);
        pbf = anynum2Pympf(b, 0);
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if(mpf_sgn(pbf->f)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpf division by zero");
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        if (!(rf = Pympf_new(0))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(rf->f, paf->f, pbf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        return (PyObject *) rf;
    }

    if(isRational(a) && isRational(b)) {
        if(options.debug) fprintf(stderr, "True divide (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if(mpq_sgn(pbq->q)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpq division by zero");
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        return (PyObject *) rq;
    }

    if(isNumber(a) && isNumber(b)) {
        if(options.debug) fprintf(stderr, "True divide (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if(mpf_sgn(pbf->f)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpf division by zero");
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        if (!(rf = Pympf_new(bits))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(rf->f, paf->f, pbf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

#if PY_MAJOR_VERSION < 3
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
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Dividing (mpz,small_int)\n");
        if((temp=Py2or3Int_AsLong(b))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_cdiv_q_ui(rz->z, ((PympzObject*)a)->z, -temp);
                mpz_neg(rz->z, rz->z);
                return (PyObject *) rz;
            }
        } else {
            if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                return NULL;
            }
            if (!(rz = Pympz_new())) return NULL;
            mpz_fdiv_q_ui(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    /* Use floordiv for integer types. */

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Floor divide (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if(mpz_sgn(pbz->z)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_fdiv_q(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return (PyObject *) rz;
    }

    /* Use truediv for rational types. */

    if(isRational(a) && isRational(b)) {
        if(options.debug) fprintf(stderr, "True divide (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if(mpq_sgn(pbq->q)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpq division by zero");
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new())) {
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        return (PyObject *) rq;
    }

    /* Use truediv for floating-point types. */

    if(isNumber(a) && isNumber(b)) {
        if(options.debug) fprintf(stderr, "True divide (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if(mpf_sgn(pbf->f)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpf division by zero");
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        if (!(rf = Pympf_new(bits))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(rf->f, paf->f, pbf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
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
    PyObject *r = 0;
    PympzObject *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Modulo (mpz,small_int)\n");
        if((temp=Py2or3Int_AsLong(b))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new())) return NULL;
                mpz_cdiv_r_ui(rz->z, ((PympzObject*)a)->z, -temp);
                mpz_neg(rz->z, rz->z);
                return (PyObject *) rz;
            }
        } else {
            if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                return NULL;
            }
            if (!(rz = Pympz_new())) return NULL;
            mpz_fdiv_r_ui(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Modulo (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if(mpz_sgn(pbz->z)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = Pympz_new())) {
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_fdiv_r(rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return (PyObject *) rz;
    }

    if(isRational(a) && isRational(b)) {
        if(options.debug) fprintf(stderr, "Modulo (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if(mpq_sgn(pbq->q)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpq modulo by zero");
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new()) || !(rz = Pympz_new())) {
            Py_XDECREF((PyObject*)rq); Py_XDECREF((PyObject*)rz);
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(rz->z, mpq_numref(rq->q), mpq_denref(rq->q));
        /* Need to calculate paq - rz * pbq */
        mpq_set_z(rq->q, rz->z);
        mpq_mul(rq->q, rq->q, pbq->q);
        mpq_sub(rq->q, paq->q, rq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        Py_DECREF((PyObject*)rz);
        return (PyObject *) rq;
    }

    if(isNumber(a) && isNumber(b)) {
        if(options.debug) fprintf(stderr, "Modulo (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if(mpf_sgn(pbf->f)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpf modulo by zero");
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        /* To prevent rounding errors, the working precision is increased. */
        temp = (paf->f->_mp_exp - pbf->f->_mp_exp) * GMP_NUMB_BITS + bits;
        if(options.debug) {
            fprintf(stderr, "Working precision %ld\n", temp);
        }
        if (!(rf = Pympf_new(temp))) {
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(rf->f, paf->f, pbf->f);
        mpf_floor(rf->f, rf->f);
        mpf_mul(rf->f, pbf->f, rf->f);
        mpf_sub(rf->f, paf->f, rf->f);
        mpf_set_prec(rf->f, bits);
        rf->rebits = bits;
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
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
    PympzObject *qz = 0, *rz = 0, *paz = 0, *pbz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *qf = 0, *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if (options.debug) fprintf(stderr, "Divmod (mpz,small_int)\n");
        if((temp=Py2or3Int_AsLong(b))<0) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                if (!(rz = Pympz_new()) || !(qz = Pympz_new())) {
                    Py_XDECREF((PyObject*)rz); Py_XDECREF((PyObject*)qz);
                    return NULL;
                }
                mpz_cdiv_qr_ui(qz->z, rz->z, ((PympzObject*)a)->z, -temp);
                mpz_neg(rz->z, rz->z);
                return Py_BuildValue("(NN)", qz, rz);
            }
        } else {
            if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
                return NULL;
            }
            if (!(rz = Pympz_new()) || !(qz = Pympz_new())) {
                Py_XDECREF((PyObject*)rz); Py_XDECREF((PyObject*)qz);
                return NULL;
            }
            mpz_fdiv_qr_ui(qz->z, rz->z, ((PympzObject*)a)->z, temp);
            return Py_BuildValue("(NN)", qz, rz);
        }
    }

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "Divmod (integer,integer)\n");
        paz = anyint2Pympz(a);
        pbz = anyint2Pympz(b);
        if(!paz || !pbz) {
            Py_XDECREF((PyObject*)paz); Py_XDECREF((PyObject*)pbz);
            PyErr_SetString(PyExc_SystemError, "Can not convert integer to mpz");
            return NULL;
        }
        if(mpz_sgn(pbz->z)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        if (!(rz = Pympz_new()) || !(qz = Pympz_new())) {
            Py_XDECREF((PyObject*)rz); Py_XDECREF((PyObject*)qz);
            Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
            return NULL;
        }
        mpz_fdiv_qr(qz->z, rz->z, paz->z, pbz->z);
        Py_DECREF((PyObject*)paz); Py_DECREF((PyObject*)pbz);
        return Py_BuildValue("(NN)", qz, rz);
    }

    if(isRational(a) && isRational(b)) {
        if(options.debug) fprintf(stderr, "Divmod (rational,rational)\n");
        paq = anyrational2Pympq(a);
        pbq = anyrational2Pympq(b);
        if(!paq || !pbq) {
            PyErr_SetString(PyExc_SystemError, "Can not convert rational to mpq");
            Py_XDECREF((PyObject*)paq); Py_XDECREF((PyObject*)pbq);
            return NULL;
        }
        if(mpq_sgn(pbq->q)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpq divmod by zero");
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        if (!(rq = Pympq_new()) || !(qz = Pympz_new())) {
            Py_XDECREF((PyObject*)rq); Py_XDECREF((PyObject*)qz);
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(qz->z, mpq_numref(rq->q), mpq_denref(rq->q));
        /* Need to calculate paq - rz * pbq */
        mpq_set_z(rq->q, qz->z);
        mpq_mul(rq->q, rq->q, pbq->q);
        mpq_sub(rq->q, paq->q, rq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        return Py_BuildValue("(NN)", qz, rq);
    }

    if(isNumber(a) && isNumber(b)) {
        if(options.debug) fprintf(stderr, "Divmod (number,number)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2Pympf(a, 0);
            pbf = anynum2Pympf(b, paf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        }
        if(!paf || !pbf) {
            PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
            Py_XDECREF((PyObject*)paf); Py_XDECREF((PyObject*)pbf);
            return NULL;
        }
        if(mpf_sgn(pbf->f)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpf divmod by zero");
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        bits = paf->rebits;
        if(pbf->rebits<bits) bits=pbf->rebits;
        /* To prevent rounding errors, the working precision is increased. */
        temp = (paf->f->_mp_exp - pbf->f->_mp_exp) * GMP_NUMB_BITS + bits;
        if(options.debug) {
            fprintf(stderr, "Working precision %ld\n", temp);
        }
        if (!(qf = Pympf_new(temp)) || !(rf = Pympf_new(temp))) {
            Py_XDECREF((PyObject*)qf); Py_XDECREF((PyObject*)rf);
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(qf->f, paf->f, pbf->f);
        mpf_floor(qf->f, qf->f);
        mpf_mul(rf->f, pbf->f, qf->f);
        mpf_sub(rf->f, paf->f, rf->f);
        mpf_set_prec(rf->f, bits);
        rf->rebits = bits;
        mpf_set_prec(qf->f, bits);
        qf->rebits = bits;
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        return Py_BuildValue("(NN)", qf, rf);
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}
