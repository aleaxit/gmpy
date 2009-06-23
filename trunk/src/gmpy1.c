/* gmpy1.c
 *
 * Generic methods for gmpy types.
 *
 * This file should be considered part of gmpy.c
 */

/* Generic addition
 *
 * Support addition for gmpy types with automatic conversion of Python types.
 * The following conversion logic is used:
 *  1) 'mpz' combined with an integer type returns an 'mpz'
 *  2) 'mpz' combined with an integer or rational type returns an 'mpq'
 *  3) 'mpz' combined with a floating-point type returns an 'mpf'
 *  4) 'mpq' combined with an integer or rational type returns an 'mpq'
 *  5) 'mpq' combines with a floating-point type returns an 'mpf'
 *
 * The most common inputs are processed as efficiently as possible.
 */

#include <math.h>

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
        paz = anyint2mpz(a);
        pbz = anyint2mpz(b);
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
        paq = anyrational2mpq(a);
        pbq = anyrational2mpq(b);
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
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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
        mpf_normalize(rf);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Generic Subtraction */

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
        paz = anyint2mpz(a);
        pbz = anyint2mpz(b);
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
        paq = anyrational2mpq(a);
        pbq = anyrational2mpq(b);
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
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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
        mpf_normalize(rf);
        return (PyObject *) rf;
    }

    r= Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

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
        paz = anyint2mpz(a);
        pbz = anyint2mpz(b);
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
        paq = anyrational2mpq(a);
        pbq = anyrational2mpq(b);
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
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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
        mpf_normalize(rf);
        return (PyObject *) rf;
    }

    r= Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Pympany_floordiv follows the // semantics from Python 3.x */

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
        paz = anyint2mpz(a);
        pbz = anyint2mpz(b);
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
        paq = anyrational2mpq(a);
        pbq = anyrational2mpq(b);
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
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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
        if (!(rf = Pympf_new(bits)) || !(rz = Pympz_new())) {
            Py_XDECREF((PyObject*)rf); Py_XDECREF((PyObject*)rz);
            Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
            return NULL;
        }
        mpf_div(rf->f, paf->f, pbf->f);
        mpf_floor(rf->f, rf->f);
        mpz_set_f(rz->z, rf->f);
        Py_DECREF((PyObject*)paf); Py_DECREF((PyObject*)pbf);
        Py_DECREF((PyObject*)rf);
        return (PyObject *) rz;
    }

    r= Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

static PyObject *
Pympany_truediv(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    unsigned int bits;

    if(isInteger(a) && isInteger(b)) {
        if(options.debug) fprintf(stderr, "True divide (integer,integer)\n");
        if(Pympf_Check(a) && Pympf_Check(b)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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

    if(isRational(a) && isRational(b)) {
        if(options.debug) fprintf(stderr, "True divide (rational,rational)\n");
        paq = anyrational2mpq(a);
        pbq = anyrational2mpq(b);
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
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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

    r= Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Pympz_div2 follows the conversions rules for Python 2.x */

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
        paz = anyint2mpz(a);
        pbz = anyint2mpz(b);
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
        paq = anyrational2mpq(a);
        pbq = anyrational2mpq(b);
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
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, 0);
        } else if(Pympf_Check(a)) {
            paf = anynum2mpf(a, 0);
            pbf = anynum2mpf(b, paf->rebits);
        } else {
            pbf = anynum2mpf(b, 0);
            paf = anynum2mpf(a, pbf->rebits);
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

    r= Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

/* Pympz_rem2 follows the conversions rules for Python 2.x */

static PyObject *
Pympz_rem2(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    PympzObject *rz = 0, *pa = 0, *pb = 0;
    long temp;
    double tempa, tempb;

    if(Pympz_Check(a) && Pympz_Check(b)) {
        if(options.debug) fprintf(stderr, "Modulo (mpz,mpz)\n");
        if(mpz_sgn(((PympzObject*)a)->z)==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
            return NULL;
        }
        if (!(rz = Pympz_new())) return NULL;
        mpz_fdiv_r(rz->z, ((PympzObject*)a)->z, ((PympzObject*)b)->z);
        return (PyObject *) rz;
    }

    if(Pympz_Check(a) && Py2or3Int_Check(b)) {
        if((temp=Py2or3Int_AsLong(b))<0) {
            PyErr_Clear();
        } else {
            if (options.debug) fprintf(stderr, "Modulo (mpz,small_int)\n");
            if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                return NULL;
            }
            if (!(rz = Pympz_new())) return NULL;
            mpz_fdiv_r_ui(rz->z, ((PympzObject*)a)->z, temp);
            return (PyObject *) rz;
        }
    }

    if(Pympz_Check(a) && PyFloat_Check(b)) {
        if (options.debug) fprintf(stderr, "Modulo (mpz,float)\n");
        tempa = mpz_get_d(((PympzObject*)a)->z);
        tempb = PyFloat_AsDouble(b);
        if(tempb==0.0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
            return NULL;
        }
        if((tempa<0 && tempb<0) || (tempa>=0 && tempb>0)) {
            return PyFloat_FromDouble(fmod(tempa,tempb));
        } else {
            return PyFloat_FromDouble(fmod(tempa,tempb)+tempb);
        }
    }

    if(PyFloat_Check(a) && Pympz_Check(b)) {
        if (options.debug) fprintf(stderr, "Modulo (float,mpz)\n");
        tempa = PyFloat_AsDouble(a);
        tempb = mpz_get_d(((PympzObject*)b)->z);
        if(tempb==0.0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
            return NULL;
        }
        if((tempa<0 && tempb<0) || (tempa>=0 && tempb>0)) {
            return PyFloat_FromDouble(fmod(tempa,tempb));
        } else {
            return PyFloat_FromDouble(fmod(tempa,tempb)+tempb);
        }
    }

    pa = anyint2mpz(a);
    pb = anyint2mpz(b);
    if(!pa || !pb) {
        r = Py_NotImplemented;
        Py_XDECREF((PyObject*)pa); Py_XDECREF((PyObject*)pb);
        Py_INCREF(r);
        return r;
    }
    if (options.debug) fprintf(stderr, "Modulo (integer, integer)\n");
    if  (mpz_sgn(pb->z) == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
        Py_DECREF((PyObject*)pa); Py_DECREF((PyObject*)pb);
        return NULL;
    }
    if (!(rz = Pympz_new())) {
        Py_DECREF((PyObject*)pa); Py_DECREF((PyObject*)pb);
        return NULL;
    }
    mpz_fdiv_r(rz->z, pa->z, pb->z);
    Py_DECREF((PyObject*)pa); Py_DECREF((PyObject*)pb);
    return (PyObject *) rz;
}
