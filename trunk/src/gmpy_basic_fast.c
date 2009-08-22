/* gmpy_basic_fast.c
 *
 * Generic methods for gmpy types. The routines are highly optimized.
 *
 * This file should be considered part of gmpy.c
 */

#include <math.h>

#define LONG_FROMONE(L) \
    ((PyLongObject*)L)->ob_digit[0]
#define LONG_FROMTWO(L) \
    ((long)((PyLongObject*)L)->ob_digit[1] << Py2or3Long_SHIFT) \
    + ((PyLongObject*)L)->ob_digit[0]

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
    mpz_t rzz, tempz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
#if PY_MAJOR_VERSION == 2
    long temp;
#endif
    unsigned int bits;

    /* Try to make mpz + small_int faster */

    if(Pympz_Check(a)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Adding (mpz,small_int)\n");
            if((temp = PyInt_AS_LONG(b))>=0) {
                mpz_add_ui(rzz, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_sub_ui(rzz, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif

        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Adding (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    mpz_set(rzz, Pympz_AS_MPZ(a));
                    break;
                case -1:
                    mpz_sub_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case -2:
                    mpz_sub_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_add(rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "Adding (mpz,mpz)\n");
            mpz_add(rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
    }

    if(Pympz_Check(b)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Adding (small_int,mpz)\n");
            if((temp = PyInt_AS_LONG(a))>=0) {
                mpz_add_ui(rzz, Pympz_AS_MPZ(b), temp);
            } else {
                mpz_sub_ui(rzz, Pympz_AS_MPZ(b), -temp);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif

        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Adding (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMTWO(a));
                    break;
                case 1:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMONE(a));
                    break;
                case 0:
                    mpz_set(rzz, Pympz_AS_MPZ(b));
                    break;
                case -1:
                    mpz_sub_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMONE(a));
                    break;
                case -2:
                    mpz_sub_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMTWO(a));
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_add(rzz, Pympz_AS_MPZ(b), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
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
        mpf_normalize(rf->f);
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
    mpz_t rzz, tempz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
#if PY_MAJOR_VERSION == 2
    long temp;
#endif
    unsigned int bits;

    if(Pympz_Check(a)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Subtracting (mpz,small_int)\n");
            if((temp = PyInt_AS_LONG(b))>=0) {
                mpz_sub_ui(rzz, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_add_ui(rzz, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif

        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Subtracting (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_sub_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_sub_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    mpz_set(rzz, Pympz_AS_MPZ(a));
                    break;
                case -1:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case -2:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_sub(rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "Subtracting (mpz,mpz)\n");
            mpz_sub(rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
    }

    if(Pympz_Check(b)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Subtracting (small_int,mpz)\n");
            if((temp = PyInt_AS_LONG(a))>=0) {
                mpz_ui_sub(rzz, temp, Pympz_AS_MPZ(b));
            } else {
                mpz_add_ui(rzz, Pympz_AS_MPZ(b), -temp);
                mpz_neg(rzz, rzz);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif

        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Subtracting (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_ui_sub(rzz, LONG_FROMTWO(a), Pympz_AS_MPZ(b));
                    break;
                case 1:
                    mpz_ui_sub(rzz, LONG_FROMONE(a), Pympz_AS_MPZ(b));
                    break;
                case 0:
                    mpz_set(rzz, Pympz_AS_MPZ(b));
                    mpz_neg(rzz, rzz);
                    break;
                case -1:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMONE(a));
                    mpz_neg(rzz, rzz);
                    break;
                case -2:
                    mpz_add_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMTWO(a));
                    mpz_neg(rzz, rzz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_sub(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
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
        mpf_normalize(rf->f);
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
    mpz_t rzz, tempz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    unsigned int bits;

    if(Pympz_Check(a)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Multiplying (mpz,small_int)\n");
            mpz_mul_si(rzz, Pympz_AS_MPZ(a), PyInt_AS_LONG(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif

        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Multiplying (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    mpz_set_si(rzz, 0);
                    break;
                case -1:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    mpz_neg(rzz, rzz);
                    break;
                case -2:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    mpz_neg(rzz, rzz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_mul(rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "Multiplying (mpz,mpz)\n");
            mpz_mul(rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
    }

    if(Pympz_Check(b)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Multiplying (small_int,mpz)\n");
            mpz_mul_si(rzz, Pympz_AS_MPZ(b), PyInt_AS_LONG(a));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif

        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Multiplying (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMTWO(a));
                    break;
                case 1:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMONE(a));
                    break;
                case 0:
                    mpz_set_si(rzz, 0);
                    break;
                case -1:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMONE(a));
                    mpz_neg(rzz, rzz);
                    break;
                case -2:
                    mpz_mul_ui(rzz, Pympz_AS_MPZ(b), LONG_FROMTWO(a));
                    mpz_neg(rzz, rzz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_mul(rzz, Pympz_AS_MPZ(b), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
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
        mpf_normalize(rf->f);
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
    mpz_t rzz, tempz;
    PympzObject *rz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
#if PY_MAJOR_VERSION == 2
    long temp;
#endif
    unsigned int bits;

    if(Pympz_Check(a)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Floor divide (mpz,small_int)\n");
            if((temp=PyInt_AS_LONG(b))>0) {
                mpz_fdiv_q_ui(rzz, Pympz_AS_MPZ(a), temp);
            } else if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                mpz_cloc(rzz);
                return NULL;
            } else {
                mpz_cdiv_q_ui(rzz, Pympz_AS_MPZ(a), -temp);
                mpz_neg(rzz, rzz);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif
        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Floor divide (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_fdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_fdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                    mpz_cloc(rzz);
                    return NULL;
                    break;
                case -1:
                    mpz_cdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    mpz_neg(rzz, rzz);
                    break;
                case -2:
                    mpz_cdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    mpz_neg(rzz, rzz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_fdiv_q(rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "Floor divide (integer,integer)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                mpz_cloc(rzz);
                return NULL;
            }
            mpz_fdiv_q(rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Floor divide (small_int,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif
        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Floor divide (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMTWO(a));
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 1:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMONE(a));
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 0:
                    mpz_set_si(rzz, 0);
                    break;
                case -1:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMONE(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case -2:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMTWO(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
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
        mpf_normalize(rf->f);
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
        mpf_normalize(rf->f);
        return (PyObject *) rf;
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}

#if PY_MAJOR_VERSION == 2
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
    mpz_t rzz, tempz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    /* Use floordiv for integer types. */

    if(Pympz_Check(a)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "True divide (mpz,small_int)\n");
            if((temp=PyInt_AS_LONG(b))>0) {
                mpz_fdiv_q_ui(rzz, Pympz_AS_MPZ(a), temp);
            } else if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                mpz_cloc(rzz);
                return NULL;
            } else {
                mpz_cdiv_q_ui(rzz, Pympz_AS_MPZ(a), -temp);
                mpz_neg(rzz, rzz);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif
        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "True divide (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_fdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_fdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                    mpz_cloc(rzz);
                    return NULL;
                    break;
                case -1:
                    mpz_cdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    mpz_neg(rzz, rzz);
                    break;
                case -2:
                    mpz_cdiv_q_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    mpz_neg(rzz, rzz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_fdiv_q(rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "True divide (integer,integer)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                mpz_cloc(rzz);
                return NULL;
            }
            mpz_fdiv_q(rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "True divide (small_int,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif
        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "True divide (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMTWO(a));
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 1:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMONE(a));
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 0:
                    mpz_set_si(rzz, 0);
                    break;
                case -1:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMONE(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case -2:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMTWO(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_fdiv_q(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
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
        mpf_normalize(rf->f);
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
    mpz_t rzz, tempz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a)) {
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Modulo (mpz,small_int)\n");
            if((temp=PyInt_AS_LONG(b))==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                mpz_cloc(rzz);
                return NULL;
            } else if(temp>0) {
                mpz_fdiv_r_ui(rzz, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_cdiv_r_ui(rzz, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif
        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "Modulo (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_fdiv_r_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_fdiv_r_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                    mpz_cloc(rzz);
                    return NULL;
                case -1:
                    mpz_cdiv_r_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case -2:
                    mpz_cdiv_r_ui(rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_fdiv_r(rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "Modulo (integer,integer)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                mpz_cloc(rzz);
                return NULL;
            }
            mpz_fdiv_r(rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        mpz_inoc(rzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Modulo (small_int,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_r(rzz, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
#endif
        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "Modulo (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_inoc(tempz);
                    mpz_set_si(tempz, LONG_FROMTWO(a));
                    mpz_fdiv_r(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 1:
                    mpz_inoc(tempz);
                    mpz_set_si(tempz, LONG_FROMONE(a));
                    mpz_fdiv_r(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 0:
                    mpz_set_si(rzz, 0);
                    break;
                case -1:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMONE(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_r(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case -2:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMTWO(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_r(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_fdiv_r(rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
            }
            return (PyObject *)Pympz_FROM_MPZ(rzz);
        }
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
        if (!(rq = Pympq_new())) {
            Py_XDECREF((PyObject*)rq);
            Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
            return NULL;
        }
        mpz_inoc(rzz);
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(rzz, mpq_numref(rq->q), mpq_denref(rq->q));
        /* Need to calculate paq - rz * pbq */
        mpq_set_z(rq->q, rzz);
        mpq_mul(rq->q, rq->q, pbq->q);
        mpq_sub(rq->q, paq->q, rq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        mpz_cloc(rzz);
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
        mpf_normalize(rf->f);
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
    mpz_t qzz, rzz, tempz;
    PympzObject *qz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *qf = 0, *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;

    if(Pympz_Check(a)) {
        if(!(r=PyTuple_New(2))) {
            return NULL;
        }
        mpz_inoc(rzz);
        mpz_inoc(qzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "divmod (mpz,small_int)\n");
            if((temp=PyInt_AS_LONG(b))>0) {
                mpz_fdiv_qr_ui(qzz, rzz, Pympz_AS_MPZ(a), temp);
            } else if(temp==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
                mpz_cloc(rzz);
                mpz_cloc(qzz);
                Py_DECREF(r);
                return NULL;
            } else {
                mpz_cdiv_qr_ui(qzz, rzz, Pympz_AS_MPZ(a), -temp);
                mpz_neg(qzz, qzz);
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)Pympz_FROM_MPZ(qzz));
            PyTuple_SET_ITEM(r, 1, (PyObject*)Pympz_FROM_MPZ(rzz));
            return r;
        }
#endif
        if(PyLong_CheckExact(b)) {
            if (options.debug) fprintf(stderr, "divmod (mpz,long)\n");
            switch(Py_SIZE(b)) {
                case 2:
                    mpz_fdiv_qr_ui(qzz, rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    break;
                case 1:
                    mpz_fdiv_qr_ui(qzz, rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    break;
                case 0:
                    PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
                    mpz_cloc(rzz);
                    mpz_cloc(qzz);
                    return NULL;
                    break;
                case -1:
                    mpz_cdiv_qr_ui(qzz, rzz, Pympz_AS_MPZ(a), LONG_FROMONE(b));
                    mpz_neg(qzz, qzz);
                    break;
                case -2:
                    mpz_cdiv_qr_ui(qzz, rzz, Pympz_AS_MPZ(a), LONG_FROMTWO(b));
                    mpz_neg(qzz, qzz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, b);
                    mpz_fdiv_qr(qzz, rzz, Pympz_AS_MPZ(a), tempz);
                    mpz_cloc(tempz);
                    break;
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)Pympz_FROM_MPZ(qzz));
            PyTuple_SET_ITEM(r, 1, (PyObject*)Pympz_FROM_MPZ(rzz));
            return r;
        }

        if(Pympz_Check(b)) {
            if(options.debug) fprintf(stderr, "divmod (integer,integer)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
                mpz_cloc(rzz);
                mpz_cloc(qzz);
                Py_DECREF(r);
                return NULL;
            }
            mpz_fdiv_qr(qzz, rzz, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)Pympz_FROM_MPZ(qzz));
            PyTuple_SET_ITEM(r, 1, (PyObject*)Pympz_FROM_MPZ(rzz));
            return r;
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
            return NULL;
        }

        if(!(r=PyTuple_New(2))) {
            return NULL;
        }
        mpz_inoc(rzz);
        mpz_inoc(qzz);
#if PY_MAJOR_VERSION == 2
        if(PyInt_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "divmod (small_int,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_qr(qzz, rzz, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)Pympz_FROM_MPZ(qzz));
            PyTuple_SET_ITEM(r, 1, (PyObject*)Pympz_FROM_MPZ(rzz));
            return r;
        }
#endif
        if(PyLong_CheckExact(a)) {
            if (options.debug) fprintf(stderr, "divmod (long,mpz)\n");
            switch(Py_SIZE(a)) {
                case 2:
                    mpz_inoc(tempz);
                    mpz_set_si(tempz, LONG_FROMTWO(a));
                    mpz_fdiv_qr(qzz, rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 1:
                    mpz_inoc(tempz);
                    mpz_set_si(tempz, LONG_FROMONE(a));
                    mpz_fdiv_qr(qzz, rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case 0:
                    mpz_set_si(rzz, 0);
                    mpz_set_si(qzz, 0);
                    break;
                case -1:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMONE(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_qr(qzz, rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                case -2:
                    mpz_inoc(tempz);
                    mpz_set_ui(tempz, LONG_FROMTWO(a));
                    mpz_neg(tempz, tempz);
                    mpz_fdiv_qr(qzz, rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
                default:
                    mpz_inoc(tempz);
                    mpz_set_PyLong(tempz, a);
                    mpz_fdiv_qr(qzz, rzz, tempz, Pympz_AS_MPZ(b));
                    mpz_cloc(tempz);
                    break;
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)Pympz_FROM_MPZ(qzz));
            PyTuple_SET_ITEM(r, 1, (PyObject*)Pympz_FROM_MPZ(rzz));
            return r;
        }
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
        mpf_normalize(qf->f);
        mpf_normalize(rf->f);
        return Py_BuildValue("(NN)", qf, rf);
    }

    r = Py_NotImplemented;
    Py_INCREF(r);
    return r;
}
