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
    PyObject *r = 0;
    mpz_t tempz;
    PympzObject *rz = 0;
    PyxmpzObject *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    int overflow;
    unsigned int bits;

    if(Pympz_Check(a)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* mpz + integer */
        if(PyIntOrLong_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_add(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *) rz;
        }

        /* mpz + mpz */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (mpz,mpz)\n");
            mpz_add(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }

        /* mpz + xmpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (mpz,xmpz)\n");
            mpz_add(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pympz_Check(b)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* integer + mpz */
        if(PyIntOrLong_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Adding (long,mpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_add(rz->z, Pympz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            } else if(temp >=0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(b), temp);
            } else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(b), -temp);
            }
            return (PyObject *) rz;
        }

        /* xmpz * mpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Adding (xmpz,mpz)\n");
            mpz_add(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz + integer */
        if(PyIntOrLong_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_add(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_add_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            } else {
                mpz_sub_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
            }
            return (PyObject *) rxz;
        }

        /* xmpz + xmpz */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (xmpz,xmpz)\n");
            mpz_add(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
        }

        /* xmpz + mpz */
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (xmpz,mpz)\n");
            mpz_add(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rxz;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer + xmpz */
        if(PyIntOrLong_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Adding (long,xmpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_add(rz->z, Pyxmpz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            } else if(temp >=0) {
                mpz_add_ui(rxz->z, Pyxmpz_AS_MPZ(b), temp);
            } else {
                mpz_sub_ui(rxz->z, Pyxmpz_AS_MPZ(b), -temp);
            }
            return (PyObject *) rxz;
        }

        /* mpz + xmpz */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Adding (mpz,xmpz)\n");
            mpz_add(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_INFINITY(d) || Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_INFINITY(d) || Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    Py_RETURN_NOTIMPLEMENTED;
}

/* Generic Subtraction
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pympany_sub(PyObject *a, PyObject *b)
{
    PyObject *r = 0;
    mpz_t tempz;
    PympzObject *rz;
    PyxmpzObject *rxz;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;
    int overflow;

    if(Pympz_Check(a)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* mpz - integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Subtracting (mpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_sub(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rz;
        }

        /* mpz - mpz */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Subtracting (mpz,mpz)\n");
            mpz_sub(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }

        /* mpz - xmpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Subtracting (mpz,xmpz)\n");
            mpz_sub(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pympz_Check(b)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* integer - mpz */
        if(PyIntOrLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Subtracting (long,mpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_sub(rz->z, tempz, Pympz_AS_MPZ(b));
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_ui_sub(rz->z, temp, Pympz_AS_MPZ(b));
            } else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(b), -temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject *)rz;
        }

        /* xmpz - mpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Subtracting (xmpz,mpz)\n");
            mpz_sub(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz - integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Subtracting (xmpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_sub(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_sub_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            } else {
                mpz_add_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rxz;
        }

        /* xmpz - xmpz */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Subtracting (xmpz,xmpz)\n");
            mpz_sub(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
        }

        /* xmpz - mpz */
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Subtracting (xmpz,mpz)\n");
            mpz_sub(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rxz;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer - xmpz */
        if(PyIntOrLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Subtracting (long,xmpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_sub(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_ui_sub(rxz->z, temp, Pyxmpz_AS_MPZ(b));
            } else {
                mpz_add_ui(rxz->z, Pyxmpz_AS_MPZ(b), -temp);
                mpz_neg(rxz->z, rxz->z);
            }
            return (PyObject *)rxz;
        }

        /* mpz - xmpz */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Subtracting (mpz,xmpz)\n");
            mpz_sub(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_INFINITY(d) || Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_INFINITY(d) || Py_IS_NAN(d)) {
                    if(Py_IS_INFINITY(d))
                        r = PyFloat_FromDouble(-d);
                    else
                        r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    Py_RETURN_NOTIMPLEMENTED;
}

/* Generic Multiplication
 *
 * Follows the same conversion rules as Pympany_add.
 */

static PyObject *
Pympany_mul(PyObject *a, PyObject *b)
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

    if(Pympz_Check(a)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* mpz * integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Multiplying (mpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_mul(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(a), temp);
            }
            return (PyObject *)rz;
        }

        /* mpz * mpz */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Multiplying (mpz,mpz)\n");
            mpz_mul(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }

        /* mpz * xmpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Multiplying (mpz,xmpz)\n");
            mpz_mul(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pympz_Check(b)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* integer * mpz */
        if(PyIntOrLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Multiplying (long,mpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_mul(rz->z, Pympz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            } else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(b), temp);
            }
            return (PyObject *)rz;
        }

        /* xmpz * mpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Multiplying (mpz,xmpz)\n");
            mpz_mul(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz * integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Multiplying (xmpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_mul(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else {
                mpz_mul_si(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            }
            return (PyObject *)rxz;
        }

        /* xmpz * xmpz */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Multiplying (xmpz,xmpz)\n");
            mpz_mul(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
        }

        /* xmpz * mpz */
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Multiplying (xmpz,mpz)\n");
            mpz_mul(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rxz;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer * xmpz */
        if(PyIntOrLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Multiplying (long,xmpz)\n");
            temp = PyLong_AsLongAndOverflow(a, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, a);
                mpz_mul(rxz->z, Pyxmpz_AS_MPZ(b), tempz);
                mpz_cloc(tempz);
            } else {
                mpz_mul_si(rxz->z, Pyxmpz_AS_MPZ(b), temp);
            }
            return (PyObject *)rxz;
        }

        /* mpz * xmpz */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Multiplying (mpz,xmpz)\n");
            mpz_mul(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)pbf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        r = PyFloat_FromDouble(Py_NAN);
                    } else if(mpf_sgn(pbf->f) < 0) {
                        r = PyFloat_FromDouble(-d);
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(paf->f) == 0) {
                        r = PyFloat_FromDouble(Py_NAN);
                    } else if(mpf_sgn(paf->f) < 0) {
                        r = PyFloat_FromDouble(-d);
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)paf);
                    return r;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
Pympany_floordiv(PyObject *a, PyObject *b)
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

    if(Pympz_Check(a)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* mpz // integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (mpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject *)rz);
                return NULL;
            } else if(temp > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), -temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject *)rz;
        }

        /* mpz // mpz */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Floor divide (mpz,mpz)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }

        /* mpz // xmpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Floor divide (mpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        if(!(rz = Pympz_new()))
            return NULL;

        /* integer // mpz */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rz;
        }
#endif

        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rz;
        }

        /* xmpz // mpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Floor divide (xmpz,mpz)\n");
            mpz_fdiv_q(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz // integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
                Py_DECREF((PyObject *)rxz);
                return NULL;
            } else if(temp > 0) {
                mpz_fdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            } else {
                mpz_cdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
                mpz_neg(rxz->z, rxz->z);
            }
            return (PyObject *)rxz;
        }

        /* xmpz // xmpz */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Floor divide (xmpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
        }

        /* xmpz // mpz */
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Floor divide (xmpz,mpz)\n");
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rxz;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
            return NULL;
        }
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer // xmpz */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rxz;
        }
#endif

        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rxz;
        }

        /* mpz // xmpz */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Floor divide (mpz,xmpz)\n");
            mpz_fdiv_q(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_NAN(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else if(mpf_sgn(pbf->f) < 0) {
                        r = PyFloat_FromDouble(-d);
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    mpf_set_d(paf->f, 0.0);
                    return (PyObject*)paf;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    PyObject *r = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    unsigned int bits;

    if(Pympz_Check(b) && (mpz_sgn(Pympz_AS_MPZ(b)) == 0)) {
        PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
        return NULL;
    }

    if(Pyxmpz_Check(b) && (mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0)) {
        PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
        return NULL;
    }

    if(Pympq_Check(b) && (mpq_sgn(Pympq_AS_MPQ(b)) == 0)) {
        PyErr_SetString(PyExc_ZeroDivisionError, "mpq division by zero");
        return NULL;
    }

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
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_NAN(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else if(mpf_sgn(pbf->f) < 0) {
                        r = PyFloat_FromDouble(-d);
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    mpf_set_d(paf->f, 0.0);
                    return (PyObject*)paf;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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

    if(Pympz_Check(a)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* mpz / integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "True divide (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            } else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }

        /* mpz / mpz */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "True divide (mpz,mpz)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }

        /* mpz / xmpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "True divide (mpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        if(!(rz = Pympz_new()))
            return NULL;

        /* integer / mpz */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "True divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rz;
        }
#endif
        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "True divide (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rz;
        }

        /* xmpz / mpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "True divide (xmpz,mpz)\n");
            mpz_fdiv_q(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz / integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "True divide (xmpz,long)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            } else {
                mpz_cdiv_q_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
                mpz_neg(rxz->z, rxz->z);
            }
            return (PyObject*)rxz;
        }

        /* xmpz / xmpz */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "True divide (xmpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
        }

        /* xmpz / mpz */
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "True divide (xmpz,mpz)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_q(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rxz;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
            return NULL;
        }
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer / xmpz */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "True divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rxz;
        }
#endif
        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "True divide (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_q(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject*)rxz;
        }

        /* mpz / xmpz */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "True divide (mpz,xmpz)\n");
            mpz_fdiv_q(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject*)rxz;
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_NAN(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else if(mpf_sgn(pbf->f) < 0) {
                        r = PyFloat_FromDouble(-d);
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    mpf_set_d(paf->f, 0.0);
                    return (PyObject*)paf;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    PyObject *r = 0;
    mpz_t tempz;
    PympzObject *rz = 0;
    PyxmpzObject *rxz = 0;
    PympqObject *rq = 0, *paq = 0, *pbq = 0;
    PympfObject *rf = 0, *paf = 0, *pbf = 0;
    long temp;
    unsigned int bits;
    int overflow;

    if(Pympz_Check(a)) {
        if(!(rz = Pympz_new()))
            return NULL;

        /* mpz % integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Modulo (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_r_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            } else {
                mpz_cdiv_r_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rz;
        }

        /* mpz % mpz */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Modulo (integer,integer)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }

        /* mpz % xmpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Modulo (mpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        if(!(rz = Pympz_new()))
            return NULL;

        /* integer % mpz */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Modulo (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_r(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rz;
        }
#endif
        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Modulo (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_r(rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rz;
        }

        /* xmpz % mpz */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Modulo (xmpz,mpz)\n");
            mpz_fdiv_r(rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rz;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* xmpz % integer */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Modulo (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_r(rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_r_ui(rxz->z, Pyxmpz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz modulo by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            } else {
                mpz_cdiv_r_ui(rxz->z, Pyxmpz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rxz;
        }

        /* xmpz % xmpz */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Modulo (xmpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz modulo by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_r(rxz->z, Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
        }

        /* xmpz % mpz */
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Modulo (xmpz,mpz)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz modulo by zero");
                Py_DECREF((PyObject*)rxz);
                return NULL;
            }
            mpz_fdiv_r(rxz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject *)rxz;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
            return NULL;
        }
        if(!(rxz = Pyxmpz_new()))
            return NULL;

        /* integer % xmpz */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Modulo (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_r(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rxz;
        }
#endif
        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "Modulo (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_r(rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            return (PyObject *)rxz;
        }

        /* mpz % xmpz */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "Modulo (mpz,xmpz)\n");
            mpz_fdiv_r(rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            return (PyObject *)rxz;
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
        mpz_inoc(tempz);
        mpq_div(rq->q, paq->q, pbq->q);
        mpz_fdiv_q(tempz, mpq_numref(rq->q), mpq_denref(rq->q));
        /* Need to calculate paq - rz * pbq */
        mpq_set_z(rq->q, tempz);
        mpq_mul(rq->q, rq->q, pbq->q);
        mpq_sub(rq->q, paq->q, rq->q);
        Py_DECREF((PyObject*)paq); Py_DECREF((PyObject*)pbq);
        mpz_cloc(tempz);
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_NAN(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(pbf->f) == 0) {
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                    } else if(mpf_sgn(pbf->f) < 0) {
                        r = PyFloat_FromDouble(-d);
                    } else {
                        r = PyFloat_FromDouble(d);
                    }
                    Py_DECREF((PyObject*)pbf);
                    return r;
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_NAN(d)) {
                    r = PyFloat_FromDouble(d);
                    Py_DECREF((PyObject*)paf);
                    return r;
                } else if(Py_IS_INFINITY(d)) {
                    mpf_set_d(paf->f, 0.0);
                    return (PyObject*)paf;
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    unsigned int bits;
    int overflow;

    if(Pympz_Check(a)) {
        if(!(r=PyTuple_New(2)) || !(rz=Pympz_new()) || !(qz=Pympz_new())) {
            Py_XDECREF((PyObject*)rz);
            Py_XDECREF((PyObject*)qz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(mpz,integer) */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "divmod (mpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_qr(qz->z, rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_qr_ui(qz->z, rz->z, Pympz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
                Py_DECREF((PyObject*)rz);
                Py_DECREF((PyObject*)qz);
                return NULL;
            } else {
                mpz_cdiv_qr_ui(qz->z, rz->z, Pympz_AS_MPZ(a), -temp);
                mpz_neg(qz->z, qz->z);
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }

        /* divmod(mpz,mpz) */
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "divmod (mpz,mpz)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
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
        if((!options.prefer_mutable) && Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "divmod (mpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz divmod by zero");
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

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b))==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
            return NULL;
        }
        if(!(r=PyTuple_New(2)) || !(rz=Pympz_new()) || !(qz=Pympz_new())) {
            Py_XDECREF((PyObject*)rz);
            Py_XDECREF((PyObject*)qz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(integer,mpz) */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "divmod (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_qr(qz->z, rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
#endif
        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "divmod (integer,mpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_qr(qz->z, rz->z, tempz, Pympz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }

        /* divmod(xmpz,mpz) */
        if((!options.prefer_mutable) && Pyxmpz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "divmod (xmpz,mpz)\n");
            mpz_fdiv_qr(qz->z, rz->z, Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rz);
            return r;
        }
    }

    if(Pyxmpz_Check(a)) {
        if(!(r=PyTuple_New(2)) || !(rxz=Pyxmpz_new()) || !(qxz=Pyxmpz_new())) {
            Py_XDECREF((PyObject*)rxz);
            Py_XDECREF((PyObject*)qxz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(xmpz,integer) */
        if(PyIntOrLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "divmod (xmpz,integer)\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_qr(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_qr_ui(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz divmod by zero");
                Py_DECREF((PyObject*)rxz);
                Py_DECREF((PyObject*)qxz);
                return NULL;
            } else {
                mpz_cdiv_qr_ui(qxz->z, rxz->z, Pyxmpz_AS_MPZ(a), -temp);
                mpz_neg(qxz->z, qxz->z);
            }
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }

        /* divmod(xmpz,xmpz) */
        if(Pyxmpz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "divmod (xmpz,xmpz)\n");
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz divmod by zero");
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
        if((options.prefer_mutable) && Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "divmod (xmpz,mpz)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "xmpz divmod by zero");
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

    if(Pyxmpz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b))==0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz divmod by zero");
            return NULL;
        }
        if(!(r=PyTuple_New(2)) || !(rxz=Pyxmpz_new()) || !(qxz=Pyxmpz_new())) {
            Py_XDECREF((PyObject*)rxz);
            Py_XDECREF((PyObject*)qxz);
            Py_XDECREF(r);
            return NULL;
        }

        /* divmod(integer,xmpz) */
#ifdef PY2
        if(PyInt_Check(a)) {
            if (options.debug)
                fprintf(stderr, "divmod (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_si(tempz, PyInt_AS_LONG(a));
            mpz_fdiv_qr(qxz->z, rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }
#endif
        if(PyLong_Check(a)) {
            if (options.debug)
                fprintf(stderr, "divmod (integer,xmpz)\n");
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, a);
            mpz_fdiv_qr(qxz->z, rxz->z, tempz, Pyxmpz_AS_MPZ(b));
            mpz_cloc(tempz);
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
            return r;
        }

        /* divmod(mpz,xmpz) */
        if((options.prefer_mutable) && Pympz_Check(a)) {
            if(options.debug)
                fprintf(stderr, "divmod (xmpz,mpz)\n");
            mpz_fdiv_qr(qxz->z, rxz->z, Pympz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
            PyTuple_SET_ITEM(r, 0, (PyObject*)qxz);
            PyTuple_SET_ITEM(r, 1, (PyObject*)rxz);
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
        } else if(Pympf_Check(b)) {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, pbf->rebits);
        } else {
            pbf = anynum2Pympf(b, 0);
            paf = anynum2Pympf(a, 0);
        }
        if(!paf || !pbf) {
            if(PyErr_Occurred()) {
                PyErr_Clear();
            } else {
                PyErr_SetString(PyExc_SystemError,
                    "Internal error status is confused.");
                return NULL;
            }
            /* Need to handle special float values. */
            if(pbf && !paf && PyFloat_Check(a)) {
                double d = PyFloat_AS_DOUBLE(a);
                if(Py_IS_INFINITY(d) || Py_IS_NAN(d)) {
                    /* divmod(inf|nan, number) */
                    if(mpf_sgn(pbf->f) == 0) {
                        /* if number == 0, raise ZeroDivisionError */
                        PyErr_SetString(PyExc_ZeroDivisionError,
                            "mpf division by zero");
                        r = NULL;
                        Py_DECREF((PyObject*)pbf);
                        return r;
                    } else {
                        /* if number != 0, return (nan,nan) */
                        Py_DECREF((PyObject*)pbf);
                        paf = (PympfObject*)PyFloat_FromDouble(Py_NAN);
                        pbf = (PympfObject*)PyFloat_FromDouble(Py_NAN);
                        return Py_BuildValue("(NN)", paf, pbf);
                    }
                }
            } else if(paf && !pbf && PyFloat_Check(b)) {
                /* divmod(number, inf|nan) */
                double d = PyFloat_AS_DOUBLE(b);
                if(Py_IS_NAN(d)) {
                    Py_DECREF((PyObject*)paf);
                    paf = (PympfObject*)PyFloat_FromDouble(d);
                    pbf = (PympfObject*)PyFloat_FromDouble(d);
                    return Py_BuildValue("(NN)", paf, pbf);
                } else if(Py_IS_INFINITY(d)) {
                    if(mpf_sgn(paf->f) == 0) {
                        /* if number == 0, return (0.0, 0.0) */
                        pbf = Pympf_new(paf->rebits);
                        mpf_set_d(pbf->f, 0.0);
                        mpf_set_d(paf->f, 0.0);
                        return Py_BuildValue("(NN)", paf, pbf);
                    } else if(mpf_sgn(paf->f) < 0) {
                        if(d < 0) {
                            /* return (0, number) */
                            /* XXX Not checking out-of-memory! */
                            pbf = Pympf_new(paf->rebits);
                            mpf_set_d(pbf->f, 0.0);
                            return Py_BuildValue("(NN)", pbf, paf);
                        } else {
                            /* return (-1.0, inf) */
                            pbf = Pympf_new(paf->rebits);
                            mpf_set_d(pbf->f, -1.0);
                            Py_DECREF((PyObject*)paf);
                            paf = (PympfObject*)PyFloat_FromDouble(d);
                            return Py_BuildValue("(NN)", pbf, paf);
                        }
                    } else {
                        if(d > 0) {
                            /* return (0, number) */
                            pbf = Pympf_new(paf->rebits);
                            mpf_set_d(pbf->f, 0.0);
                            return Py_BuildValue("(NN)", pbf, paf);
                        } else {
                            /* return (-1.0, inf) */
                            pbf = Pympf_new(paf->rebits);
                            mpf_set_d(pbf->f, -1.0);
                            Py_DECREF((PyObject*)paf);
                            paf = (PympfObject*)PyFloat_FromDouble(d);
                            return Py_BuildValue("(NN)", pbf, paf);
                        }
                    }
                }
            } else {
                PyErr_SetString(PyExc_SystemError, "Can not convert number to mpf");
                Py_XDECREF((PyObject*)paf);
                Py_XDECREF((PyObject*)pbf);
                return NULL;
            }
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
    Py_RETURN_NOTIMPLEMENTED;
}
