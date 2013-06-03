/* gmpy_mpz_inplace.c
 *
 * Provides inplace operations for mpz
 *
 * NOTE: These functions do NOT mutate the mpz.
 *
 * This file should be considered part of gmpy.c.
 */

#include <math.h>

#define Py_RETURN_NOTIMPLEMENTED\
    return Py_INCREF(Py_NotImplemented), Py_NotImplemented

/* Inplace mpz addition. Does NOT mutate!
 */

static PyObject *
Pympz_inplace_add(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    mpz_t tempz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    /* Try to make mpz + small_int faster */
    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {
#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (mpz,small_int)\n");
            if((temp = PyInt_AS_LONG(b)) >= 0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rz;
        }
#endif
        if(PyLong_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Adding (mpz,long)\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
                PyErr_Clear();
#endif
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_add(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rz;
        }

        if(Pympz_Check(b)) {
            mpz_add(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }
    if(options.debug)
        fprintf(stderr, "Pympz_inplace_add returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz subtraction. Does NOT mutate!
 */

static PyObject *
Pympz_inplace_sub(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    mpz_t tempz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {

#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Subtracting (mpz,small_int)\n");
            if((temp = PyInt_AS_LONG(b)) >= 0) {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rz;
        }
#endif

        if(PyLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Subtracting (mpz,long)\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
                PyErr_Clear();
#endif
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_sub(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp >= 0) {
                mpz_sub_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_add_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rz;
        }

        if(Pympz_Check(b)) {
            mpz_sub(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }
    if(!options.debug)
        fprintf(stderr, "Pympz_inplace_sub returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz multiplication. Does NOT mutate!
 */

static PyObject *
Pympz_inplace_mul(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    mpz_t tempz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {

#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Multiplying (mpz,small_int)\n");
            mpz_mul_si(rz->z, Pympz_AS_MPZ(a), PyInt_AS_LONG(b));
            return (PyObject *)rz;
        }
#endif

        if(PyLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Multiplying (mpz,long)\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
                PyErr_Clear();
#endif
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_mul(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else {
                mpz_mul_si(rz->z, Pympz_AS_MPZ(a), temp);
            }
            return (PyObject*)rz;
        }
        if(Pympz_Check(b)) {
            mpz_mul(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }
    if(!options.debug)
        fprintf(stderr, "Pympz_inplace_mul returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
Pympz_inplace_floordiv(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    mpz_t tempz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {

#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (mpz,small_int)\n");
            if((temp=PyInt_AS_LONG(b)) > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                Py_DECREF((PyObject *)rz);
                return NULL;
            } else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), -temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject *)rz;
        }
#endif

        if(PyLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Floor divide (mpz,long)\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
                PyErr_Clear();
#endif
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                return NULL;
            } else if(temp > 0) {
                mpz_fdiv_q_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else {
                mpz_cdiv_q_ui(rz->z, Pympz_AS_MPZ(a), -temp);
                mpz_neg(rz->z, rz->z);
            }
            return (PyObject*)rz;
        }
        if(Pympz_Check(b)) {
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz division by zero");
                return NULL;
            }
            mpz_fdiv_q(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }
    if(options.debug)
        fprintf(stderr, "Pympz_inplace_floordiv returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz remainder. Does NOT mutate!
 */

static PyObject *
Pympz_inplace_rem(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    mpz_t tempz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {

#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Modulo (mpz,small_int)\n");
            if((temp=PyInt_AS_LONG(b)) > 0) {
                mpz_fdiv_r_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                Py_DECREF((PyObject *)rz);
                return NULL;
            } else {
                mpz_cdiv_r_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject *)rz;
        }
#endif

        if(PyLong_Check(b)) {
            if (options.debug)
                fprintf(stderr, "Modulo (mpz,long)\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
                PyErr_Clear();
#endif
                mpz_inoc(tempz);
                mpz_set_PyLong(tempz, b);
                mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), tempz);
                mpz_cloc(tempz);
            } else if(temp > 0) {
                mpz_fdiv_r_ui(rz->z, Pympz_AS_MPZ(a), temp);
            } else if(temp == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                return NULL;
            } else {
                mpz_cdiv_r_ui(rz->z, Pympz_AS_MPZ(a), -temp);
            }
            return (PyObject*)rz;
        }
        if(Pympz_Check(b)) {
            if(options.debug)
                fprintf(stderr, "Modulo (integer,integer)\n");
            if(mpz_sgn(Pympz_AS_MPZ(b)) == 0) {
                PyErr_SetString(PyExc_ZeroDivisionError, "mpz modulo by zero");
                return NULL;
            }
            mpz_fdiv_r(rz->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
            return (PyObject*)rz;
        }
    }
    if(options.debug)
        fprintf(stderr, "Pympz_inplace_rem returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz rshift. Does NOT mutate!
 */

static PyObject *
Pympz_inplace_rshift(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    /* Try to make mpz + small_int faster */
    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {
#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if(options.debug)
                fprintf(stderr, "right shift\n");
            if((temp = PyInt_AS_LONG(b)) >= 0) {
                mpz_fdiv_q_2exp(rz->z, Pympz_AS_MPZ(a), temp);
                return (PyObject *)rz;
            } else {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
        }
#endif
        if(PyLong_Check(b)) {
            if(options.debug)
                fprintf(stderr, "right shift\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
#endif
                PyErr_SetString(PyExc_ValueError, "outrageous shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            } else if(temp >= 0) {
                mpz_fdiv_q_2exp(rz->z, Pympz_AS_MPZ(a), temp);
                return (PyObject *)rz;
            } else {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
        }

        if(Pympz_Check(b)) {
            if(mpz_sgn(Pympz_AS_MPZ(b)) < 0) {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            if(!mpz_fits_slong_p(Pympz_AS_MPZ(b))) {
                PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            temp = mpz_get_si(Pympz_AS_MPZ(b));
            mpz_fdiv_q_2exp(rz->z, Pympz_AS_MPZ(a), temp);
            return (PyObject*)rz;
        }
    }
    if(options.debug)
        fprintf(stderr, "Pympz_inplace_rshift returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz lshift. Does NOT mutate!
 */

static PyObject *
Pympz_inplace_lshift(PyObject *a, PyObject *b)
{
    PympzObject *rz;
    long temp;
#if PY_MAJOR_VERSION == 3
    int overflow;
#endif

    /* Try to make mpz + small_int faster */
    if(!(rz = Pympz_new()))
        return NULL;
    if(Pympz_Check(a)) {
#if PY_MAJOR_VERSION == 2
        if(PyInt_Check(b)) {
            if(options.debug)
                fprintf(stderr, "left shift\n");
            if((temp = PyInt_AS_LONG(b)) >= 0) {
                mpz_mul_2exp(rz->z, Pympz_AS_MPZ(a), temp);
                return (PyObject *)rz;
            } else {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
        }
#endif
        if(PyLong_Check(b)) {
            if(options.debug)
                fprintf(stderr, "left shift\n");
#if PY_MAJOR_VERSION == 3
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
#else
            temp = PyLong_AsLong(b);
            if(PyErr_Occurred()) {
#endif
                PyErr_SetString(PyExc_ValueError, "outrageous shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            } else if(temp >= 0) {
                mpz_mul_2exp(rz->z, Pympz_AS_MPZ(a), temp);
                return (PyObject *)rz;
            } else {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
        }

        if(Pympz_Check(b)) {
            if(mpz_sgn(Pympz_AS_MPZ(b)) < 0) {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            if(!mpz_fits_slong_p(Pympz_AS_MPZ(b))) {
                PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
                Py_DECREF((PyObject*)rz);
                return NULL;
            }
            temp = mpz_get_si(Pympz_AS_MPZ(b));
            mpz_mul_2exp(rz->z, Pympz_AS_MPZ(a), temp);
            return (PyObject*)rz;
        }
    }
    if(options.debug)
        fprintf(stderr, "Pympz_inplace_lshift returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz_pow. Does NOT mutate.
 */
static PyObject *
Pympany_pow(PyObject *in_b, PyObject *in_e, PyObject *in_m);

static PyObject *
Pympz_inplace_pow(PyObject *in_b, PyObject *in_e, PyObject *in_m)
{
    PympzObject *r, *e = 0;
    unsigned long el;

    if(options.debug)
        fprintf(stderr, "Pympz_inplace_pow\n");

    if(!Pympz_Check(in_b)) {
        PyErr_SetString(PyExc_TypeError, "bogus base type");
        return NULL;
    }

    e = Pympz_From_Integer(in_e);
    if(!e || (in_m != Py_None)) {
        PyErr_Clear();
        Py_XDECREF((PyObject*)e);
        return Pympany_pow(in_b, in_e, in_m);
    }
    if(mpz_sgn(e->z) < 0) {
        PyErr_SetString(PyExc_ValueError, "mpz.pow with negative power");
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    if(!mpz_fits_ulong_p(e->z)) {
        PyErr_SetString(PyExc_ValueError, "mpz.pow outrageous exponent");
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    if(!(r = Pympz_new())) {
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    el = mpz_get_ui(e->z);
    mpz_pow_ui(r->z, Pympz_AS_MPZ(in_b), el);
    Py_DECREF((PyObject*)e);
    return (PyObject*)r;
}

