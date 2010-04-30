/* gmpy_xmpz_inplace.c
 *
 * Provides inplace mutating operations for xmpz
 *
 * This file should be considered part of gmpy.c
 */

#include <math.h>

/* Inplace xmpz addition.
 */

static PyObject *
Pyxmpz_inplace_add(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    long temp;
    int overflow;

    /* Try to make mpz + small_int faster */
    if(PyIntOrLong_Check(b)) {
        if(options.debug)
            fprintf(stderr, "Adding (xmpz,long)\n");
        temp = PyLong_AsLongAndOverflow(b, &overflow);
        if(overflow) {
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, b);
            mpz_add(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), tempz);
            mpz_cloc(tempz);
        } else if(temp >= 0) {
            mpz_add_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        } else {
            mpz_sub_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), -temp);
        }
        Py_INCREF(a);
        return a;
    }

    if(Pyxmpz_Check(b)) {
        mpz_add(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(Pympz_Check(b)) {
        mpz_add(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(options.debug)
        fprintf(stderr, "Pyxmpz_inplace_add returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace mpz subtraction.
 */

static PyObject *
Pyxmpz_inplace_sub(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    long temp;
    int overflow;

    if(PyIntOrLong_Check(b)) {
        if (options.debug)
            fprintf(stderr, "Subtracting (xmpz,long)\n");
        temp = PyLong_AsLongAndOverflow(b, &overflow);
        if(overflow) {
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, b);
            mpz_sub(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), tempz);
            mpz_cloc(tempz);
        } else if(temp >= 0) {
            mpz_sub_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        } else {
            mpz_add_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), -temp);
        }
        Py_INCREF(a);
        return a;
    }

    if(Pyxmpz_Check(b)) {
        mpz_sub(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(Pympz_Check(b)) {
        mpz_sub(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(!options.debug)
        fprintf(stderr, "Pyxmpz_inplace_sub returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz multiplication.
 */

static PyObject *
Pyxmpz_inplace_mul(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    long temp;
    int overflow;

    if(PyIntOrLong_Check(b)) {
        if (options.debug)
            fprintf(stderr, "Multiplying (xmpz,long)\n");
        temp = PyLong_AsLongAndOverflow(b, &overflow);
        if(overflow) {
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, b);
            mpz_mul(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), tempz);
            mpz_cloc(tempz);
        } else {
            mpz_mul_si(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        }
        Py_INCREF(a);
        return a;
    }

    if(Pyxmpz_Check(b)) {
        mpz_mul(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(Pympz_Check(b)) {
        mpz_mul(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(!options.debug)
        fprintf(stderr, "Pyxmpz_inplace_mul returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Pympany_floordiv follows the // semantics from Python 3.x. The result is
 * an mpz when the arguments are mpz or mpq, but the result is an mpf when
 * the arguments are mpf.
 */

static PyObject *
Pyxmpz_inplace_floordiv(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    long temp;
    int overflow;

    if(PyIntOrLong_Check(b)) {
        if (options.debug)
            fprintf(stderr, "Floor divide (xmpz,long)\n");
        temp = PyLong_AsLongAndOverflow(b, &overflow);
        if(overflow) {
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, b);
            mpz_fdiv_q(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), tempz);
            mpz_cloc(tempz);
        } else if(temp == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
            return NULL;
        } else if(temp > 0) {
            mpz_fdiv_q_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        } else {
            mpz_cdiv_q_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), -temp);
            mpz_neg(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a));
        }
        Py_INCREF(a);
        return a;
    }

    if(Pyxmpz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
            return NULL;
        }
        mpz_fdiv_q(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz division by zero");
            return NULL;
        }
        mpz_fdiv_q(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(options.debug)
        fprintf(stderr, "Pyxmpz_inplace_floordiv returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz remainder.
 */

static PyObject *
Pyxmpz_inplace_rem(PyObject *a, PyObject *b)
{
    mpz_t tempz;
    long temp;
    int overflow;

    if(PyIntOrLong_Check(b)) {
        if (options.debug)
            fprintf(stderr, "Modulo (xmpz,long)\n");
        temp = PyLong_AsLongAndOverflow(b, &overflow);
        if(overflow) {
            mpz_inoc(tempz);
            mpz_set_PyLong(tempz, b);
            mpz_fdiv_r(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), tempz);
            mpz_cloc(tempz);
        } else if(temp > 0) {
            mpz_fdiv_r_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        } else if(temp == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz modulo by zero");
            return NULL;
        } else {
            mpz_cdiv_r_ui(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), -temp);
        }
        Py_INCREF(a);
        return a;
    }

    if(Pyxmpz_Check(b)) {
        if(options.debug)
            fprintf(stderr, "Modulo (integer,integer)\n");
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz modulo by zero");
            return NULL;
        }
        mpz_fdiv_r(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(Pympz_Check(b)) {
        if(options.debug)
            fprintf(stderr, "Modulo (integer,integer)\n");
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) == 0) {
            PyErr_SetString(PyExc_ZeroDivisionError, "xmpz modulo by zero");
            return NULL;
        }
        mpz_fdiv_r(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), Pympz_AS_MPZ(b));
        Py_INCREF(a);
        return a;
    }

    if(options.debug)
        fprintf(stderr, "Pyxmpz_inplace_rem returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz rshift.
 */

static PyObject *
Pyxmpz_inplace_rshift(PyObject *a, PyObject *b)
{
    long temp;
    int overflow;

    /* Try to make mpz + small_int faster */
        if(PyIntOrLong_Check(b)) {
            if(options.debug)
                fprintf(stderr, "right shift\n");
            temp = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                PyErr_SetString(PyExc_ValueError, "outrageous shift count");
                return NULL;
            } else if(temp >= 0) {
                mpz_fdiv_q_2exp(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
                Py_INCREF(a);
                return a;
            } else {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                return NULL;
            }
        }

        if(Pyxmpz_Check(b)) {
            if(mpz_sgn(Pyxmpz_AS_MPZ(b)) < 0) {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                return NULL;
            }
            if(!mpz_fits_slong_p(Pyxmpz_AS_MPZ(b))) {
                PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
                return NULL;
            }
            temp = mpz_get_si(Pyxmpz_AS_MPZ(b));
            mpz_fdiv_q_2exp(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
            Py_INCREF(a);
            return a;
        }

        if(Pympz_Check(b)) {
            if(mpz_sgn(Pympz_AS_MPZ(b)) < 0) {
                PyErr_SetString(PyExc_ValueError, "negative shift count");
                return NULL;
            }
            if(!mpz_fits_slong_p(Pympz_AS_MPZ(b))) {
                PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
                return NULL;
            }
            temp = mpz_get_si(Pympz_AS_MPZ(b));
            mpz_fdiv_q_2exp(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
            Py_INCREF(a);
            return a;
        }

    if(options.debug)
        fprintf(stderr, "Pyxmpz_inplace_rshift returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz lshift.
 */

static PyObject *
Pyxmpz_inplace_lshift(PyObject *a, PyObject *b)
{
    long temp;
    int overflow;

    /* Try to make mpz + small_int faster */
    if(PyIntOrLong_Check(b)) {
        if(options.debug)
            fprintf(stderr, "left shift\n");
        temp = PyLong_AsLongAndOverflow(b, &overflow);
        if(overflow) {
            PyErr_SetString(PyExc_ValueError, "outrageous shift count");
            return NULL;
        } else if(temp >= 0) {
            mpz_mul_2exp(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        } else {
            PyErr_SetString(PyExc_ValueError, "negative shift count");
            return NULL;
        }
    }

    if(Pyxmpz_Check(b)) {
        if(mpz_sgn(Pyxmpz_AS_MPZ(b)) < 0) {
            PyErr_SetString(PyExc_ValueError, "negative shift count");
            return NULL;
        }
        if(!mpz_fits_slong_p(Pyxmpz_AS_MPZ(b))) {
            PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
            return NULL;
        }
        temp = mpz_get_si(Pyxmpz_AS_MPZ(b));
        mpz_mul_2exp(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        Py_INCREF(a);
        return a;
    }

    if(Pympz_Check(b)) {
        if(mpz_sgn(Pympz_AS_MPZ(b)) < 0) {
            PyErr_SetString(PyExc_ValueError, "negative shift count");
            return NULL;
        }
        if(!mpz_fits_slong_p(Pympz_AS_MPZ(b))) {
            PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
            return NULL;
        }
        temp = mpz_get_si(Pympz_AS_MPZ(b));
        mpz_mul_2exp(Pyxmpz_AS_MPZ(a), Pyxmpz_AS_MPZ(a), temp);
        Py_INCREF(a);
        return a;
    }

    if(options.debug)
        fprintf(stderr, "Pyxmpz_inplace_lshift returned NotImplemented\n");
    Py_RETURN_NOTIMPLEMENTED;
}

/* Inplace xmpz_pow.
 */

static PyObject *
Pyxmpz_inplace_pow(PyObject *in_b, PyObject *in_e, PyObject *in_m)
{
    PympzObject *e = 0;
    unsigned long el;

    if(options.debug)
        fprintf(stderr, "Pyxmpz_inplace_pow\n");

    if(!Pyxmpz_Check(in_b)) {
        PyErr_SetString(PyExc_TypeError, "bogus base type");
        return NULL;
    }
    if(in_m != Py_None) {
        PyErr_SetString(PyExc_SystemError, "modulo not expected");
        return NULL;
    }
    e = Pympz_From_Integer(in_e);
    if(!e) {
        PyErr_SetString(PyExc_TypeError, "expected an integer exponent");
        return NULL;
    }
    if(mpz_sgn(e->z) < 0) {
        PyErr_SetString(PyExc_ValueError, "xmpz.pow with negative power");
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    if(!mpz_fits_ulong_p(e->z)) {
        PyErr_SetString(PyExc_ValueError, "xmpz.pow outrageous exponent");
        Py_DECREF((PyObject*)e);
        return NULL;
    }
    el = mpz_get_ui(e->z);
    mpz_pow_ui(Pyxmpz_AS_MPZ(in_b), Pyxmpz_AS_MPZ(in_b), el);
    Py_DECREF((PyObject*)e);
    Py_INCREF((PyObject*)in_b);
    return (PyObject*)in_b;
}

