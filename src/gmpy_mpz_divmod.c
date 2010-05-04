/* gmpy_mpz_divmod.c
 *
 * This file should be considered part of gmpy.c
 *
 * This file contains functions related to division and remainder.
 * Functions are optimized by writing distinct functions for
 * gmpy2.function versus mpz.function.
 */

/*
 **************************************************************************
 * Floor division and remainder.
 **************************************************************************
 */

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_cdivmodg[]="\
cdivmod(x,y): returns the quotient of x divided by y. The quotient\n\
is rounded towards +Inf and the remainder will have the opposite\n\
sign to y. x and y must be mpz, or else get coerced to mpz.\n\
";
static PyObject *
Pympz_cdivmod(PyObject *self, PyObject *args)
{
    PyObject *other, *result;
    PympzObject *quot, *rem;

    PARSE_TWO_MPZ(other, "cdivmod() requires 'mpz','mpz' arguments");

    quot = Pympz_new();
    rem = Pympz_new();
    result = PyTuple_New(2);
    if(!quot || !rem || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)quot);
        Py_XDECREF((PyObject*)rem);
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpz_cdiv_qr(quot->z, rem->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    PyTuple_SET_ITEM(result, 0, (PyObject*)quot);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

static char doc_fdivmodg[]="\
fdivmod(x,y): returns the quotient of x divided by y. The quotient\n\
is rounded towards -Inf and the remainder will have the same sign\n\
as y. x and y must be mpz, or else get coerced to mpz.\n\
";
static PyObject *
Pympz_fdivmod(PyObject *self, PyObject *args)
{
    PyObject *other, *result;
    PympzObject *quot, *rem;

    PARSE_TWO_MPZ(other, "fdivmod() requires 'mpz','mpz' arguments");

    quot = Pympz_new();
    rem = Pympz_new();
    result = PyTuple_New(2);
    if(!quot || !rem || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)quot);
        Py_XDECREF((PyObject*)rem);
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpz_fdiv_qr(quot->z, rem->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    PyTuple_SET_ITEM(result, 0, (PyObject*)quot);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

static char doc_tdivmodg[]="\
tdivmod(x,y): returns the quotient of x divided by y. The quotient\n\
is rounded towards zero and the remaider will have the same sign\n\
as x. x and y must be mpz, or else get coerced to mpz.\n\
";
static PyObject *
Pympz_tdivmod(PyObject *self, PyObject *args)
{
    PyObject *other, *result;
    PympzObject *quot, *rem;

    PARSE_TWO_MPZ(other, "tdivmod() requires 'mpz','mpz' arguments");

    quot = Pympz_new();
    rem = Pympz_new();
    result = PyTuple_New(2);
    if(!quot || !rem || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)quot);
        Py_XDECREF((PyObject*)rem);
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpz_tdiv_qr(quot->z, rem->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    PyTuple_SET_ITEM(result, 0, (PyObject*)quot);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

