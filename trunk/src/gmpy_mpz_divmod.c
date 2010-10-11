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
 * Ceiling division and remainder.
 **************************************************************************
 */

static char doc_cdivmodg[]="\
cdivmod(x,y): returns the quotient and remainder of x divided by\n\
y. The quotient is rounded towards +Inf and the remainder will\n\
have the opposite sign of y. x and y must be integers.\n\
";

static PyObject *
Pygmpy_cdivmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *q, *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cdivmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE_TWO_MPZ_TUPLE(q, r, result);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cdivmod() division by 0");
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_qr(q->z, r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("cdivmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cdivmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_qr(q->z, r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_cdivg[]="\
cdiv(x,y): returns the quotient of x divided by y. The quotient\n\
is rounded towards +Inf. x and y must be integers.\n\
";

static PyObject *
Pygmpy_cdiv(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *q, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cdiv() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(q = Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_cdiv_q(q->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("cdiv() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_cdiv_q(q->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)q;
}

static char doc_cmodg[]="\
cmod(x,y): returns the remainder of x divided by y. The remainder\n\
will have the opposite sign of y. x and y must be integers.\n\
";

static PyObject *
Pygmpy_cmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(r = Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_cdiv_r(r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_cdiv_r(r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)r;
}

static char doc_cdivm[]="\
x.cdiv(y): returns the quotient of x divided by y. The quotient is\n\
rounded towards +Inf. x and y must be integers.\n\
";

static PyObject *
Pympz_cdiv(PyObject *self, PyObject *other)
{
    PympzObject *result, *tempx;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            return NULL;
        }
        if (!(result = Pympz_new()))
            return NULL;
        mpz_cdiv_q(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("cdiv() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if (!(result = Pympz_new())) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_cdiv_q(result->z, Pympz_AS_MPZ(self), tempx->z);
        Py_DECREF((PyObject*)tempx);
    }
    return (PyObject*)result;
}

static char doc_cmodm[]="\
x.cmod(y): returns the remainder of x divided by y. The remainder will\n\
have the opposite sign of y. x and y must be integers.\n\
";

static PyObject *
Pympz_cmod(PyObject *self, PyObject *other)
{
    PympzObject *result, *tempx;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            return NULL;
        }
        if (!(result = Pympz_new()))
            return NULL;
        mpz_cdiv_r(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(!(result = Pympz_new())) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_cdiv_r(result->z, Pympz_AS_MPZ(self), tempx->z);
        Py_DECREF((PyObject*)tempx);
    }
    return (PyObject*)result;
}

/*
 **************************************************************************
 * Floor division and remainder.
 **************************************************************************
 */

static char doc_fdivmodg[]="\
fdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards -Inf and the remainder will have the same\n\
sign as y. x and y must be integers.\n\
";
static PyObject *
Pygmpy_fdivmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *q, *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fdivmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE_TWO_MPZ_TUPLE(q, r, result);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fdivmod() division by 0");
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(q->z, r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("fdivmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fdivmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(q->z, r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_fdivg[]="\
fdiv(x,y): returns the quotient of x divided by y. The quotient is\n\
rounded towards -Inf. x and y must be integers.\n\
";

static PyObject *
Pygmpy_fdiv(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *q, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fdiv() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(q = Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_fdiv_q(q->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("fdiv() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_fdiv_q(q->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)q;
}

static char doc_fmodg[]="\
fmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as y. x and y must be integers.\n\
";

static PyObject *
Pygmpy_fmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(r = Pympz_new()))
        return NULL;

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_fdiv_r(r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_fdiv_r(r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)r;
}

static char doc_fdivm[]="\
x.fdiv(y): returns the quotient of x divided by y. The quotient is rounded\n\
towards -Inf. x and y must be integers.\n\
";

static PyObject *
Pympz_fdiv(PyObject *self, PyObject *other)
{
    PympzObject *result, *tempx;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            return NULL;
        }
        if (!(result = Pympz_new()))
            return NULL;
        mpz_fdiv_q(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("fdiv() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if (!(result = Pympz_new())) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_fdiv_q(result->z, Pympz_AS_MPZ(self), tempx->z);
        Py_DECREF((PyObject*)tempx);
    }
    return (PyObject*)result;
}

static char doc_fmodm[]="\
x.fmod(y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as y. x and y must be integers. Will mutate x if it\n\
is an 'xmpz'.\n\
";

static PyObject *
Pympz_fmod(PyObject *self, PyObject *other)
{
    PympzObject *result, *tempx;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            return NULL;
        }
        if(!(result = Pympz_new()))
            return NULL;
        mpz_fdiv_r(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(!(result = Pympz_new())) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_fdiv_r(result->z, Pympz_AS_MPZ(self), tempx->z);
        Py_DECREF((PyObject*)tempx);
    }
    return (PyObject*)result;
}

/*
 **************************************************************************
 * Truncating division and remainder.
 **************************************************************************
 */

static char doc_tdivmodg[]="\
tdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards zero and the remainder will have the same\n\
sign as x. x and y must be integers.\n\
";
static PyObject *
Pygmpy_tdivmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *q, *r, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tdivmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE_TWO_MPZ_TUPLE(q, r, result);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tdivmod() division by 0");
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_qr(q->z, r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("tdivmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tdivmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            Py_DECREF((PyObject*)r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_qr(q->z, r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_tdivg[]="\
tdiv(x,y): returns the quotient of x divided by y. The quotient is \n\
rounded towards 0. x and y must be integers.\n\
";

static PyObject *
Pygmpy_tdiv(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *q, *tempx, *tempy;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tdiv() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(q = Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_tdiv_q(q->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("tdiv() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)q);
            return NULL;
        }
        mpz_tdiv_q(q->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)q;
}

static char doc_tmodg[]="\
tmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as x. x and y must be integers..\n\
";

static PyObject *
Pygmpy_tmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *r, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if (!(r = Pympz_new()))
        return NULL;

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_tdiv_r(r->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)r);
            return NULL;
        }
        mpz_tdiv_r(r->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)r;
}

static char doc_tdivm[]="\
x.tdiv(y): returns the quotient of x divided by y. The quotient is rounded\n\
towards 0. x and y must be integers.\n\
";

static PyObject *
Pympz_tdiv(PyObject *self, PyObject *other)
{
    PympzObject *result, *tempx;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            return NULL;
        }
        if(!(result = Pympz_new()))
            return NULL;
        mpz_tdiv_q(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("tdiv() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if (!(result = Pympz_new())) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_tdiv_q(result->z, Pympz_AS_MPZ(self), tempx->z);
        Py_DECREF((PyObject*)tempx);
    }
    return (PyObject*)result;
}

static char doc_tmodm[]="\
x.tmod(y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as x. x and y must be integers.\n\
";

static PyObject *
Pympz_tmod(PyObject *self, PyObject *other)
{
    PympzObject *result, *tempx;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            return NULL;
        }
        if (!(result = Pympz_new()))
            return NULL;
        mpz_tdiv_r(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(!(result = Pympz_new())) {
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_tdiv_r(result->z, Pympz_AS_MPZ(self), tempx->z);
        Py_DECREF((PyObject*)tempx);
    }
    return (PyObject*)result;
}

