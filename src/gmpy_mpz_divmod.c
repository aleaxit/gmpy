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
cdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards +Inf and the remainder will have the opposite\n\
sign of y. x and y must be integers. If x is an 'xmpz', the results will\n\
be 'xmpz' and neither x or y will be mutated. Otherwise, the results will\n\
be 'mpz'.\n\
";

static PyObject *
Pygmpy_cdivmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *q, *r, *result;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cdivmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        q = (PyObject*)Pyxmpz_new();
        r = (PyObject*)Pyxmpz_new();
    } else {
        q = (PyObject*)Pympz_new();
        r = (PyObject*)Pympz_new();
    }
    result = PyTuple_New(2);
    if(!q || !r || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)q);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cdivmod() division by 0");
            return NULL;
        }
        mpz_cdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("cdivmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cdivmod() division by 0");
            return NULL;
        }
        mpz_cdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_cdivg[]="\
cdiv(x,y): returns the quotient of x divided by y. The quotient is rounded\n\
to +Inf. x and y must be integers. If x is an 'xmpz', the result will be an\n\
'xmpz' and x will NOT be mutated. Otherwise, the result will be an 'mpz'.\n\
";

static PyObject *
Pygmpy_cdiv(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *q;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cdiv() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        q = (PyObject*)Pyxmpz_new();
    } else {
        q = (PyObject*)Pympz_new();
    }
    if(!q) {
        Py_XDECREF((PyObject*)q);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            return NULL;
        }
        mpz_cdiv_q(Pympz_AS_MPZ(q), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("cdiv() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            return NULL;
        }
        mpz_cdiv_q(Pympz_AS_MPZ(q), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    return q;
}

static char doc_cmodg[]="\
cmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the opposite sign of y. x and y must be integers. If x is an 'xmpz',\n\
the result will be an 'xmpz' and x will NOT be mutated. Otherwise, the\n\
result will be 'mpz'.\n\
";

static PyObject *
Pygmpy_cmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *r;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        r = (PyObject*)Pyxmpz_new();
    } else {
        r = (PyObject*)Pympz_new();
    }
    if(!r) {
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            return NULL;
        }
        mpz_cdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            return NULL;
        }
        mpz_cdiv_r(Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    return r;
}

static char doc_cmodm[]="\
x.cmod(y): returns the remainder of x divided by y. The remainder will\n\
have the opposite sign of y. x and y must be integers. If x is an 'xmpz',\n\
the result will replace x and None will be returned. Otherwise, the result\n\
will be an 'mpz'.\n\
";

static PyObject *
Pympz_cmod(PyObject *self, PyObject *other)
{
    PyObject *r;
    PympzObject *temp;

    if(Pyxmpz_Check(self)) {
        r = self;
    } else {
        r = (PyObject*)Pympz_new();
        if(!r) {
            Py_XDECREF((PyObject*)r);
            return NULL;
        }
    }

    if((Pympz_Check(other) || Pyxmpz_Check(other))) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            return NULL;
        }
        mpz_cdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    } else {
        temp = Pympz_From_Integer(other);
        if(!temp) {
            TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)temp);
            Py_XDECREF((PyObject*)r);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(temp)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            return NULL;
        }
        mpz_cdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), temp->z);
        Py_DECREF(temp);
    }
    if(Pyxmpz_Check(self))
        Py_RETURN_NONE;
    else
        return r;
}

static char doc_fdivmodg[]="\
fdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards -Inf and the remainder will have the same sign\n\
as y. x and y must be integers. If x is an 'xmpz', the results will be\n\
'xmpz' and neither x or y will be mutated. Otherwise, the results will be\n\
'mpz'\n\
";
static PyObject *
Pygmpy_fdivmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *q, *r, *result;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fdivmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        q = (PyObject*)Pyxmpz_new();
        r = (PyObject*)Pyxmpz_new();
    } else {
        q = (PyObject*)Pympz_new();
        r = (PyObject*)Pympz_new();
    }
    result = PyTuple_New(2);
    if(!q || !r || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)q);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fdivmod() division by 0");
            return NULL;
        }
        mpz_fdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
        PyTuple_SET_ITEM(result, 0, (PyObject*)q);
        PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("fdivmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fdivmod() division by 0");
            return NULL;
        }
        mpz_fdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_fdivg[]="\
fdiv(x,y): returns the quotient of x divided by y. The quotient is rounded\n\
to -Inf. x and y must be integers. If x is an 'xmpz', the result will be an\n\
'xmpz' and x will NOT be mutated. Otherwise, the result will be an 'mpz'.\n\
";

static PyObject *
Pygmpy_fdiv(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *q;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fdiv() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        q = (PyObject*)Pyxmpz_new();
    } else {
        q = (PyObject*)Pympz_new();
    }
    if(!q) {
        Py_XDECREF((PyObject*)q);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            return NULL;
        }
        mpz_fdiv_q(Pympz_AS_MPZ(q), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("fdiv() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            return NULL;
        }
        mpz_fdiv_q(Pympz_AS_MPZ(q), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    return q;
}

static char doc_fmodg[]="\
fmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as y. x and y must be integers. If x is an 'xmpz',\n\
the result will be an 'xmpz' and x will NOT be mutated. Otherwise, the\n\
result will be 'mpz'.\n\
";

static PyObject *
Pygmpy_fmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *r;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        r = (PyObject*)Pyxmpz_new();
    } else {
        r = (PyObject*)Pympz_new();
    }
    if(!r) {
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            return NULL;
        }
        mpz_fdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            return NULL;
        }
        mpz_fdiv_r(Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    return r;
}

static char doc_fmodm[]="\
x.fmod(y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as y. x and y must be integers. If x is an 'xmpz', the\n\
result will replace x and None will be returned. Otherwise, the result will\n\
be an 'mpz'.\n\
";

static PyObject *
Pympz_fmod(PyObject *self, PyObject *other)
{
    PyObject *r;
    PympzObject *temp;

    if(Pyxmpz_Check(self)) {
        r = self;
    } else {
        r = (PyObject*)Pympz_new();
        if(!r) {
            Py_XDECREF((PyObject*)r);
            return NULL;
        }
    }

    if((Pympz_Check(other) || Pyxmpz_Check(other))) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            return NULL;
        }
        mpz_fdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    } else {
        temp = Pympz_From_Integer(other);
        if(!temp) {
            TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)temp);
            Py_XDECREF((PyObject*)r);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(temp)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            return NULL;
        }
        mpz_fdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), temp->z);
        Py_DECREF(temp);
    }
    if(Pyxmpz_Check(self))
        Py_RETURN_NONE;
    else
        return r;
}

static char doc_tdivmodg[]="\
tdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards zero and the remainder will have the same sign\n\
as x. x and y must be integers. If x is an 'xmpz', the results will be\n\
'xmpz' and neither x or y will be mutated. Otherwise, the results will be\n\
'mpz'.\n\
";
static PyObject *
Pygmpy_tdivmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *q, *r, *result;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tdivmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        q = (PyObject*)Pyxmpz_new();
        r = (PyObject*)Pyxmpz_new();
    } else {
        q = (PyObject*)Pympz_new();
        r = (PyObject*)Pympz_new();
    }
    result = PyTuple_New(2);
    if(!q || !r || !result) {
        Py_XDECREF((PyObject*)result);
        Py_XDECREF((PyObject*)q);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tdivmod() division by 0");
            return NULL;
        }
        mpz_tdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("tdivmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tdivmod() division by 0");
            return NULL;
        }
        mpz_tdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_tdivg[]="\
tdiv(x,y): returns the quotient of x divided by y. The quotient is rounded\n\
towards 0. x and y must be integers. If x is an 'xmpz', the result will be an\n\
'xmpz' and x will NOT be mutated. Otherwise, the result will be an 'mpz'.\n\
";

static PyObject *
Pygmpy_tdiv(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *q;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tdiv() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        q = (PyObject*)Pyxmpz_new();
    } else {
        q = (PyObject*)Pympz_new();
    }
    if(!q) {
        Py_XDECREF((PyObject*)q);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            return NULL;
        }
        mpz_tdiv_q(Pympz_AS_MPZ(q), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("tdiv() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            return NULL;
        }
        mpz_tdiv_q(Pympz_AS_MPZ(q), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    return q;
}

static char doc_tmodg[]="\
tmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as x. x and y must be integers. If x is an 'xmpz',\n\
the result will be an 'xmpz' and x will NOT be mutated. Otherwise, the\n\
result will be 'mpz'.\n\
";

static PyObject *
Pygmpy_tmod(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *r;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    if(Pyxmpz_Check(x)) {
        r = (PyObject*)Pyxmpz_new();
    } else {
        r = (PyObject*)Pympz_new();
    }
    if(!r) {
        Py_XDECREF((PyObject*)r);
        return NULL;
    }

    if((Pympz_Check(x) || Pyxmpz_Check(x)) && ((Pympz_Check(y) || Pyxmpz_Check(y)))) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            return NULL;
        }
        mpz_tdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            return NULL;
        }
        mpz_tdiv_r(Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF(tempx);
        Py_DECREF(tempy);
    }
    return r;
}

static char doc_tmodm[]="\
x.tmod(y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as x. x and y must be integers. If x is an 'xmpz', the\n\
result will replace x and None will be returned. Otherwise, the result will\n\
be an 'mpz'.\n\
";

static PyObject *
Pympz_tmod(PyObject *self, PyObject *other)
{
    PyObject *r;
    PympzObject *temp;

    if(Pyxmpz_Check(self)) {
        r = self;
    } else {
        r = (PyObject*)Pympz_new();
        if(!r) {
            Py_XDECREF((PyObject*)r);
            return NULL;
        }
    }

    if((Pympz_Check(other) || Pyxmpz_Check(other))) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            return NULL;
        }
        mpz_tdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    } else {
        temp = Pympz_From_Integer(other);
        if(!temp) {
            TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)temp);
            Py_XDECREF((PyObject*)r);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(temp)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            return NULL;
        }
        mpz_tdiv_r(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), temp->z);
        Py_DECREF(temp);
    }
    if(Pyxmpz_Check(self))
        Py_RETURN_NONE;
    else
        return r;
}

