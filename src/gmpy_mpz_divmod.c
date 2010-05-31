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
cdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards +Inf and the remainder will have the opposite\n\
sign of y. x and y must be integers.\n\
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
    CREATE1_TWO_MPZANY_TUPLE(x, q, r, result);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cdivmod() division by 0");
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
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
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cdivmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_cdivg[]="\
cdiv(x,y): returns the quotient of x divided by y. The quotient is rounded\n\
towards +Inf. x and y must be integers.\n\
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
    CREATE1_ONE_MPZANY(x, q);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            Py_DECREF(q);
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
            Py_DECREF(q);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(q);
            return NULL;
        }
        mpz_cdiv_q(Pympz_AS_MPZ(q), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return q;
}

static char doc_cmodg[]="\
cmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the opposite sign of y. x and y must be integers.\n\
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
    CREATE1_ONE_MPZANY(x, r);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            Py_DECREF(r);
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
            Py_DECREF(r);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(r);
            return NULL;
        }
        mpz_cdiv_r(Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return r;
}

static char doc_cdivm[]="\
x.cdiv(y): returns the quotient of x divided by y. The quotient is rounded\n\
towards +Inf. x and y must be integers. Will mutate x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_cdiv(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_cdiv_q(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_cdiv_q(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("cdiv() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("cdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_cdiv_q(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_cdiv_q(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

static char doc_cmodm[]="\
x.cmod(y): returns the remainder of x divided by y. The remainder will\n\
have the opposite sign of y. x and y must be integers. Will mutate x if\n\
if is an 'xmpz'.\n\
";

static PyObject *
Pympz_cmod(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_cdiv_r(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_cdiv_r(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("cmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("cmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_cdiv_r(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_cdiv_r(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

/*
 **************************************************************************
 * Floor division and remainder.
 **************************************************************************
 */

static char doc_fdivmodg[]="\
fdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards -Inf and the remainder will have the same sign\n\
as y. x and y must be integers.\n\
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
    CREATE1_TWO_MPZANY_TUPLE(x, q, r, result);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fdivmod() division by 0");
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("fdivmod() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fdivmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_fdivg[]="\
fdiv(x,y): returns the quotient of x divided by y. The quotient is rounded\n\
towards -Inf. x and y must be integers.\n\
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
    CREATE1_ONE_MPZANY(x, q);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            Py_DECREF(q);
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
            Py_DECREF(q);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(q);
            return NULL;
        }
        mpz_fdiv_q(Pympz_AS_MPZ(q), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return q;
}

static char doc_fmodg[]="\
fmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as y. x and y must be integers.\n\
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
    CREATE1_ONE_MPZANY(x, r);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            Py_DECREF(r);
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
            Py_DECREF(r);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(r);
            return NULL;
        }
        mpz_fdiv_r(Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return r;
}
static char doc_fdivm[]="\
x.fdiv(y): returns the quotient of x divided by y. The quotient is rounded\n\
towards -Inf. x and y must be integers. Will mutate x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_fdiv(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_fdiv_q(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_fdiv_q(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("fdiv() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("fdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_fdiv_q(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_fdiv_q(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

static char doc_fmodm[]="\
x.fmod(y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as y. x and y must be integers. Will mutate x if it\n\
is an 'xmpz'.\n\
";

static PyObject *
Pympz_fmod(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_fdiv_r(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_fdiv_r(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("fmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("fmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_fdiv_r(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_fdiv_r(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

/*
 **************************************************************************
 * Truncating division and remainder.
 **************************************************************************
 */

static char doc_tdivmodg[]="\
tdivmod(x,y): returns the quotient and remainder of x divided by y. The\n\
quotient is rounded towards zero and the remainder will have the same sign\n\
as x. x and y must be integers.\n\
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
    CREATE1_TWO_MPZANY_TUPLE(x, q, r, result);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tdivmod() division by 0");
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
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
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tdivmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_qr(Pympz_AS_MPZ(q), Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_tdivg[]="\
tdiv(x,y): returns the quotient of x divided by y. The quotient is rounded\n\
towards 0. x and y must be integers.\n\
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
    CREATE1_ONE_MPZANY(x, q);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            Py_DECREF(q);
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
            Py_DECREF(q);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(q);
            return NULL;
        }
        mpz_tdiv_q(Pympz_AS_MPZ(q), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return q;
}

static char doc_tmodg[]="\
tmod(x,y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as x. x and y must be integers. Will mutate x if it\n\
is an 'xmpz'.\n\
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
    CREATE1_ONE_MPZANY(x, r);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            Py_DECREF(r);
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
            Py_DECREF(r);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(r);
            return NULL;
        }
        mpz_tdiv_r(Pympz_AS_MPZ(r), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return r;
}

static char doc_tdivm[]="\
x.tdiv(y): returns the quotient of x divided by y. The quotient is rounded\n\
towards 0. x and y must be integers. Will mutate x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_tdiv(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_tdiv_q(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_tdiv_q(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("tdiv() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("tdiv() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_tdiv_q(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_tdiv_q(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

static char doc_tmodm[]="\
x.tmod(y): returns the remainder of x divided by y. The remainder will\n\
have the same sign as x. x and y must be integers. Will mutate x if it\n\
is an 'xmpz'.\n\
";

static PyObject *
Pympz_tmod(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_tdiv_r(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_tdiv_r(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("tmod() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("tmod() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_tdiv_r(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_tdiv_r(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

