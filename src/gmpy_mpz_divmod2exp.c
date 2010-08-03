/* gmpy_mpz_divmod2exp.c
 *
 * This file should be considered part of gmpy.c
 *
 * This file contains functions related to division and remainder by a power
 * of two. Functions are optimized by writing distinct functions for
 * gmpy2.function versus mpz.function.
 */

/*
 **************************************************************************
 * Ceiling division and remainder by power of two.
 **************************************************************************
 */

static char doc_cdivmod2expg[]="\
cdivmod2exp(x,n): returns the quotient and remainder of x divided by 2**n.\n\
The quotient is rounded towards +Inf and the remainder will be negative.\n\
x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_cdivmod2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *q, *r, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cdivmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("cdivmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("cdivmod2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_TWO_MPZANY_TUPLE(x, q, r, result);

    if(CHECK_MPZANY(x)) {
        mpz_cdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(x), nbits);
        mpz_cdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("cdivmod2exp() requires 'mpz','int' arguments");
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_q_2exp(Pympz_AS_MPZ(q), tempx->z, nbits);
        mpz_cdiv_r_2exp(Pympz_AS_MPZ(r), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_cdiv2expg[]="\
cdiv2exp(x,n): returns the quotient of x divided by 2**n. The quotient is\n\
rounded towards +inf. x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_cdiv2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cdiv2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("cdiv2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("cdiv2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x)) {
        mpz_cdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("cdiv2exp() requires 'mpz','int' arguments");
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_q_2exp(Pympz_AS_MPZ(result), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    return result;
}

static char doc_cmod2expg[]="\
cmod2exp(x,n): returns the remainder of x divided by 2**n. The remainder\n\
will be negative. x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_cmod2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("cmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("cmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("cmod2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x)) {
        mpz_cdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("cmod2exp() requires expects 'mpz','int' arguments");
            Py_DECREF(result);
            return NULL;
        }
        mpz_cdiv_r_2exp(Pympz_AS_MPZ(result), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    return result;
}

static char doc_cdiv2expm[]="\
x.cdiv2exp(n): returns the quotient of x divided by 2**n. The quotient\n\
is rounded towards +Inf. x must be an integer. n must be > 0. Will\n\
mutate x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_cdiv2exp(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("cdiv2exp() requires 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("cdiv2exp() requires n > 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_cdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    } else {
        if(!(result = (PyObject*)Pympz_new())) {
            return NULL;
        }
        mpz_cdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    }
}

static char doc_cmod2expm[]="\
x.cmod2exp(n): returns the remainder of x divided by 2**n. The remainder\n\
will be negative. x must be an integer. n must be > 0. Will mutate x if it\n\
is an 'xmpz'.\n\
";

static PyObject *
Pympz_cmod2exp(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("cmod2exp() requires 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("cmod2exp() requires n > 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_cdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    } else {
        if(!(result = (PyObject*)Pympz_new())) {
            return NULL;
        }
        mpz_cdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    }
}

/*
 **************************************************************************
 * Floor division and remainder by power of two.
 **************************************************************************
 */

static char doc_fdivmod2expg[]="\
fdivmod2exp(x,n): returns quotient and remainder after dividing x by 2**n.\n\
The quotient is rounded towards -Inf and the remainder will be positive.\n\
x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_fdivmod2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *q, *r, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fdivmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("fdivmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("fdivmod2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_TWO_MPZANY_TUPLE(x, q, r, result);

    if(CHECK_MPZANY(x)) {
        mpz_fdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(x), nbits);
        mpz_fdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("fdivmod2exp() requires 'mpz','int' arguments");
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_q_2exp(Pympz_AS_MPZ(q), tempx->z, nbits);
        mpz_fdiv_r_2exp(Pympz_AS_MPZ(r), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_fdiv2expg[]="\
fdiv2exp(x,n): returns the quotient of x divided by 2**n. The quotient\n\
is rounded towards -Inf. x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_fdiv2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fdiv2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("fdiv2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("fdiv2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x)) {
        mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("fdiv2exp() requires 'mpz','int' arguments");
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    return result;
}

static char doc_fmod2expg[]="\
fmod2exp(x,n): returns remainder of x divided by 2**n. The remainder will\n\
be positive. x must be an integer. n must be greater than 0.\n\
";

static PyObject *
Pygmpy_fmod2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("fmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("fmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("fmod2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x)) {
        mpz_fdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("fmod2exp() requires 'mpz','int' arguments");
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(tempx), nbits);
        Py_DECREF((PyObject*)tempx);
    }
    return result;
}

static char doc_fdiv2expm[]="\
x.fdiv2exp(n): returns the quotient of x divided by 2**n. The quotient\n\
is rounded towards -Inf. x must be an integer. n must be > 0. Will mutate\n\
x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_fdiv2exp(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("fdiv2exp() requires 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("fdiv2exp() requires n > 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_fdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    } else {
        if(!(result = (PyObject*)Pympz_new())) {
            return NULL;
        }
        mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    }
}

static char doc_fmod2expm[]="\
x.fmod2exp(n): returns the remainder of x divided by 2**n. The remainder\n\
will be positive. x must be an integer. n must be > 0. Will mutate x if\n\
it is an 'xmpz'.\n\
";

static PyObject *
Pympz_fmod2exp(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("fmod2exp() requires 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("fmod2exp() requires n > 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_fdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    } else {
        if(!(result = (PyObject*)Pympz_new())) {
            return NULL;
        }
        mpz_fdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    }
}

/*
 **************************************************************************
 * Truncating division and remainder by power of two.
 **************************************************************************
 */

static char doc_tdivmod2expg[]="\
tdivmod2exp(x,n): returns the quotient and remainder of x divided by 2**n.\n\
The quotient is rounded towards zero and the remainder will have the same\n\
sign as x. x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_tdivmod2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *q, *r, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tdivmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("tdivmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("tdivmod2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_TWO_MPZANY_TUPLE(x, q, r, result);

    if(CHECK_MPZANY(x)) {
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(x), nbits);
        mpz_tdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("fdivmod2exp() requires 'mpz','int' arguments");
            Py_DECREF(q);
            Py_DECREF(r);
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(q), tempx->z, nbits);
        mpz_tdiv_r_2exp(Pympz_AS_MPZ(r), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static char doc_tdiv2expg[]="\
tdiv2exp(x,n): returns the quotient of x divided by 2**n. Uses 'truncate'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pygmpy_tdiv2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tdiv2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("tdiv2exp() requires 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("tdiv2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x)) {
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("tdiv2exp() requires 'mpz','int' arguments");
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(result), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    return result;
}

static char doc_tmod2expg[]="\
tmod2exp(x,n): returns the remainder of x divided by 2**n. The remainder\n\
will have the same sign as x. x must be an integer. n must be > 0.\n\
";

static PyObject *
Pygmpy_tmod2exp(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("tmod2exp() requires 'mpz','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("tmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("tmod2exp() requires n > 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x)) {
        mpz_tdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), nbits);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("tmod2exp() requires 'mpz','int' arguments");
            Py_DECREF(result);
            return NULL;
        }
        mpz_tdiv_r_2exp(Pympz_AS_MPZ(result), tempx->z, nbits);
        Py_DECREF((PyObject*)tempx);
    }
    return result;
}

static char doc_tdiv2expm[]="\
x.tdiv2exp(n): returns the quotient of x divided by 2**n. The quotient\n\
is rounded towards 0. x must be an integer. n must be > 0. Will mutate\n\
x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_tdiv2exp(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("tdiv2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("tdiv2exp() requires n > 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_tdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    } else {
        if(!(result = (PyObject*)Pympz_new())) {
            return NULL;
        }
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    }
}

static char doc_tmod2expm[]="\
x.tmod2exp(n): returns the remainder of x divided by 2**n. The remainder\n\
will have the same sign as x. x must be an integer. n must be > 0. Will\n\
mutate x if it is an 'xmpz'.\n\
";

static PyObject *
Pympz_tmod2exp(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("tmod2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("tmod2exp() requires n > 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_tdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    } else {
        if(!(result = (PyObject*)Pympz_new())) {
            return NULL;
        }
        mpz_tdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    }
}
