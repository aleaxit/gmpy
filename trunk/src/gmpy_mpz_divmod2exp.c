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
cdivmod2exp(x,n): returns the quotient and remainder of x divided\n\
by 2**n. The quotient is rounded towards +Inf and the remainder\n\
will be negative. x must be an integer. n must be > 0.\n\
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
cdiv2exp(x,n): returns the quotient of x divided by 2**n. The\n\
quotient is rounded towards +inf. x must be an integer. n must\n\
be > 0.\n\
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
cmod2exp(x,n): returns the remainder of x divided by 2**n. The\n\
remainder will be negative. x must be an integer. n must be > 0.\n\
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
x.cdiv2exp(n): returns the quotient of x divided by 2**n. The\n\
quotient is rounded towards +Inf. x must be an integer. n must\n\
be > 0. If x is an 'xmpz', it will be modified in-place.\n\
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
x.cmod2exp(n): returns the remainder of x divided by 2**n. The\n\
remainder will be negative. x must be an integer. n must be > 0.\n\
If x is an 'xmpz', it will be modified in-place.\n\
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
fdivmod2exp(x,n): returns quotient and remainder after dividing x by\n\
2**n. The quotient is rounded towards -Inf and the remainder will be\n\
positive. x must be an integer. n must be > 0.\n\
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
fmod2exp(x,n): returns remainder of x divided by 2**n. The remainder\n\
will be positive. x must be an integer. n must be greater than 0.\n\
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
tdivmod2exp(x,n): returns the quotient and remainder of x divided by\n\
2**n. The quotient is rounded towards zero and the remainder will have\n\
the same sign as x. x must be an integer. n must be > 0.\n\
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
            TYPE_ERROR("tdivmod2exp() requires 'mpz','int' arguments");
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

/*
 **************************************************************************
 * pack and unpack methods
 *
 * Both pack and unpack use a devious trick when stuffing values into the
 * internal structure of an mpz_t. By setting a bit beyond the range of
 * interest, we are guaranteed that memory will remain available. When
 * the bit is cleared, it also triggers normalization of the value by
 * accounting for leading bits that are zero.
 **************************************************************************
 */

static char doc_packg[]="\
pack(l,n): packs a list of integers 'l' into a single 'mpz' (or 'xmpz')\n\
by concatenating each integer after padding to length n bits. Raises an\n\
error if any integer is negative or greater than n bits in length.\n\
";

static PyObject *
Pygmpy_pack(PyObject *self, PyObject *args)
{
    Py_ssize_t nbits, total_bits, index, lst_count, i, temp_bits, limb_count, extra_bits;
    PyObject *lst, *result;
    mpz_t temp;
    PympzObject *tempx = 0;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("pack() requires 'list','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("pack() requires 'list','int' arguments");
        return NULL;
    }

    if(nbits <= 0) {
        VALUE_ERROR("pack() requires n > 0");
        return NULL;
    }

    if(!PyList_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("pack() requires 'list','int' arguments");
        return NULL;
    }

    lst = PyTuple_GET_ITEM(args, 0);
    lst_count = PyList_GET_SIZE(lst);
    total_bits = nbits * lst_count;

    CREATE0_ONE_MPZANY(result);
    mpz_setbit(Pympz_AS_MPZ(result), total_bits+1);

    mpz_inoc(temp);
    mpz_set_ui(temp, 0);
    limb_count = 0;
    extra_bits = 0;

    for(index = 0; index < lst_count; index++) {
        if(!(tempx = Pympz_From_Integer(PyList_GetItem(lst, index)))
            || (mpz_sgn(Pympz_AS_MPZ(tempx)) < 0)
            || (mpz_sizeinbase(Pympz_AS_MPZ(tempx),2) > nbits)) {
            TYPE_ERROR("pack() requires list elements be positive integers < 2^nbits");
            mpz_cloc(temp);
            Py_XDECREF((PyObject*)tempx);
            Py_DECREF(result);
            return NULL;
        }
        mpz_mul_2exp(Pympz_AS_MPZ(tempx), Pympz_AS_MPZ(tempx), extra_bits);
        mpz_add(temp, temp, Pympz_AS_MPZ(tempx));
        extra_bits += nbits;
        i = 0;
        temp_bits = mpz_sizeinbase(temp, 2) * mpz_sgn(temp);
        while(extra_bits >= mp_bits_per_limb) {
            if(temp_bits > 0)
                Pympz_AS_MPZ(result)->_mp_d[limb_count] = mpz_getlimbn(temp, i);
            i += 1;
            extra_bits -= mp_bits_per_limb;
            limb_count += 1;
            temp_bits -= mp_bits_per_limb;
        }
        if(temp_bits > 0) {
            mpz_tdiv_q_2exp(temp, temp, mp_bits_per_limb * i);
        } else {
            mpz_set_ui(temp, 0);
        }
        Py_DECREF((PyObject*)tempx);
    }
    Pympz_AS_MPZ(result)->_mp_d[limb_count] = mpz_getlimbn(temp, 0);
    mpz_clrbit(Pympz_AS_MPZ(result), total_bits+1);
    mpz_cloc(temp);
    return result;
}

static char doc_unpackg[]="\
unpack(x,n): unpacks an integer 'x' into a list of n-bit values. Raises\n\
error if 'x' is negative.\n\
";

static PyObject *
Pygmpy_unpack(PyObject *self, PyObject *args)
{
    Py_ssize_t nbits, total_bits, index = 0, lst_count, i, temp_bits = 0, extra_bits = 0;
    Py_ssize_t guard_bit, lst_ptr = 0;
    PyObject *item, *result;
    mpz_t temp;
    mp_limb_t extra = 0;
    PympzObject *tempx = 0;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("unpack() requires 'int','int' arguments");
        return NULL;
    }

    nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("unpack() requires 'int','int' arguments");
        return NULL;
    }

    if(nbits <= 0) {
        VALUE_ERROR("unpack() requires n > 0");
        return NULL;
    }

    if(!(tempx = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
        TYPE_ERROR("unpack() requires 'int','int' arguments");
        return NULL;
    }

    total_bits = mpz_sizeinbase(Pympz_AS_MPZ(tempx), 2) * mpz_sgn(Pympz_AS_MPZ(tempx));
    lst_count = total_bits / nbits;
    if((total_bits % nbits) || !lst_count)
        lst_count += 1;

    if(!(result = PyList_New(lst_count))) {
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    if(mpz_sgn(Pympz_AS_MPZ(tempx))==0) {
        if(options.prefer_mutable) {
            item = (PyObject*)Pyxmpz_new();
        } else {
            item = (PyObject*)Pympz_new();
        }
        if(!item) {
            Py_DECREF((PyObject*)tempx);
            Py_DECREF(result);
            return NULL;
        }
        mpz_set_ui(Pympz_AS_MPZ(item), 0);
        PyList_SET_ITEM(result, 0, item);
        Py_DECREF((PyObject*)tempx);
        return result;
    }

    mpz_inoc(temp);
    guard_bit = nbits + (2 * mp_bits_per_limb);

    while(lst_ptr < lst_count) {
        i = 0;
        temp_bits = 0;
        mpz_set_ui(temp, 0);
        mpz_setbit(temp, guard_bit);
        while(temp_bits + extra_bits < nbits) {
            temp->_mp_d[i++] = mpz_getlimbn(Pympz_AS_MPZ(tempx), index++);
            temp_bits += mp_bits_per_limb;
        }
        mpz_clrbit(temp, guard_bit);
        mpz_mul_2exp(temp, temp, extra_bits);
        mpz_add_ui(temp, temp, extra);
        temp_bits += extra_bits;

        while((lst_ptr <lst_count) && (temp_bits >= nbits)) {
            if(options.prefer_mutable) {
                item = (PyObject*)Pyxmpz_new();
            } else {
                item = (PyObject*)Pympz_new();
            }
            if(!item) {
                mpz_cloc(temp);
                Py_DECREF((PyObject*)tempx);
                Py_DECREF(result);
                return NULL;
            }
            mpz_tdiv_r_2exp(Pympz_AS_MPZ(item), temp, nbits);
            PyList_SET_ITEM(result, lst_ptr++, item);
            mpz_tdiv_q_2exp(temp, temp, nbits);
            temp_bits -= nbits;
        }
        extra = mpz_get_ui(temp);
        extra_bits = temp_bits;
    }
    Py_DECREF((PyObject*)tempx);
    mpz_cloc(temp);
    return result;
}


