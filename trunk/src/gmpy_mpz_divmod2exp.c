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
 * Floor division and remainder by power of two.
 **************************************************************************
 */

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_fmod2expg[]="\
fmod2exp(x,n): returns remainder after dividing x by 2**n. Uses 'floor'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_fmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.fmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.fmod2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);
        if(Pympz_Check(self)) {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_fdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            return result;
        } else if(Pyxmpz_Check(self)) {
            mpz_fdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
            Py_RETURN_NONE;
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.fmod2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_fdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
            return result;
        }
    }
}

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static char doc_fmod2expm[]="\
x.fmod2exp(n): returns remainder after dividing x by 2**n. Uses 'floor'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_fmod2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.fmod2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
    if(Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new())) {
            Py_DECREF(self);
            return NULL;
        }
        mpz_fdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    } else {
        mpz_fdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    }
}


/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_fdiv2expg[]="\
fdiv2exp(x,n): returns quotient after dividing x by 2**n. Uses 'floor'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_fdiv2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.fdiv2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.fdiv2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);
        if(Pympz_Check(self)) {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            return result;
        } else if(Pyxmpz_Check(self)) {
            mpz_fdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
            Py_RETURN_NONE;
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.fdiv2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
            return result;
        }
    }
}

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static char doc_fdiv2expm[]="\
x.fdiv2exp(n): returns quotient after dividing x by 2**n. Uses 'floor'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_fdiv2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.fdiv2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
    if(Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new())) {
            Py_DECREF(self);
            return NULL;
        }
        mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    } else {
        mpz_fdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    }
}

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_fdivmod2expg[]="\
fdivmod2exp(x,n): returns quotient and remainder after dividing x by 2**n.\n\
Uses 'floor' rounding. Both quotient and remainder are new objects. The result\n\
types will match the type of x. n must be > 0.\n\
";

static PyObject *
Pympz_fdivmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *q, *r, *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.fdivmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.fdivmod2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);

        if(Pyxmpz_Check(self)) {
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

        if(Pympz_Check(self) || Pyxmpz_Check(self)) {
            mpz_fdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
            mpz_fdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.fdivmod2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            mpz_fdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
            mpz_fdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
        }
        PyTuple_SET_ITEM(result, 0, (PyObject*)q);
        PyTuple_SET_ITEM(result, 1, (PyObject*)r);
        return result;
    }
}

/*
 **************************************************************************
 * Ceiling division and remainder by power of two.
 **************************************************************************
 */

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_cmod2expg[]="\
cmod2exp(x,n): returns remainder after dividing x by 2**n. Uses 'ceiling'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_cmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.cmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.cmod2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);
        if(Pympz_Check(self)) {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_cdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            return result;
        } else if(Pyxmpz_Check(self)) {
            mpz_cdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
            Py_RETURN_NONE;
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.cmod2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_cdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
            return result;
        }
    }
}

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static char doc_cmod2expm[]="\
x.cmod2exp(n): returns remainder after dividing x by 2**n. Uses 'ceiling'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_cmod2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.cmod2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
    if(Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new())) {
            Py_DECREF(self);
            return NULL;
        }
        mpz_cdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    } else {
        mpz_cdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    }
}

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_cdiv2expg[]="\
cdiv2exp(x,n): returns quotient after dividing x by 2**n. Uses 'ceiling'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_cdiv2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.cdiv2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.cdiv2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);
        if(Pympz_Check(self)) {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_cdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            return result;
        } else if(Pyxmpz_Check(self)) {
            mpz_cdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
            Py_RETURN_NONE;
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.cdiv2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_cdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
            return result;
        }
    }
}

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static char doc_cdiv2expm[]="\
x.cdiv2exp(n): returns quotient after dividing x by 2**n. Uses 'ceiling'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_cdiv2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.cdiv2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
    if(Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new())) {
            Py_DECREF(self);
            return NULL;
        }
        mpz_cdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    } else {
        mpz_cdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    }
}

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_cdivmod2expg[]="\
cdivmod2exp(x,n): returns quotient and remainder after dividing x by 2**n.\n\
Uses 'ceiling' rounding. Both quotient and remainder are new objects. The result\n\
types will match the type of x. n must be > 0.\n\
";

static PyObject *
Pympz_cdivmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *q, *r, *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.cdivmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.cdivmod2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);

        if(Pyxmpz_Check(self)) {
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

        if(Pympz_Check(self) || Pyxmpz_Check(self)) {
            mpz_cdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
            mpz_cdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.cdivmod2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            mpz_cdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
            mpz_cdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
        }
        PyTuple_SET_ITEM(result, 0, (PyObject*)q);
        PyTuple_SET_ITEM(result, 1, (PyObject*)r);
        return result;
    }
}

/*
 **************************************************************************
 *Truncating division and remainder by power of two.
 **************************************************************************
 */

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_tmod2expg[]="\
tmod2exp(x,n): returns remainder after dividing x by 2**n. Uses 'truncate'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_tmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.tmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.tmod2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);
        if(Pympz_Check(self)) {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_tdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            return result;
        } else if(Pyxmpz_Check(self)) {
            mpz_tdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
            Py_RETURN_NONE;
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.tmod2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_tdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
            return result;
        }
    }
}

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static char doc_tmod2expm[]="\
x.tmod2exp(n): returns remainder after dividing x by 2**n. Uses 'truncate'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_tmod2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.tmod2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
    if(Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new())) {
            Py_DECREF(self);
            return NULL;
        }
        mpz_tdiv_r_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    } else {
        mpz_tdiv_r_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    }
}

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_tdiv2expg[]="\
tdiv2exp(x,n): returns quotient after dividing x by 2**n. Uses 'truncate'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_tdiv2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.tdiv2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.tdiv2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);
        if(Pympz_Check(self)) {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_tdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            return result;
        } else if(Pyxmpz_Check(self)) {
            mpz_tdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
            Py_RETURN_NONE;
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.tdiv2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_tdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
            return result;
        }
    }
}

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static char doc_tdiv2expm[]="\
x.tdiv2exp(n): returns quotient after dividing x by 2**n. Uses 'truncate'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

static PyObject *
Pympz_tdiv2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.tdiv2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
    if(Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new())) {
            Py_DECREF(self);
            return NULL;
        }
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), nbits);
        return result;
    } else {
        mpz_tdiv_q_2exp(Pyxmpz_AS_MPZ(self), Pyxmpz_AS_MPZ(self), nbits);
        Py_RETURN_NONE;
    }
}

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_tdivmod2expg[]="\
tdivmod2exp(x,n): returns quotient and remainder after dividing x by 2**n.\n\
Uses 'truncate' rounding. Both quotient and remainder are new objects. The\n\
result types will match the type of x. n must be > 0.\n\
";

static PyObject *
Pympz_tdivmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *q, *r, *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.tdivmod2exp() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.tdivmod2exp() requires expects 'mpz','int' arguments");
            return NULL;
        }
        if(nbits <= 0) {
            VALUE_ERROR("nbits must be > 0");
            return NULL;
        }

        self = PyTuple_GET_ITEM(args, 0);

        if(Pyxmpz_Check(self)) {
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

        if(Pympz_Check(self) || Pyxmpz_Check(self)) {
            mpz_tdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
            mpz_tdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
        } else {
            if(!(self = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0)))) {
                TYPE_ERROR("gmpy2.tdivmod2exp() requires expects 'mpz','int' arguments");
                return NULL;
            }
            mpz_tdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
            mpz_tdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
            Py_DECREF(self);
        }
        PyTuple_SET_ITEM(result, 0, (PyObject*)q);
        PyTuple_SET_ITEM(result, 1, (PyObject*)r);
        return result;
    }
}
