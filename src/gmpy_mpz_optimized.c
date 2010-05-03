/* gmpy_mpz_optimized.c
 *
 * Optimized functions that operate strictly on mpz or xmpz.
 *
 * This file should be considered part of gmpy.c
 *
 * This file contains functions that are optimized for performance. Functions
 * are optimized by writing distinct functions for gmpy2.function versus
 * mpz.function.
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

static char doc_fmod2expm[]="\
x.fmod2exp(n): returns remainder after dividing x by 2**n. Uses 'floor'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

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

static char doc_fdiv2expm[]="\
x.fdiv2exp(n): returns quotient after dividing x by 2**n. Uses 'floor'\n\
rounding. Will mutate x if it is an 'xmpz'. n must be > 0.\n\
";

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

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

static char doc_fdivmod2expm[]="\
x.fdivmod2exp(n): returns quotient and remainder after dividing x by 2**n.\n\
Uses 'floor' rounding. Both quotient and remainder are new objects. The result\n\
types will match the type of x. n must be > 0.\n\
";

/* This function can only be used as a bound method of mpz or xmpz. The
 * function assumes 'self' contains an instance of 'mpz' or 'xmpz'. 'other'
 * containts the numbers of bits.
 *
 * This function will mutate an xmpz.
 */

static PyObject *
Pympz_fdivmod2expm(PyObject *self, PyObject *other)
{
    long nbits;
    PyObject *q, *r, *result;

    nbits = clong_From_Integer(other);
    if(nbits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz.fdivmod2exp() requires expects 'int' argument");
        return NULL;
    }
    if(nbits <= 0) {
        VALUE_ERROR("nbits must be > 0");
        return NULL;
    }
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

    mpz_fdiv_q_2exp(Pympz_AS_MPZ(q), Pympz_AS_MPZ(self), nbits);
    mpz_fdiv_r_2exp(Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), nbits);
    PyTuple_SET_ITEM(result, 0, (PyObject*)q);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

