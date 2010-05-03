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

/* return N lowest bits from an mpz or xmpz*/

/* This function can only be used as an unbound method of gmpy2. The function
 * ignores 'self' and assumes all arguments are passed in the tuple 'args'.
 */

static char doc_fmod2expg[]="\
lowbits(x,n): returns the n lowest bits of x; n must be an\n\
ordinary Python int, >0; x must be an mpz, or else gets\n\
coerced to one.\n\
";

static PyObject *
Pympz_fmod2expg(PyObject *self, PyObject *args)
{
    long nbits;
    PyObject *result;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gmpy2.lowbits() requires expects 'mpz','int' arguments");
        return NULL;
    } else {
        nbits = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
        if(nbits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("gmpy2.lowbits() requires expects 'mpz','int' arguments");
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
                TYPE_ERROR("gmpy2.lowbits() requires expects 'mpz','int' arguments");
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
x.lowbits(n): returns the n lowest bits of x; n must be an\n\
ordinary Python int, >0.\n\
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
        TYPE_ERROR("mpz.lowbits() requires expects 'int' argument");
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

