/* gmpy_mpany.c
 *
 * Generic module-level methods for gmpy types.
 *
 * These methods are designed to accept any number type as input and call
 * the appropriate type-specific method. For example, gmpy2.digits(n) will
 * call gmpy2.mpz(n).digits() if n is an integer, gmpy2.mpq(n).digits() if
 * n is a rational, or gmpy2.mpf(n).digits() is n is a float.
 *
 * This file should be considered part of gmpy2.c
 */

/* gmpy_square is only intended to be used at the module level!
 * gmpy_square uses the METH_O/METH_NOARGS calling convention!
 * gmpy_square assumes mpX_square also use the METH_O/METH_NOARGS convention!
 */

PyDoc_STRVAR(doc_gmpy_square,
"square(x) -> number\n\n"
"Return x * x. If x is an integer, then the result is an 'mpz'.\n"
"If x is a rational, then the result is an 'mpq'. If x is a float,\n"
"then the result is an 'mpf'.");

static PyObject *
Pygmpy_square(PyObject *self, PyObject *other)
{
    if (isInteger(other)) {
        TYPE_ERROR("square() not supported for integers");
        return NULL;
    }
    else if (isRational(other)) {
        TYPE_ERROR("square() not supported for rationals");
        return NULL;
    }
    else if (isFloat(other)) {
        return Pympf_sqr(self, other);
    }

    TYPE_ERROR("square() not supported");
    return NULL;
}

/* gmpy_digits is only intended to be used at the module level!
 * gmpy_digits uses the METH_VARARGS calling convention!
 * gmpy_digits assumes mpX_digits also use the METH_VARARGS convention!
 */

PyDoc_STRVAR(doc_gmpy_digits,
"digits(x,[base]) -> string\n\n"
"Return string representing x. Calls mpz.digits, mpq.digits, or\n"
"mpf.digits as appropriate.");

static PyObject *
Pygmpy_digits(PyObject *self, PyObject *args)
{
    PyObject *temp;

    if (PyTuple_GET_SIZE(args) == 0) {
        TYPE_ERROR("digits() requires at least one argument");
        return NULL;
    }

    temp = PyTuple_GET_ITEM(args, 0);
    if (isInteger(temp))
        return Pympz_digits(self, args);
    else if (isRational(temp))
        return Pympq_digits(self, args);
    else if (isFloat(temp))
        return Pympf_digits(self, args);

    TYPE_ERROR("digits() not supported");
    return NULL;
}

/* gmpy_sign is only intended to be used at the module level!
 * gmpy_sign uses the METH_O/METH_NOARGS calling convention!
 * gmpy_sign assumes mpX_sign also use the METH_O/METH_NOARGS convention!
 */

PyDoc_STRVAR(doc_gmpy_sign,
"sign(x) -> number\n\n"
"Return -1 if x < 0, 0 if x == 0, or +1 if x >0.");

static PyObject *
Pygmpy_sign(PyObject *self, PyObject *other)
{
    if (isInteger(other))
        return Pympz_sign(self, other);
    else if (isRational(other))
        return Pympq_sign(self, other);
    else if (isFloat(other))
        return Pympf_sign(self, other);

    TYPE_ERROR("sign() not supported");
    return NULL;
}

