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

PyDoc_STRVAR(doc_g_mpany_square,
"square(x) -> number\n\n"
"Return x * x. If x is an integer, then the result is an 'mpz'.\n"
"If x is a rational, then the result is an 'mpq'. If x is a float,\n"
"then the result is an 'mpf'.");

static PyObject *
Pympany_square(PyObject *self, PyObject *other)
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
        return Pympfr_sqr(self, other);
    }

    TYPE_ERROR("square() not supported");
    return NULL;
}

/* gmpy_sqrt is only intended to be used at the module level!
 * gmpy_sqrt uses the METH_O/METH_NOARGS calling convention!
 * gmpy_sqrt assumes mpX_square also use the METH_O/METH_NOARGS convention!
 */

PyDoc_STRVAR(doc_g_mpany_sqrt,
"sqrt(x) -> number\n\n"
"Return square root of x. If x is an integer, then the result is the\n"
"integer portion of the square root. If x is a rational or a float,\n"
"then the result is an 'mpf'.");

static PyObject *
Pympany_sqrt(PyObject *self, PyObject *other)
{
    if (isInteger(other))
        return Pympz_sqrt(self, other);
    else if (isRational(other) || isFloat(other)) {
        return Pympfr_sqrt(self, other);
    }

    TYPE_ERROR("sqrt() not supported");
    return NULL;
}

/* gmpy_root is only intended to be used at the module level!
 * gmpy_root uses the METH_VARARGS calling convention!
 * gmpy_root assumes mpX_square also use the METH_VARARGS convention!
 */

PyDoc_STRVAR(doc_g_mpany_root,
"root(x,n) -> number\n\n"
"Return n-th root of x. If x is an integer, then the result is a\n"
"tuple containing the integer portion of the root and True if the\n"
"result is exact. If x is a rational or a float, then the result\n"
"is an 'mpf'.");

static PyObject *
Pympany_root(PyObject *self, PyObject *args)
{
    PyObject *temp;

    if (PyTuple_GET_SIZE(args) == 0) {
        TYPE_ERROR("root() requires at least one argument");
        return NULL;
    }

    temp = PyTuple_GET_ITEM(args, 0);
    if (isInteger(temp))
        return Pympz_root(self, args);
    else if (isRational(temp) || isFloat(temp)) {
        return Pympfr_root(self, args);
    }

    TYPE_ERROR("root() not supported");
    return NULL;
}

/* gmpy_digits is only intended to be used at the module level!
 * gmpy_digits uses the METH_VARARGS calling convention!
 * gmpy_digits assumes mpX_digits also use the METH_VARARGS convention!
 */

PyDoc_STRVAR(doc_g_mpany_digits,
"digits(x,[base,[prec]]) -> string\n\n"
"Return string representing x. Calls mpz.digits, mpq.digits, or\n"
"mpfr.digits as appropriate.");

static PyObject *
Pympany_digits(PyObject *self, PyObject *args)
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
        return Pympfr_digits(self, args);

    TYPE_ERROR("digits() not supported");
    return NULL;
}

/* gmpy_sign is only intended to be used at the module level!
 * gmpy_sign uses the METH_O/METH_NOARGS calling convention!
 * gmpy_sign assumes mpX_sign also use the METH_O/METH_NOARGS convention!
 */

PyDoc_STRVAR(doc_g_mpany_sign,
"sign(x) -> number\n\n"
"Return -1 if x < 0, 0 if x == 0, or +1 if x >0.");

static PyObject *
Pympany_sign(PyObject *self, PyObject *other)
{
    if (isInteger(other))
        return Pympz_sign(self, other);
    else if (isRational(other))
        return Pympq_sign(self, other);
    else if (isFloat(other))
        return Pympfr_sign(self, other);

    TYPE_ERROR("sign() not supported");
    return NULL;
}

