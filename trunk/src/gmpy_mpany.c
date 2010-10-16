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
 * gmpy_square uses the METH_O calling convention!
 * gmpy_square assumes mpX_square also use the METH_O convention!
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
