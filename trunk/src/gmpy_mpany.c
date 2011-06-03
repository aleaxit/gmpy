/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpany.c                                                            *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008 - 2011 Case Van Horsen                              *
 *                                                                         *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2.1 of the License, or  *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This library is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * Lesser General Public License for more details.                         *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with this library; if not, write to the Free Software     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA           *
 * 02110-1301  USA                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


/* Generic module-level methods for gmpy types.
 *
 * These methods are designed to accept any number type as input and call
 * the appropriate type-specific method. For example, gmpy2.digits(n) will
 * call gmpy2.mpz(n).digits() if n is an integer, gmpy2.mpq(n).digits() if
 * n is a rational, or gmpy2.mpf(n).digits() is n is a float.
 */

/* gmpy_square is only intended to be used at the module level!
 * gmpy_square uses the METH_O/METH_NOARGS calling convention!
 * gmpy_square assumes mpX_square also use the METH_O/METH_NOARGS convention!
 */

PyDoc_STRVAR(doc_mpany_square,
"square(x) -> number\n\n"
"Return x * x. If x is an integer, then the result is an 'mpz'.\n"
"If x is a rational, then the result is an 'mpq'. If x is a float,\n"
"then the result is an 'mpf'.");

static PyObject *
Pympany_square(PyObject *self, PyObject *other)
{
    if (isInteger(other)) {
        return Pympz_square(self, other);
    }
    else if (isRational(other)) {
        return Pympq_square(self, other);
    }
#ifdef WITHMPFR
    else if (isReal(other)) {
        return Pympfr_sqr(self, other);
    }
#endif

    TYPE_ERROR("square() not supported");
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
#ifdef WITHMPFR
    else if (isReal(temp))
        return Pympfr_digits(self, args);
#endif

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
#ifdef WITHMPFR
    else if (isReal(other))
        return Pympfr_sign(self, other);
#endif

    TYPE_ERROR("sign() not supported");
    return NULL;
}

/* create a copy of a gmpy2 object */
PyDoc_STRVAR(doc_copyg,
"_copy(x): -> gmpy2_object\n\n"
"Return a copy of x. Raises TypeError if x is not a gmpy2 object.");
static PyObject *
Pympany_copy(PyObject *self, PyObject *other)
{
    if (self && Pympz_Check(self))
        return (PyObject*)Pympz2Pympz(self);
    else if (self && Pyxmpz_Check(self))
        return (PyObject*)Pyxmpz2Pyxmpz(self);
    else if (self && Pympq_Check(self))
        return (PyObject*)Pympq2Pympq(self);
#ifdef WITHMPFR
    else if (self && Pympfr_Check(self))
        return (PyObject*)Pympfr2Pympfr(self, 0);
#endif
    else if (Pympz_Check(other))
        return (PyObject*)Pympz2Pympz(other);
    else if (Pyxmpz_Check(other))
        return (PyObject*)Pyxmpz2Pyxmpz(other);
    else if (Pympq_Check(other))
        return (PyObject*)Pympq2Pympq(other);
#ifdef WITHMPFR
    else if (Pympfr_Check(other))
        return (PyObject*)Pympfr2Pympfr(other, 0);
#endif
    TYPE_ERROR("_copy() requires a gmpy2 object as argument");
    return NULL;
}

PyDoc_STRVAR(doc_binarym,
"x.binary() -> binary string\n\n"
"Return a Python string (or bytes for Python 3+) that is a portable\n"
"binary representation of a gmpy2 object x. The binary string can\n"
"later be passed to the appropriate constructor function to obtain\n"
"an exact copy of x's value.");
PyDoc_STRVAR(doc_binaryg,
"binary(x) -> binary string\n\n"
"Return a Python string (or bytes for Python 3+) that is a portable\n"
"binary representation of a gmpy2 object x. The binary string can\n"
"later be passed to the appropriate constructor function to obtain\n"
"an exact copy of x's value. Raises TypeError if x is not a gmpy2\n"
"object.");

static PyObject *
Pympany_binary(PyObject *self, PyObject *other)
{
    if (self && Pympz_Check(self))
        return Pympz2binary((PympzObject*)self);
    else if(self && Pyxmpz_Check(self))
        return Pyxmpz2binary((PyxmpzObject*)self);
    else if(self && Pympq_Check(self))
        return Pympq2binary((PympqObject*)self);
#ifdef WITHMPFR
    else if(self && Pympfr_Check(self))
        return Pympfr2binary((PympfrObject*)self);
#endif
    else if(Pympz_Check(other))
        return Pympz2binary((PympzObject*)other);
    else if(Pyxmpz_Check(other))
        return Pyxmpz2binary((PyxmpzObject*)other);
    else if(Pympq_Check(other))
        return Pympq2binary((PympqObject*)other);
#ifdef WITHMPFR
    else if(Pympfr_Check(other))
        return Pympfr2binary((PympfrObject*)other);
#endif
    TYPE_ERROR("binary() requires a gmpy2 object as argument");
    return NULL;
}

static PyObject *
Pympany_pow(PyObject *base, PyObject *exp, PyObject *mod)
{
#ifndef WITHMPFR
    PyObject *result = 0, *temp;
#endif

    if (isInteger(base) && isInteger(exp))
        return Pympz_pow(base, exp, mod);
    else if (isRational(base) && isRational(exp))
        return Pympq_pow(base, exp, mod);
#ifdef WITHMPFR
    else if (isReal(base) && isReal(exp));
        return Pympfr2_pow(base, exp, mod);
#else
    /* Support mpz**float and float**mpz. */
    if (CHECK_MPZANY(base) && PyFloat_Check(exp)) {
        temp = PyFloat_FromDouble(mpz_get_d(Pympz_AS_MPZ(base)));
        if (temp) {
            result = PyNumber_Power(temp, exp, mod);
            Py_DECREF(temp);
        }
        return result;
    }
    if (CHECK_MPZANY(exp) && PyFloat_Check(base)) {
        temp = PyFloat_FromDouble(mpz_get_d(Pympz_AS_MPZ(exp)));
        if (temp) {
            result = PyNumber_Power(base, temp, mod);
            Py_DECREF(temp);
        }
        return result;
    }
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

PyDoc_STRVAR(doc_printf,
"printf(fmt, x) -> string\n\n"
"Return a Python string by formatting 'x' using the format string\n"
"'fmt'. Note: invalid format strings will cause a crash. Please\n"
"see the GMP/MPFR/MPC manuals for details on the format code.");

static PyObject *
Pympany_printf(PyObject *self, PyObject *args)
{
    PyObject *result = 0, *x = 0;
    char *buffer = 0, *fmtcode = 0;
    int buflen;
    void *generic;

    if (!PyArg_ParseTuple(args, "sO", &fmtcode, &x))
        return NULL;

    if (CHECK_MPZANY(x) || Pympq_Check(x)) {
        if (CHECK_MPZANY(x))
            generic = Pympz_AS_MPZ(x);
        else
            generic = Pympq_AS_MPQ(x);
        buflen = gmp_asprintf(&buffer, fmtcode, generic);
        result = Py_BuildValue("s", buffer);
        PyMem_Free(buffer);
        return result;
    }
#ifdef WITHMPFR
    else if(Pympfr_Check(x)) {
        generic = Pympfr_AS_MPFR(x);
        buflen = mpfr_asprintf(&buffer, fmtcode, generic);
        result = Py_BuildValue("s", buffer);
        PyMem_Free(buffer);
        return result;
    }
#endif
    else {
        TYPE_ERROR("printf() requires a gmpy2 object as argument");
        return NULL;
    }
}

