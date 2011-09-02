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
"then the result is an 'mpf'. If x is a comples number, then the\n"
"result is an 'mpc'.");

static PyObject *
Pympany_square(PyObject *self, PyObject *other)
{
    if (isInteger(other))
        return Pympz_square(self, other);
    else if (isRational(other))
        return Pympq_square(self, other);
#ifdef WITHMPFR
    else if (isReal(other))
        return Pympfr_sqr(self, other);
#endif
#ifdef WITHMPC
    else if (isComplex(other))
        return Pympc_sqr(self, other);
#endif

    TYPE_ERROR("square() argument type not supported");
    return NULL;
}

/* gmpy_digits is only intended to be used at the module level!
 * gmpy_digits uses the METH_VARARGS calling convention!
 * gmpy_digits assumes mpX_digits also use the METH_VARARGS convention!
 */

PyDoc_STRVAR(doc_g_mpany_digits,
"digits(x,[base,[prec]]) -> string\n\n"
"Return string representing x. Calls mpz.digits, mpq.digits,\n"
"mpfr.digits, or mpc.digits as appropriate.");

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
#ifdef WITHMPC
    else if (isComplex(temp))
        return Pympc_digits(self, args);
#endif

    TYPE_ERROR("digits() argument type not supported");
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

    TYPE_ERROR("sign() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_add,
"add(x, y) -> number\n\n"
"Return x + y.");

static PyObject *
Pympany_add(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("add() requires 2 arguments.");
        return NULL;
    }

    if (isInteger(PyTuple_GET_ITEM(args, 0)) &&
        isInteger(PyTuple_GET_ITEM(args, 1)))
        return Pympz_add(self, args);

    if (isRational(PyTuple_GET_ITEM(args, 0)) &&
        isRational(PyTuple_GET_ITEM(args, 1)))
        return Pympq_add(self, args);

#ifdef WITHMPFR
    if (isReal(PyTuple_GET_ITEM(args, 0)) &&
        isReal(PyTuple_GET_ITEM(args, 1)))
        return Pympfr_add(self, args);

#endif
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)) &&
        isComplex(PyTuple_GET_ITEM(args, 1)))
        return Pympc_add(self, args);

#endif
    TYPE_ERROR("add() argument types not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_sub,
"sub(x, y) -> number\n\n"
"Return x - y.");

static PyObject *
Pympany_sub(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("sub() requires 2 arguments.");
        return NULL;
    }

    if (isInteger(PyTuple_GET_ITEM(args, 0)) &&
        isInteger(PyTuple_GET_ITEM(args, 1)))
        return Pympz_sub(self, args);

    if (isRational(PyTuple_GET_ITEM(args, 0)) &&
        isRational(PyTuple_GET_ITEM(args, 1)))
        return Pympq_sub(self, args);

#ifdef WITHMPFR
    if (isReal(PyTuple_GET_ITEM(args, 0)) &&
        isReal(PyTuple_GET_ITEM(args, 1)))
        return Pympfr_sub(self, args);

#endif
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)) &&
        isComplex(PyTuple_GET_ITEM(args, 1)))
        return Pympc_sub(self, args);

#endif
    TYPE_ERROR("sub() argument types not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_mul,
"mul(x, y) -> number\n\n"
"Return x * y.");

static PyObject *
Pympany_mul(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mul() requires 2 arguments.");
        return NULL;
    }

    if (isInteger(PyTuple_GET_ITEM(args, 0)) &&
        isInteger(PyTuple_GET_ITEM(args, 1)))
        return Pympz_mul(self, args);

    if (isRational(PyTuple_GET_ITEM(args, 0)) &&
        isRational(PyTuple_GET_ITEM(args, 1)))
        return Pympq_mul(self, args);

#ifdef WITHMPFR
    if (isReal(PyTuple_GET_ITEM(args, 0)) &&
        isReal(PyTuple_GET_ITEM(args, 1)))
        return Pympfr_mul(self, args);

#endif
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)) &&
        isComplex(PyTuple_GET_ITEM(args, 1)))
        return Pympc_mul(self, args);

#endif
    TYPE_ERROR("mul() argument types not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_div,
"div(x, y) -> number\n\n"
"Return x / y.");

static PyObject *
Pympany_div(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("div() requires 2 arguments.");
        return NULL;
    }

    if (isInteger(PyTuple_GET_ITEM(args, 0)) &&
        isInteger(PyTuple_GET_ITEM(args, 1)))
        return Pympz_div(self, args);

    if (isRational(PyTuple_GET_ITEM(args, 0)) &&
        isRational(PyTuple_GET_ITEM(args, 1)))
        return Pympq_div(self, args);

#ifdef WITHMPFR
    if (isReal(PyTuple_GET_ITEM(args, 0)) &&
        isReal(PyTuple_GET_ITEM(args, 1)))
        return Pympfr_div(self, args);

#endif
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)) &&
        isComplex(PyTuple_GET_ITEM(args, 1)))
        return Pympc_div(self, args);

#endif
    TYPE_ERROR("div() argument types not supported");
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
    TYPE_ERROR("binary() argument type not supported");
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
    else if (isReal(base) && isReal(exp))
        return Pympfr_pow(base, exp, mod);
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
#ifdef WITHMPC
    else if (isComplex(base) && isComplex(exp))
        return Pympc_pow(base, exp, mod);
#endif

    Py_RETURN_NOTIMPLEMENTED;
}

PyDoc_STRVAR(doc_printf,
"printf(fmt, x) -> string\n\n"
"Return a Python string by formatting 'x' using the format string\n"
"'fmt'. Note: invalid format strings will cause a crash. Please\n"
"see the GMP and MPFR manuals for details on the format code. 'mpc'\n"
"objects are not supported.");

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
        GMPY_FREE(buffer);
        return result;
    }
#ifdef WITHMPFR
    else if(Pympfr_Check(x)) {
        generic = Pympfr_AS_MPFR(x);
        buflen = mpfr_asprintf(&buffer, fmtcode, generic);
        if (buflen < 0) {
            VALUE_ERROR("printf() could not format the 'mpfr' object");
        }
        else {
            result = Py_BuildValue("s", buffer);
            mpfr_free_str(buffer);
        }
        return result;
    }
#endif
#ifdef WITHMPC
    else if(Pympc_Check(x)) {
        TYPE_ERROR("printf() does not support 'mpc'");
        return NULL;
    }
#endif
    else {
        TYPE_ERROR("printf() argument type not supported");
        return NULL;
    }
}

#ifdef WITHMPC
#define MPANY_MPFR_MPC(NAME) \
static PyObject * \
Pympany_##NAME(PyObject *self, PyObject *other) \
{ \
    if (isReal(other)) \
        return Pympfr_##NAME(self, other); \
    else if (isComplex(other)) \
        return Pympc_##NAME(self, other); \
    TYPE_ERROR(#NAME"() argument type not supported"); \
    return NULL; \
}
#else
#define MPANY_MPFR_MPC(NAME) \
static PyObject * \
Pympany_##NAME(PyObject *self, PyObject *other) \
{ \
    if (isReal(other)) \
        return Pympfr_##NAME(self, other); \
    TYPE_ERROR(#NAME"() argument type not supported"); \
    return NULL; \
}
#endif

#ifdef WITHMPFR
PyDoc_STRVAR(doc_mpany_is_nan,
"is_nan(x) -> boolean\n\n"
"Return True if x is 'nan' (Not-A-Number). Calls mpfr.is_nan or\n"
"mpc.is_nan as appropriate.");

static PyObject *
Pympany_is_nan(PyObject *self, PyObject *other)
{
    if (isReal(other))
        return Pympfr_is_nan(self, other);
#ifdef WITHMPC
    else if (isComplex(other))
        return Pympc_is_NAN(self, other);
#endif
    TYPE_ERROR("is_nan() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_is_inf,
"is_inf(x) -> boolean\n\n"
"Return True if x is '+Infinity' or '-Infinity'. Calls mpfr.is_inf\n"
"or mpc.is_inf as appropriate.");

static PyObject *
Pympany_is_inf(PyObject *self, PyObject *other)
{
    if (isReal(other))
        return Pympfr_is_inf(self, other);
#ifdef WITHMPC
    else if (isComplex(other))
        return Pympc_is_INF(self, other);
#endif
    TYPE_ERROR("is_inf() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_is_zero,
"is_zero(x) -> boolean\n\n"
"Return True if x is zero. Calls mpfr.is_zero or mpc.is_zero as\n"
"appropriate.");

static PyObject *
Pympany_is_zero(PyObject *self, PyObject *other)
{
    if (isReal(other))
        return Pympfr_is_zero(self, other);
#ifdef WITHMPC
    else if (isComplex(other))
        return Pympc_is_ZERO(self, other);
#endif
    TYPE_ERROR("is_zero() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_log,
"log(x) -> number\n\n"
"Return the natural logarithm of x.\n");

MPANY_MPFR_MPC(log)

PyDoc_STRVAR(doc_mpany_exp,
"exp(x) -> number\n\n"
"Return the exponential of x.\n");

MPANY_MPFR_MPC(exp)

PyDoc_STRVAR(doc_mpany_cos,
"cos(x) -> number\n\n"
"Return the cosine of x; x in radians.\n");

MPANY_MPFR_MPC(cos)

PyDoc_STRVAR(doc_mpany_sin,
"sin(x) -> number\n\n"
"Return the sine of x; x in radians.\n");

MPANY_MPFR_MPC(sin)

PyDoc_STRVAR(doc_mpany_tan,
"tan(x) -> number\n\n"
"Return the tangent of x; x in radians.\n");

MPANY_MPFR_MPC(tan)

PyDoc_STRVAR(doc_mpany_acos,
"acos(x) -> number\n\n"
"Return the arc-cosine of x; x in radians.\n");

MPANY_MPFR_MPC(acos)

PyDoc_STRVAR(doc_mpany_asin,
"asin(x) -> number\n\n"
"Return the arc-sine of x; x in radians.\n");

MPANY_MPFR_MPC(asin)

PyDoc_STRVAR(doc_mpany_atan,
"atan(x) -> number\n\n"
"Return the arc-tangent of x; x in radians.\n");

MPANY_MPFR_MPC(atan)

PyDoc_STRVAR(doc_mpany_cosh,
"cosh(x) -> number\n\n"
"Return the hyperbolic cosine of x.\n");

MPANY_MPFR_MPC(cosh)

PyDoc_STRVAR(doc_mpany_sinh,
"sinh(x) -> number\n\n"
"Return the hyberbolic sine of x.\n");

MPANY_MPFR_MPC(sinh)

PyDoc_STRVAR(doc_mpany_tanh,
"tanh(x) -> number\n\n"
"Return the hyperbolic tangent of x.\n");

MPANY_MPFR_MPC(tanh)

PyDoc_STRVAR(doc_mpany_acosh,
"acosh(x) -> number\n\n"
"Return the inverse hyperbolic cosine of x.\n");

MPANY_MPFR_MPC(acosh)

PyDoc_STRVAR(doc_mpany_asinh,
"asinh(x) -> number\n\n"
"Return the inverse hyperbolic sine of x.\n");

MPANY_MPFR_MPC(asinh)

PyDoc_STRVAR(doc_mpany_atanh,
"atanh(x) -> number\n\n"
"Return the inverse hyperbolic tangent of x.\n");

MPANY_MPFR_MPC(atanh)

PyDoc_STRVAR(doc_mpany_sqrt,
"sqrt(x) -> number\n\n"
"Return the square root of x. If x is integer, rational, or real,\n"
"then an 'mpf' will be returned. If x is complex, then an 'mpc' will\n"
"be returned. If context.allow_complex is True, negative values of x\n"
"will return an 'mpc'.\n");

MPANY_MPFR_MPC(sqrt)

PyDoc_STRVAR(doc_mpany_sin_cos,
"sin_cos(x) -> (number, number)\n\n"
"Return a tuple containing the sine and cosine of x; x in radians.\n");

MPANY_MPFR_MPC(sin_cos)

PyDoc_STRVAR(doc_mpany_fma,
"fma(x,y,z) -> number\n\n"
"Return correctly rounded result of (x * y) + z.");

static PyObject *
Pympany_fma(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fma() requires 3 arguments.");
        return NULL;
    }

    if (isReal(PyTuple_GET_ITEM(args, 0)) &&
        isReal(PyTuple_GET_ITEM(args, 1)) &&
        isReal(PyTuple_GET_ITEM(args, 2)))
        return Pympfr_fma(self, args);
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)) &&
        isComplex(PyTuple_GET_ITEM(args, 1)) &&
        isComplex(PyTuple_GET_ITEM(args, 2)))
        return Pympc_fma(self, args);
#endif
    TYPE_ERROR("fma() argument types not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_fms,
"fms(x,y,z) -> number\n\n"
"Return correctly rounded result of (x * y) - z.");

static PyObject *
Pympany_fms(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fms() requires 3 arguments.");
        return NULL;
    }

    if (isReal(PyTuple_GET_ITEM(args, 0)) &&
        isReal(PyTuple_GET_ITEM(args, 1)) &&
        isReal(PyTuple_GET_ITEM(args, 2)))
        return Pympfr_fms(self, args);
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)) &&
        isComplex(PyTuple_GET_ITEM(args, 1)) &&
        isComplex(PyTuple_GET_ITEM(args, 2)))
        return Pympc_fms(self, args);
#endif
    TYPE_ERROR("fms() argument types not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_div_2exp,
"div_2exp(x, n) -> number\n\n"
"Return 'mpfr' (or 'mpc') divided by 2**n.");

static PyObject *
Pympany_div_2exp(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("div_2exp() requires 2 arguments.");
        return NULL;
    }

    if (isReal(PyTuple_GET_ITEM(args, 0)))
        return Pympfr_div_2exp(self, args);
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)))
        return Pympc_div_2exp(self, args);
#endif
    TYPE_ERROR("div_2exp() argument types not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_mul_2exp,
"mul_2exp(x, n) -> number\n\n"
"Return 'mpfr' (or 'mpc') multiplied by 2**n.");

static PyObject *
Pympany_mul_2exp(PyObject *self, PyObject *args)
{
    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mul_2exp() requires 2 arguments.");
        return NULL;
    }

    if (isReal(PyTuple_GET_ITEM(args, 0)))
        return Pympfr_mul_2exp(self, args);
#ifdef WITHMPC
    if (isComplex(PyTuple_GET_ITEM(args, 0)))
        return Pympc_mul_2exp(self, args);
#endif
    TYPE_ERROR("mul_2exp() argument types not supported");
    return NULL;
}

#endif
