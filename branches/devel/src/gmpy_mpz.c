/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Case Van Horsen      *
 *                                                                         *
 * This file is part of GMPY2.                                             *
 *                                                                         *
 * GMPY2 is free software: you can redistribute it and/or modify it under  *
 * the terms of the GNU Lesser General Public License as published by the  *
 * Free Software Foundation, either version 3 of the License, or (at your  *
 * option) any later version.                                              *
 *                                                                         *
 * GMPY2 is distributed in the hope that it will be useful, but WITHOUT    *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or   *
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public    *
 * License for more details.                                               *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with GMPY2; if not, see <http://www.gnu.org/licenses/>    *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

PyDoc_STRVAR(doc_mpz,
"mpz() -> mpz(0)\n\n"
"     If no argument is given, return mpz(0).\n\n"
"mpz(n) -> mpz\n\n"
"     Return an 'mpz' object with a numeric value 'n' (truncating n\n"
"     to its integer part if it's a Fraction, 'mpq', Decimal, float\n"
"     or 'mpfr').\n\n"
"mpz(s[, base=0]):\n\n"
"     Return an 'mpz' object from a string 's' made of digits in the\n"
"     given base.  If base=0, binary, octal, or hex Python strings\n"
"     are recognized by leading 0b, 0o, or 0x characters, otherwise\n"
"     the string is assumed to be decimal. Values for base can range\n"
"     between 2 and 62.");

static PyObject *
Pygmpy_mpz(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPZ_Object *result = 0;
    PyObject *n = 0;
    int base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"n", "base", NULL };
    CTXT_Object *context = NULL;

    /* Optimize the most common use case */
    argc = PyTuple_Size(args);
    if (argc == 0) {
        if ((result = GMPy_MPZ_New(context))) {
            mpz_set_ui(result->z, 0);
        }
        return (PyObject*)result;
    }
    if (argc == 1) {
        n = PyTuple_GetItem(args, 0);
        if (IS_REAL(n) && !keywds) {
            result = GMPy_MPZ_From_Number(n, context);
            return (PyObject*)result;
        }
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist,
                                     &n, &base))
        return NULL;

    if ((base!=0) && ((base<2)||(base>62))) {
        VALUE_ERROR("base for mpz() must be 0 or in the "
                    "interval 2 ... 62");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        /* build-from-string (ascii or unicode) */
        result = GMPy_MPZ_From_PyStr(n, base, context);
    }
    else {
        if (argc==2 || (argc == 1 && keywds))
            TYPE_ERROR("mpz() with non-string argument needs exactly "
                       "1 argument");
        else {
            result = GMPy_MPZ_From_Number(n, context);
            if (!result)
                TYPE_ERROR("mpz() requires numeric or string argument");
        }
    }
    return (PyObject*)result;
}

/* Functions that operate strictly on mpz or xmpz. */

/* return number-of-digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(doc_num_digitsm,
"x.num_digits([base=10]) -> int\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

PyDoc_STRVAR(doc_num_digitsg,
"num_digits(x[, base=10]) -> int\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

static PyObject *
Pympz_num_digits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;
    CTXT_Object *context = NULL;

    PARSE_ONE_MPZ_OPT_CLONG(&base, "num_digits() requires 'mpz',['int'] arguments");

    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        Py_DECREF(self);
        return NULL;
    }
    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(MPZ(self), (int)base));
    Py_DECREF(self);
    return result;
}


PyDoc_STRVAR(doc_mpz_iroot,
"iroot(x,n) -> (number, boolean)\n\n"
"Return the integer n-th root of x and boolean value that is True\n"
"iff the root is exact. x >= 0. n > 0.");

static PyObject *
Pympz_iroot(PyObject *self, PyObject *args)
{
    mpir_si n;
    int exact;
    MPZ_Object *s = 0;
    PyObject *result = 0;
    CTXT_Object *context = NULL;

    PARSE_ONE_MPZ_REQ_SI(&n, "iroot() requires 'mpz','int' arguments");

    if (n <= 0) {
        VALUE_ERROR("n must be > 0");
        Py_DECREF(self);
        return NULL;
    }
    else if (n>1) {
        if (mpz_sgn(MPZ(self))<0) {
            VALUE_ERROR("iroot() of negative number");
            Py_DECREF(self);
            return NULL;
        }
    }
    if (!(s = GMPy_MPZ_New(context)) || !(result = PyTuple_New(2))) {
        Py_DECREF(self);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF(result);
        return NULL;
    }
    exact = mpz_root(s->z, MPZ(self), n);
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, (PyObject*)s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)PyBool_FromLong(exact));
    return result;
}

PyDoc_STRVAR(doc_mpz_iroot_rem,
"iroot_rem(x,n) -> (number, number)\n\n"
"Return a 2-element tuple (y,r), such that y is the integer n-th\n"
"root of x and x=y**n + r. x >= 0. n > 0.");

static PyObject *
Pympz_iroot_rem(PyObject *self, PyObject *args)
{
    mpir_si n;
    MPZ_Object *r = 0, *y = 0;
    PyObject *result = 0;
    CTXT_Object *context = NULL;

    PARSE_ONE_MPZ_REQ_SI(&n,
            "iroot_rem() requires 'mpz','int' arguments");

    if (n <= 0) {
        VALUE_ERROR("n must be > 0");
        Py_DECREF(self);
        return NULL;
    }
    else if (n>1) {
        if (mpz_sgn(MPZ(self))<0) {
            VALUE_ERROR("iroot_rem() of negative number");
            Py_DECREF(self);
            return NULL;
        }
    }
    y = GMPy_MPZ_New(context);
    r = GMPy_MPZ_New(context);
    result = PyTuple_New(2);
    if (!y || !r || !result) {
        Py_DECREF(self);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)y);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }
    mpz_rootrem(y->z, r->z, MPZ(self), n);
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, (PyObject*)y);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static PyObject *
Pympz_sign(PyObject *self, PyObject *other)
{
    long res;
    MPZ_Object* tempx;
    CTXT_Object *context = NULL;

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_sgn(MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_sgn(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("sign() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_sgn(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }
    return PyIntOrLong_FromLong(res);
}

PyDoc_STRVAR(doc_mpz_ceil, "Ceiling of an mpz returns itself.");

static PyObject *
Pympz_ceil(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(doc_mpz_floor, "Floor of an mpz returns itself.");

static PyObject *
Pympz_floor(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(doc_mpz_trunc, "Truncating an mpz returns itself.");

static PyObject *
Pympz_trunc(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(doc_mpz_round, "Round an mpz to power of 10.");

static PyObject *
Pympz_round(PyObject *self, PyObject *args)
{
    Py_ssize_t round_digits;
    MPZ_Object *result;
    mpz_t temp, rem;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) == 0) {
        Py_INCREF(self);
        return self;
    }

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("Too many arguments for __round__().");
        return NULL;
    }

    round_digits = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 0));
    if (round_digits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("__round__() requires 'int' argument");
        return NULL;
    }

    if (round_digits >= 0) {
        Py_INCREF(self);
        return self;
    }
    round_digits = -round_digits;

    if ((result = GMPy_MPZ_New(context))) {
        if (round_digits >= mpz_sizeinbase(MPZ(self), 10)) {
            mpz_set_ui(result->z, 0);
        }
        else {
            mpz_inoc(temp);
            mpz_inoc(rem);
            mpz_ui_pow_ui(temp, 10, round_digits);
            mpz_fdiv_qr(result->z, rem, MPZ(self), temp);
            mpz_mul_2exp(rem, rem, 1);
            if (mpz_cmp(rem, temp) > 0) {
                mpz_add_ui(result->z, result->z, 1);
            }
            else if (mpz_cmp(rem, temp) == 0) {
                if (mpz_odd_p(result->z)) {
                    mpz_add_ui(result->z, result->z, 1);
                }
            }
            mpz_mul(result->z, result->z, temp);
            mpz_cloc(rem);
            mpz_cloc(temp);
        }
    }

    return (PyObject*)result;
}

static int
Pympz_nonzero(MPZ_Object *self)
{
    return mpz_sgn(self->z) != 0;
}

#if PY_MAJOR_VERSION < 3
/* hex/oct formatting (mpz-only) */
static PyObject *
Pympz_oct(MPZ_Object *self)
{
    CTXT_Object *context = NULL;

    return GMPy_PyStr_From_MPZ(self, 8, 0, context);
}

static PyObject *
Pympz_hex(MPZ_Object *self)
{
    CTXT_Object *context = NULL;

    return GMPy_PyStr_From_MPZ(self, 16, 0, context);
}
#endif

/* Miscellaneous gmpy functions */
PyDoc_STRVAR(doc_gcd,
"gcd(a, b) -> mpz\n\n"
"Return the greatest common denominator of integers a and b.");

static PyObject *
Pygmpy_gcd(PyObject *self, PyObject *args)
{
    PyObject *a, *b;
    MPZ_Object *result, *tempa, *tempb;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcd() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_gcd(result->z, MPZ(a), MPZ(b));
    }
    else {
        tempa = GMPy_MPZ_From_Integer(a, context);
        tempb = GMPy_MPZ_From_Integer(b, context);
        if (!tempa || !tempb) {
            TYPE_ERROR("gcd() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_gcd(result->z, tempa->z, tempb->z);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_lcm,
"lcm(a, b) -> mpz\n\n"
"Return the lowest common multiple of integers a and b.");

static PyObject *
Pygmpy_lcm(PyObject *self, PyObject *args)
{
    PyObject *a, *b;
    MPZ_Object *result, *tempa, *tempb;
    CTXT_Object *context = NULL;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("lcm() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_lcm(result->z, MPZ(a), MPZ(b));
    }
    else {
        tempa = GMPy_MPZ_From_Integer(a, context);
        tempb = GMPy_MPZ_From_Integer(b, context);
        if (!tempa || !tempb) {
            TYPE_ERROR("lcm() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_lcm(result->z, tempa->z, tempb->z);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_gcdext,
"gcdext(a, b) - > tuple\n\n"
"Return a 3-element tuple (g,s,t) such that\n"
"    g == gcd(a,b) and g == a*s + b*t");

static PyObject *
Pygmpy_gcdext(PyObject *self, PyObject *args)
{
    PyObject *a, *b, *result = 0;
    MPZ_Object *g = 0, *s = 0, *t = 0, *tempa, *tempb;
    CTXT_Object *context = NULL;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
        return NULL;
    }

    g = GMPy_MPZ_New(context);
    s = GMPy_MPZ_New(context);
    t = GMPy_MPZ_New(context);
    result = PyTuple_New(3);
    if (!g || !s || !t || !result) {
        Py_XDECREF((PyObject*)g);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)t);
        Py_XDECREF(result);
        return NULL;
    }

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_gcdext(g->z, s->z, t->z, MPZ(a), MPZ(b));
    }
    else {
        tempa = GMPy_MPZ_From_Integer(a, context);
        tempb = GMPy_MPZ_From_Integer(b, context);
        if(!tempa || !tempb) {
            TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)g);
            Py_DECREF((PyObject*)s);
            Py_DECREF((PyObject*)t);
            Py_DECREF(result);
            return NULL;
        }
        mpz_gcdext(g->z, s->z, t->z, tempa->z, tempb->z);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)g);
    PyTuple_SET_ITEM(result, 1, (PyObject*)s);
    PyTuple_SET_ITEM(result, 2, (PyObject*)t);
    return result;
}

PyDoc_STRVAR(doc_divm,
"divm(a, b, m) -> mpz\n\n"
"Return x such that b*x == a mod m. Raises a ZeroDivisionError\n"
"exception if no such value x exists.");

static PyObject *
Pygmpy_divm(PyObject *self, PyObject *args)
{
    MPZ_Object *result, *num, *den, *mod;
    mpz_t gcdz;
    int ok;
    CTXT_Object *context = NULL;

    if(PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    num = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context);
    den = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), context);
    mod = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), context);

    if(!num || !den || !mod) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        Py_XDECREF((PyObject*)mod);
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    if (mpz_invert(result->z, den->z, mod->z)) { /* inverse exists */
        ok = 1;
    }
    else {
        /* last-ditch attempt: do num, den AND mod have a gcd>1 ? */
        mpz_inoc(gcdz);
        mpz_gcd(gcdz, num->z, den->z);
        mpz_gcd(gcdz, gcdz, mod->z);
        mpz_divexact(num->z, num->z, gcdz);
        mpz_divexact(den->z, den->z, gcdz);
        mpz_divexact(mod->z, mod->z, gcdz);
        mpz_cloc(gcdz);
        ok = mpz_invert(result->z, den->z, mod->z);
    }

    if (ok) {
        mpz_mul(result->z, result->z, num->z);
        mpz_mod(result->z, result->z, mod->z);
        Py_DECREF((PyObject*)num);
        Py_DECREF((PyObject*)den);
        Py_DECREF((PyObject*)mod);
        return (PyObject*)result;
    }
    else {
        ZERO_ERROR("not invertible");
        Py_DECREF((PyObject*)num);
        Py_DECREF((PyObject*)den);
        Py_DECREF((PyObject*)mod);
        Py_DECREF((PyObject*)result);
        return NULL;
    }
}

PyDoc_STRVAR(doc_fac,
"fac(n) -> mpz\n\n"
"Return the exact factorial of n.\n\n"
"See factorial(n) to get the floating-point approximation.");

static PyObject *
Pygmpy_fac(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    mpir_si n;
    CTXT_Object *context = NULL;

    n = SI_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("fac() requires 'int' argument");
        return NULL;
    }
    else if (n < 0) {
        VALUE_ERROR("fac() of negative number");
        return NULL;
    }
    else {
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_fac_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_fib,
"fib(n) -> mpz\n\n"
"Return the n-th Fibonacci number.");

static PyObject *
Pygmpy_fib(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    mpir_si n;
    CTXT_Object *context = NULL;

    n = SI_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("fib() requires 'int' argument");
        return NULL;
    }
    else if (n < 0) {
        VALUE_ERROR("Fibonacci of negative number");
        return NULL;
    }
    else {
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_fib_ui(MPZ(result), n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_fib2,
"fib2(n) -> tuple\n\n"
"Return a 2-tuple with the (n-1)-th and n-th Fibonacci numbers.");

static PyObject *
Pygmpy_fib2(PyObject *self, PyObject *other)
{
    PyObject *result;
    MPZ_Object *fib1 = 0, *fib2 = 0;
    mpir_si n;
    CTXT_Object *context = NULL;

    n = SI_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("fib2() requires 'int' argument");
        return NULL;
    }
    else if (n < 0) {
        VALUE_ERROR("Fibonacci of negative number");
        return NULL;
    }
    else {
        CREATE_TWO_MPZ_TUPLE(fib1, fib2, result);
        mpz_fib2_ui(fib1->z, fib2->z, n);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)fib1);
    PyTuple_SET_ITEM(result, 1, (PyObject*)fib2);
    return result;
}

PyDoc_STRVAR(doc_lucas,
"lucas(n) -> mpz\n\n"
"Return the n-th Lucas number.");

static PyObject *
Pygmpy_lucas(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    mpir_si n;
    CTXT_Object *context = NULL;

    n = SI_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("luc() requires 'int' argument");
        return NULL;
    }
    else if (n < 0) {
        VALUE_ERROR("Lucas of negative number");
        return NULL;
    }
    else {
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_lucnum_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_lucas2,
"lucas2(n) -> tuple\n\n"
"Return a 2-tuple with the (n-1)-th and n-th Lucas numbers.");

static PyObject *
Pygmpy_lucas2(PyObject *self, PyObject *other)
{
    PyObject *result;
    MPZ_Object *luc1, *luc2;
    mpir_si n;
    CTXT_Object *context = NULL;

    n = SI_From_Integer(other);
    if ((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("luc2() requires 'int' argument");
        return NULL;
    }
    else if (n < 0) {
        VALUE_ERROR("Lucas of negative number");
        return NULL;
    }
    else {
        CREATE_TWO_MPZ_TUPLE(luc1, luc2, result);
        mpz_lucnum2_ui(luc1->z, luc2->z, n);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)luc1);
    PyTuple_SET_ITEM(result, 1, (PyObject*)luc2);
    return result;
}

PyDoc_STRVAR(doc_bincoefg,
"bincoef(x, n) -> mpz\n\n"
"Return the binomial coefficient ('x over n'). n >= 0.");

PyDoc_STRVAR(doc_combg,
"comb(x, n) -> mpz\n\n"
"Return the number of combinations of 'x things, taking n at a\n"
"time'. n >= 0.");

static PyObject *
Pympz_bincoef(PyObject *self, PyObject *args)
{
    MPZ_Object *result;
    mpir_si k;
    CTXT_Object *context = NULL;

    PARSE_ONE_MPZ_REQ_SI(&k, "bincoef() requires 'mpz','int' arguments");

    if (k < 0) {
        VALUE_ERROR("binomial coefficient with negative k");
        Py_DECREF(self);
        return NULL;
    }

    if(!(result = GMPy_MPZ_New(context))) {
        Py_DECREF(self);
        return NULL;
    }
    mpz_bin_ui(result->z, MPZ(self), k);
    Py_DECREF(self);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_isqrt,
"isqrt(x) -> mpz\n\n"
"Return the integer square root of an integer x. x >= 0.");

static PyObject *
Pympz_isqrt(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if (self && (CHECK_MPZANY(self))) {
        if (mpz_sgn(MPZ(self)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_sqrt(result->z, MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_sqrt(result->z, MPZ(other));
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("isqrt() requires 'mpz' argument");
            return NULL;
        }
        if (mpz_sgn(result->z) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_sqrt(result->z, result->z);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_isqrt_rem,
"isqrt_rem(x) -> tuple\n\n"
"Return a 2-element tuple (s,t) such that s=isqrt(x) and t=x-s*s.\n"
"x >=0.");

static PyObject *
Pympz_isqrt_rem(PyObject *self, PyObject *args)
{
    MPZ_Object *root = 0, *rem = 0;
    PyObject *result = 0;
    CTXT_Object *context = NULL;

    PARSE_ONE_MPZ("isqrt_rem() requires 'mpz' argument");

    if (mpz_sgn(MPZ(self)) < 0) {
        VALUE_ERROR("isqrt_rem() of negative number");
        Py_DECREF(self);
        return NULL;
    }

    root = GMPy_MPZ_New(context);
    rem = GMPy_MPZ_New(context);
    result = PyTuple_New(2);
    if (!root || !rem || !result) {
        Py_DECREF(self);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF((PyObject*)rem);
        return NULL;
    }

    mpz_sqrtrem(root->z, rem->z, MPZ(self));
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

PyDoc_STRVAR(doc_removeg,
"remove(x, f) -> tuple\n\n"
"Return a 2-element tuple (y,m) such that x=y*(f**m) and f does\n"
"not divide y. Remove the factor f from x as many times as\n"
"possible. m is the multiplicity f in x. f > 1.");

static PyObject *
Pympz_remove(PyObject *self, PyObject *args)
{
    MPZ_Object *result;
    PyObject *factor;
    size_t multiplicity;
    CTXT_Object *context = NULL;

    PARSE_TWO_MPZ(factor, "remove() requires 'mpz','mpz' arguments");

    if (mpz_cmp_si(MPZ(factor), 2) < 0) {
        VALUE_ERROR("factor must be > 1");
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context))) {
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }
    multiplicity = mpz_remove(result->z, MPZ(self), MPZ(factor));
    Py_DECREF(self);
    Py_DECREF(factor);
    return Py_BuildValue("(Nk)", result, multiplicity);
}

PyDoc_STRVAR(doc_invertg,
"invert(x, m) -> mpz\n\n"
"Return y such that x*y == 1 modulo m. Raises ZeroDivisionError i no \n"
"inverse exists.");

static PyObject *
Pygmpy_invert(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    MPZ_Object *result, *tempx, *tempy;
    int success;
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;
    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        success = mpz_invert(result->z, MPZ(x), MPZ(y));
        if (!success) {
            ZERO_ERROR("invert() no inverse exists");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }
    else {
        tempx = GMPy_MPZ_From_Integer(x, context);
        tempy = GMPy_MPZ_From_Integer(y, context);
        if (!tempx || !tempy) {
            TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        success = mpz_invert(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        if (!success) {
            ZERO_ERROR("invert() no inverse exists");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_divexactg,
"divexact(x, y) -> mpz\n\n"
"Return the quotient of x divided by y. Faster than standard\n"
"division but requires the remainder is zero!");

static PyObject *
Pygmpy_divexact(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    MPZ_Object *result, *tempx, *tempy;
    CTXT_Object *context = NULL;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;
    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_divexact(result->z, MPZ(x), MPZ(y));
    }
    else {
        tempx = GMPy_MPZ_From_Integer(x, context);
        tempy = GMPy_MPZ_From_Integer(y, context);
        if (!tempx || !tempy) {
            TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(MPZ(tempy)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_divexact(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_is_squareg,
"is_square(x) -> bool\n\n"
"Returns True if x is a perfect square, else return False.");

static PyObject *
Pympz_is_square(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;
    CTXT_Object *context = NULL;

    if (CHECK_MPZANY(other)) {
        res = mpz_perfect_square_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("is_square() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_perfect_square_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_powerg,
"is_power(x) -> bool\n\n"
"Return True if x is a perfect power (there exists a y and an\n"
"n > 1, such that x=y**n), else return False.");

static PyObject *
Pympz_is_power(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object* tempx;
    CTXT_Object *context = NULL;

    if (CHECK_MPZANY(other)) {
        res = mpz_perfect_power_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("is_power() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_perfect_power_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }
    if(res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_primeg,
"is_prime(x[, n=25]) -> bool\n\n"
"Return True if x is _probably_ prime, else False if x is\n"
"definately composite. x is checked for small divisors and up\n"
"to n Miller-Rabin tests are performed.");

static PyObject *
Pympz_is_prime(PyObject *self, PyObject *args)
{
    int i, reps = 25;
    CTXT_Object *context = NULL;

    PARSE_ONE_MPZ_OPT_CLONG(&reps,
            "is_prime() requires 'mpz'[,'int'] arguments");

    if (reps <= 0) {
        VALUE_ERROR("repetition count for is_prime() must be positive");
        Py_DECREF(self);
        return NULL;
    }
    i = mpz_probab_prime_p(MPZ(self), reps);
    Py_DECREF(self);
    if (i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_next_primeg,
"next_prime(x) -> mpz\n\n"
"Return the next _probable_ prime number > x.");

static PyObject *
Pympz_next_prime(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if(CHECK_MPZANY(other)) {
        if(!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_nextprime(result->z, MPZ(other));
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("next_prime() requires 'mpz' argument");
            return NULL;
        }
        else {
            mpz_nextprime(result->z, result->z);
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_jacobig,
"jacobi(x, y) -> mpz\n\n"
"Return the Jacobi symbol (x|y). y must be odd and >0.");

static PyObject *
Pympz_jacobi(PyObject *self, PyObject *args)
{
    PyObject *other;
    long i;
    CTXT_Object *context = NULL;

    PARSE_TWO_MPZ(other, "jacobi() requires 'mpz','mpz' arguments");

    if (mpz_sgn(MPZ(other)) <= 0 ||
        mpz_even_p(MPZ(other))) {
        VALUE_ERROR("y must be odd and >0");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    i = (long)(mpz_jacobi(MPZ(self), MPZ(other)));
    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(i);
}

PyDoc_STRVAR(doc_legendreg,
"legendre(x, y) -> mpz\n\n"
"Return the Legendre symbol (x|y). y is assumed to be an odd prime.");

static PyObject *
Pympz_legendre(PyObject *self, PyObject *args)
{
    PyObject *other;
    long i;
    CTXT_Object *context = NULL;

    PARSE_TWO_MPZ(other, "legendre() requires 'mpz','mpz' arguments");

    if (mpz_sgn(MPZ(other)) <= 0 ||
        mpz_even_p(MPZ(other))) {
        VALUE_ERROR("y must be odd and >0");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    i = (long) mpz_legendre(MPZ(self), MPZ(other));
    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(i);
}

PyDoc_STRVAR(doc_kroneckerg,
"kronecker(x, y) -> mpz\n\n"
"Return the Kronecker-Jacobi symbol (x|y).");

static PyObject *
Pympz_kronecker(PyObject *self, PyObject *args)
{
    PyObject *other;
    long ires;
    CTXT_Object *context = NULL;

    PARSE_TWO_MPZ(other, "kronecker() requires 'mpz','mpz' arguments");

    ires = (long) mpz_kronecker(MPZ(self), (MPZ(other)));

    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(ires);
}

PyDoc_STRVAR(doc_is_eveng,
"is_even(x) -> bool\n\n"
"Return True if x is even, False otherwise.");

static PyObject *
Pympz_is_even(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;
    CTXT_Object *context = NULL;

    if (CHECK_MPZANY(other)) {
        res = mpz_even_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("is_even() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_even_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_is_oddg,
"is_odd(x) -> bool\n\n"
"Return True if x is odd, False otherwise.");

static PyObject *
Pympz_is_odd(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;
    CTXT_Object *context = NULL;

    if (CHECK_MPZANY(other)) {
        res = mpz_odd_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, context))) {
            TYPE_ERROR("is_odd() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_odd_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/*
 * Add mapping support to mpz objects.
 */

static Py_ssize_t
Pympz_nbits(MPZ_Object *self)
{
    return mpz_sizeinbase(self->z, 2);
}

static PyObject *
Pympz_subscript(MPZ_Object *self, PyObject *item)
{
    CTXT_Object *context = NULL;

    if (PyIndex_Check(item)) {
        Py_ssize_t i;

        i = PyNumber_AsSsize_t(item, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred())
            return NULL;
        if (i < 0)
            i += mpz_sizeinbase(self->z, 2);
        return PyIntOrLong_FromLong(mpz_tstbit(self->z, i));
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step, slicelength, cur, i;
        MPZ_Object *result;

#if PY_VERSION_HEX > 0x030200A4
        if (PySlice_GetIndicesEx(item,
#else
        if (PySlice_GetIndicesEx((PySliceObject*)item,
#endif
                        mpz_sizeinbase(self->z, 2),
                        &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }

        if ((step < 0 && start < stop) ||
            (step > 0 && start > stop))
            stop = start;

        if (!(result = GMPy_MPZ_New(context)))
            return NULL;

        mpz_set_ui(result->z, 0);
        if (slicelength > 0) {
            for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                if(mpz_tstbit(self->z, cur)) {
                    mpz_setbit(result->z, i);
                }
            }
        }
        return (PyObject*)result;
    }
    else {
        TYPE_ERROR("bit positions must be integers");
        return NULL;
    }
}

static PyObject *
Pympz_getnumer(MPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
Pympz_getdenom(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set_ui(result->z, 1);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_sizeof,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x. Note: deleted mpz objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
Pympz_sizeof(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPZ_Object) + \
        (MPZ(self)->_mp_alloc * sizeof(mp_limb_t)));
}

#ifdef PY3
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) GMPy_MPZ_Add_Slot,        /* nb_add                  */
    (binaryfunc) GMPy_MPZ_Sub_Slot,        /* nb_subtract             */
    (binaryfunc) GMPy_MPZ_Mul_Slot,        /* nb_multiply             */
    (binaryfunc) GMPy_MPZ_Mod_Slot,        /* nb_remainder            */
    (binaryfunc) GMPy_MPZ_DivMod_Slot,     /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,     /* nb_power                */
    (unaryfunc) GMPy_MPZ_Minus_Slot,       /* nb_negative             */
    (unaryfunc) GMPy_MPZ_Plus_Slot,        /* nb_positive             */
    (unaryfunc) GMPy_MPZ_Abs_Slot,         /* nb_absolute             */
    (inquiry) Pympz_nonzero,               /* nb_bool                 */
    (unaryfunc) GMPy_MPZ_Invert_Slot,      /* nb_invert               */
    (binaryfunc) GMPy_MPZ_Lshift_Slot,     /* nb_lshift               */
    (binaryfunc) GMPy_MPZ_Rshift_Slot,     /* nb_rshift               */
    (binaryfunc) GMPy_MPZ_And_Slot,        /* nb_and                  */
    (binaryfunc) GMPy_MPZ_Xor_Slot,        /* nb_xor                  */
    (binaryfunc) GMPy_MPZ_Ior_Slot,        /* nb_or                   */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_int                  */
        0,                                 /* nb_reserved             */
    (unaryfunc) GMPy_MPZ_Float_Slot,       /* nb_float                */
    (binaryfunc) GMPy_MPZ_IAdd_Slot,       /* nb_inplace_add          */
    (binaryfunc) GMPy_MPZ_ISub_Slot,       /* nb_inplace_subtract     */
    (binaryfunc) GMPy_MPZ_IMul_Slot,       /* nb_inplace_multiply     */
    (binaryfunc) GMPy_MPZ_IRem_Slot,       /* nb_inplace_remainder    */
    (ternaryfunc) GMPy_MPZ_IPow_Slot,      /* nb_inplace_power        */
    (binaryfunc) GMPy_MPZ_ILshift_Slot,    /* nb_inplace_lshift       */
    (binaryfunc) GMPy_MPZ_IRshift_Slot,    /* nb_inplace_rshift       */
        0,                                 /* nb_inplace_and          */
        0,                                 /* nb_inplace_xor          */
        0,                                 /* nb_inplace_or           */
    (binaryfunc) GMPy_MPZ_FloorDiv_Slot,   /* nb_floor_divide         */
    (binaryfunc) GMPy_MPZ_TrueDiv_Slot,    /* nb_true_divide          */
    (binaryfunc) GMPy_MPZ_IFloorDiv_Slot,  /* nb_inplace_floor_divide */
        0,                                 /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_index                */
};

#else
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) GMPy_MPZ_Add_Slot,        /* nb_add                  */
    (binaryfunc) GMPy_MPZ_Sub_Slot,        /* nb_subtract             */
    (binaryfunc) GMPy_MPZ_Mul_Slot,        /* nb_multiply             */
    (binaryfunc) GMPy_MPZ_Div2_Slot,       /* nb_divide               */
    (binaryfunc) GMPy_MPZ_Mod_Slot,        /* nb_remainder            */
    (binaryfunc) GMPy_MPZ_DivMod_Slot,     /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,     /* nb_power                */
    (unaryfunc) GMPy_MPZ_Minus_Slot,       /* nb_negative             */
    (unaryfunc) GMPy_MPZ_Plus_Slot,        /* nb_positive             */
    (unaryfunc) GMPy_MPZ_Abs_Slot,         /* nb_absolute             */
    (inquiry) Pympz_nonzero,               /* nb_bool                 */
    (unaryfunc) GMPy_MPZ_Invert_Slot,      /* nb_invert               */
    (binaryfunc) GMPy_MPZ_Lshift_Slot,     /* nb_lshift               */
    (binaryfunc) GMPy_MPZ_Rshift_Slot,     /* nb_rshift               */
    (binaryfunc) GMPy_MPZ_And_Slot,        /* nb_and                  */
    (binaryfunc) GMPy_MPZ_Xor_Slot,        /* nb_xor                  */
    (binaryfunc) GMPy_MPZ_Ior_Slot,        /* nb_or                   */
        0,                                 /* nb_coerce               */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_int                  */
    (unaryfunc) GMPy_MPZ_Long_Slot,        /* nb_long                 */
    (unaryfunc) GMPy_MPZ_Float_Slot,       /* nb_float                */
    (unaryfunc) Pympz_oct,                 /* nb_oct                  */
    (unaryfunc) Pympz_hex,                 /* nb_hex                  */
    (binaryfunc) GMPy_MPZ_IAdd_Slot,       /* nb_inplace_add          */
    (binaryfunc) GMPy_MPZ_ISub_Slot,       /* nb_inplace_subtract     */
    (binaryfunc) GMPy_MPZ_IMul_Slot,       /* nb_inplace_multiply     */
        0,                                 /* nb_inplace_divide       */
    (binaryfunc) GMPy_MPZ_IRem_Slot,       /* nb_inplace_remainder    */
    (ternaryfunc) GMPy_MPZ_IPow_Slot,      /* nb_inplace_power        */
    (binaryfunc) GMPy_MPZ_ILshift_Slot,    /* nb_inplace_lshift       */
    (binaryfunc) GMPy_MPZ_IRshift_Slot,    /* nb_inplace_rshift       */
        0,                                 /* nb_inplace_and          */
        0,                                 /* nb_inplace_xor          */
        0,                                 /* nb_inplace_or           */
    (binaryfunc) GMPy_MPZ_FloorDiv_Slot,   /* nb_floor_divide         */
    (binaryfunc) GMPy_MPZ_TrueDiv_Slot,    /* nb_true_divide          */
    (binaryfunc) GMPy_MPZ_IFloorDiv_Slot,  /* nb_inplace_floor_divide */
        0,                                 /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_index                */
};
#endif

static PyMappingMethods mpz_mapping_methods = {
    (lenfunc)Pympz_nbits,
    (binaryfunc)Pympz_subscript,
    NULL
};

static PyGetSetDef Pympz_getseters[] =
{
    { "numerator", (getter)Pympz_getnumer, NULL, "numerator", NULL },
    { "denominator", (getter)Pympz_getdenom, NULL, "denominator", NULL },
    {NULL}
};

static PyMethodDef Pympz_methods [] =
{
    { "__format__", GMPy_MPZ_Format, METH_VARARGS, GMPy_doc_mpz_format },
    { "__ceil__", Pympz_ceil, METH_NOARGS, doc_mpz_ceil },
    { "__floor__", Pympz_floor, METH_NOARGS, doc_mpz_floor },
    { "__round__", Pympz_round, METH_VARARGS, doc_mpz_round },
    { "__sizeof__", Pympz_sizeof, METH_NOARGS, doc_mpz_sizeof },
    { "__trunc__", Pympz_trunc, METH_NOARGS, doc_mpz_trunc },
    { "bit_clear", GMPy_MPZ_bit_clear_method, METH_O, doc_bit_clear_method },
    { "bit_flip", GMPy_MPZ_bit_flip_method, METH_O, doc_bit_flip_method },
    { "bit_length", GMPy_MPZ_bit_length_method, METH_NOARGS, doc_bit_length_method },
    { "bit_scan0", GMPy_MPZ_bit_scan0_method, METH_VARARGS, doc_bit_scan0_method },
    { "bit_scan1", GMPy_MPZ_bit_scan1_method, METH_VARARGS, doc_bit_scan1_method },
    { "bit_set", GMPy_MPZ_bit_set_method, METH_O, doc_bit_set_method },
    { "bit_test", GMPy_MPZ_bit_test_method, METH_O, doc_bit_test_method },
    { "digits", GMPy_MPZ_Digits_Method, METH_VARARGS, GMPy_doc_mpz_digits_method },
    { "num_digits", Pympz_num_digits, METH_VARARGS, doc_num_digitsm },
    { NULL, NULL, 1 }
};

static PyTypeObject MPZ_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpz",                                  /* tp_name          */
    sizeof(MPZ_Object),                     /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPZ_Dealloc,          /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPZ_Repr_Slot,          /* tp_repr          */
    &mpz_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &mpz_mapping_methods,                   /* tp_as_mapping    */
    (hashfunc) GMPy_MPZ_Hash_Slot,          /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPZ_Str_Slot,           /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_INDEX|Py_TPFLAGS_HAVE_RICHCOMPARE| \
    Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_HAVE_CLASS| \
    Py_TPFLAGS_HAVE_INPLACEOPS,
#endif
    "Multiple precision integer",           /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympz_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympz_getseters,                        /* tp_getset        */
};

