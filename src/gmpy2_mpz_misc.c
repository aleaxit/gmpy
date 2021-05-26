/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_misc.c                                                        *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
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

/* return number-of-digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(GMPy_doc_mpz_method_num_digits,
"x.num_digits([base]) -> int\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

PyDoc_STRVAR(GMPy_doc_mpz_function_num_digits,
"num_digits(x[, base]) -> int\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

static PyObject *
GMPy_MPZ_Method_NumDigits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;

    if (PyTuple_GET_SIZE(args) == 1) {
        base = PyIntOrLong_AsLong(PyTuple_GET_ITEM(args, 0));
        if (base == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }

    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval [2, 62]");
        return NULL;
    }

    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(MPZ(self), (int)base));
    return result;
}

static PyObject *
GMPy_MPZ_Function_NumDigits(PyObject *self, PyObject *args)
{
    long base = 10;
    Py_ssize_t argc;
    MPZ_Object *temp;
    PyObject *result;

    argc = PyTuple_GET_SIZE(args);
    if (argc == 0 || argc > 2) {
        TYPE_ERROR("num_digits() requires 'mpz',['int'] arguments");
        return NULL;
    }

    if (argc == 2) {
        base = PyIntOrLong_AsLong(PyTuple_GET_ITEM(args, 1));
        if (base == -1 && PyErr_Occurred()) {
            return NULL;
        }
    }

    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval [2, 62]");
        return NULL;
    }

    if (!(temp = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL))) {
        return NULL;
    }

    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(temp->z, (int)base));
    Py_DECREF((PyObject*)temp);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_iroot,
"iroot(x,n) -> (number, boolean)\n\n"
"Return the integer n-th root of x and boolean value that is True\n"
"iff the root is exact. x >= 0. n > 0.");

static PyObject *
GMPy_MPZ_Function_Iroot(PyObject *self, PyObject *args)
{
    unsigned long n;
    int exact;
    MPZ_Object *root = NULL, *tempx = NULL;
    PyObject *result = NULL;

    if ((PyTuple_GET_SIZE(args) != 2) ||
        ((!IS_INTEGER(PyTuple_GET_ITEM(args, 0))) ||
         (!IS_INTEGER(PyTuple_GET_ITEM(args, 1))))) {

        TYPE_ERROR("iroot() requires 'int','int' arguments");
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 1));
    if ((n == 0) || ((n == (unsigned long)(-1)) && PyErr_Occurred())) {
        VALUE_ERROR("n must be > 0");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_sgn(tempx->z) < 0) {
        VALUE_ERROR("iroot() of negative number");
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(root = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF(result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    exact = mpz_root(root->z, tempx->z, n);
    Py_DECREF((PyObject*)tempx);

    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)PyBool_FromLong(exact));
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_iroot_rem,
"iroot_rem(x,n) -> (number, number)\n\n"
"Return a 2-element tuple (y,r), such that y is the integer n-th\n"
"root of x and x=y**n + r. x >= 0. n > 0.");

static PyObject *
GMPy_MPZ_Function_IrootRem(PyObject *self, PyObject *args)
{
    unsigned long n;
    MPZ_Object *root = NULL, *rem = NULL, *tempx = NULL;
    PyObject *result = NULL;

    if ((PyTuple_GET_SIZE(args) != 2) ||
        ((!IS_INTEGER(PyTuple_GET_ITEM(args, 0))) ||
         (!IS_INTEGER(PyTuple_GET_ITEM(args, 1))))) {

        TYPE_ERROR("iroot_rem() requires 'int','int' arguments");
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 1));
    if ((n == 0) || ((n == (unsigned long)(-1)) && PyErr_Occurred())) {
        VALUE_ERROR("n must be > 0");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (mpz_sgn(tempx->z) < 0) {
        VALUE_ERROR("iroot_rem() of negative number");
        Py_DECREF((PyObject*)tempx);
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(root = GMPy_MPZ_New(NULL)) ||
        !(rem = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)tempx);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF((PyObject*)rem);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_rootrem(root->z, rem->z, tempx->z, n);
    Py_DECREF((PyObject*)tempx);

    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_ceil, "Ceiling of an mpz returns itself.");

static PyObject *
GMPy_MPZ_Method_Ceil(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_floor, "Floor of an mpz returns itself.");

static PyObject *
GMPy_MPZ_Method_Floor(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_trunc, "Truncating an mpz returns itself.");

static PyObject *
GMPy_MPZ_Method_Trunc(PyObject *self, PyObject *other)
{
    Py_INCREF(self);
    return self;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_round, "Round an mpz to power of 10.");

static PyObject *
GMPy_MPZ_Method_Round(PyObject *self, PyObject *args)
{
    Py_ssize_t round_digits;
    MPZ_Object *result;
    mpz_t temp, rem;

    if (PyTuple_GET_SIZE(args) == 0) {
        Py_INCREF(self);
        return self;
    }

    round_digits = GMPy_Integer_AsSsize_t(PyTuple_GET_ITEM(args, 0));
    if (round_digits == -1 && PyErr_Occurred()) {
        TYPE_ERROR("__round__() requires 'int' argument");
        return NULL;
    }

    if (round_digits >= 0) {
        Py_INCREF(self);
        return self;
    }

    /* We can now assume round_digits > 0. */
    round_digits = -round_digits;

    if ((result = GMPy_MPZ_New(NULL))) {
        if ((unsigned)round_digits >= mpz_sizeinbase(MPZ(self), 10)) {
            mpz_set_ui(result->z, 0);
        }
        else {
            mpz_init(temp);
            mpz_init(rem);
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
            mpz_clear(rem);
            mpz_clear(temp);
        }
    }

    return (PyObject*)result;
}

static int
GMPy_MPZ_NonZero_Slot(MPZ_Object *self)
{
    return mpz_sgn(self->z) != 0;
}

#if PY_MAJOR_VERSION < 3

/* hex/oct formatting (mpz-only) */

static PyObject *
GMPy_MPZ_Oct_Slot(MPZ_Object *self)
{
    return GMPy_PyStr_From_MPZ(self, 8, 0, NULL);
}

static PyObject *
GMPy_MPZ_Hex_Slot(MPZ_Object *self)
{
    return GMPy_PyStr_From_MPZ(self, 16, 0, NULL);
}
#endif

/* Miscellaneous gmpy functions */

PyDoc_STRVAR(GMPy_doc_mpz_function_gcd,
"gcd(*integers) -> mpz\n\n"
"Return the greatest common divisor of integers.");

static PyObject *
GMPy_MPZ_Function_GCD(PyObject *self, PyObject *args)
{
    PyObject *arg;
    MPZ_Object *result = NULL, *tempa = NULL;
    CTXT_Object *context = NULL;
    Py_ssize_t i, nargs;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    nargs = PyTuple_Size(args);

    for (i = 0; i < nargs; i++) {
        arg = PyTuple_GET_ITEM(args, i);

        if MPZ_Check(arg) {
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_gcd(result->z, MPZ(arg), result->z);
            GMPY_MAYBE_END_ALLOW_THREADS(context);\
        } else {
            if (!(tempa = GMPy_MPZ_From_Integer(arg, NULL))) {
                TYPE_ERROR("gcd() requires 'mpz' arguments");
                Py_XDECREF((PyObject*)tempa);
                Py_DECREF((PyObject*)result);
                return NULL;
            }

            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_gcd(result->z, tempa->z, result->z);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            Py_DECREF((PyObject*)tempa);
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_lcm,
"lcm(*integers) -> mpz\n\n"
"Return the lowest common multiple of integers.");

static PyObject *
GMPy_MPZ_Function_LCM(PyObject *self, PyObject *args)
{
    PyObject *arg;
    MPZ_Object *result = NULL, *tempa = NULL;
    CTXT_Object *context = NULL;
    Py_ssize_t i, nargs;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_set_ui(result->z, 1);
    nargs = PyTuple_Size(args);

    for (i = 0; i < nargs; i++) {
        arg = PyTuple_GET_ITEM(args, i);

        if MPZ_Check(arg) {
            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_lcm(result->z, MPZ(arg), result->z);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
        } else {
            if (!(tempa = GMPy_MPZ_From_Integer(arg, NULL))) {
                TYPE_ERROR("lcm() requires 'mpz' arguments");
                Py_XDECREF((PyObject*)tempa);
                Py_DECREF((PyObject*)result);
                return NULL;
            }

            GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
            mpz_lcm(result->z, tempa->z, result->z);
            GMPY_MAYBE_END_ALLOW_THREADS(context);
            Py_DECREF((PyObject*)tempa);
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_gcdext,
"gcdext(a, b) - > tuple\n\n"
"Return a 3-element tuple (g,s,t) such that\n"
"    g == gcd(a,b) and g == a*s + b*t");

static PyObject *
GMPy_MPZ_Function_GCDext(PyObject *self, PyObject *args)
{
    PyObject *arg0, *arg1, *result = NULL;
    MPZ_Object *g = NULL, *s = NULL, *t = NULL, *tempa = NULL, *tempb = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = PyTuple_New(3)) ||
        !(g = GMPy_MPZ_New(NULL)) ||
        !(s = GMPy_MPZ_New(NULL)) ||
        !(t = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_XDECREF((PyObject*)g);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)t);
        Py_XDECREF(result);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    arg0 = PyTuple_GET_ITEM(args, 0);
    arg1 = PyTuple_GET_ITEM(args, 1);

    if (MPZ_Check(arg0) && MPZ_Check(arg1)) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_gcdext(g->z, s->z, t->z, MPZ(arg0), MPZ(arg1));
        GMPY_MAYBE_END_ALLOW_THREADS(context);
    }
    else {
        if(!(tempa = GMPy_MPZ_From_Integer(arg0, NULL)) ||
           !(tempb = GMPy_MPZ_From_Integer(arg1, NULL))) {

            TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)g);
            Py_DECREF((PyObject*)s);
            Py_DECREF((PyObject*)t);
            Py_DECREF(result);
            return NULL;
        }
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_gcdext(g->z, s->z, t->z, tempa->z, tempb->z);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    PyTuple_SET_ITEM(result, 0, (PyObject*)g);
    PyTuple_SET_ITEM(result, 1, (PyObject*)s);
    PyTuple_SET_ITEM(result, 2, (PyObject*)t);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_divm,
"divm(a, b, m) -> mpz\n\n"
"Return x such that b*x == a mod m. Raises a ZeroDivisionError\n"
"exception if no such value x exists.");

static PyObject *
GMPy_MPZ_Function_Divm(PyObject *self, PyObject *args)
{
    MPZ_Object *result = NULL, *num = NULL, *den = NULL, *mod = NULL;
    mpz_t numz, denz, modz, gcdz;
    int ok = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    if (!(num = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(den = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL)) ||
        !(mod = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL))) {

        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        Py_XDECREF((PyObject*)mod);
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    /* Make copies so we don't destroy the input. */
    mpz_init(numz);
    mpz_init(denz);
    mpz_init(modz);
    mpz_set(numz, num->z);
    mpz_set(denz, den->z);
    mpz_set(modz, mod->z);
    Py_DECREF((PyObject*)num);
    Py_DECREF((PyObject*)den);
    Py_DECREF((PyObject*)mod);


    GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
    ok = mpz_invert(result->z, denz, modz);
    GMPY_MAYBE_END_ALLOW_THREADS(context);

    if (!ok) {
        /* last-ditch attempt: do num, den AND mod have a gcd>1 ? */
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_init(gcdz);
        mpz_gcd(gcdz, numz, denz);
        mpz_gcd(gcdz, gcdz, modz);
        mpz_divexact(numz, numz, gcdz);
        mpz_divexact(denz, denz, gcdz);
        mpz_divexact(modz, modz, gcdz);
        mpz_clear(gcdz);
        ok = mpz_invert(result->z, denz, modz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
    }

    if (ok) {
        GMPY_MAYBE_BEGIN_ALLOW_THREADS(context);
        mpz_mul(result->z, result->z, numz);
        mpz_mod(result->z, result->z, modz);
        mpz_clear(numz);
        mpz_clear(denz);
        mpz_clear(modz);
        GMPY_MAYBE_END_ALLOW_THREADS(context);
        return (PyObject*)result;
    }
    else {
        ZERO_ERROR("not invertible");
        mpz_clear(numz);
        mpz_clear(denz);
        mpz_clear(modz);
        Py_DECREF((PyObject*)result);
        return NULL;
    }
}

PyDoc_STRVAR(GMPy_doc_mpz_function_fac,
"fac(n) -> mpz\n\n"
"Return the exact factorial of n.\n\n"
"See factorial(n) to get the floating-point approximation.");

static PyObject *
GMPy_MPZ_Function_Fac(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_fac_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_double_fac,
"double_fac(n) -> mpz\n\n"
"Return the exact double factorial (n!!) of n. The double\n"
"factorial is defined as n*(n-2)*(n-4)...");

static PyObject *
GMPy_MPZ_Function_DoubleFac(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_2fac_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_primorial,
"primorial(n) -> mpz\n\n"
"Return the product of all positive prime numbers less than or"
"equal to n.");

static PyObject *
GMPy_MPZ_Function_Primorial(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_primorial_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_multi_fac,
"multi_fac(n,m) -> mpz\n\n"
"Return the exact m-multi factorial of n. The m-multi"
"factorial is defined as n*(n-m)*(n-2m)...");

static PyObject *
GMPy_MPZ_Function_MultiFac(PyObject *self, PyObject *args)
{
    MPZ_Object *result = NULL;
    unsigned long n, m;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("multi_fac() requires 2 integer arguments");
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 0));
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    m = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 1));
    if (m == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_mfac_uiui(result->z, n, m);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_fib,
"fib(n) -> mpz\n\n"
"Return the n-th Fibonacci number.");

static PyObject *
GMPy_MPZ_Function_Fib(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long  n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }
    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_fib_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_fib2,
"fib2(n) -> tuple\n\n"
"Return a 2-tuple with the (n-1)-th and n-th Fibonacci numbers.");

static PyObject *
GMPy_MPZ_Function_Fib2(PyObject *self, PyObject *other)
{
    PyObject *result = NULL;
    MPZ_Object *fib1 = NULL, *fib2 = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(fib1 = GMPy_MPZ_New(NULL)) ||
        !(fib2 = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)fib1);
        Py_XDECREF((PyObject*)fib2);
        result = NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_fib2_ui(fib1->z, fib2->z, n);
    PyTuple_SET_ITEM(result, 0, (PyObject*)fib1);
    PyTuple_SET_ITEM(result, 1, (PyObject*)fib2);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_lucas,
"lucas(n) -> mpz\n\n"
"Return the n-th Lucas number.");

static PyObject *
GMPy_MPZ_Function_Lucas(PyObject *self, PyObject *other)
{
    MPZ_Object *result = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_lucnum_ui(result->z, n);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_lucas2,
"lucas2(n) -> tuple\n\n"
"Return a 2-tuple with the (n-1)-th and n-th Lucas numbers.");

static PyObject *
GMPy_MPZ_Function_Lucas2(PyObject *self, PyObject *other)
{
    PyObject *result = NULL;
    MPZ_Object *luc1 = NULL, *luc2 = NULL;
    unsigned long n;

    n = GMPy_Integer_AsUnsignedLong(other);
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(luc1 = GMPy_MPZ_New(NULL)) ||
        !(luc2 = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)luc1);
        Py_XDECREF((PyObject*)luc2);
        result = NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_lucnum2_ui(luc1->z, luc2->z, n);
    PyTuple_SET_ITEM(result, 0, (PyObject*)luc1);
    PyTuple_SET_ITEM(result, 1, (PyObject*)luc2);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_bincoef,
"bincoef(n, k) -> mpz\n\n"
"Return the binomial coefficient ('n choose k'). k >= 0.");

PyDoc_STRVAR(GMPy_doc_mpz_function_comb,
"comb(n, k) -> mpz\n\n"
"Return the number of combinations of 'n things, taking k at a\n"
"time'. k >= 0. Same as bincoef(n, k)");

static PyObject *
GMPy_MPZ_Function_Bincoef(PyObject *self, PyObject *args)
{
    MPZ_Object *result = NULL, *tempx;
    unsigned long n, k;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bincoef() requires two integer arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    k = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 1));
    if (k == (unsigned long)(-1) && PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    n = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 0));
    if (n == (unsigned long)(-1) && PyErr_Occurred()) {
        /* Since we plan to skip the else clause and continue,
         * we need to clear the error since we aren't acting on it.
         */
        PyErr_Clear();
    }
    else {
        /* Use mpz_bin_uiui which should be faster. */
        mpz_bin_uiui(result->z, n, k);
        return (PyObject*)result;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL))) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    mpz_bin_ui(result->z, tempx->z, k);
    Py_DECREF((PyObject*)tempx);
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_isqrt,
"isqrt(x) -> mpz\n\n"
"Return the integer square root of an integer x. x >= 0.");

static PyObject *
GMPy_MPZ_Function_Isqrt(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if ((result = GMPy_MPZ_New(NULL))) {
            mpz_sqrt(result->z, MPZ(other));
        }
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(other, NULL))) {
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

PyDoc_STRVAR(GMPy_doc_mpz_function_isqrt_rem,
"isqrt_rem(x) -> tuple\n\n"
"Return a 2-element tuple (s,t) such that s=isqrt(x) and t=x-s*s.\n"
"x >=0.");

static PyObject *
GMPy_MPZ_Function_IsqrtRem(PyObject *self, PyObject *other)
{
    MPZ_Object *root = NULL, *rem = NULL, *temp = NULL;
    PyObject *result;

    if (!(temp = GMPy_MPZ_From_Integer(other, NULL))) {
        TYPE_ERROR("isqrt_rem() requires 'mpz' argument");
        return NULL;
    }

    if (mpz_sgn(temp->z) < 0) {
        VALUE_ERROR("isqrt_rem() of negative number");
        Py_DECREF((PyObject*)temp);
        return NULL;
    }

    if (!(result = PyTuple_New(2)) ||
        !(root = GMPy_MPZ_New(NULL)) ||
        !(rem = GMPy_MPZ_New(NULL))) {

        /* LCOV_EXCL_START */
        Py_DECREF((PyObject*)temp);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF((PyObject*)rem);
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    mpz_sqrtrem(root->z, rem->z, temp->z);
    Py_DECREF((PyObject*)temp);
    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_remove,
"remove(x, f) -> tuple\n\n"
"Return a 2-element tuple (y,m) such that x=y*(f**m) and f does\n"
"not divide y. Remove the factor f from x as many times as\n"
"possible. m is the multiplicity f in x. f > 1.");

static PyObject *
GMPy_MPZ_Function_Remove(PyObject *self, PyObject *args)
{
    MPZ_Object *result = NULL, *tempx = NULL, *tempf = NULL;
    PyObject *x, *f;
    size_t multiplicity;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("remove() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    x = PyTuple_GET_ITEM(args, 0);
    f = PyTuple_GET_ITEM(args, 1);

    if (MPZ_Check(x) && MPZ_Check(f)) {
        if (mpz_cmp_si(MPZ(f), 2) < 0) {
            VALUE_ERROR("factor must be > 1");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        multiplicity = mpz_remove(result->z, MPZ(x), MPZ(f));
        return Py_BuildValue("(Nk)", result, multiplicity);
    }
    else {


        if (!(tempx = GMPy_MPZ_From_Integer(x, NULL)) ||
            !(tempf = GMPy_MPZ_From_Integer(f, NULL))) {

            TYPE_ERROR("remove() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempf);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_cmp_si(MPZ(tempf), 2) < 0) {
            VALUE_ERROR("factor must be > 1");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempf);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        multiplicity = mpz_remove(result->z, tempx->z, tempf->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempf);
        return Py_BuildValue("(Nk)", result, multiplicity);
    }
}

PyDoc_STRVAR(GMPy_doc_mpz_function_invert,
"invert(x, m) -> mpz\n\n"
"Return y such that x*y == 1 modulo m. Raises ZeroDivisionError if no\n"
"inverse exists.");

static PyObject *
GMPy_MPZ_Function_Invert(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    MPZ_Object *result = NULL, *tempx = NULL, *tempy = NULL;
    int success;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    if (MPZ_Check(x) && MPZ_Check(y)) {
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


        if (!(tempx = GMPy_MPZ_From_Integer(x, NULL)) ||
            !(tempy = GMPy_MPZ_From_Integer(y, NULL))) {

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

PyDoc_STRVAR(GMPy_doc_mpz_function_divexact,
"divexact(x, y) -> mpz\n\n"
"Return the quotient of x divided by y. Faster than standard\n"
"division but requires the remainder is zero!");

static PyObject *
GMPy_MPZ_Function_Divexact(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    MPZ_Object *result, *tempx= NULL, *tempy = NULL;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(NULL))) {
        /* LCOV_EXCL_START */
        return NULL;
        /* LCOV_EXCL_STOP */
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    if (MPZ_Check(x) && MPZ_Check(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_divexact(result->z, MPZ(x), MPZ(y));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(x, NULL)) ||
            !(tempy = GMPy_MPZ_From_Integer(y, NULL))) {

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

PyDoc_STRVAR(GMPy_doc_mpz_function_is_square,
"is_square(x) -> bool\n\n"
"Returns True if x is a perfect square, else return False.");

static PyObject *
GMPy_MPZ_Function_IsSquare(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;

    if (MPZ_Check(other)) {
        res = mpz_perfect_square_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
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

PyDoc_STRVAR(GMPy_doc_mpz_method_is_square,
"x.is_square() -> bool\n\n"
"Returns True if x is a perfect square, else return False.");

static PyObject *
GMPy_MPZ_Method_IsSquare(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_perfect_square_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_divisible,
"is_divisible(x, d) -> bool\n\n"
"Returns True if x is divisible by d, else return False.");

static PyObject *
GMPy_MPZ_Function_IsDivisible(PyObject *self, PyObject *args)
{
    unsigned long temp;
    int res = 0;
    MPZ_Object *tempx, *tempd;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("is_divisible() requires 2 integer arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL))) {
        return NULL;
    }

    temp = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 1));
    if (temp == (unsigned long)-1 && PyErr_Occurred()) {
        PyErr_Clear();
        /* Implement mpz_divisible_p here. */

        if (!(tempd = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL))) {
            TYPE_ERROR("is_divisible() requires 2 integer arguments");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }

        res = mpz_divisible_p(tempx->z, tempd->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempd);
    }
    else {
        /* Implement mpz_divisible_ui_p here. */

        res = mpz_divisible_ui_p(tempx->z, temp);
        Py_DECREF((PyObject*)tempx);
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_divisible,
"x.is_divisible(d) -> bool\n\n"
"Returns True if x is divisible by d, else return False.");

static PyObject *
GMPy_MPZ_Method_IsDivisible(PyObject *self, PyObject *other)
{
    unsigned long temp;
    int res;
    MPZ_Object *tempd;

    temp = GMPy_Integer_AsUnsignedLong(other);
    if (temp == (unsigned long)-1 && PyErr_Occurred()) {
        PyErr_Clear();
        /* Implement mpz_divisible_p here. */

        if (!(tempd = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_divisible() requires 2 integer arguments");
            return NULL;
        }

        res = mpz_divisible_p(MPZ(self), tempd->z);
        Py_DECREF((PyObject*)tempd);
    }
    else {
        /* Implement mpz_divisible_ui_p here. */

        res = mpz_divisible_ui_p(MPZ(self), temp);
    }
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_congruent,
"is_congruent(x, y, m) -> bool\n\n"
"Returns True if x is congruent to y modulo m, else return False.");

static PyObject *
GMPy_MPZ_Function_IsCongruent(PyObject *self, PyObject *args)
{
    int res;
    MPZ_Object *tempx = NULL, *tempy = NULL, *tempm = NULL;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("is_congruent() requires 3 integer arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL)) ||
        !(tempm = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 2), NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)tempm);
        TYPE_ERROR("is_congruent() requires 3 integer arguments");
        return NULL;
    }

    res = mpz_congruent_p(tempx->z, tempy->z, tempm->z);
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempm);
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_congruent,
"x.is_congruent(y, m) -> bool\n\n"
"Returns True if x is congruent to y modulo m, else return False.");

static PyObject *
GMPy_MPZ_Method_IsCongruent(PyObject *self, PyObject *args)
{
    int res;
    MPZ_Object *tempy = NULL, *tempm = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("is_congruent() requires 2 integer arguments");
        return NULL;
    }

    if (!(tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(tempm = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL))) {

        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)tempm);
        TYPE_ERROR("is_congruent() requires 2 integer arguments");
        return NULL;
    }

    res = mpz_congruent_p(MPZ(self), tempy->z, tempm->z);
    Py_DECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)tempm);
    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_power,
"is_power(x) -> bool\n\n"
"Return True if x is a perfect power (there exists a y and an\n"
"n > 1, such that x=y**n), else return False.");

static PyObject *
GMPy_MPZ_Function_IsPower(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object* tempx;

    if (MPZ_Check(other)) {
        res = mpz_perfect_power_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("is_power() requires 'mpz' argument");
            return NULL;
        }
        else {
            res = mpz_perfect_power_p(tempx->z);
            Py_DECREF((PyObject*)tempx);
        }
    }

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_power,
"x.is_power() -> bool\n\n"
"Return True if x is a perfect power (there exists a y and an\n"
"n > 1, such that x=y**n), else return False.");

static PyObject *
GMPy_MPZ_Method_IsPower(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_perfect_power_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_prime,
"is_prime(x[, n=25]) -> bool\n\n"
"Return True if x is _probably_ prime, else False if x is\n"
"definitely composite. x is checked for small divisors and up\n"
"to n Miller-Rabin tests are performed.");

static PyObject *
GMPy_MPZ_Function_IsPrime(PyObject *self, PyObject *args)
{
    int i;
    unsigned long reps = 25;
    MPZ_Object* tempx;
    Py_ssize_t argc;

    argc = PyTuple_GET_SIZE(args);

    if (argc == 0 || argc > 2) {
        TYPE_ERROR("is_prime() requires 'mpz'[,'int'] arguments");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 2) {
        reps = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 1));
        if (reps == (unsigned long)(-1) && PyErr_Occurred()) {
            return NULL;
        }
        /* Silently limit n to a reasonable value. */
        if (reps > 1000) {
            reps = 1000;
        }
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL))) {
        return NULL;
    }

    i = mpz_probab_prime_p(tempx->z, (int)reps);
    Py_DECREF((PyObject*)tempx);

    if (i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_is_prime,
"x.is_prime([n=25]) -> bool\n\n"
"Return True if x is _probably_ prime, else False if x is\n"
"definitely composite. x is checked for small divisors and up\n"
"to n Miller-Rabin tests are performed.");

static PyObject *
GMPy_MPZ_Method_IsPrime(PyObject *self, PyObject *args)
{
    int i;
    unsigned long reps = 25;
    Py_ssize_t argc;

    argc = PyTuple_GET_SIZE(args);

    if (argc > 1) {
        TYPE_ERROR("is_prime() takes at most 1 argument");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        reps = GMPy_Integer_AsUnsignedLong(PyTuple_GET_ITEM(args, 0));
        if (reps == (unsigned long)(-1) && PyErr_Occurred()) {
            return NULL;
        }
        /* Silently limit n to a reasonable value. */
        if (reps > 1000) {
            reps = 1000;
        }
    }

    i = mpz_probab_prime_p(MPZ(self), (int)reps);

    if (i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_next_prime,
"next_prime(x) -> mpz\n\n"
"Return the next _probable_ prime number > x.");

static PyObject *
GMPy_MPZ_Function_NextPrime(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if(MPZ_Check(other)) {
        if(!(result = GMPy_MPZ_New(NULL))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpz_nextprime(result->z, MPZ(other));
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(other, NULL))) {
            TYPE_ERROR("next_prime() requires 'mpz' argument");
            return NULL;
        }
        else {
            mpz_nextprime(result->z, result->z);
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_jacobi,
"jacobi(x, y) -> mpz\n\n"
"Return the Jacobi symbol (x|y). y must be odd and >0.");

static PyObject *
GMPy_MPZ_Function_Jacobi(PyObject *self, PyObject *args)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    long res;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("jacobi() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    if (mpz_sgn(tempy->z) <= 0 || mpz_even_p(tempy->z)) {
        VALUE_ERROR("y must be odd and >0");
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    res = (long)(mpz_jacobi(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return PyIntOrLong_FromLong(res);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_legendre,
"legendre(x, y) -> mpz\n\n"
"Return the Legendre symbol (x|y). y is assumed to be an odd prime.");

static PyObject *
GMPy_MPZ_Function_Legendre(PyObject *self, PyObject *args)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    long res;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("legendre() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    if (mpz_sgn(tempy->z) <= 0 || mpz_even_p(tempy->z)) {
        VALUE_ERROR("y must be odd, prime, and >0");
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    res = (long)(mpz_legendre(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return PyIntOrLong_FromLong(res);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_kronecker,
"kronecker(x, y) -> mpz\n\n"
"Return the Kronecker-Jacobi symbol (x|y).");

static PyObject *
GMPy_MPZ_Function_Kronecker(PyObject *self, PyObject *args)
{
    MPZ_Object *tempx = NULL, *tempy = NULL;
    long res;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("kronecker() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(tempx = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), NULL)) ||
        !(tempy = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), NULL))) {

        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        return NULL;
    }

    res = (long)(mpz_kronecker(tempx->z, tempy->z));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);
    return PyIntOrLong_FromLong(res);
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_even,
"is_even(x) -> bool\n\n"
"Return True if x is even, False otherwise.");

static PyObject *
GMPy_MPZ_Function_IsEven(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;

    if (MPZ_Check(other)) {
        res = mpz_even_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
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

PyDoc_STRVAR(GMPy_doc_mpz_method_is_even,
"x.is_even() -> bool\n\n"
"Return True if x is even, False otherwise.");

static PyObject *
GMPy_MPZ_Method_IsEven(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_even_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(GMPy_doc_mpz_function_is_odd,
"is_odd(x) -> bool\n\n"
"Return True if x is odd, False otherwise.");

static PyObject *
GMPy_MPZ_Function_IsOdd(PyObject *self, PyObject *other)
{
    int res;
    MPZ_Object *tempx;

    if (CHECK_MPZANY(other)) {
        res = mpz_odd_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
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

PyDoc_STRVAR(GMPy_doc_mpz_method_is_odd,
"x.is_odd() -> bool\n\n"
"Return True if x is odd, False otherwise.");

static PyObject *
GMPy_MPZ_Method_IsOdd(PyObject *self, PyObject *other)
{
    int res;

    res = mpz_odd_p(MPZ(self));

    if (res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/*
 * Add mapping support to mpz objects.
 */

static Py_ssize_t
GMPy_MPZ_Method_Length(MPZ_Object *self)
{
    return mpz_sizeinbase(self->z, 2);
}

static PyObject *
GMPy_MPZ_Method_SubScript(MPZ_Object *self, PyObject *item)
{
    if (PyIndex_Check(item)) {
        Py_ssize_t i;

        i = PyIntOrLong_AsSsize_t(item);
        if (i == -1 && PyErr_Occurred()) {
            INDEX_ERROR("argument too large to convert to an index");
            return NULL;
        }
        if (i < 0) {
            i += mpz_sizeinbase(self->z, 2);
        }
        return PyIntOrLong_FromLong(mpz_tstbit(self->z, i));
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step, slicelength, cur, i;
        MPZ_Object *result;

#if PY_VERSION_HEX > 0x030200A4
        if (PySlice_GetIndicesEx(item,
                        mpz_sizeinbase(self->z, 2),
                        &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }
#else
        if (PySlice_GetIndicesEx((PySliceObject*)item,
                        mpz_sizeinbase(self->z, 2),
                        &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }
#endif

        if ((step < 0 && start < stop) || (step > 0 && start > stop)) {
            stop = start;
        }

        if (!(result = GMPy_MPZ_New(NULL))) {
            return NULL;
        }

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
GMPy_MPZ_Attrib_GetNumer(MPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPZ_Attrib_GetReal(MPZ_Object *self, void *closure)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
GMPy_MPZ_Attrib_GetDenom(MPZ_Object *self, void *closure)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_set_ui(result->z, 1);
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Attrib_GetImag(MPZ_Object *self, void *closure)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(NULL))) {
        mpz_set_ui(result->z, 0);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(GMPy_doc_mpz_method_sizeof,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x. Note: deleted mpz objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
GMPy_MPZ_Method_SizeOf(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPZ_Object) + \
        (MPZ(self)->_mp_alloc * sizeof(mp_limb_t)));
}

/* Note: this particular function is also used for xmpz, mpq, and mpfr. Only
 * mpc.conjugate() does more that just return another reference to the original
 * object.
 */

PyDoc_STRVAR(GMPy_doc_mp_method_conjugate,
"x.conjugate() -> number\n\n"
"Return the conjugate of x (which is just a new reference to x since x is\n"
"not a complex number).");

static PyObject *
GMPy_MP_Method_Conjugate(PyObject *self, PyObject *args)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}



