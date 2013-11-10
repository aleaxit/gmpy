/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013 Case Van Horsen            *
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

    /* Optimize the most common use case */
    argc = PyTuple_Size(args);
    if (argc == 0) {
        if ((result = GMPy_MPZ_New())) {
            mpz_set_ui(result->z, 0);
        }
        return (PyObject*)result;
    }
    if (argc == 1) {
        n = PyTuple_GetItem(args, 0);
        if (IS_REAL(n) && !keywds) {
            /* _Temp can be used here since another reference to an existing
             * mpz can be returned.
             */
            result = GMPy_MPZ_From_Number_Temp(n);
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
        result = GMPy_MPZ_From_PyStr(n, base);
    }
    else {
        if (argc==2 || (argc == 1 && keywds))
            TYPE_ERROR("mpz() with non-string argument needs exactly "
                       "1 argument");
        else {
            result = GMPy_MPZ_From_Number_Temp(n);
            if (!result)
                TYPE_ERROR("mpz() requires numeric or string argument");
        }
    }
    return (PyObject*)result;
}

/* Functions that operate strictly on mpz or xmpz. */

/* produce digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(doc_mpz_digits,
"x.digits([base=10]) -> string\n\n"
"Return Python string representing x in the given base. Values for\n"
"base can range between 2 to 62. A leading '-' is present if x<0\n"
"but no leading '+' is present if x>=0.");

static PyObject *
Pympz_digits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "digits() requires 'int' argument for base");
    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        Py_DECREF(self);
        return NULL;
    }
    result = GMPy_PyStr_From_MPZ((MPZ_Object*)self, (int)base, 16);
    Py_DECREF(self);
    return result;
}

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

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "num_digits() requires 'mpz',['int'] arguments");
    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        Py_DECREF(self);
        return NULL;
    }
    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(MPZ(self),
                                    (int)base));
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_bit_lengthm,
"x.bit_length() -> int\n\n"
"Return the number of significant bits in the radix-2\n"
"representation of x. Note: bit_length(0) returns 0.");

PyDoc_STRVAR(doc_bit_lengthg,
"x.bit_length() -> int\n\n"
"Return the number of significant bits in the radix-2\n"
"representation of x. Note: mpz(0).bit_length() returns 0.");

static PyObject *
Pympz_bit_length(PyObject *self, PyObject *other)
{
    size_t i = 0;
    MPZ_Object* tempx;

    if (self && (CHECK_MPZANY(self))) {
        if (mpz_size(MPZ(self)))
            i = mpz_sizeinbase(MPZ(self), 2);
    }
    else if(CHECK_MPZANY(other)) {
        if (mpz_size(MPZ(other)))
            i = mpz_sizeinbase(MPZ(other), 2);
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
            TYPE_ERROR("bit_length() requires 'mpz' argument");
            return NULL;
        }
        else {
            if (mpz_size(MPZ(tempx)))
                i = mpz_sizeinbase(tempx->z, 2);
            Py_DECREF((PyObject*)tempx);
        }
    }
    return PyIntOrLong_FromSize_t(i);
}

PyDoc_STRVAR(doc_bit_maskg,
"bit_mask(n) -> mpz\n\n"
"Return an 'mpz' exactly n bits in length with all bits set.\n");

static PyObject *
Pympz_bit_mask(PyObject *self, PyObject *other)
{
    Py_ssize_t i = 0;
    MPZ_Object* result;

    i = ssize_t_From_Integer(other);

    if (i == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_mask() requires 'int' argument");
        return NULL;
    }

    if (i < 0) {
        VALUE_ERROR("mask length must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    mpz_set_ui(result->z, 1);
    mpz_mul_2exp(result->z, result->z, i);
    mpz_sub_ui(result->z, result->z, 1);

    return (PyObject*)result;
}

/* return scan0/scan1 for an mpz */
PyDoc_STRVAR(doc_bit_scan0m,
"x.bit_scan0(n=0) -> int\n\n"
"Return the index of the first 0-bit of x with index >= n. n >= 0.\n"
"If there are no more 0-bits in x at or above index n (which can\n"
"only happen for x<0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

PyDoc_STRVAR(doc_bit_scan0g,
"bit_scan0(x, n=0) -> int\n\n"
"Return the index of the first 0-bit of x with index >= n. n >= 0.\n"
"If there are no more 0-bits in x at or above index n (which can\n"
"only happen for x<0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

static PyObject *
Pympz_bit_scan0(PyObject *self, PyObject *args)
{
    Py_ssize_t maxbit, starting_bit = 0;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_SSIZE_T(&starting_bit,
            "bit_scan0() requires 'mpz',['int'] arguments");

    if (starting_bit < 0) {
        VALUE_ERROR("starting bit must be >= 0");
        Py_DECREF(self);
        return NULL;
    }
    maxbit = mpz_sizeinbase(MPZ(self), 2);
    if (starting_bit > maxbit) {
        if (mpz_sgn(MPZ(self))<0) {
            Py_DECREF(self);
            Py_RETURN_NONE;
        }
        else {
            result = PyIntOrLong_FromSsize_t(starting_bit);
        }
    }
    else {
        result = PyIntOrLong_FromSsize_t(mpz_scan0(MPZ(self), starting_bit));
    }
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_bit_scan1m,
"x.bit_scan1(n=0) -> int\n\n"
"Return the index of the first 1-bit of x with index >= n. n >= 0.\n"
"If there are no more 1-bits in x at or above index n (which can\n"
"only happen for x>=0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

PyDoc_STRVAR(doc_bit_scan1g,
"bit_scan1(x, n=0) -> int\n\n"
"Return the index of the first 1-bit of x with index >= n. n >= 0.\n"
"If there are no more 1-bits in x at or above index n (which can\n"
"only happen for x>=0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

static PyObject *
Pympz_bit_scan1(PyObject *self, PyObject *args)
{
    Py_ssize_t maxbit, starting_bit = 0;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_SSIZE_T(&starting_bit,
            "bit_scan1() requires 'mpz',['int'] arguments");

    if (starting_bit < 0) {
        VALUE_ERROR("starting bit must be >= 0");
        Py_DECREF(self);
        return NULL;
    }
    maxbit = mpz_sizeinbase(MPZ(self), 2);
    if (starting_bit >= maxbit) {
        if (mpz_sgn(MPZ(self))>=0) {
            Py_DECREF(self);
            Py_RETURN_NONE;
        }
        else {
            result = PyIntOrLong_FromSsize_t(starting_bit);
        }
    }
    else {
        result = PyIntOrLong_FromSsize_t(mpz_scan1(MPZ(self), starting_bit));
    }
    Py_DECREF(self);
    return result;
}

/* return population-count (# of 1-bits) for an mpz */

PyDoc_STRVAR(doc_popcountg,
"popcount(x) -> int\n\n"
"Return the number of 1-bits set in x. If x<0, the number of\n"
"1-bits is infinite so -1 is returned in that case.");

static PyObject *
Pympz_popcount(PyObject *self, PyObject *other)
{
    Py_ssize_t temp;
    MPZ_Object *tempx;

    if (self && (CHECK_MPZANY(self)))
        return PyIntOrLong_FromSsize_t(mpz_popcount(MPZ(self)));
    else if(CHECK_MPZANY(other))
        return PyIntOrLong_FromSsize_t(mpz_popcount(MPZ(other)));
    else {
        if ((tempx = GMPy_MPZ_From_Integer_Temp(other))) {
            temp = mpz_popcount(tempx->z);
            Py_DECREF((PyObject*)tempx);
            return PyIntOrLong_FromSsize_t(temp);
        }
        else {
            TYPE_ERROR("popcount() requires 'mpz' argument");
            return NULL;
        }
    }
}

/* get & return one bit from an mpz */
PyDoc_STRVAR(doc_bit_testg,
"bit_test(x, n) -> bool\n\n"
"Return the value of the n-th bit of x.");

static PyObject *
Pygmpy_bit_test(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    int temp;
    PyObject *x;
    MPZ_Object *tempx;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    if (CHECK_MPZANY(x)) {
        temp = mpz_tstbit(MPZ(x), bit_index);
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(x))) {
            TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
            return NULL;
        }
        temp = mpz_tstbit(tempx->z, bit_index);
        Py_DECREF((PyObject*)tempx);
    }
    if (temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_bit_testm,
"x.bit_test(n) -> bool\n\n"
"Return the value of the n-th bit of x.");

static PyObject *
Pympz_bit_test(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (mpz_tstbit(MPZ(self), bit_index))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyDoc_STRVAR(doc_bit_clearg,
"bit_clear(x, n) -> mpz\n\n"
"Return a copy of x with the n-th bit cleared.");

static PyObject *
Pygmpy_bit_clear(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    PyObject *x;
    MPZ_Object *result;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    if (CHECK_MPZANY(x)) {
        if (!(result = GMPy_MPZ_New()))
            return NULL;
        mpz_set(result->z, MPZ(x));
        mpz_clrbit(result->z, bit_index);
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer_Temp(x))) {
            TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_clrbit(result->z, bit_index);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_clearm,
"x.bit_clear(n) -> mpz\n\n"
"Return a copy of x with the n-th bit cleared.");

static PyObject *
Pympz_bit_clear(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;
    MPZ_Object *result;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;
    mpz_set(result->z, MPZ(self));
    mpz_clrbit(result->z, bit_index);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_setg,
"bit_set(x, n) -> mpz\n\n"
"Return a copy of x with the n-th bit set.");

static PyObject *
Pygmpy_bit_set(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    PyObject *x;
    MPZ_Object *result;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    if (CHECK_MPZANY(x)) {
        if (!(result = GMPy_MPZ_New()))
            return NULL;
        mpz_set(result->z, MPZ(x));
        mpz_setbit(result->z, bit_index);
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer_Temp(x))) {
            TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_setbit(result->z, bit_index);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_setm,
"x.bit_set(n) -> mpz\n\n"
"Return a copy of x with the n-th bit set.");

static PyObject *
Pympz_bit_set(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;
    MPZ_Object *result;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;
    mpz_set(result->z, MPZ(self));
    mpz_setbit(result->z, bit_index);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_flipg,
"bit_flip(x, n) -> mpz\n\n"
"Return a copy of x with the n-th bit inverted.");

static PyObject *
Pygmpy_bit_flip(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    PyObject *x;
    MPZ_Object *result;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 1));
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    if (CHECK_MPZANY(x)) {
        if (!(result = GMPy_MPZ_New()))
            return NULL;
        mpz_set(result->z, MPZ(x));
        mpz_combit(result->z, bit_index);
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer_Temp(x))) {
            TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_combit(result->z, bit_index);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_flipm,
"x.bit_flip(n) -> mpz\n\n"
"Return a copy of x with the n-th bit inverted.");

static PyObject *
Pympz_bit_flip(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;
    MPZ_Object *result;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;
    mpz_set(result->z, MPZ(self));
    mpz_combit(result->z, bit_index);
    return (PyObject*)result;
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

    PARSE_ONE_MPZ_REQ_SI(&n,
                         "iroot() requires 'mpz','int' arguments");

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
    if (!(s = GMPy_MPZ_New()) || !(result = PyTuple_New(2))) {
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
    y = GMPy_MPZ_New();
    r = GMPy_MPZ_New();
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

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_sgn(MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_sgn(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
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

static PyObject *
Pympz_neg(MPZ_Object *self)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New()))
        mpz_neg(result->z, self->z);

    return (PyObject*)result;
}

static PyObject *
Pympz_pos(MPZ_Object *self)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
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

    if ((result = GMPy_MPZ_New())) {
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

static PyObject *
Pympz_square(PyObject *self, PyObject *other)
{
    MPZ_Object *tempx, *result;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    if (self && (CHECK_MPZANY(self))) {
        mpz_mul(result->z, MPZ(self), MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        mpz_mul(result->z, MPZ(other), MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
            TYPE_ERROR("square() requires 'mpz' argument");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        else {
            mpz_mul(result->z, MPZ(tempx), MPZ(tempx));
            Py_DECREF((PyObject*)tempx);
        }
    }
    return (PyObject*)result;
}

static int
Pympz_nonzero(MPZ_Object *self)
{
    return mpz_sgn(self->z) != 0;
}

/* BIT OPERATIONS */

static PyObject *
Pympz_com(MPZ_Object *self)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New()))
        mpz_com(result->z, MPZ(self));

    return (PyObject*)result;
}

#define MPZ_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *self, PyObject *other) \
{ \
    MPZ_Object *result = 0; \
    if (CHECK_MPZANY(self)) { \
        if (CHECK_MPZANY(other)) { \
            if (!(result = GMPy_MPZ_New())) \
                return NULL; \
            NAME(result->z, MPZ(self), MPZ(other)); \
        } \
        else { \
            if (!(result = GMPy_MPZ_From_Integer_Temp(other))) \
                return NULL; \
            NAME(result->z, MPZ(self), result->z); \
        } \
    } \
    else if (CHECK_MPZANY(other)) { \
        if (!(result = GMPy_MPZ_From_Integer_Temp(self))) \
            return NULL; \
        NAME(result->z, result->z, MPZ(other)); \
    } \
    else { \
        Py_RETURN_NOTIMPLEMENTED; \
    } \
    return (PyObject*)result; \
}

MPZ_BINOP(mpz_and)
MPZ_BINOP(mpz_ior)
MPZ_BINOP(mpz_xor)

static PyObject *
Pympz_rshift(PyObject *self, PyObject *other)
{
    mpir_si count_si;
    int overflow;
    MPZ_Object *result, *tempa, *tempb;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    /* Try to make mpz >> Python int/long as fast as possible. */
    if (CHECK_MPZANY(self)) {
        if (PyIntOrLong_Check(other)) {
            count_si = PyLong_AsSIAndOverflow(other, &overflow);
            if (overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else if (count_si >= 0) {
                mpz_fdiv_q_2exp(result->z, MPZ(self), count_si);
                return (PyObject*)result;
            }
            else {
                VALUE_ERROR("negative shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }
    }

    tempa = GMPy_MPZ_From_Integer_Temp(self);
    tempb = GMPy_MPZ_From_Integer_Temp(other);
    if (!tempb || !tempa) {
        TYPE_ERROR("Pympz_rshift() expects integer arguments");
        goto err;
    }
    if (mpz_sgn(MPZ(tempb)) < 0) {
        VALUE_ERROR("negative shift count");
        goto err;
    }
    if(!mpz_fits_si_p(MPZ(tempb))) {
        OVERFLOW_ERROR("outrageous shift count");
        goto err;
    }
    count_si = mpz_get_si(tempb->z);
    mpz_fdiv_q_2exp(result->z, tempa->z, count_si);
    Py_DECREF((PyObject*)tempa);
    Py_DECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_DECREF((PyObject*)result);
    Py_DECREF((PyObject*)tempa);
    Py_DECREF((PyObject*)tempb);
    return NULL;
}

static PyObject *
Pympz_lshift(PyObject *self, PyObject *other)
{
    mpir_si count_si;
    int overflow;
    MPZ_Object *result, *tempa, *tempb;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    /* Try to make mpz >> Python int/long as fast as possible. */
    if (CHECK_MPZANY(self)) {
        if (PyIntOrLong_Check(other)) {
            count_si = PyLong_AsSIAndOverflow(other, &overflow);
            if (overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else if (count_si >= 0) {
                mpz_mul_2exp(result->z, MPZ(self), count_si);
                return (PyObject*)result;
            }
            else {
                VALUE_ERROR("negative shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }
    }

    tempa = GMPy_MPZ_From_Integer_Temp(self);
    tempb = GMPy_MPZ_From_Integer_Temp(other);
    if (!tempb || !tempa) {
        TYPE_ERROR("Pympz_lshift() expects integer arguments");
        goto err;
        }
    if (mpz_sgn(MPZ(tempb)) < 0) {
        VALUE_ERROR("negative shift count");
        goto err;
    }
    if(!mpz_fits_si_p(MPZ(tempb))) {
        OVERFLOW_ERROR("outrageous shift count");
        goto err;
    }
    count_si = mpz_get_si(tempb->z);
    mpz_mul_2exp(result->z, tempa->z, count_si);
    Py_DECREF((PyObject*)tempa);
    Py_DECREF((PyObject*)tempb);
    return (PyObject*)result;

  err:
    Py_DECREF((PyObject*)result);
    Py_DECREF((PyObject*)tempa);
    Py_DECREF((PyObject*)tempb);
    return NULL;
}

#if PY_MAJOR_VERSION < 3
/* hex/oct formatting (mpz-only) */
static PyObject *
Pympz_oct(MPZ_Object *self)
{
    return GMPy_PyStr_From_MPZ(self, 8, 0);
}

static PyObject *
Pympz_hex(MPZ_Object *self)
{
    return GMPy_PyStr_From_MPZ(self, 16, 0);
}
#endif

static Py_hash_t
Pympz_hash(MPZ_Object *self)
{
#ifdef _PyHASH_MODULUS
    Py_hash_t hash;

    if (self->hash_cache != -1)
        return self->hash_cache;

    hash = (Py_hash_t)mpn_mod_1(self->z->_mp_d, mpz_size(self->z), _PyHASH_MODULUS);
    if (mpz_sgn(self->z)<0)
        hash = -hash;
    if (hash == -1)
        hash = -2;
    return (self->hash_cache = hash);
#else
    if (self->hash_cache != -1)
        return self->hash_cache;
    else
        return (self->hash_cache = mpz_pythonhash(self->z));
#endif
}

/* Miscellaneous gmpy functions */
PyDoc_STRVAR(doc_gcd,
"gcd(a, b) -> mpz\n\n"
"Return the greatest common denominator of integers a and b.");

static PyObject *
Pygmpy_gcd(PyObject *self, PyObject *args)
{
    PyObject *a, *b;
    MPZ_Object *result, *tempa, *tempb;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcd() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_gcd(result->z, MPZ(a), MPZ(b));
    }
    else {
        tempa = GMPy_MPZ_From_Integer_Temp(a);
        tempb = GMPy_MPZ_From_Integer_Temp(b);
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

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("lcm() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_lcm(result->z, MPZ(a), MPZ(b));
    }
    else {
        tempa = GMPy_MPZ_From_Integer_Temp(a);
        tempb = GMPy_MPZ_From_Integer_Temp(b);
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

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
        return NULL;
    }

    g = GMPy_MPZ_New();
    s = GMPy_MPZ_New();
    t = GMPy_MPZ_New();
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
        tempa = GMPy_MPZ_From_Integer_Temp(a);
        tempb = GMPy_MPZ_From_Integer_Temp(b);
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

    if(PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    num = GMPy_MPZ_From_Integer_Temp(PyTuple_GET_ITEM(args, 0));
    den = GMPy_MPZ_From_Integer_Temp(PyTuple_GET_ITEM(args, 1));
    mod = GMPy_MPZ_From_Integer_Temp(PyTuple_GET_ITEM(args, 2));

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
        if (!(result = GMPy_MPZ_New()))
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
        if (!(result = GMPy_MPZ_New()))
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
        if (!(result = GMPy_MPZ_New()))
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
        mpz_fib2_ui(luc1->z, luc2->z, n);
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

    PARSE_ONE_MPZ_REQ_SI(&k, "bincoef() requires 'mpz','int' arguments");

    if (k < 0) {
        VALUE_ERROR("binomial coefficient with negative k");
        Py_DECREF(self);
        return NULL;
    }

    if(!(result = GMPy_MPZ_New())) {
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

    if (self && (CHECK_MPZANY(self))) {
        if (mpz_sgn(MPZ(self)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if (!(result = GMPy_MPZ_New()))
            return NULL;
        mpz_sqrt(result->z, MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        if (mpz_sgn(MPZ(other)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if (!(result = GMPy_MPZ_New()))
            return NULL;
        mpz_sqrt(result->z, MPZ(other));
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer_Temp(other))) {
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

    PARSE_ONE_MPZ("isqrt_rem() requires 'mpz' argument");

    if (mpz_sgn(MPZ(self)) < 0) {
        VALUE_ERROR("isqrt_rem() of negative number");
        Py_DECREF(self);
        return NULL;
    }

    root = GMPy_MPZ_New();
    rem = GMPy_MPZ_New();
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

    PARSE_TWO_MPZ(factor, "remove() requires 'mpz','mpz' arguments");

    if (mpz_cmp_si(MPZ(factor), 2) < 0) {
        VALUE_ERROR("factor must be > 1");
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }

    if (!(result = GMPy_MPZ_New())) {
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

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
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
        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
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

PyDoc_STRVAR(doc_hamdistg,
"hamdist(x, y) -> int\n\n"
"Return the Hamming distance (number of bit-positions where the\n"
"bits differ) between integers x and y.");

static PyObject *
Pympz_hamdist(PyObject *self, PyObject *args)
{
    PyObject *result, *other;

    PARSE_TWO_MPZ(other, "hamdist() requires 'mpz','mpz' arguments");

    result = PyIntOrLong_FromSize_t(
            mpz_hamdist(MPZ(self),MPZ(other)));
    Py_DECREF(self);
    Py_DECREF(other);
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

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New()))
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
        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
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

    if (CHECK_MPZANY(other)) {
        res = mpz_perfect_square_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
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

    if (CHECK_MPZANY(other)) {
        res = mpz_perfect_power_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
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

    if(CHECK_MPZANY(other)) {
        if(!(result = GMPy_MPZ_New()))
            return NULL;
        mpz_nextprime(result->z, MPZ(other));
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer_Temp(other))) {
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

    if (CHECK_MPZANY(other)) {
        res = mpz_even_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
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

    if (CHECK_MPZANY(other)) {
        res = mpz_odd_p(MPZ(other));
    }
    else {
        if (!(tempx = GMPy_MPZ_From_Integer_Temp(other))) {
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

        if (!(result = GMPy_MPZ_New()))
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

PyDoc_STRVAR(doc_mpz_format,
"x.__format__(fmt) -> string\n\n"
"Return a Python string by formatting mpz 'x' using the format string\n"
"'fmt'. A valid format string consists of:\n"
"     optional alignment code:\n"
"        '<' -> left shifted in field\n"
"        '>' -> right shifted in field\n"
"        '^' -> centered in field\n"
"     optional leading sign code:\n"
"        '+' -> always display leading sign\n"
"        '-' -> only display minus sign\n"
"        ' ' -> minus for negative values, space for positive values\n"
"     optional base indicator\n"
"        '#' -> precede binary, octal, or hex with 0b, 0o or 0x\n"
"     optional width\n"
"     optional conversion code:\n"
"        'd' -> decimal format\n"
"        'b' -> binary format\n"
"        'o' -> octal format\n"
"        'x' -> hex format\n"
"The default format is 'd'.");

/* Formatting occurs in two phases. Pympz_ascii() is used to create a string
 * with the appropriate binary/octal/decimal/hex formatting, including the
 * leading sign character (+ , -, or space) and base encoding (0b, 0o, or 0x).
 * Left/right/centering using the specified width is done by creating a
 * format string and calling the __format__() method of the string object
 * returned by Pympz_ascii().
 */

static PyObject *
Pympz_format(PyObject *self, PyObject *args)
{
    PyObject *result = 0, *mpzstr = 0;
    char *fmtcode = 0, *p1, *p2;
    char fmt[30];
    int base = 10, option = 16;
    int seensign = 0, seenindicator = 0, seenalign = 0, seendigits = 0;

    if (!CHECK_MPZANY(self)) {
        TYPE_ERROR("requires mpz type");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "s", &fmtcode))
        return NULL;

    p2 = fmt;
    for (p1 = fmtcode; *p1 != '\00'; p1++) {
        if (*p1 == '<' || *p1 == '>' || *p1 == '^') {
            if (seenalign || seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seenalign = 1;
                continue;
            }
        }
        if (*p1 == '+') {
            if (seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                option |= 2;
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '-') {
            if (seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                seensign = 1;
                continue;
            }
        }
        if (*p1 == ' ') {
            if (seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                option |= 4;
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '#') {
            if (seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                option |= 8;
                seenindicator = 1;
                continue;
            }
        }
        if (isdigit(*p1)) {
            if (!seenalign) {
                *(p2++) = '>';
                seenalign = 1;
            }
            *(p2++) = *p1;
            seendigits = 1;
            continue;
        }
        if (*p1 == 'b') {
            base = 2;
            break;
        }
        if (*p1 == 'o') {
            base = 8;
            break;
        }
        if (*p1 == 'x') {
            base = 16;
            break;
        }
        if (*p1 == 'd') {
            base = 10;
            break;
        }
        if (*p1 == 'X') {
            base = -16;
            break;
        }
        VALUE_ERROR("Invalid conversion specification");
        return NULL;
    }
    *(p2++) = '\00';

    if (!(mpzstr = mpz_ascii(MPZ(self), base, option)))
        return NULL;

    result = PyObject_CallMethod(mpzstr, "__format__", "(s)", fmt);
    Py_DECREF(mpzstr);
    return result;
}

/* Divide two Integer objects (see convert.c/isInteger) using floor (//)
 * division. If an error occurs, NULL is returned and an exception is set.
 * If either x or y can't be converted into an mpz, Py_NotImplemented is
 * returned. */

static PyObject *
Pympz_FloorDiv_Integer(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPZ_Object *result;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_fdiv_q(result->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_q_ui(result->z, MPZ(x), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else {
                mpz_cdiv_q_ui(result->z, MPZ(x), -temp_si);
                mpz_neg(result->z, result->z);
            }
            return (PyObject*)result;
        }

        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpz_fdiv_q(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
    }

    if (CHECK_MPZANY(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (PyIntOrLong_Check(x)) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, x);
            mpz_fdiv_q(result->z, tempz, MPZ(y));
            mpz_cloc(tempz);
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        MPZ_Object *tempx, *tempy;

        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Integer to mpz.");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }

        mpz_fdiv_q(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

/* Implement floordiv for Pympz. On entry, one of the two arguments must
 * be a Pympz. If the other object is an Integer, add and return a Pympz.
 * If the other object isn't a Pympz, call the appropriate function. If
 * no appropriate function can be found, return NotImplemented. */

static PyObject *
Pympz_floordiv_fast(PyObject *x, PyObject *y)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return Pympz_FloorDiv_Integer(x, y, context);
    else if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return Pympq_FloorDiv_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        return Pympfr_FloorDiv_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return Pympc_FloorDiv_Complex(x, y, context);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Divide two Integer objects (see convert.c/isInteger) using true division.
 * If an error occurs, NULL is returned and an exception is set. If either x
 * or y can't be converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
Pympz_TrueDiv_Integer(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPZ_Object *tempx, *tempy;
    mpq_t tempq;
    MPFR_Object *result;

    if (!(result = (MPFR_Object*)Pympfr_new_context(context)))
        return NULL;

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Integer to mpz.");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpq_init(tempq);
        mpq_set_num(tempq, tempx->z);
        mpq_set_den(tempq, tempy->z);
        mpq_canonicalize(tempq);
        mpfr_clear_flags();
        result->rc = mpfr_set_q(result->f, tempq, GET_MPFR_ROUND(context));
        mpq_clear(tempq);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        MPFR_CLEANUP_RESULT("division");
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_truediv_fast(PyObject *x, PyObject *y)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return Pympz_TrueDiv_Integer(x, y, context);
    else if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return Pympq_TrueDiv_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        return Pympfr_TrueDiv_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return Pympc_TrueDiv_Complex(x, y, context);

    Py_RETURN_NOTIMPLEMENTED;
}

#ifdef PY2
static PyObject *
Pympz_div2_fast(PyObject *x, PyObject *y)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return Pympz_FloorDiv_Integer(x, y, context);
    else if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return Pympq_TrueDiv_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        return Pympfr_TrueDiv_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return Pympc_TrueDiv_Complex(x, y, context);

    Py_RETURN_NOTIMPLEMENTED;
}
#endif

/* Divide two Integer objects (see convert.c/isInteger) and return remainder.
 * If an error occurs, NULL is returned and an exception is set. If either x
 * or y can't be converted into an mpz, Py_NotImplemented is returned. */

static PyObject *
Pympz_Mod_Integer(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    MPZ_Object *tempx, *tempy, *result;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    if (!(result = GMPy_MPZ_New()))
        return NULL;

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_fdiv_r(result->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_r_ui(result->z, MPZ(x), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else {
                mpz_cdiv_r_ui(result->z, MPZ(x), -temp_si);
            }
            return (PyObject*)result;
        }
        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            mpz_fdiv_r(result->z, MPZ(x), MPZ(y));
            return (PyObject*)result;
        }
    }

    if (CHECK_MPZANY(y)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (PyIntOrLong_Check(x)) {
            mpz_inoc(tempz);
            mpz_set_PyIntOrLong(tempz, x);
            mpz_fdiv_r(result->z, tempz, MPZ(y));
            mpz_cloc(tempz);
            return (PyObject*)result;
        }
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Integer to mpz.");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_fdiv_r(result->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_mod_fast(PyObject *x, PyObject *y)
{
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        return Pympz_Mod_Integer(x, y, context);
    else if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return Pympq_Mod_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        return Pympfr_Mod_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return Pympc_Mod_Complex(x, y, context);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Divide two Integer objects (see convert.c/isInteger) and return quotient
 * and remainder. If an error occurs, NULL is returned and an exception is set.
 * If either x or y can't be converted into an mpz, Py_NotImplemented is
 * returned. */

static PyObject *
Pympz_DivMod_Integer(PyObject *x, PyObject *y, GMPyContextObject *context)
{
    PyObject *result;
    MPZ_Object *tempx, *tempy, *rem, *quo;
    mpz_t tempz;
    mpir_si temp_si;
    int overflow;

    result = PyTuple_New(2);
    rem = GMPy_MPZ_New();
    quo = GMPy_MPZ_New();
    if (!result || !rem || !quo) {
        Py_XDECREF((PyObject*)rem);
        Py_XDECREF((PyObject*)quo);
        Py_XDECREF(result);
        return NULL;
    }

    if (CHECK_MPZANY(x)) {
        if (PyIntOrLong_Check(y)) {
            temp_si = PyLong_AsSIAndOverflow(y, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, y);
                mpz_fdiv_qr(quo->z, rem->z, MPZ(x), tempz);
                mpz_cloc(tempz);
            }
            else if (temp_si > 0) {
                mpz_fdiv_qr_ui(quo->z, rem->z, MPZ(x), temp_si);
            }
            else if (temp_si == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rem);
                Py_DECREF((PyObject*)quo);
                Py_DECREF(result);
                return NULL;
            }
            else {
                mpz_cdiv_qr_ui(quo->z, rem->z, MPZ(x), -temp_si);
                mpz_neg(quo->z, quo->z);
            }
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }
        if (CHECK_MPZANY(y)) {
            if (mpz_sgn(MPZ(y)) == 0) {
                ZERO_ERROR("division or modulo by zero");
                Py_DECREF((PyObject*)rem);
                Py_DECREF((PyObject*)quo);
                Py_DECREF(result);
                return NULL;
            }
            mpz_fdiv_qr(quo->z, rem->z, MPZ(x), MPZ(y));
            PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
            PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
            return result;
        }
    }

    if (CHECK_MPZANY(y) && PyIntOrLong_Check(x)) {
        if (mpz_sgn(MPZ(y)) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_DECREF((PyObject*)rem);
            Py_DECREF((PyObject*)quo);
            Py_DECREF(result);
            return NULL;
        }
        mpz_inoc(tempz);
        mpz_set_PyIntOrLong(tempz, x);
        mpz_fdiv_qr(quo->z, rem->z, tempz, MPZ(y));
        mpz_cloc(tempz);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return (PyObject*)result;
    }

    if (IS_INTEGER(x) && IS_INTEGER(y)) {
        tempx = GMPy_MPZ_From_Integer_Temp(x);
        tempy = GMPy_MPZ_From_Integer_Temp(y);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Integer to mpz.");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)rem);
            Py_DECREF((PyObject*)quo);
            Py_DECREF(result);
            return NULL;
        }
        if (mpz_sgn(tempy->z) == 0) {
            ZERO_ERROR("division or modulo by zero");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)rem);
            Py_DECREF((PyObject*)quo);
            Py_DECREF(result);
            return NULL;
        }
        mpz_fdiv_qr(quo->z, rem->z, tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_divmod_fast(PyObject *x, PyObject *y)
{
    PyObject *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (IS_INTEGER(x) && IS_INTEGER(y))
        result = Pympz_DivMod_Integer(x, y, context);
    else if (IS_RATIONAL(x) && IS_RATIONAL(y))
        result = Pympq_DivMod_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        result = Pympfr_DivMod_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        result = Pympc_DivMod_Complex(x, y, context);
    else {
        Py_INCREF(Py_NotImplemented);
        result = Py_NotImplemented;
    }
    return result;
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

    if ((result = GMPy_MPZ_New()))
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
    (binaryfunc) GMPy_mpz_add_fast,        /* nb_add                  */
    (binaryfunc) GMPy_mpz_sub_fast,        /* nb_subtract             */
    (binaryfunc) GMPy_mpz_mul_fast,        /* nb_multiply             */
    (binaryfunc) Pympz_mod_fast,           /* nb_remainder            */
    (binaryfunc) Pympz_divmod_fast,        /* nb_divmod               */
    (ternaryfunc) GMPy_mpany_pow_fast,     /* nb_power                */
    (unaryfunc) Pympz_neg,                 /* nb_negative             */
    (unaryfunc) Pympz_pos,                 /* nb_positive             */
    (unaryfunc) GMPy_mpz_abs_fast,         /* nb_absolute             */
    (inquiry) Pympz_nonzero,               /* nb_bool                 */
    (unaryfunc) Pympz_com,                 /* nb_invert               */
    (binaryfunc) Pympz_lshift,             /* nb_lshift               */
    (binaryfunc) Pympz_rshift,             /* nb_rshift               */
    (binaryfunc) Pympz_and,                /* nb_and                  */
    (binaryfunc) Pympz_xor,                /* nb_xor                  */
    (binaryfunc) Pympz_ior,                /* nb_or                   */
    (unaryfunc) GMPy_PyLong_From_MPZ,      /* nb_int                  */
        0,                                 /* nb_reserved             */
    (unaryfunc) GMPy_PyFloat_From_MPZ,     /* nb_float                */
    (binaryfunc) Pympz_inplace_add,        /* nb_inplace_add          */
    (binaryfunc) Pympz_inplace_sub,        /* nb_inplace_subtract     */
    (binaryfunc) Pympz_inplace_mul,        /* nb_inplace_multiply     */
    (binaryfunc) Pympz_inplace_rem,        /* nb_inplace_remainder    */
    (ternaryfunc) Pympz_inplace_pow,       /* nb_inplace_power        */
    (binaryfunc) Pympz_inplace_lshift,     /* nb_inplace_lshift       */
    (binaryfunc) Pympz_inplace_rshift,     /* nb_inplace_rshift       */
        0,                                 /* nb_inplace_and          */
        0,                                 /* nb_inplace_xor          */
        0,                                 /* nb_inplace_or           */
    (binaryfunc) Pympz_floordiv_fast,      /* nb_floor_divide         */
    (binaryfunc) Pympz_truediv_fast,       /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv,   /* nb_inplace_floor_divide */
        0,                                 /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_PyIntOrLong_From_MPZ, /* nb_index                */
};

#else
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) GMPy_mpz_add_fast,        /* nb_add                  */
    (binaryfunc) GMPy_mpz_sub_fast,        /* nb_subtract             */
    (binaryfunc) GMPy_mpz_mul_fast,        /* nb_multiply             */
    (binaryfunc) Pympz_div2_fast,          /* nb_divide               */
    (binaryfunc) Pympz_mod_fast,           /* nb_remainder            */
    (binaryfunc) Pympz_divmod_fast,        /* nb_divmod               */
    (ternaryfunc) GMPy_mpany_pow_fast,     /* nb_power                */
    (unaryfunc) Pympz_neg,                 /* nb_negative             */
    (unaryfunc) Pympz_pos,                 /* nb_positive             */
    (unaryfunc) GMPy_mpz_abs_fast,         /* nb_absolute             */
    (inquiry) Pympz_nonzero,               /* nb_bool                 */
    (unaryfunc) Pympz_com,                 /* nb_invert               */
    (binaryfunc) Pympz_lshift,             /* nb_lshift               */
    (binaryfunc) Pympz_rshift,             /* nb_rshift               */
    (binaryfunc) Pympz_and,                /* nb_and                  */
    (binaryfunc) Pympz_xor,                /* nb_xor                  */
    (binaryfunc) Pympz_ior,                /* nb_or                   */
        0,                                 /* nb_coerce               */
    (unaryfunc) GMPy_PyIntOrLong_From_MPZ, /* nb_int                  */
    (unaryfunc) GMPy_PyLong_From_MPZ,      /* nb_long                 */
    (unaryfunc) GMPy_PyFloat_From_MPZ,     /* nb_float                */
    (unaryfunc) Pympz_oct,                 /* nb_oct                  */
    (unaryfunc) Pympz_hex,                 /* nb_hex                  */
    (binaryfunc) Pympz_inplace_add,        /* nb_inplace_add          */
    (binaryfunc) Pympz_inplace_sub,        /* nb_inplace_subtract     */
    (binaryfunc) Pympz_inplace_mul,        /* nb_inplace_multiply     */
        0,                                 /* nb_inplace_divide       */
    (binaryfunc) Pympz_inplace_rem,        /* nb_inplace_remainder    */
    (ternaryfunc) Pympz_inplace_pow,       /* nb_inplace_power        */
    (binaryfunc) Pympz_inplace_lshift,     /* nb_inplace_lshift       */
    (binaryfunc) Pympz_inplace_rshift,     /* nb_inplace_rshift       */
        0,                                 /* nb_inplace_and          */
        0,                                 /* nb_inplace_xor          */
        0,                                 /* nb_inplace_or           */
    (binaryfunc) Pympz_floordiv_fast,      /* nb_floor_divide         */
    (binaryfunc) Pympz_truediv_fast,       /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv,   /* nb_inplace_floor_divide */
        0,                                 /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_PyIntOrLong_From_MPZ, /* nb_index                */
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
    { "__format__", Pympz_format, METH_VARARGS, doc_mpz_format },
    { "__ceil__", Pympz_ceil, METH_NOARGS, doc_mpz_ceil },
    { "__floor__", Pympz_floor, METH_NOARGS, doc_mpz_floor },
    { "__round__", Pympz_round, METH_VARARGS, doc_mpz_round },
    { "__sizeof__", Pympz_sizeof, METH_NOARGS, doc_mpz_sizeof },
    { "__trunc__", Pympz_trunc, METH_NOARGS, doc_mpz_trunc },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "digits", Pympz_digits, METH_VARARGS, doc_mpz_digits },
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
    (hashfunc) Pympz_hash,                  /* tp_hash          */
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

