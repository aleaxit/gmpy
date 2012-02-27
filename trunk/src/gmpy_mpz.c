/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen                  *
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
    PympzObject *result = 0;
    PyObject *n = 0;
    long base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"n", "base", NULL };

    TRACE("Pygmpy_mpz() called...\n");

    /* Optimize the most common use case */
    argc = PyTuple_Size(args);
    if (argc == 1) {
        n = PyTuple_GetItem(args, 0);
#ifdef WITHMPFR
        if (isReal(n) && !keywds) {
#else
        if ((isRational(n) || PyFloat_Check(n)) && !keywds) {
#endif
            result = anynum2Pympz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("mpz() requires numeric or string argument");
            return (PyObject*)result;
        }
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist,
                                     &n, &base))
        return NULL;

    if ((base!=0) && ((base<2)||(base>62))) {
        VALUE_ERROR("base for mpz() must be 0 or in the "
                    "interval 2 ... 62");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        /* build-from-string (ascii or unicode) */
        result = PyStr2Pympz(n, base);
    }
    else {
        if (argc==2 || (argc == 1 && keywds))
            TYPE_ERROR("mpz() with non-string argument needs exactly "
                       "1 argument");
        else {
            result = anynum2Pympz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("mpz() requires numeric or string argument");
        }
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_xmpz,
"xmpz(n) -> xmpz\n\n"
"     Return an 'xmpz' object with a numeric value 'n' (truncating n\n"
"     to its integer part if it's a Fraction, 'mpq', Decimal, float\n"
"     or 'mpfr').\n\n"
"xmpz(s[, base=0]):\n\n"
"     Return an 'xmpz' object from a string 's' made of digits in the\n"
"     given base.  If base=0, binary, octal, or hex Python strings\n"
"     are recognized by leading 0b, 0o, or 0x characters, otherwise\n"
"     the string is assumed to be decimal. Values for base can range\n"
"     between 2 and 62.\n\n"
"     Note: 'xmpz' is a mutable integer. It can be faster for when\n"
"     used for augmented assignment (+=, *=, etc.). 'xmpz' objects\n"
"     cannot be used as dictionary keys. The use of 'mpz' objects is\n"
"     recommended in most cases.");

static PyObject *
Pygmpy_xmpz(PyObject *self, PyObject *args, PyObject *keywds)
{
    PyxmpzObject *result = 0;
    PyObject *n = 0;
    long base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"n", "base", NULL };

    TRACE("Pygmpy_xmpz() called...\n");

    /* Optimize the most common use case */
    argc = PyTuple_Size(args);
    if (argc == 1) {
        n = PyTuple_GetItem(args, 0);
#ifdef WITHMPFR
        if (isReal(n) && !keywds) {
#else
        if ((isRational(n) || PyFloat_Check(n)) && !keywds) {
#endif
            result = anynum2Pyxmpz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("xmpz() requires numeric or string argument");
            return (PyObject*)result;
        }
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist,
                                     &n, &base))
        return NULL;

    if ((base!=0) && ((base<2)||(base>62))) {
        VALUE_ERROR("base for xmpz() must be 0 or in the "
                    "interval 2 ... 62");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        /* build-from-string (ascii or unicode) */
        result = PyStr2Pyxmpz(n, base);
    }
    else {
        if (argc==2 || (argc == 1 && keywds))
            TYPE_ERROR("xmpz() with non-string argument needs exactly "
                       "1 argument");
        else {
            result = anynum2Pyxmpz(n);
            if (!result && !PyErr_Occurred())
                TYPE_ERROR("xmpz() requires numeric or string argument");
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
    result = Pympz_ascii((PympzObject*)self, base, 16);
    Py_DECREF(self);
    return result;
}

static PyObject *
Pyxmpz_digits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "digits() requires 'int' argument for base");
    result = Pyxmpz_ascii((PyxmpzObject*)self, base, 0);
    Py_DECREF(self);
    return result;
}

/* return number-of-digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(doc_numdigitsm,
"x.numdigits([base=10]) -> int\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

PyDoc_STRVAR(doc_numdigitsg,
"numdigits(x[, base=10]) -> int\n\n"
"Return length of string representing the absolute value of x in\n"
"the given base. Values  for base can range between 2 and 62. The\n"
"value returned may be 1 too large.");

static PyObject *
Pympz_numdigits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "numdigits() requires 'mpz',['int'] arguments");
    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        Py_DECREF(self);
        return NULL;
    }
    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(Pympz_AS_MPZ(self),
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
    PympzObject* tempx;

    if (self && (CHECK_MPZANY(self))) {
        if (mpz_size(Pympz_AS_MPZ(self)))
            i = mpz_sizeinbase(Pympz_AS_MPZ(self), 2);
    }
    else if(CHECK_MPZANY(other)) {
        if (mpz_size(Pympz_AS_MPZ(other)))
            i = mpz_sizeinbase(Pympz_AS_MPZ(other), 2);
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("bit_length() requires 'mpz' argument");
            return NULL;
        }
        else {
            if (mpz_size(Pympz_AS_MPZ(tempx)))
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
    PympzObject* result;

    i = ssize_t_From_Integer(other);

    if (i == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_mask() requires 'int' argument");
        return NULL;
    }

    if (i < 0) {
        VALUE_ERROR("mask length must be >= 0");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;

    mpz_set_ui(result->z, 1);
    mpz_mul_2exp(result->z, result->z, i);
    mpz_sub_ui(result->z, result->z, 1);

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_xbit_maskg,
"xbit_mask(n) -> xmpz\n\n"
"Return an 'xmpz' exactly n bits in length with all bits set.\n");

static PyObject *
Pyxmpz_xbit_mask(PyObject *self, PyObject *other)
{
    Py_ssize_t i = 0;
    PyxmpzObject* result;

    i = ssize_t_From_Integer(other);
    if (i == -1 && PyErr_Occurred()) {
        TYPE_ERROR("xbit_mask() requires 'int' argument");
        return NULL;
    }

    if (i < 0) {
        VALUE_ERROR("mask length must be >= 0");
        return NULL;
    }

    if (!(result = Pyxmpz_new()))
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
    maxbit = mpz_sizeinbase(Pympz_AS_MPZ(self), 2);
    if (starting_bit > maxbit) {
        if (mpz_sgn(Pympz_AS_MPZ(self))<0) {
            Py_DECREF(self);
            Py_RETURN_NONE;
        }
        else {
            result = PyIntOrLong_FromSsize_t(starting_bit);
        }
    }
    else {
        result = PyIntOrLong_FromSsize_t(mpz_scan0(Pympz_AS_MPZ(self), starting_bit));
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
    maxbit = mpz_sizeinbase(Pympz_AS_MPZ(self), 2);
    if (starting_bit >= maxbit) {
        if (mpz_sgn(Pympz_AS_MPZ(self))>=0) {
            Py_DECREF(self);
            Py_RETURN_NONE;
        }
        else {
            result = PyIntOrLong_FromSsize_t(starting_bit);
        }
    }
    else {
        result = PyIntOrLong_FromSsize_t(mpz_scan1(Pympz_AS_MPZ(self), starting_bit));
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
    PympzObject *tempx;

    if (self && (CHECK_MPZANY(self)))
        return PyIntOrLong_FromSsize_t(mpz_popcount(Pympz_AS_MPZ(self)));
    else if(CHECK_MPZANY(other))
        return PyIntOrLong_FromSsize_t(mpz_popcount(Pympz_AS_MPZ(other)));
    else {
        if ((tempx = Pympz_From_Integer(other))) {
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
    PympzObject *tempx;

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
        temp = mpz_tstbit(Pympz_AS_MPZ(x), bit_index);
    }
    else {
        if (!(tempx = Pympz_From_Integer(x))) {
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

    if (mpz_tstbit(Pympz_AS_MPZ(self), bit_index))
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
    PympzObject *result;

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
        if (!(result = Pympz_new()))
            return NULL;
        mpz_set(result->z, Pympz_AS_MPZ(x));
        mpz_clrbit(result->z, bit_index);
    }
    else {
        if (!(result = Pympz_From_Integer(x))) {
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
    PympzObject *result;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;
    mpz_set(result->z, Pympz_AS_MPZ(self));
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
    PympzObject *result;

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
        if (!(result = Pympz_new()))
            return NULL;
        mpz_set(result->z, Pympz_AS_MPZ(x));
        mpz_setbit(result->z, bit_index);
    }
    else {
        if (!(result = Pympz_From_Integer(x))) {
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
    PympzObject *result;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;
    mpz_set(result->z, Pympz_AS_MPZ(self));
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
    PympzObject *result;

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
        if (!(result = Pympz_new()))
            return NULL;
        mpz_set(result->z, Pympz_AS_MPZ(x));
        mpz_combit(result->z, bit_index);
    }
    else {
        if (!(result = Pympz_From_Integer(x))) {
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
    PympzObject *result;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;
    mpz_set(result->z, Pympz_AS_MPZ(self));
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
    PympzObject *s = 0;
    PyObject *result = 0;

    PARSE_ONE_MPZ_REQ_SI(&n,
                         "iroot() requires 'mpz','int' arguments");

    if (n <= 0) {
        VALUE_ERROR("n must be > 0");
        Py_DECREF(self);
        return NULL;
    }
    else if (n>1) {
        if (mpz_sgn(Pympz_AS_MPZ(self))<0) {
            VALUE_ERROR("iroot() of negative number");
            Py_DECREF(self);
            return NULL;
        }
    }
    if (!(s = Pympz_new()) || !(result = PyTuple_New(2))) {
        Py_DECREF(self);
        Py_XDECREF((PyObject*)s);
        Py_XDECREF(result);
        return NULL;
    }
    exact = mpz_root(s->z, Pympz_AS_MPZ(self), n);
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
    PympzObject *r = 0, *y = 0;
    PyObject *result = 0;

    PARSE_ONE_MPZ_REQ_SI(&n,
            "iroot_rem() requires 'mpz','int' arguments");

    if (n <= 0) {
        VALUE_ERROR("n must be > 0");
        Py_DECREF(self);
        return NULL;
    }
    else if (n>1) {
        if (mpz_sgn(Pympz_AS_MPZ(self))<0) {
            VALUE_ERROR("iroot_rem() of negative number");
            Py_DECREF(self);
            return NULL;
        }
    }
    y = Pympz_new();
    r = Pympz_new();
    result = PyTuple_New(2);
    if (!y || !r || !result) {
        Py_DECREF(self);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)y);
        Py_XDECREF((PyObject*)r);
        return NULL;
    }
    mpz_rootrem(y->z, r->z, Pympz_AS_MPZ(self), n);
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, (PyObject*)y);
    PyTuple_SET_ITEM(result, 1, (PyObject*)r);
    return result;
}

static PyObject *
Pympz_sign(PyObject *self, PyObject *other)
{
    long res;
    PympzObject* tempx;

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_sgn(Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_sgn(Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
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
Pympz_abs(PympzObject *x)
{
    PympzObject *result;

    if (!(result = Pympz_new()))
        return NULL;
    mpz_abs(result->z, x->z);
    return (PyObject *)result;
}

static PyObject *
Pyxmpz_abs(PyxmpzObject *x)
{
    mpz_abs(x->z, x->z);
    Py_RETURN_NONE;
}

static PyObject *
Pympz_neg(PympzObject *x)
{
    PympzObject *result;

    if (!(result = Pympz_new()))
        return NULL;
    mpz_neg(result->z, x->z);
    return (PyObject*)result;
}

static PyObject *
Pyxmpz_neg(PyxmpzObject *x)
{
    mpz_neg(x->z, x->z);
    Py_RETURN_NONE;
}

static PyObject *
Pympz_pos(PympzObject *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject*)x;
}

static PyObject *
Pyxmpz_pos(PyxmpzObject *x)
{
    Py_RETURN_NONE;
}

static PyObject *
Pympz_square(PyObject *self, PyObject *other)
{
    PympzObject *tempx, *result;

    if (!(result = Pympz_new()))
        return NULL;

    if (self && (CHECK_MPZANY(self))) {
        mpz_mul(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        mpz_mul(result->z, Pympz_AS_MPZ(other), Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("square() requires 'mpz' argument");
            return NULL;
        }
        else {
            mpz_mul(result->z, Pympz_AS_MPZ(tempx), Pympz_AS_MPZ(tempx));
            Py_DECREF((PyObject*)tempx);
        }
    }
    return (PyObject*)result;
}

/* Pympz_pow is called by Pympany_pow after verifying that all the
 * arguments are integers, but not necessarily mpz.
 */

static PyObject *
Pympz_pow(PyObject *b, PyObject *e, PyObject *m)
{
    PympzObject *result, *tempb = 0, *tempe = 0, *tempm = 0;

    if (!(result = Pympz_new()))
        return NULL;

    tempb = Pympz_From_Integer(b);
    tempe = Pympz_From_Integer(e);

    /* m will either be a number or Py_None. */
    if (m != Py_None) {
        tempm = Pympz_From_Integer(m);
    }

    if (!tempb || !tempe || (!tempm && (m != Py_None))) {
        TYPE_ERROR("Unsupported operand in mpz.pow()");
        goto err;
    }

    if (m == Py_None) {
        /* When no modulo is present, the exponent must fit in mpir_ui
         * the exponent must be positive.
         */
        mpir_ui el;
        if (mpz_sgn(tempe->z) < 0) {
            VALUE_ERROR("pow() exponent cannot be negative");
            goto err;
        }
        if (!mpz_fits_ui_p(tempe->z)) {
            VALUE_ERROR("pow() outrageous exponent");
            goto err;
        }
        el = mpz_get_ui(tempe->z);
        mpz_pow_ui(result->z, tempb->z, el);
    }
    else { /* Modulo exponentiation */
        int sign;
        mpz_t mm, base, exp;

        sign = mpz_sgn(tempm->z);
        if (sign == 0) {
            VALUE_ERROR("pow() 3rd argument cannot be 0");
            goto err;
        }
        mpz_inoc(mm);
        mpz_abs(mm, tempm->z);
        /* A negative exponent is allowed if inverse exists. */
        if (mpz_sgn(tempe->z) < 0) {
            mpz_inoc(base);
            if (!mpz_invert(base, tempb->z, mm)) {
                VALUE_ERROR("pow() base not invertible");
                mpz_cloc(base);
                mpz_cloc(mm);
                goto err;
            }
            else {
                mpz_inoc(exp);
                mpz_abs(exp, tempe->z);
            }
            mpz_powm(result->z, base, exp, mm);
            mpz_cloc(base);
            mpz_cloc(exp);
        }
        else {
            mpz_powm(result->z, tempb->z, tempe->z, mm);
        }
        mpz_cloc(mm);

        /* Python uses a rather peculiar convention for negative modulos
         * If the modulo is negative, result should be in the interval
         * m < r <= 0 .
         */
        if ((sign<0) && (mpz_sgn(Pympz_AS_MPZ(result)) > 0)) {
            mpz_add(result->z, result->z, tempm->z);
        }
    }
    Py_XDECREF((PyObject*)tempb);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempm);
    return (PyObject*)result;

  err:
    Py_XDECREF((PyObject*)tempb);
    Py_XDECREF((PyObject*)tempe);
    Py_XDECREF((PyObject*)tempm);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static int
Pympz_nonzero(PympzObject *x)
{
    return mpz_sgn(x->z) != 0;
}

static int
Pyxmpz_nonzero(PyxmpzObject *x)
{
    return mpz_sgn(x->z) != 0;
}

/* BIT OPERATIONS */

static PyObject *
Pympz_com(PympzObject *x)
{
    PympzObject *result;

    if (!(result = Pympz_new()))
        return NULL;
    mpz_com(result->z, Pympz_AS_MPZ(x));
    return (PyObject*)result;
}

static PyObject *
Pyxmpz_com(PyxmpzObject *x)
{
    mpz_com(x->z, x->z);
    Py_RETURN_NONE;
}

#define MPZ_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *a, PyObject *b) \
{ \
    PympzObject *result = 0; \
    if (CHECK_MPZANY(a)) { \
        if (CHECK_MPZANY(b)) { \
            if (!(result = Pympz_new())) \
                return NULL; \
            NAME(result->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b)); \
        } \
        else { \
            if (!(result = Pympz_From_Integer(b))) \
                return NULL; \
            NAME(result->z, Pympz_AS_MPZ(a), result->z); \
        } \
    } \
    else if (CHECK_MPZANY(b)) { \
        if (CHECK_MPZANY(a)) { \
            if (!(result = Pympz_new())) \
                return NULL; \
            NAME(result->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b)); \
        } \
        else { \
            if (!(result = Pympz_From_Integer(a))) \
                return NULL; \
            NAME(result->z, result->z, Pympz_AS_MPZ(b)); \
        } \
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
Pympz_rshift(PyObject *a, PyObject *b)
{
    mpir_si count_si;
    int overflow;
    PympzObject *result, *tempa, *tempb;

    if (!(result = Pympz_new()))
        return NULL;

    /* Try to make mpz >> Python int/long as fast as possible. */
    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            count_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else if (count_si >= 0) {
                mpz_fdiv_q_2exp(result->z, Pympz_AS_MPZ(a), count_si);
                return (PyObject*)result;
            }
            else {
                VALUE_ERROR("negative shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }
    }

    tempa = Pympz_From_Integer(a);
    tempb = Pympz_From_Integer(b);
    if (!tempb || !tempa) {
        TYPE_ERROR("Pympz_rshift() expects integer arguments");
        goto err;
    }
    if (mpz_sgn(Pympz_AS_MPZ(tempb)) < 0) {
        VALUE_ERROR("negative shift count");
        goto err;
    }
    if(!mpz_fits_si_p(Pympz_AS_MPZ(tempb))) {
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
Pympz_lshift(PyObject *a, PyObject *b)
{
    mpir_si count_si;
    int overflow;
    PympzObject *result, *tempa, *tempb;

    if (!(result = Pympz_new()))
        return NULL;

    /* Try to make mpz >> Python int/long as fast as possible. */
    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            count_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else if (count_si >= 0) {
                mpz_mul_2exp(result->z, Pympz_AS_MPZ(a), count_si);
                return (PyObject*)result;
            }
            else {
                VALUE_ERROR("negative shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
        }
    }

    tempa = Pympz_From_Integer(a);
    tempb = Pympz_From_Integer(b);
    if (!tempb || !tempa) {
        TYPE_ERROR("Pympz_lshift() expects integer arguments");
        goto err;
        }
    if (mpz_sgn(Pympz_AS_MPZ(tempb)) < 0) {
        VALUE_ERROR("negative shift count");
        goto err;
    }
    if(!mpz_fits_si_p(Pympz_AS_MPZ(tempb))) {
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
Pympz_oct(PympzObject *self)
{
    return Pympz_ascii(self, 8, 0);
}

static PyObject *
Pyxmpz_oct(PyxmpzObject *self)
{
    return Pyxmpz_ascii(self, 8, 0);
}

static PyObject *
Pympz_hex(PympzObject *self)
{
    return Pympz_ascii(self, 16, 0);
}

static PyObject *
Pyxmpz_hex(PyxmpzObject *self)
{
    return Pyxmpz_ascii(self, 16, 0);
}
#endif

static Py_hash_t
Pympz_hash(PympzObject *self)
{
#ifdef _PyHASH_MODULUS
    Py_hash_t hash;

    if (self->hash_cache != -1)
        return self->hash_cache;

    hash = (Py_hash_t)mpn_mod_1(self->z->_mp_d, mpz_size(self->z), _PyHASH_MODULUS);
    if (mpz_sgn(self->z)<0)
        hash = -hash;
    if (hash==-1)
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
    PympzObject *result, *tempa, *tempb;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcd() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_gcd(result->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
    }
    else {
        tempa = Pympz_From_Integer(a);
        tempb = Pympz_From_Integer(b);
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
    PympzObject *result, *tempa, *tempb;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("lcm() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_lcm(result->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
    }
    else {
        tempa = Pympz_From_Integer(a);
        tempb = Pympz_From_Integer(b);
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
    PympzObject *g = 0, *s = 0, *t = 0, *tempa, *tempb;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
        return NULL;
    }

    g = Pympz_new();
    s = Pympz_new();
    t = Pympz_new();
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
        mpz_gcdext(g->z, s->z, t->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
    }
    else {
        tempa = Pympz_From_Integer(a);
        tempb = Pympz_From_Integer(b);
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
        mpz_gcdext(g->z, s->z, t->z, Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
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
"Return x such that b*x==a modulo m. Raises a ZeroDivisionError\n"
"exception if no such value x exists.");

static PyObject *
Pygmpy_divm(PyObject *self, PyObject *args)
{
    PympzObject *result, *num, *den, *mod;
    mpz_t gcdz;
    int ok;

    if(PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;

    num = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    den = Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
    mod = Pympz_From_Integer(PyTuple_GET_ITEM(args, 2));
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
    PympzObject *result;
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
        if (!(result = Pympz_new()))
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
    PympzObject *result;
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
        if (!(result = Pympz_new()))
            return NULL;
        mpz_fib_ui(Pympz_AS_MPZ(result), n);
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
    PympzObject *fib1 = 0, *fib2 = 0;
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
    PympzObject *result;
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
        if (!(result = Pympz_new()))
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
    PympzObject *luc1, *luc2;
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
    PympzObject *result;
    mpir_si k;

    PARSE_ONE_MPZ_REQ_SI(&k, "bincoef() requires 'mpz','int' arguments");

    if (k < 0) {
        VALUE_ERROR("binomial coefficient with negative k");
        Py_DECREF(self);
        return NULL;
    }

    if(!(result = Pympz_new())) {
        Py_DECREF(self);
        return NULL;
    }
    mpz_bin_ui(result->z, Pympz_AS_MPZ(self), k);
    Py_DECREF(self);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_isqrt,
"isqrt(x) -> mpz\n\n"
"Return the integer square root of an integer x. x >= 0.");

static PyObject *
Pympz_isqrt(PyObject *self, PyObject *other)
{
    PympzObject *result;

    if (self && (CHECK_MPZANY(self))) {
        if (mpz_sgn(Pympz_AS_MPZ(self)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if (!(result = Pympz_new()))
            return NULL;
        mpz_sqrt(result->z, Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) < 0) {
            VALUE_ERROR("isqrt() of negative number");
            return NULL;
        }
        if (!(result = Pympz_new()))
            return NULL;
        mpz_sqrt(result->z, Pympz_AS_MPZ(other));
    }
    else {
        if (!(result = Pympz_From_Integer(other))) {
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
    PympzObject *root = 0, *rem = 0;
    PyObject *result = 0;

    PARSE_ONE_MPZ("isqrt_rem() requires 'mpz' argument");

    if (mpz_sgn(Pympz_AS_MPZ(self)) < 0) {
        VALUE_ERROR("isqrt_rem() of negative number");
        Py_DECREF(self);
        return NULL;
    }

    root = Pympz_new();
    rem = Pympz_new();
    result = PyTuple_New(2);
    if (!root || !rem || !result) {
        Py_DECREF(self);
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)root);
        Py_XDECREF((PyObject*)rem);
        return NULL;
    }

    mpz_sqrtrem(root->z, rem->z, Pympz_AS_MPZ(self));
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
    PympzObject *result;
    PyObject *factor;
    size_t multiplicity;

    PARSE_TWO_MPZ(factor, "remove() requires 'mpz','mpz' arguments");

    if (mpz_cmp_si(Pympz_AS_MPZ(factor), 2) < 0) {
        VALUE_ERROR("factor must be > 1");
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }

    if (!(result = Pympz_new())) {
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }
    multiplicity = mpz_remove(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(factor));
    Py_DECREF(self);
    Py_DECREF(factor);
    return Py_BuildValue("(Nk)", result, multiplicity);
}

PyDoc_STRVAR(doc_invertg,
"invert(x, m) -> mpz\n\n"
"Return the y such that x*y==1 modulo m, or 0 if no such y exists.");

static PyObject *
Pygmpy_invert(PyObject *self, PyObject *args)
{
    PyObject *x, *y;
    PympzObject *result, *tempx, *tempy;
    int success;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;
    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        success = mpz_invert(result->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
        if (!success)
            mpz_set_ui(result->z, 0);
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
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
        if (!success)
            mpz_set_ui(result->z, 0);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
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
            mpz_hamdist(Pympz_AS_MPZ(self),Pympz_AS_MPZ(other)));
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
    PympzObject *result, *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
        return NULL;
    }

    if (!(result = Pympz_new()))
        return NULL;
    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);

    if (CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if (mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        mpz_divexact(result->z, Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    }
    else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if (!tempx || !tempy) {
            TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        if (mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
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

PyDoc_STRVAR(doc_is_squarem,
"x.is_square() -> bool\n\n"
"Returns True if x is a perfect square, else return False.");

PyDoc_STRVAR(doc_is_squareg,
"is_square(x) -> bool\n\n"
"Returns True if x is a perfect square, else return False.");

static PyObject *
Pympz_is_square(PyObject *self, PyObject *other)
{
    int res;
    PympzObject *tempx;

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_perfect_square_p(Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_perfect_square_p(Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
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

PyDoc_STRVAR(doc_is_powerm,
"x.is_power() -> bool\n\n"
"Return True if x is a perfect power (there exists a y and an\n"
"n > 1, such that x=y**n), else return False.");

PyDoc_STRVAR(doc_is_powerg,
"is_power(x) -> bool\n\n"
"Return True if x is a perfect power (there exists a y and an\n"
"n > 1, such that x=y**n), else return False.");

static PyObject *
Pympz_is_power(PyObject *self, PyObject *other)
{
    int res;
    PympzObject* tempx;

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_perfect_power_p(Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_perfect_power_p(Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
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

PyDoc_STRVAR(doc_is_primem,
"x.is_prime([n=25]) -> bool\n\n"
"Return True if x is _probably_ prime, else False if x is\n"
"definately composite. x is checked for small divisors and up\n"
"to n Miller-Rabin tests are performed.");

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
            "is_prime() requires 'mpz',['int'] arguments");

    if (reps <= 0) {
        VALUE_ERROR("repetition count for is_prime must be positive");
        Py_DECREF(self);
        return NULL;
    }
    i = mpz_probab_prime_p(Pympz_AS_MPZ(self), reps);
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
    PympzObject *result;

    if(self && Pympz_Check(self)) {
        if(!(result = Pympz_new()))
            return NULL;
        mpz_nextprime(result->z, Pympz_AS_MPZ(self));
    }
    else if(CHECK_MPZANY(other)) {
        if(!(result = Pympz_new()))
            return NULL;
        mpz_nextprime(result->z, Pympz_AS_MPZ(other));
    }
    else {
        if (!(result = Pympz_From_Integer(other))) {
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

    if (mpz_sgn(Pympz_AS_MPZ(other)) <= 0 ||
        mpz_even_p(Pympz_AS_MPZ(other))) {
        VALUE_ERROR("y must be odd and >0");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    i = (long)(mpz_jacobi(Pympz_AS_MPZ(self), Pympz_AS_MPZ(other)));
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

    if (mpz_sgn(Pympz_AS_MPZ(other)) <= 0 ||
        mpz_even_p(Pympz_AS_MPZ(other))) {
        VALUE_ERROR("y must be odd and >0");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    i = (long) mpz_legendre(Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
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

    ires = (long) mpz_kronecker(Pympz_AS_MPZ(self), (Pympz_AS_MPZ(other)));

    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(ires);
}

PyDoc_STRVAR(doc_is_evenm,
"x.is_even() -> bool\n\n"
"Return True if x is even, False otherwise.");

PyDoc_STRVAR(doc_is_eveng,
"is_even(x) -> bool\n\n"
"Return True if x is even, False otherwise.");

static PyObject *
Pympz_is_even(PyObject *self, PyObject *other)
{
    int res;
    PympzObject *tempx;

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_even_p(Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_even_p(Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
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

PyDoc_STRVAR(doc_is_oddm,
"x.is_odd() -> bool\n\n"
"Return True if x is odd, False otherwise.");

PyDoc_STRVAR(doc_is_oddg,
"is_odd(x) -> bool\n\n"
"Return True if x is odd, False otherwise.");

static PyObject *
Pympz_is_odd(PyObject *self, PyObject *other)
{
    int res;
    PympzObject *tempx;

    if (self && (CHECK_MPZANY(self))) {
        res = mpz_odd_p(Pympz_AS_MPZ(self));
    }
    else if (CHECK_MPZANY(other)) {
        res = mpz_odd_p(Pympz_AS_MPZ(other));
    }
    else {
        if (!(tempx = Pympz_From_Integer(other))) {
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

PyDoc_STRVAR(doc_make_mpzm,
"xmpz.make_mpz() -> mpz\n\n"
"Return an mpz by converting an 'xmpz' to an 'mpz' as quickly as\n"
"possible.\n\n"
"NOTE: Optimized for speed so the original xmpz is set to 0!.");

static PyObject *
Pyxmpz_make_mpz(PyObject *self, PyObject *other)
{
    PympzObject* result;

    if (!(result = Pympz_new()))
        return NULL;
    mpz_swap(result->z, Pympz_AS_MPZ(self));
    mpz_set_ui(Pympz_AS_MPZ(self), 0);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_xmpz_copy,
"xmpz.copy() -> xmpz\n\n"
"Return a copy of an xmpz.");

static PyObject *
Pyxmpz_copy(PyObject *self, PyObject *other)
{
    return (PyObject*)Pyxmpz2Pyxmpz(self);
}

/*
 * Add mapping support to xmpz objects.
 */

static Py_ssize_t
Pyxmpz_nbits(PyxmpzObject *obj)
{
    return mpz_sizeinbase(obj->z, 2);
}

static PyObject *
Pyxmpz_subscript(PyxmpzObject* self, PyObject* item)
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
        PyObject* result;

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

        if (!(result = (PyObject*)Pyxmpz_new()))
            return NULL;
        mpz_set_ui(Pympz_AS_MPZ(result), 0);
        if (slicelength > 0) {
            for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                if (mpz_tstbit(self->z, cur)) {
                    mpz_setbit(Pympz_AS_MPZ(result), i);
                }
            }
        }
        return result;
    }
    else {
        TYPE_ERROR("bit positions must be integers");
        return NULL;
    }
}

static int
Pyxmpz_assign_subscript(PyxmpzObject* self, PyObject* item, PyObject* value)
{
    if (PyIndex_Check(item)) {
        Py_ssize_t bit_value, i;

        i = PyNumber_AsSsize_t(item, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred())
            return -1;
        if (i < 0)
            i += mpz_sizeinbase(self->z, 2);

        bit_value = PyNumber_AsSsize_t(value, PyExc_ValueError);
        if (bit_value == -1 && PyErr_Occurred()) {
            VALUE_ERROR("bit value must be 0 or 1");
            return -1;
        }
        if (bit_value == 1) {
            mpz_setbit(self->z, i);
            return 0;
        }
        else if (bit_value == 0) {
            mpz_clrbit(self->z, i);
            return 0;
        }
        else {
            VALUE_ERROR("bit value must be 0 or 1");
            return -1;
        }
    }
    else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step, slicelength;

#if PY_VERSION_HEX > 0x030200A4
        if (PySlice_GetIndicesEx(item,
#else
        if (PySlice_GetIndicesEx((PySliceObject*)item,
#endif
                        mpz_sizeinbase(self->z, 2),
                        &start, &stop, &step, &slicelength) < 0) {
            return -1;
        }

        if ((step < 0 && start < stop) ||
            (step > 0 && start > stop))
            stop = start;

        if (value == NULL) {
            TYPE_ERROR("deleting bits not supported");
            return -1;
        }
        else {
            Py_ssize_t cur, i;
            int bit;
            PympzObject *tempx;

            if (!(tempx = Pympz_From_Integer(value))) {
                VALUE_ERROR("must specify bit sequence as an integer");
                return -1;
            }
            if (mpz_sgn(tempx->z) == 0) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_clrbit(self->z, cur);
                }
            }
            else if (!(mpz_cmp_si(tempx->z, -1))) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_setbit(self->z, cur);
                }
            }
            else {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    bit = mpz_tstbit(tempx->z, i);
                    if (bit)
                        mpz_setbit(self->z, cur);
                    else
                        mpz_clrbit(self->z, cur);
                }
            }
            Py_DECREF((PyObject*)tempx);
            return 0;
        }
    }
    else {
        TYPE_ERROR("bit positions must be integers");
        return -1;
    }
    return -1;
}

/*
 * Add mapping support to mpz objects.
 */

static Py_ssize_t
Pympz_nbits(PyxmpzObject *obj)
{
    return mpz_sizeinbase(obj->z, 2);
}

static PyObject *
Pympz_subscript(PyxmpzObject* self, PyObject* item)
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
        PyObject* result;

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

        if (!(result = (PyObject*)Pympz_new()))
            return NULL;
        mpz_set_ui(Pympz_AS_MPZ(result), 0);
        if (slicelength > 0) {
            for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                if(mpz_tstbit(self->z, cur)) {
                    mpz_setbit(Pympz_AS_MPZ(result), i);
                }
            }
        }
        return result;
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
    int base = 10, option = 0;
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
                seensign = 1;
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
        if (*p1 == 'X') {
            base = -16;
            break;
        }
        VALUE_ERROR("Invalid conversion specification");
        return NULL;
    }
    *(p2++) = '\00';

    mpzstr = mpz_ascii(Pympz_AS_MPZ(self), base, option);
    if (!mpzstr)
        return NULL;
    result = PyObject_CallMethod(mpzstr, "__format__", "(s)", fmt);
    if (!result) {
        Py_DECREF(mpzstr);
        return NULL;
    }
    Py_DECREF(mpzstr);
    return result;
}

static PyObject *
Pympz_add(PyObject *self, PyObject *args)
{
    PympzObject *result;
    PyObject *other;

    PARSE_TWO_MPZ(other, "add() requires 'mpz','mpz' arguments");

    if ((result = Pympz_new()))
        mpz_add(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympz_sub(PyObject *self, PyObject *args)
{
    PympzObject *result;
    PyObject *other;

    PARSE_TWO_MPZ(other, "sub() requires 'mpz','mpz' arguments");

    if ((result = Pympz_new()))
        mpz_sub(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympz_mul(PyObject *self, PyObject *args)
{
    PympzObject *result;
    PyObject *other;

    PARSE_TWO_MPZ(other, "mul() requires 'mpz','mpz' arguments");

    if ((result = Pympz_new()))
        mpz_mul(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static PyObject *
Pympz_div(PyObject *self, PyObject *args)
{
    PympzObject *result;
    PyObject *other;

    PARSE_TWO_MPZ(other, "div() requires 'mpz','mpz' arguments");

    if ((result = Pympz_new())) {
        if (mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("mpz division by zero");
            Py_DECREF((PyObject*)result);
            result = 0;
        }
        else {
            mpz_div(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
        }
    }

    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

#ifdef PY3
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympz_neg,               /* nb_negative             */
    (unaryfunc) Pympz_pos,               /* nb_positive             */
    (unaryfunc) Pympz_abs,               /* nb_absolute             */
    (inquiry) Pympz_nonzero,             /* nb_bool                 */
    (unaryfunc) Pympz_com,               /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
    (unaryfunc) Pympz2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (binaryfunc) Pympz_inplace_add,      /* nb_inplace_add          */
    (binaryfunc) Pympz_inplace_sub,      /* nb_inplace_subtract     */
    (binaryfunc) Pympz_inplace_mul,      /* nb_inplace_multiply     */
    (binaryfunc) Pympz_inplace_rem,      /* nb_inplace_remainder    */
    (ternaryfunc) Pympz_inplace_pow,     /* nb_inplace_power        */
    (binaryfunc) Pympz_inplace_lshift,   /* nb_inplace_lshift       */
    (binaryfunc) Pympz_inplace_rshift,   /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv, /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc)  Pympz_To_Integer,       /* nb_index                */
};

#else
static PyNumberMethods mpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympz_neg,               /* nb_negative             */
    (unaryfunc) Pympz_pos,               /* nb_positive             */
    (unaryfunc) Pympz_abs,               /* nb_absolute             */
    (inquiry) Pympz_nonzero,             /* nb_bool                 */
    (unaryfunc) Pympz_com,               /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympz_To_Integer,        /* nb_int                  */
    (unaryfunc) Pympz2PyLong,            /* nb_long                 */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (unaryfunc) Pympz_oct,               /* nb_oct                  */
    (unaryfunc) Pympz_hex,               /* nb_hex                  */
    (binaryfunc) Pympz_inplace_add,      /* nb_inplace_add          */
    (binaryfunc) Pympz_inplace_sub,      /* nb_inplace_subtract     */
    (binaryfunc) Pympz_inplace_mul,      /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
    (binaryfunc) Pympz_inplace_rem,      /* nb_inplace_remainder    */
    (ternaryfunc) Pympz_inplace_pow,     /* nb_inplace_power        */
    (binaryfunc) Pympz_inplace_lshift,   /* nb_inplace_lshift       */
    (binaryfunc) Pympz_inplace_rshift,   /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pympz_inplace_floordiv, /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc) Pympz_To_Integer,        /* nb_index                */
};
#endif

static PyMappingMethods mpz_mapping_methods = {
    (lenfunc)Pympz_nbits,
    (binaryfunc)Pympz_subscript,
    NULL
};

#ifdef PY3
static PyNumberMethods xmpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pyxmpz_neg,              /* nb_negative             */
    (unaryfunc) Pyxmpz_pos,              /* nb_positive             */
    (unaryfunc) Pyxmpz_abs,              /* nb_absolute             */
    (inquiry) Pyxmpz_nonzero,            /* nb_bool                 */
    (unaryfunc) Pyxmpz_com,              /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
    (unaryfunc) Pympz2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (binaryfunc) Pyxmpz_inplace_add,     /* nb_inplace_add          */
    (binaryfunc) Pyxmpz_inplace_sub,     /* nb_inplace_subtract     */
    (binaryfunc) Pyxmpz_inplace_mul,     /* nb_inplace_multiply     */
    (binaryfunc) Pyxmpz_inplace_rem,     /* nb_inplace_remainder    */
    (ternaryfunc) Pyxmpz_inplace_pow,    /* nb_inplace_power        */
    (binaryfunc) Pyxmpz_inplace_lshift,  /* nb_inplace_lshift       */
    (binaryfunc) Pyxmpz_inplace_rshift,  /* nb_inplace_rshift       */
    (binaryfunc) Pyxmpz_inplace_and,     /* nb_inplace_and          */
    (binaryfunc) Pyxmpz_inplace_xor,     /* nb_inplace_xor          */
    (binaryfunc) Pyxmpz_inplace_ior,     /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pyxmpz_inplace_floordiv,/* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc)  Pyxmpz_To_Integer,      /* nb_index                */
};

#else
static PyNumberMethods xmpz_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pyxmpz_neg,              /* nb_negative             */
    (unaryfunc) Pyxmpz_pos,              /* nb_positive             */
    (unaryfunc) Pyxmpz_abs,              /* nb_absolute             */
    (inquiry) Pyxmpz_nonzero,            /* nb_bool                 */
    (unaryfunc) Pyxmpz_com,              /* nb_invert               */
    (binaryfunc) Pympz_lshift,           /* nb_lshift               */
    (binaryfunc) Pympz_rshift,           /* nb_rshift               */
    (binaryfunc) Pympz_and,              /* nb_and                  */
    (binaryfunc) Pympz_xor,              /* nb_xor                  */
    (binaryfunc) Pympz_ior,              /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympz_To_Integer,        /* nb_int                  */
    (unaryfunc) Pympz2PyLong,            /* nb_long                 */
    (unaryfunc) Pympz2PyFloat,           /* nb_float                */
    (unaryfunc) Pyxmpz_oct,              /* nb_oct                  */
    (unaryfunc) Pyxmpz_hex,              /* nb_hex                  */
    (binaryfunc) Pyxmpz_inplace_add,     /* nb_inplace_add          */
    (binaryfunc) Pyxmpz_inplace_sub,     /* nb_inplace_subtract     */
    (binaryfunc) Pyxmpz_inplace_mul,     /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
    (binaryfunc) Pyxmpz_inplace_rem,     /* nb_inplace_remainder    */
    (ternaryfunc) Pyxmpz_inplace_pow,    /* nb_inplace_power        */
    (binaryfunc) Pyxmpz_inplace_lshift,  /* nb_inplace_lshift       */
    (binaryfunc) Pyxmpz_inplace_rshift,  /* nb_inplace_rshift       */
    (binaryfunc) Pyxmpz_inplace_and,     /* nb_inplace_and          */
    (binaryfunc) Pyxmpz_inplace_xor,     /* nb_inplace_xor          */
    (binaryfunc) Pyxmpz_inplace_ior,     /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
    (binaryfunc) Pyxmpz_inplace_floordiv,/* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
    (unaryfunc) Pyxmpz_To_Integer,       /* nb_index                */
};
#endif

static PyMappingMethods xmpz_mapping_methods = {
    (lenfunc)Pyxmpz_nbits,
    (binaryfunc)Pyxmpz_subscript,
    (objobjargproc)Pyxmpz_assign_subscript
};

static PyMethodDef Pympz_methods [] =
{
    { "__format__", Pympz_format, METH_VARARGS, doc_mpz_format },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "digits", Pympz_digits, METH_VARARGS, doc_mpz_digits },
    { "is_even", Pympz_is_even, METH_NOARGS, doc_is_evenm },
    { "is_odd", Pympz_is_odd, METH_NOARGS, doc_is_oddm },
    { "is_square", Pympz_is_square, METH_NOARGS, doc_is_squarem },
    { "is_power", Pympz_is_power, METH_NOARGS, doc_is_powerm },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primem },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsm },
    { NULL, NULL, 1 }
};

static PyMethodDef Pyxmpz_methods [] =
{
    { "__format__", Pympz_format, METH_VARARGS, doc_mpz_format },
    { "bit_clear", Pympz_bit_clear, METH_O, doc_bit_clearm },
    { "bit_flip", Pympz_bit_flip, METH_O, doc_bit_flipm },
    { "bit_length", Pympz_bit_length, METH_NOARGS, doc_bit_lengthm },
    { "bit_scan0", Pympz_bit_scan0, METH_VARARGS, doc_bit_scan0m },
    { "bit_scan1", Pympz_bit_scan1, METH_VARARGS, doc_bit_scan1m },
    { "bit_set", Pympz_bit_set, METH_O, doc_bit_setm },
    { "bit_test", Pympz_bit_test, METH_O, doc_bit_testm },
    { "copy", Pyxmpz_copy, METH_NOARGS, doc_xmpz_copy },
    { "digits", Pyxmpz_digits, METH_VARARGS, doc_mpz_digits },
    { "is_even", Pympz_is_even, METH_NOARGS, doc_is_evenm },
    { "is_odd", Pympz_is_odd, METH_NOARGS, doc_is_oddm },
    { "is_square", Pympz_is_square, METH_VARARGS, doc_is_squarem },
    { "is_power", Pympz_is_power, METH_VARARGS, doc_is_powerm },
    { "is_prime", Pympz_is_prime, METH_VARARGS, doc_is_primem },
    { "make_mpz", Pyxmpz_make_mpz, METH_NOARGS, doc_make_mpzm },
    { "numdigits", Pympz_numdigits, METH_VARARGS, doc_numdigitsm },
    { NULL, NULL, 1 }
};

static PyTypeObject Pympz_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpz",                                  /* tp_name          */
    sizeof(PympzObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympz_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympz2repr,                  /* tp_repr          */
    &mpz_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &mpz_mapping_methods,                   /* tp_as_mapping    */
    (hashfunc) Pympz_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympz2str,                   /* tp_str           */
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
};

static PyTypeObject Pyxmpz_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "xmpz",                                 /* tp_name          */
    sizeof(PyxmpzObject),                   /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pyxmpz_dealloc,            /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pyxmpz2repr,                 /* tp_repr          */
    &xmpz_number_methods,                   /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &xmpz_mapping_methods,                  /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pyxmpz2str,                  /* tp_str           */
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
    Pyxmpz_methods,                         /* tp_methods       */
};



