/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz.c                                                              *
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


/* Functions that operate strictly on mpz or xmpz. */

/* produce digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(doc_mpz_digits,
"x.digits([base]) -> string\n\n"
"Return Python string representing x in the given base (2 to 62,\n"
"default 10 if omitted or 0); leading '-' is present if x<0, but\n"
"no leading '+' if x>=0.");
static PyObject *
Pympz_digits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "digits() requires 'mpz',['int'] arguments");
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
            "digits() requires 'xmpz',['int'] arguments");
    result = Pyxmpz_ascii((PyxmpzObject*)self, base, 0);
    Py_DECREF(self);
    return result;
}

/* return number-of-digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(doc_numdigitsm,
"x.numdigits([base]) -> int\n\n"
"Return length of string representing x in the given base (2 to 62,\n"
"default 10 if omitted or 0); the value returned may sometimes be 1\n"
"more than necessary; no provision for any 'sign' character, nor\n"
"leading '0' or '0x' decoration, is made in the returned length.");

PyDoc_STRVAR(doc_numdigitsg,
"numdigits(x[,base]) -> int\n\n"
"Return length of string representing x in the given base (2 to 62,\n"
"default 10 if omitted or 0); the value returned may sometimes be 1\n"
"more than necessary; no provision for any 'sign' character, nor\n"
"leading '0' or '0x' decoration, is made in the returned length.");
static PyObject *
Pympz_numdigits(PyObject *self, PyObject *args)
{
    long base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "numdigits() requires 'mpz',['int'] arguments");
    if (base == 0)
        base = 10;
    if ((base < 2) || (base > 62)) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 62");
        Py_DECREF(self);
        return NULL;
    }
    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(Pympz_AS_MPZ(self),
                                    (int)base));
    Py_DECREF(self);
    return result;
}

static char doc_bit_lengthm[]="\
x.bit_length(): returns the number of significant bits in the\n\
radix-2 representation of x. Note: bit_length(0) returns 0.\n\
";
static char doc_bit_lengthg[]="\
bit_length(x): returns the number of significant bits in the\n\
radix-2 representation of x. Note: bit_length(0) returns 0.\n\
";
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

static char doc_bit_maskg[]="\
bit_mask(n): create an 'mpz' exactly n bits in length with\n\
all bits set.\n\
";
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

static char doc_xbit_maskg[]="\
xbit_mask(n): create an 'xmpz' exactly n bits in length with\n\
all bits set.\n\
";
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
static char doc_bit_scan0m[]="\
x.bit_scan0(n=0): returns the index of the first 0-bit of x that\n\
is >= n. n must be >= 0.  If no more 0-bits are in x at or above\n\
index n (which can only happen for x<0, notionally extended with\n\
infinite 1-bits), None is returned.\n\
";
static char doc_bit_scan0g[]="\
bit_scan0(x, n=0): returns the index of the first 0-bit of x that\n\
is >= n. n must be >= 0.  If no more 0-bits are in x at or above\n\
index n (which can only happen for x<0, notionally extended with\n\
infinite 1-bits), None is returned.\n\
";
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

static char doc_bit_scan1m[]="\
x.bit_scan1(n=0): returns the index of the first 1-bit of x that\n\
is >= n. n must be >= 0.  If no more 1-bits are in x at or above\n\
index n (which can only happen for x>=0, notionally extended with\n\
infinite 0-bits), None is returned.\n\
";
static char doc_bit_scan1g[]="\
bit_scan1(x, n=0): returns the index of the first 1-bit of x that\n\
is >= n. n must be >= 0.  If no more 1-bits are in x at or above\n\
index n (which can only happen for x>=0, notionally extended with\n\
infinite 0-bits), None is returned.\n\
";
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
static char doc_popcountg[]="\
popcount(x): returns the number of 1-bits set in x; note that\n\
this is 'infinite' if x<0, and in that case, -1 is returned.\n\
x must be an mpz, or else gets coerced to one.\n\
";
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
static char doc_bit_testm[]="\
x.bit_test(n): return the value of the n-th bit of x.\n\
";
static char doc_bit_testg[]="\
bit_test(x,n): return the value of the n-th bit of x.\n\
";
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

static char doc_bit_clearm[]="\
x.bit_clear(n): clear the n-th bit of x. \n\
";
static char doc_bit_clearg[]="\
bit_clear(x,n): clear the n-th bit of x.\n\
";

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

static char doc_bit_setm[]="\
x.bit_set(n): set the n-th bit of x. \n\
";
static char doc_bit_setg[]="\
bit_set(x,n): set the n-th bit of x.\n\
";

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

static char doc_bit_flipm[]="\
x.bit_flip(n): complements the n-th bit of x.";
static char doc_bit_flipg[]="\
bit_flip(x,n): complements the n-th bit of x.";

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
"iff the root is exact. x must be a non-negative integer.");

static PyObject *
Pympz_iroot(PyObject *self, PyObject *args)
{
    long n;
    int exact;
    PympzObject *s = 0;
    PyObject *result = 0;

    PARSE_ONE_MPZ_REQ_CLONG(&n,
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
    if(!(s = Pympz_new()) || !(result = PyTuple_New(2))) {
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
"root of x and x=y**n + r. x must be a non-negative integer. n must\n"
"be a positive integer.");

static PyObject *
Pympz_iroot_rem(PyObject *self, PyObject *args)
{
    long n;
    PympzObject *r = 0, *y = 0;
    PyObject *result = 0;

    PARSE_ONE_MPZ_REQ_CLONG(&n,
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
        /* When no modulo is present, the exponent must fit in C ulong
         * the exponent must be positive.
         */
        unsigned long el;
        if (mpz_sgn(tempe->z) < 0) {
            VALUE_ERROR("pow() exponent cannot be negative");
            goto err;
        }
        if (!mpz_fits_ulong_p(tempe->z)) {
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

/* BIT OPERATIONS (mpz-only) */

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
    long count;
    int overflow;
    PympzObject *result, *tempa, *tempb;

    if (!(result = Pympz_new()))
        return NULL;

    /* Try to make mpz >> Python int/long as fast as possible. */
    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            count = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else if (count >= 0) {
                mpz_fdiv_q_2exp(result->z, Pympz_AS_MPZ(a), count);
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
    if(!mpz_fits_slong_p(Pympz_AS_MPZ(tempb))) {
        PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
        goto err;
    }
    count = mpz_get_si(tempb->z);
    mpz_fdiv_q_2exp(result->z, tempa->z, count);
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
    long count;
    int overflow;
    PympzObject *result, *tempa, *tempb;

    if (!(result = Pympz_new()))
        return NULL;

    /* Try to make mpz >> Python int/long as fast as possible. */
    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            count = PyLong_AsLongAndOverflow(b, &overflow);
            if (overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF((PyObject*)result);
                return NULL;
            }
            else if (count >= 0) {
                mpz_mul_2exp(result->z, Pympz_AS_MPZ(a), count);
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
    if(!mpz_fits_slong_p(Pympz_AS_MPZ(tempb))) {
        PyErr_SetString(PyExc_OverflowError, "outrageous shift count");
        goto err;
    }
    count = mpz_get_si(tempb->z);
    mpz_mul_2exp(result->z, tempa->z, count);
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
static char doc_gcd[]="\
gcd(a,b): returns the greatest common denominator of numbers a and b\n\
(which must be mpz objects, or else get coerced to mpz)\n\
";
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

static char doc_lcm[]="\
lcm(a,b): returns the lowest common multiple of numbers a and b\n\
(which must be mpz objects, or else get coerced to mpz)\n\
";
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

static char doc_gcdext[]="\
gcdext(a,b): returns a 3-element tuple (g,s,t) such that\n\
    g==gcd(a,b) and g == a*s + b*t\n\
(a and b must be mpz objects, or else get coerced to mpz)\n\
";
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

static char doc_divm[]="\
divm(a,b,m): returns x such that b*x==a modulo m, or else raises\n\
a ZeroDivisionError exception if no such value x exists\n\
(a, b and m must be mpz objects, or else get coerced to mpz)\n\
";
static PyObject *
Pygmpy_divm(PyObject *self, PyObject *args)
{
    PyObject *a, *b, *m;
    PympzObject *result, *num, *den, *mod;
    mpz_t gcdz;
    int ok;

    if(PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);
    m = PyTuple_GET_ITEM(args, 2);

    if (!(result = Pympz_new()))
        return NULL;

    num = Pympz_From_Integer(a);
    den = Pympz_From_Integer(b);
    mod = Pympz_From_Integer(m);
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
"Return the exact factorial of n.\n"
"See factorial(n) to get the floating-point approximation.");
static PyObject *
Pygmpy_fac(PyObject *self, PyObject *other)
{
    PympzObject *result;
    long n;

    n = clong_From_Integer(other);
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

static char doc_fib[]="\
fib(n): returns the n-th Fibonacci number; takes O(n) time; n must be\n\
an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_fib(PyObject *self, PyObject *other)
{
    PympzObject *result;
    long n;

    n = clong_From_Integer(other);
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

static char doc_fib2[]="\
fib2(n): returns the n-th and n+1-th Fibonacci numbers; takes O(n) time;\n\
n must be an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_fib2(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *fib1 = 0, *fib2 = 0;
    long n;

    n = clong_From_Integer(other);
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

static char doc_lucas[]="\
lucas(n): returns the n-th Lucas number; takes O(n) time; n must be\n\
an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_lucas(PyObject *self, PyObject *other)
{
    PympzObject *result;
    long n;

    n = clong_From_Integer(other);
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

static char doc_lucas2[]="\
lucas2(n): returns the n-1 and n-th Lucas number; takes O(n) time;\n\
n must be an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_lucas2(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *luc1, *luc2;
    long n;

    n = clong_From_Integer(other);
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

static char doc_bincoefg[]="\
bincoef(x,n): returns the 'binomial coefficient' that is 'x\n\
over n'; n is an ordinary Python int, >=0; x must be an mpz,\n\
or else gets converted to one.\n\
";
static char doc_combg[]="\
comb(x,n): returns the 'number of combinations' of 'x things,\n\
taken n at a time'; n is an ordinary Python int, >=0; x must be\n\
an mpz, or else gets converted to one.\n\
";
static PyObject *
Pympz_bincoef(PyObject *self, PyObject *args)
{
    PympzObject *result;
    long k;

    PARSE_ONE_MPZ_REQ_CLONG(&k,
                            "bincoef() requires 'mpz','int' arguments");

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
"isqrt(x) -> number\n\n"
"Return the integer square root of x; x must be a non-\n"
"negative integer.");

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
"isqrt_rem(x) -> (square_root, remainder)\n\n"
"Return a 2-element tuple (s,t), such that s=isqrt(x) and x=s*s+t.\n"
"x must be a non-negative integer.");

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

static char doc_removeg[]="\
remove(x,f): returns a 2-element tuple (y,m) such that\n\
x==y*(f**m), and y%f==0; i.e., y is x with any factor f\n\
removed, and m (an ordinary Python int) is the multiplicity\n\
of the factor f in x (m=0, and y=x, unless x%f==0). x must\n\
be an mpz, or else gets coerced to one; f must be > 0.\n\
";
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

static char doc_invertg[]="\
invert(x,m): returns the inverse of x modulo m, i.e., that y\n\
such that x*y==1 modulo m, or 0 if no such y exists.\n\
m must be an ordinary Python int, !=0; x must be an mpz,\n\
or else gets converted to one.\n\
";
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

static char doc_hamdistg[]="\
hamdist(x,y): returns the Hamming distance (number of bit-positions\n\
where the bits differ) between x and y.  x and y must be mpz, or else\n\
get coerced to mpz.\n\
";
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

static char doc_divexactg[]="\
divexact(x,y): returns the quotient of x divided by y. Faster than\n\
standard division but requires the remainder is zero!  x and y must\n\
be mpz, or else get coerced to mpz.\n\
";
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

static char doc_is_squarem[]="\
x.is_square(): returns 1 if x is a perfect square, else 0.\n\
";
static char doc_is_squareg[]="\
is_square(x): returns 1 if x is a perfect square, else 0.\n\
x must be an mpz, or else gets coerced to one.\n\
";
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

static char doc_is_powerm[]="\
x.is_power(): returns 1 if x is a perfect power, i.e., there exist\n\
y, and n>1, such that x==y**n; else, 0.\n\
";
static char doc_is_powerg[]="\
is_power(x): returns 1 if x is a perfect power, i.e., there exist\n\
y, and n>1, such that x==y**n; else, 0. x must be an mpz, or else\n\
gets coerced to one.\n\
";
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

static char doc_is_primem[]="\
x.is_prime(n=25): returns 2 if x is _certainly_ prime, 1 if x is\n\
_probably_ prime (probability > 1 - 1/2**n), 0 if x is composite.\n\
If x<0, GMP considers x 'prime' iff -x is prime; gmpy reflects this\n\
GMP design choice.\n\
";
static char doc_is_primeg[]="\
is_prime(x,n=25): returns 2 if x is _certainly_ prime, 1 if x is\n\
_probably_ prime (probability > 1 - 1/2**n), 0 if x is composite.\n\
If x<0, GMP considers x 'prime' iff -x is prime; gmpy reflects this\n\
GMP design choice. x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_is_prime(PyObject *self, PyObject *args)
{
    long i;
    int reps = 25;

    PARSE_ONE_MPZ_OPT_CLONG(&reps,
            "is_prime() requires 'mpz',['int'] arguments");

    if (reps <= 0) {
        VALUE_ERROR("repetition count for is_prime must be positive");
        Py_DECREF(self);
        return NULL;
    }
    i = (long) mpz_probab_prime_p(Pympz_AS_MPZ(self), reps);
    Py_DECREF(self);
    if (i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static char doc_next_primeg[]="\
next_prime(x): returns the smallest prime number > x.  Note that\n\
GMP may use a probabilistic definition of 'prime', and also that\n\
if x<0 GMP considers x 'prime' iff -x is prime; gmpy reflects these\n\
GMP design choices. x must be an mpz, or else gets coerced to one.\n\
";
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

static char doc_jacobig[]="\
jacobi(x,y): returns the Jacobi symbol (x|y) (y should be odd and\n\
must be positive); x and y must be mpz, or else get coerced to mpz.\n\
";
static PyObject *
Pympz_jacobi(PyObject *self, PyObject *args)
{
    PyObject *other;
    long i;

    PARSE_TWO_MPZ(other, "jacobi() requires 'mpz','mpz' arguments");

    if (mpz_sgn(Pympz_AS_MPZ(other)) <= 0) {
        VALUE_ERROR("jacobi's y must be odd prime > 0");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    i = (long)(mpz_jacobi(Pympz_AS_MPZ(self), Pympz_AS_MPZ(other)));
    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(i);
}

static char doc_legendreg[]="\
legendre(x,y): returns the Legendre symbol (x|y) (y should be odd\n\
and must be positive); x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_legendre(PyObject *self, PyObject *args)
{
    PyObject *other;
    long i;

    PARSE_TWO_MPZ(other, "legendre() requires 'mpz','mpz' arguments");

    if (mpz_sgn(Pympz_AS_MPZ(other)) <= 0) {
        VALUE_ERROR("legendre's y must be odd and > 0");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    i = (long) mpz_legendre(Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(i);
}

static char doc_kroneckerg[]="\
kronecker(x,y): returns the Kronecker-Jacobi symbol (x|y).\n\
x and y must be mpz, or else get coerced to mpz (at least\n\
one of x and y, however, must also fit in a plain int).\n\
";
static PyObject *
Pympz_kronecker(PyObject *self, PyObject *args)
{
    PyObject *other;
    int ires;

    PARSE_TWO_MPZ(other, "kronecker() requires 'mpz','mpz' arguments");

    if (mpz_fits_ulong_p(Pympz_AS_MPZ(self))) {
        ires = mpz_ui_kronecker(
                mpz_get_ui(Pympz_AS_MPZ(self)),Pympz_AS_MPZ(other));
    }
    else if (mpz_fits_ulong_p(Pympz_AS_MPZ(other))) {
        ires = mpz_kronecker_ui(
                Pympz_AS_MPZ(self),mpz_get_ui(Pympz_AS_MPZ(other)));
    }
    else if (mpz_fits_slong_p(Pympz_AS_MPZ(self))) {
        ires = mpz_si_kronecker(
                mpz_get_si(Pympz_AS_MPZ(self)),Pympz_AS_MPZ(other));
    }
    else if (mpz_fits_slong_p(Pympz_AS_MPZ(other))) {
        ires = mpz_kronecker_si(
                Pympz_AS_MPZ(self),mpz_get_si(Pympz_AS_MPZ(other)));
    }
    else {
        VALUE_ERROR("Either arg in Kronecker must fit in an int");
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }
    Py_DECREF(self);
    Py_DECREF(other);
    return PyIntOrLong_FromLong(ires);
}

static char doc_is_evenm[]="\
x.is_even(): returns True if x is even, False otherwise.\n\
";
static char doc_is_eveng[]="\
is_even(x): returns True if x is even, False otherwise.\n\
";
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
    if(res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static char doc_is_oddm[]="\
x.is_odd(): returns True if x is odd, False otherwise.\n\
";
static char doc_is_oddg[]="\
is_odd(x): returns True if x is odd, False otherwise.\n\
";
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
    if(res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static char doc_make_mpzm[]="\
x.make_mpz(): returns x converted to an 'mpz'.\n\
NOTE: Optimized for speed so x is set to 0!.\n\
";

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

static char doc_xmpz_copy[]="x.copy(): returns a copy of x.\n";

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
"Return a Python string by formatting 'x' using the format string\n"
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



