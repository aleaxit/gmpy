/* gmpy_mpz.c
 *
 * Functions that operate strictly on mpz or xmpz.
 *
 * This file should be considered part of gmpy.c
 */

/* produce digits for an mpz in requested base, default 10 */
static char doc_digitsm[]="\
x.digits([base]): returns Python string representing x in the\n\
given base (2 to 36, default 10 if omitted or 0); leading '-'\n\
is present if x<0, but no leading '+' if x>=0.\n\
";
static char doc_digitsg[]="\
digits(x[,base]): returns Python string representing x in the\n\
given base (2 to 36, default 10 if omitted or 0); leading '-'\n\
present if x<0, but no leading '+' if x>=0. x must be an mpz,\n\
or else gets coerced into one.\n\
";
static PyObject *
Pympz_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "digits() requires 'mpz',['int'] arguments");
    result = Pympz_ascii((PympzObject*)self, base, 0);
    Py_DECREF(self);
    return result;
}

static PyObject *
Pyxmpz_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "digits() requires 'xmpz',['int'] arguments");
    result = Pyxmpz_ascii((PyxmpzObject*)self, base, 0);
    Py_DECREF(self);
    return result;
}

/* return number-of-digits for an mpz in requested base, default 10 */
static char doc_numdigitsm[]="\
x.numdigits([base]): returns length of string representing x in\n\
the given base (2 to 36, default 10 if omitted or 0); the value\n\
returned may sometimes be 1 more than necessary; no provision\n\
for any 'sign' character, nor leading '0' or '0x' decoration,\n\
is made in the returned length.\n\
";
static char doc_numdigitsg[]="\
numdigits(x[,base]): returns length of string representing x in\n\
the given base (2 to 36, default 10 if omitted or 0); the value\n\
returned may sometimes be 1 more than necessary; no provision\n\
for any 'sign' character, nor leading '0' or '0x' decoration,\n\
is made in the returned length.  x must be an mpz, or else gets\n\
coerced into one.\n\
";
static PyObject *
Pympz_numdigits(PyObject *self, PyObject *args)
{
    int base = 10;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&base,
            "numdigits() requires 'mpz',['int'] arguments");
    if(base==0) base=10;
    if((base < 2) || (base > 36)) {
        VALUE_ERROR("base must be either 0 or in the interval 2 ... 36");
        Py_DECREF(self);
        return NULL;
    }
    result = PyIntOrLong_FromSize_t(mpz_sizeinbase(Pympz_AS_MPZ(self), base));
    Py_DECREF(self);
    return result;
}

static char doc_bit_lengthm[]="\
x.bit_length(): returns length of string representing x in base 2\n\
";
static char doc_bit_lengthg[]="\
bit_length(x): returns length of string representing x in base 2\n\
";
static PyObject *
Pympz_bit_length(PyObject *self, PyObject *other)
{
    size_t i = 0;
    PympzObject* tempx;

    if(self && (CHECK_MPZANY(self))) {
        if (mpz_size(Pympz_AS_MPZ(self)))
            i = mpz_sizeinbase(Pympz_AS_MPZ(self), 2);
    } else if(CHECK_MPZANY(other)) {
        if (mpz_size(Pympz_AS_MPZ(other)))
            i = mpz_sizeinbase(Pympz_AS_MPZ(other), 2);
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("bit_length() requires 'mpz' argument");
            return NULL;
        } else {
            if (mpz_size(Pympz_AS_MPZ(tempx)))
                i = mpz_sizeinbase(Pympz_AS_MPZ(tempx), 2);
            Py_DECREF((PyObject*)tempx);
        }
    }
    return PyIntOrLong_FromSize_t(i);
}

static char doc_bit_maskg[]="\
bit_mask(x): create an 'mpz' or 'xmpz' exactly x bits in length, all bits\n\
are set\n\
";
static PyObject *
Pympz_bit_mask(PyObject *self, PyObject *other)
{
    long i = 0;
    PyObject* result;

    CREATE0_ONE_MPZANY(result);

    i = clong_From_Integer(other);
    if(i == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_mask() requires 'int' argument");
        return NULL;
    }

    if(i < 0) {
        VALUE_ERROR("mask length must be >= 0");
        return NULL;
    }

    mpz_set_ui(Pympz_AS_MPZ(result), 1);
    mpz_mul_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result), i);
    mpz_sub_ui(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result), 1);

    return result;
}

/* return scan0/scan1 for an mpz */
static char doc_bit_scan0m[]="\
x.bit_scan0(n=0): returns the bit-index of the first 0-bit of x (that\n\
is at least n); n must be an ordinary Python int, >=0.  If no more\n\
0-bits are in x at or above bit-index n (which can only happen for\n\
x<0, notionally extended with infinite 1-bits), None is returned.\n\
";
static char doc_bit_scan0g[]="\
bit_scan0(x, n=0): returns the bit-index of the first 0-bit of x (that\n\
is at least n); n must be an ordinary Python int, >=0.  If no more\n\
0-bits are in x at or above bit-index n (which can only happen for\n\
x<0, notionally extended with infinite 1-bits), None is returned.\n\
x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_bit_scan0(PyObject *self, PyObject *args)
{
    long starting_bit = 0;
    size_t maxbit;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&starting_bit,
            "bit_scan0() requires 'mpz',['int'] arguments");

    if(starting_bit < 0) {
        VALUE_ERROR("starting bit must be >= 0");
        Py_DECREF(self);
        return NULL;
    }
    maxbit = mpz_sizeinbase(Pympz_AS_MPZ(self), 2);
    if(starting_bit > maxbit) {
        if(mpz_sgn(Pympz_AS_MPZ(self))<0) {
            Py_DECREF(self);
            Py_RETURN_NONE;
        } else {
            result = PyIntOrLong_FromLong(starting_bit);
        }
    } else {
        result = PyIntOrLong_FromLong(mpz_scan0(Pympz_AS_MPZ(self), starting_bit));
    }
    Py_DECREF(self);
    return result;
}

static char doc_bit_scan1m[]="\
x.bit_scan1(n=0): returns the bit-index of the first 1-bit of x (that\n\
is at least n); n must be an ordinary Python int, >=0.  If no more\n\
1-bits are in x at or above bit-index n (which can only happen for\n\
x>=0, notionally extended with infinite 0-bits), None is returned.\n\
";
static char doc_bit_scan1g[]="\
bit_scan1(x, n=0): returns the bit-index of the first 1-bit of x (that\n\
is at least n); n must be an ordinary Python int, >=0.  If no more\n\
1-bits are in x at or above bit-index n (which can only happen for\n\
x>=0, notionally extended with infinite 0-bits), None is returned.\n\
x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_bit_scan1(PyObject *self, PyObject *args)
{
    long starting_bit = 0;
    long maxbit;
    PyObject *result;

    PARSE_ONE_MPZ_OPT_CLONG(&starting_bit,
            "bit_scan1() requires 'mpz',['int'] arguments");

    if(starting_bit < 0) {
        VALUE_ERROR("starting bit must be >= 0");
        Py_DECREF(self);
        return NULL;
    }
    maxbit = mpz_sizeinbase(Pympz_AS_MPZ(self), 2);
    if(starting_bit >= maxbit) {
        if(mpz_sgn(Pympz_AS_MPZ(self))>=0) {
            Py_DECREF(self);
            Py_RETURN_NONE;
        } else {
            result = PyIntOrLong_FromLong(starting_bit);
        }
    } else {
        result = PyIntOrLong_FromLong(mpz_scan1(Pympz_AS_MPZ(self), starting_bit));
    }
    Py_DECREF(self);
    return result;
}

/* return population-count (# of 1-bits) for an mpz */
static char doc_popcountm[]="\
x.popcount(): returns the number of 1-bits set in x; note that\n\
this is 'infinite' if x<0, and in that case, -1 is returned.\n\
";
static char doc_popcountg[]="\
popcount(x): returns the number of 1-bits set in x; note that\n\
this is 'infinite' if x<0, and in that case, -1 is returned.\n\
x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_popcount(PyObject *self, PyObject *other)
{
    long temp;
    PympzObject *tempx;

    if(self && (CHECK_MPZANY(self)))
        return PyIntOrLong_FromLong(mpz_popcount(Pympz_AS_MPZ(self)));
    else if(CHECK_MPZANY(other))
        return PyIntOrLong_FromLong(mpz_popcount(Pympz_AS_MPZ(other)));
    else {
        if((tempx = Pympz_From_Integer(other))) {
            temp = mpz_popcount(Pympz_AS_MPZ(tempx));
            Py_DECREF((PyObject*)tempx);
            return PyIntOrLong_FromLong(temp);
        } else {
            TYPE_ERROR("popcount() requires 'mpz' argument");
            return NULL;
        }
    }
}

/* get & return one bit from an mpz */
static char doc_bit_testm[]="\
x.bit_test(n): return the value of the nth bit of x.\n\
";
static char doc_bit_testg[]="\
bit_test(x,n): return the value of the nth bit of x.\n\
";
static PyObject *
Pygmpy_bit_test(PyObject *self, PyObject *args)
{
    long bit_index;
    int temp;
    PyObject *x;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    if(CHECK_MPZANY(x)) {
        temp = mpz_tstbit(Pympz_AS_MPZ(x), bit_index);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
            return NULL;
        }
        temp = mpz_tstbit(Pympz_AS_MPZ(tempx), bit_index);
        Py_DECREF(tempx);
    }
    if(temp)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static PyObject *
Pympz_bit_test(PyObject *self, PyObject *other)
{
    long bit_index;

    bit_index = clong_From_Integer(other);
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_test() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if(mpz_tstbit(Pympz_AS_MPZ(self), bit_index))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static char doc_bit_clearm[]="\
x.bit_clear(n): clear the nth bit of x. If x is an xmpz, x is mutated.\n\
If x is an mpz, a new object is returned.\n\
";
static char doc_bit_clearg[]="\
bit_clear(x,n): clear the nth bit of x. A new object is always returned.\n\
";

static PyObject *
Pygmpy_bit_clear(PyObject *self, PyObject *args)
{
    long bit_index;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);
    if(CHECK_MPZANY(x)) {
        mpz_set(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x));
        mpz_clrbit(Pympz_AS_MPZ(result), bit_index);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_swap(Pympz_AS_MPZ(result), Pympz_AS_MPZ(tempx));
        Py_DECREF((PyObject*)tempx);
        mpz_clrbit(Pympz_AS_MPZ(result), bit_index);
    }
    return result;
}

static PyObject *
Pympz_bit_clear(PyObject *self, PyObject *other)
{
    long bit_index;
    PyObject *result;

    bit_index = clong_From_Integer(other);
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_clrbit(Pympz_AS_MPZ(self), bit_index);
        Py_RETURN_NONE;
    } else {
        CREATE0_ONE_MPZANY(result);
        mpz_set(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self));
        mpz_clrbit(Pympz_AS_MPZ(result), bit_index);
        return result;
    }
}

static char doc_bit_setm[]="\
x.bit_set(n): set the nth bit of x. If x is an xmpz, x is mutated.\n\
If x is an mpz, a new object is returned.\n\
";
static char doc_bit_setg[]="\
bit_set(x,n): set the nth bit of x. A new object is always returned.\n\
";

static PyObject *
Pygmpy_bit_set(PyObject *self, PyObject *args)
{
    long bit_index;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);
    if(CHECK_MPZANY(x)) {
        mpz_set(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x));
        mpz_setbit(Pympz_AS_MPZ(result), bit_index);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_swap(Pympz_AS_MPZ(result), Pympz_AS_MPZ(tempx));
        Py_DECREF((PyObject*)tempx);
        mpz_setbit(Pympz_AS_MPZ(result), bit_index);
    }
    return result;
}

static PyObject *
Pympz_bit_set(PyObject *self, PyObject *other)
{
    long bit_index;
    PyObject *result;

    bit_index = clong_From_Integer(other);
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_setbit(Pympz_AS_MPZ(self), bit_index);
        Py_RETURN_NONE;
    } else {
        CREATE0_ONE_MPZANY(result);
        mpz_set(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self));
        mpz_setbit(Pympz_AS_MPZ(result), bit_index);
        return result;
    }
}

static char doc_bit_flipm[]="\
x.bit_flip(n): complements the nth bit of x. If x is an xmpz,\n\
x is mutated. If x is an mpz, a new object is returned.\n\
";
static char doc_bit_flipg[]="\
bit_flip(x,n): complements the nth bit of x. A new object is \n\
always returned.\n\
";

static PyObject *
Pygmpy_bit_flip(PyObject *self, PyObject *args)
{
    long bit_index;
    PyObject *x, *result;
    PympzObject *tempx;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    bit_index = clong_From_Integer(PyTuple_GET_ITEM(args, 1));
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    CREATE1_ONE_MPZANY(x, result);
    if(CHECK_MPZANY(x)) {
        mpz_set(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x));
        mpz_combit(Pympz_AS_MPZ(result), bit_index);
    } else {
        if(!(tempx = Pympz_From_Integer(x))) {
            TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_swap(Pympz_AS_MPZ(result), Pympz_AS_MPZ(tempx));
        Py_DECREF((PyObject*)tempx);
        mpz_combit(Pympz_AS_MPZ(result), bit_index);
    }
    return result;
}

static PyObject *
Pympz_bit_flip(PyObject *self, PyObject *other)
{
    long bit_index;
    PyObject *result;

    bit_index = clong_From_Integer(other);
    if(bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    if(bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        mpz_combit(Pympz_AS_MPZ(self), bit_index);
        Py_RETURN_NONE;
    } else {
        CREATE0_ONE_MPZANY(result);
        mpz_set(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self));
        mpz_combit(Pympz_AS_MPZ(result), bit_index);
        return result;
    }
}

/* return nth-root of an mpz (in a 2-el tuple: 2nd is int, non-0 iff exact) */
static char doc_rootm[]="\
x.root(n): returns a 2-element tuple (y,m), such that y is the\n\
(possibly truncated) n-th root of x; m, an ordinary Python int,\n\
is 1 if the root is exact (x==y**n), else 0.  n must be an ordinary\n\
Python int, >=0.\n\
";
static char doc_rootg[]="\
root(x,n): returns a 2-element tuple (y,m), such that y is the\n\
(possibly truncated) n-th root of x; m, an ordinary Python int,\n\
is 1 if the root is exact (x==y**n), else 0.  n must be an ordinary\n\
Python int, >=0. x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_root(PyObject *self, PyObject *args)
{
    long n;
    int exact;
    PyObject *s;
    PyObject *result;

    PARSE_ONE_MPZ_REQ_CLONG(&n,
                            "root() requires 'mpz','int' arguments");

    if(n <= 0) {
        VALUE_ERROR("n must be > 0");
        Py_DECREF(self);
        return NULL;
    } else if(n>1) {
        if(mpz_sgn(Pympz_AS_MPZ(self))<0) {
            VALUE_ERROR("root of negative number");
            Py_DECREF(self);
            return NULL;
        }
    }
    if(Pyxmpz_Check(self)) {
        s = (PyObject*)Pyxmpz_new();
    } else {
        s = (PyObject*)Pympz_new();
    }
    result = PyTuple_New(2);
    if(!s || !result) {
        Py_DECREF(self);
        Py_XDECREF(s);
        Py_XDECREF(result);
        return NULL;
    }
    exact = mpz_root(Pympz_AS_MPZ(s), Pympz_AS_MPZ(self), n);
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, s);
    PyTuple_SET_ITEM(result, 1, (PyObject*)PyIntOrLong_FromLong(exact));
    return result;
}

static char doc_rootremm[]="\
x.rootrem(n): returns a 2-element tuple (y,r), such that y is the\n\
(possibly truncated) n-th root of x; r is the remainder. n must be an\n\
ordinary Python int, >=0.\n\
";
static char doc_rootremg[]="\
rootrem(x,n): returns a 2-element tuple (y,r), such that y is the\n\
(possibly truncated) n-th root of x; r is the remainder. n must be an\n\
ordinary Python int, >=0. x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_rootrem(PyObject *self, PyObject *args)
{
    long n;
    PyObject *y = 0, *r = 0, *result;

    PARSE_ONE_MPZ_REQ_CLONG(&n,
            "rootrem() requires 'mpz','int' arguments");

    if(n <= 0) {
        VALUE_ERROR("n must be > 0");
        Py_DECREF(self);
        return NULL;
    } else if(n>1) {
        if(mpz_sgn(Pympz_AS_MPZ(self))<0) {
            VALUE_ERROR("root of negative number");
            Py_DECREF(self);
            return NULL;
        }
    }
    if(Pyxmpz_Check(self)) {
        y = (PyObject*)Pyxmpz_new();
        r = (PyObject*)Pyxmpz_new();
    } else {
        y = (PyObject*)Pympz_new();
        r = (PyObject*)Pympz_new();
    }
    result = PyTuple_New(2);
    if(!y || !r || !result) {
        Py_DECREF(self);
        Py_XDECREF(result);
        Py_XDECREF(y);
        Py_XDECREF(r);
        return NULL;
    }
    mpz_rootrem(Pympz_AS_MPZ(y), Pympz_AS_MPZ(r), Pympz_AS_MPZ(self), n);
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, y);
    PyTuple_SET_ITEM(result, 1, r);
    return result;
}

static char doc_signm[]="\
x.sign(): returns -1, 0, or +1, if x is negative, 0, positive.\n\
";
static char doc_signg[]="\
sign(x): returns -1, 0, or +1, if x is negative, 0, positive;\n\
x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_sign(PyObject *self, PyObject *other)
{
    long res;
    PympzObject* tempx;

    if(self && (CHECK_MPZANY(self))) {
        res = mpz_sgn(Pympz_AS_MPZ(self));
    } else if(CHECK_MPZANY(other)) {
        res = mpz_sgn(Pympz_AS_MPZ(other));
    } else {
        tempx = Pympz_From_Integer(other);
        if(tempx) {
            res = mpz_sgn(Pympz_AS_MPZ(tempx));
            Py_DECREF((PyObject*)tempx);
        } else {
            TYPE_ERROR("sign() requires 'mpz' argument");
            return NULL;
        }
    }
    return PyIntOrLong_FromLong(res);
}

static PyObject *
Pympz_abs(PympzObject *x)
{
    PympzObject *result;
    if(!(result = Pympz_new()))
        return NULL;
    mpz_abs(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x));
    return (PyObject *)result;
}

static PyObject *
Pyxmpz_abs(PyxmpzObject *x)
{
    mpz_abs(Pympz_AS_MPZ(x), Pympz_AS_MPZ(x));
    Py_RETURN_NONE;
}

static PyObject *
Pympz_neg(PympzObject *x)
{
    PympzObject *result;
    if(!(result = Pympz_new()))
        return NULL;
    mpz_neg(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x));
    return (PyObject *)result;
}

static PyObject *
Pyxmpz_neg(PyxmpzObject *x)
{
    mpz_neg(Pympz_AS_MPZ(x), Pympz_AS_MPZ(x));
    Py_RETURN_NONE;
}

static PyObject *
Pympz_pos(PympzObject *x)
{
    Py_INCREF((PyObject*)x);
    return (PyObject *) x;
}

static PyObject *
Pyxmpz_pos(PyxmpzObject *x)
{
    Py_RETURN_NONE;
}

/* Pympz_pow is called by Pympany_pow after verifying that all the
 * arguments are integers, but not necessarily mpz.
 */

static PyObject *
Pympz_pow(PyObject *b, PyObject *e, PyObject *m)
{
    PyObject *result;
    PympzObject *tempb, *tempe, *tempm = 0;

    CREATE1_ONE_MPZANY(b, result);

    tempb = Pympz_From_Integer(b);
    tempe = Pympz_From_Integer(e);

    /* m will either be a number or Py_None. */
    if(m != Py_None) {
        tempm = Pympz_From_Integer(m);
    }

    if(!tempb || !tempe || (!tempm && (m != Py_None))) {
        PyErr_Clear();
        Py_XDECREF((PyObject*)tempb);
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempm);
        Py_DECREF(result);
        Py_RETURN_NOTIMPLEMENTED;
    }

    if(m == Py_None) {
        /* When no modulo is present, the exponent must fit in C ulong
         * the exponent must be positive.
         */
        unsigned long el;
        if(mpz_sgn(tempe->z) < 0) {
            VALUE_ERROR("pow() exponent cannot be negative");
            goto err;
        }
        if(!mpz_fits_ulong_p(tempe->z)) {
            VALUE_ERROR("pow() outrageous exponent");
            goto err;
        }
        el = mpz_get_ui(tempe->z);
        mpz_pow_ui(Pympz_AS_MPZ(result), tempb->z, el);
    } else { /* Modulo exponentiation */
        int sign;
        mpz_t mm, base, exp;

        sign = mpz_sgn(tempm->z);
        if(sign == 0) {
            VALUE_ERROR("pow() 3rd argument cannot be 0");
            goto err;
        }
        mpz_inoc(mm);
        mpz_abs(mm, tempm->z);
        /* A negative exponent is allowed if inverse exists. */
        if(mpz_sgn(tempe->z) < 0) {
            mpz_inoc(base);
            if(!mpz_invert(base, tempb->z, mm)) {
                VALUE_ERROR("pow() base not invertible");
                mpz_cloc(base);
                mpz_cloc(mm);
                goto err;
            } else {
                mpz_inoc(exp);
                mpz_abs(exp, tempe->z);
            }
            mpz_powm(Pympz_AS_MPZ(result), base, exp, mm);
            mpz_cloc(base);
            mpz_cloc(exp);
        } else {
            mpz_powm(Pympz_AS_MPZ(result), tempb->z, tempe->z, mm);
        }
        mpz_cloc(mm);
        if((sign<0) && (mpz_sgn(Pympz_AS_MPZ(result)) > 0)) {
        /* Python uses a rather peculiar convention for negative modulos
         * If the modulo is negative, result should be in the interval
         * m < r <= 0 .
         */
            mpz_add(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result), tempm->z);
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
    Py_DECREF(result);
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
    if(!(result = Pympz_new()))
        return NULL;
    mpz_com(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x));
    return (PyObject *)result;
}

static PyObject *
Pyxmpz_com(PyxmpzObject *x)
{
    mpz_com(Pympz_AS_MPZ(x), Pympz_AS_MPZ(x));
    Py_RETURN_NONE;
}

#define MPZ_BINOP(NAME) \
static PyObject * \
Py##NAME(PyObject *a, PyObject *b) \
{ \
    PyObject *r = 0; \
    PympzObject *tempz = 0; \
    CREATE2_ONE_MPZANY(a, b, r);\
    if(CHECK_MPZANY(a)) {\
        if(CHECK_MPZANY(b)) {\
            NAME(Pympz_AS_MPZ(r), Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));\
        } else {\
            if(!(tempz = Pympz_From_Integer(b))) {\
                Py_DECREF(r);\
                return NULL;\
            }\
            NAME(Pympz_AS_MPZ(r), Pympz_AS_MPZ(a), tempz->z);\
            Py_DECREF((PyObject*)tempz);\
        }\
    } else if(CHECK_MPZANY(b)) {\
        if(CHECK_MPZANY(a)) {\
          NAME(Pympz_AS_MPZ(r), Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));\
        } else {\
            if(!(tempz = Pympz_From_Integer(a))) {\
                Py_DECREF(r);\
                return NULL;\
            }\
            NAME(Pympz_AS_MPZ(r), tempz->z, Pympz_AS_MPZ(b));\
            Py_DECREF((PyObject*)tempz);\
        }\
    } else {\
        Py_DECREF(r);\
        Py_RETURN_NOTIMPLEMENTED;\
    }\
    return r;\
}

MPZ_BINOP(mpz_and)
MPZ_BINOP(mpz_ior)
MPZ_BINOP(mpz_xor)

static PyObject *
Pympz_rshift(PyObject *a, PyObject *b)
{
    PyObject *result;
    long count;
    int overflow;

    CREATE1_ONE_MPZANY(a, result);

    /* Try to make mpz >> Python int/long as fast as possible. */
    if(CHECK_MPZANY(a)) {
        if(PyIntOrLong_Check(b)) {
            count = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF(result);
                return NULL;
            } else if(count >= 0) {
                mpz_fdiv_q_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(a), count);
                return result;
            } else {
                VALUE_ERROR("negative shift count");
                Py_DECREF(result);
                return NULL;
            }
        }
    }
    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
Pympz_lshift(PyObject *a, PyObject *b)
{
    PyObject *result;
    long count;
    int overflow;

    CREATE1_ONE_MPZANY(a, result);

    /* Try to make mpz >> Python int/long as fast as possible. */
    if(CHECK_MPZANY(a)) {
        if(PyIntOrLong_Check(b)) {
            count = PyLong_AsLongAndOverflow(b, &overflow);
            if(overflow) {
                VALUE_ERROR("outrageous shift count");
                Py_DECREF(result);
                return NULL;
            } else if(count >= 0) {
                mpz_mul_2exp(Pympz_AS_MPZ(result), Pympz_AS_MPZ(a), count);
                return result;
            } else {
                VALUE_ERROR("negative shift count");
                Py_DECREF(result);
                return NULL;
            }
        }
    }
    Py_RETURN_NOTIMPLEMENTED;
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

static long
Pympz_hash(PympzObject *self)
{
    if(self->hash_cache != -1)
        return self->hash_cache;
#ifdef _PyHASH_MODULUS
    long hash = (long)mpz_tdiv_ui((self->z), _PyHASH_MODULUS);
    if(mpz_sgn(self->z)<0)
        hash = -hash;
    if(hash==-1) hash = -2;
    return (self->hash_cache = hash);
#else
    return (self->hash_cache = mpz_pythonhash(Pympz_AS_MPZ(self)));
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
    PyObject *result, *a, *b;
    PympzObject *tempa, *tempb;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcd() requires 'mpz','mpz' arguments");
        return NULL;
    }

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);
    CREATE2_ONE_MPZANY(a, b, result);

    if(CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_gcd(Pympz_AS_MPZ(result), Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
    } else {
        tempa = Pympz_From_Integer(a);
        tempb = Pympz_From_Integer(b);
        if(!tempa || !tempb) {
            TYPE_ERROR("gcd() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF(result);
            return NULL;
        }
        mpz_gcd(Pympz_AS_MPZ(result), tempa->z, tempb->z);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    return result;
}

static char doc_lcm[]="\
lcm(a,b): returns the lowest common multiple of numbers a and b\n\
(which must be mpz objects, or else get coerced to mpz)\n\
";
static PyObject *
Pygmpy_lcm(PyObject *self, PyObject *args)
{
    PyObject *result, *a, *b;
    PympzObject *tempa, *tempb;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("lcm() requires 'mpz','mpz' arguments");
        return NULL;
    }

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);
    CREATE2_ONE_MPZANY(a, b, result);

    if(CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_lcm(Pympz_AS_MPZ(result), Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
    } else {
        tempa = Pympz_From_Integer(a);
        tempb = Pympz_From_Integer(b);
        if(!tempa || !tempb) {
            TYPE_ERROR("lcm() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF(result);
            return NULL;
        }
        mpz_lcm(Pympz_AS_MPZ(result), tempa->z, tempb->z);
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    return result;
}

static char doc_gcdext[]="\
gcdext(a,b): returns a 3-element tuple (g,s,t) such that\n\
    g==gcd(a,b) and g == a*s + b*t\n\
(a and b must be mpz objects, or else get coerced to mpz)\n\
";
static PyObject *
Pygmpy_gcdext(PyObject *self, PyObject *args)
{
    PyObject *g, *s, *t, *a, *b, *result;
    PympzObject *tempa, *tempb;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
        return NULL;
    }

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);
    CREATE2_THREE_MPZANY(a, b, g, s, t, result);

    if(CHECK_MPZANY(a) && CHECK_MPZANY(b)) {
        mpz_gcdext(Pympz_AS_MPZ(g), Pympz_AS_MPZ(s), Pympz_AS_MPZ(t),
                    Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
    } else {
        tempa = Pympz_From_Integer(a);
        tempb = Pympz_From_Integer(b);
        if(!tempa || !tempb) {
            TYPE_ERROR("gcdext() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempa);
            Py_XDECREF((PyObject*)tempb);
            Py_DECREF(g);
            Py_DECREF(s);
            Py_DECREF(t);
            Py_DECREF(result);
            return NULL;
        }
        mpz_gcdext(Pympz_AS_MPZ(g), Pympz_AS_MPZ(s), Pympz_AS_MPZ(t),
                    Pympz_AS_MPZ(a), Pympz_AS_MPZ(b));
        Py_DECREF((PyObject*)tempa);
        Py_DECREF((PyObject*)tempb);
    }
    PyTuple_SET_ITEM(result, 0, g);
    PyTuple_SET_ITEM(result, 1, s);
    PyTuple_SET_ITEM(result, 2, t);
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
    PyObject *a, *b, *m, *result;
    PympzObject *num, *den, *mod;
    mpz_t numz, denz, modz, gcdz;
    int ok;

    if(PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        return NULL;
    }

    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);
    m = PyTuple_GET_ITEM(args, 2);

    if(Pyxmpz_Check(a) && Pyxmpz_Check(b) && Pyxmpz_Check(m)) {
        result = (PyObject*)Pyxmpz_new();
    } else if(Pympz_Check(a) && Pympz_Check(b) && Pympz_Check(m)) {
        result = (PyObject*)Pympz_new();
    } else if(options.prefer_mutable) {
        result = (PyObject*)Pyxmpz_new();
    } else {
        result = (PyObject*)Pympz_new();
    }
    if(!result) return NULL;

    num = Pympz_From_Integer(a);
    den = Pympz_From_Integer(b);
    mod = Pympz_From_Integer(m);
    if(!num || !den || !mod) {
        TYPE_ERROR("divm() requires 'mpz','mpz','mpz' arguments");
        Py_XDECREF((PyObject*)num);
        Py_XDECREF((PyObject*)den);
        Py_XDECREF((PyObject*)mod);
        Py_DECREF(result);
        return NULL;
    }

    mpz_inoc(numz);
    mpz_inoc(denz);
    mpz_inoc(modz);
    mpz_set(numz, num->z);
    mpz_set(denz, den->z);
    mpz_set(modz, mod->z);

    if(mpz_invert(Pympz_AS_MPZ(result), denz, modz)) { /* inverse exists */
        ok = 1;
    } else {
        /* last-ditch attempt: do num, den AND mod have a gcd>1 ? */
        mpz_inoc(gcdz);
        mpz_gcd(gcdz, numz, denz);
        mpz_gcd(gcdz, gcdz, modz);
        mpz_divexact(numz, numz, gcdz);
        mpz_divexact(denz, denz, gcdz);
        mpz_divexact(modz, modz, gcdz);
        mpz_cloc(gcdz);
        ok = mpz_invert(Pympz_AS_MPZ(result), denz, modz);
    }

    if (ok) {
        mpz_mul(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result), numz);
        mpz_mod(Pympz_AS_MPZ(result), Pympz_AS_MPZ(result), modz);
        mpz_cloc(numz);
        mpz_cloc(denz);
        mpz_cloc(modz);
        Py_DECREF((PyObject*)num);
        Py_DECREF((PyObject*)den);
        Py_DECREF((PyObject*)mod);
        return result;
    } else {
        ZERO_ERROR("not invertible");
        mpz_cloc(numz);
        mpz_cloc(denz);
        mpz_cloc(modz);
        Py_DECREF((PyObject*)num);
        Py_DECREF((PyObject*)den);
        Py_DECREF((PyObject*)mod);
        Py_DECREF(result);
        return NULL;
    }
}

static char doc_fac[]="\
fac(n): returns the factorial of n; takes O(n) time; n must be\n\
an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_fac(PyObject *self, PyObject *other)
{
    PyObject *result;
    long n;

    n = clong_From_Integer(other);
    if((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("fac() requires 'int' argument");
        return NULL;
    } else if(n < 0) {
        VALUE_ERROR("factorial of negative number");
        return NULL;
    } else {
        CREATE0_ONE_MPZANY(result);
        mpz_fac_ui(Pympz_AS_MPZ(result), n);
    }
    return result;
}

static char doc_fib[]="\
fib(n): returns the n-th Fibonacci number; takes O(n) time; n must be\n\
an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_fib(PyObject *self, PyObject *other)
{
    PyObject *result;
    long n;

    n = clong_From_Integer(other);
    if((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("fib() requires 'int' argument");
        return NULL;
    } else if(n < 0) {
        VALUE_ERROR("Fibonacci of negative number");
        return NULL;
    } else {
        CREATE0_ONE_MPZANY(result);
        mpz_fib_ui(Pympz_AS_MPZ(result), n);
    }
    return result;
}

static char doc_fib2[]="\
fib2(n): returns the n-th and n+1-th Fibonacci numbers; takes O(n) time;\n\
n must be an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_fib2(PyObject *self, PyObject *other)
{
    PyObject *fib1, *fib2, *result;
    long n;

    n = clong_From_Integer(other);
    if((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("fib2() requires 'int' argument");
        return NULL;
    } else if(n < 0) {
        VALUE_ERROR("Fibonacci of negative number");
        return NULL;
    } else {
        CREATE0_TWO_MPZANY_TUPLE(fib1, fib2, result);
        mpz_fib2_ui(Pympz_AS_MPZ(fib1), Pympz_AS_MPZ(fib2), n);
    }
    PyTuple_SET_ITEM(result, 0, fib1);
    PyTuple_SET_ITEM(result, 1, fib2);
    return result;
}

static char doc_lucas[]="\
lucas(n): returns the n-th Lucas number; takes O(n) time; n must be\n\
an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_lucas(PyObject *self, PyObject *other)
{
    PyObject *result;
    long n;

    n = clong_From_Integer(other);
    if((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("luc() requires 'int' argument");
        return NULL;
    } else if(n < 0) {
        VALUE_ERROR("Lucas of negative number");
        return NULL;
    } else {
        CREATE0_ONE_MPZANY(result);
        mpz_lucnum_ui(Pympz_AS_MPZ(result), n);
    }
    return result;
}

static char doc_lucas2[]="\
lucas2(n): returns the n-1 and n-th Lucas number; takes O(n) time; n must\n\
be an ordinary Python int, >=0.\n\
";
static PyObject *
Pygmpy_lucas2(PyObject *self, PyObject *other)
{
    PyObject *luc1, *luc2, *result;
    long n;

    n = clong_From_Integer(other);
    if((n == -1) && PyErr_Occurred()) {
        TYPE_ERROR("luc2() requires 'int' argument");
        return NULL;
    } else if(n < 0) {
        VALUE_ERROR("Lucas of negative number");
        return NULL;
    } else {
        CREATE0_TWO_MPZANY_TUPLE(luc1, luc2, result);
        mpz_fib2_ui(Pympz_AS_MPZ(luc1), Pympz_AS_MPZ(luc2), n);
    }
    PyTuple_SET_ITEM(result, 0, luc1);
    PyTuple_SET_ITEM(result, 1, luc2);
    return result;
}

static char doc_bincoefm[]="\
x.bincoef(n): returns the 'binomial coefficient' that is 'x\n\
over n'; n is an ordinary Python int, >=0.\n\
";
static char doc_bincoefg[]="\
bincoef(x,n): returns the 'binomial coefficient' that is 'x\n\
over n'; n is an ordinary Python int, >=0; x must be an mpz,\n\
or else gets converted to one.\n\
";
static char doc_combm[]="\
x.comb(n): returns the 'number of combinations' of 'x things,\n\
taken n at a time'; n is an ordinary Python int, >=0.\n\
";
static char doc_combg[]="\
comb(x,n): returns the 'number of combinations' of 'x things,\n\
taken n at a time'; n is an ordinary Python int, >=0; x must be\n\
an mpz, or else gets converted to one.\n\
";
static PyObject *
Pympz_bincoef(PyObject *self, PyObject *args)
{
    PyObject *result;
    long k;

    PARSE_ONE_MPZ_REQ_CLONG(&k,
                            "bincoef() requires 'mpz','int' arguments");

    if(k < 0) {
        VALUE_ERROR("binomial coefficient with negative k");
        Py_DECREF(self);
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        result = (PyObject*)Pyxmpz_new();
    } else {
        result = (PyObject*)Pympz_new();
    }
    if(!(result)) {
        Py_DECREF(self);
        return NULL;
    }
    mpz_bin_ui(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), k);
    Py_DECREF(self);
    return result;
}

static char doc_sqrtm[]="\
x.sqrt(): returns the integer, truncated square root of x, i.e. the\n\
largest y such that x>=y*y. x must be >= 0.\n\
";
static char doc_sqrtg[]="\
sqrt(x): returns the integer, truncated square root of x, i.e. the\n\
largest y such that x>=y*y. x must be an mpz, or else gets coerced\n\
to one; further, x must be >= 0.\n\
";
static PyObject *
Pygmpy_sqrt(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(mpz_sgn(Pympz_AS_MPZ(other)) < 0) {
        VALUE_ERROR("sqrt of negative number");
        return NULL;
    }

    CREATE1_ONE_MPZANY(other, result);

    if(CHECK_MPZANY(other)) {
        mpz_sqrt(Pympz_AS_MPZ(result), Pympz_AS_MPZ(other));
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("sqrt requires 'mpz' argument");
            return NULL;
        }
        if(mpz_sgn(tempx->z) < 0) {
            VALUE_ERROR("sqrt of negative number");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        mpz_sqrt(Pympz_AS_MPZ(result), tempx->z);
        Py_DECREF(tempx);
    }
    return result;
}

static PyObject *
Pympz_sqrt(PyObject *self, PyObject *args)
{
    PyObject *result;

    if(mpz_sgn(Pympz_AS_MPZ(self)) < 0) {
        VALUE_ERROR("sqrt of negative number");
        return NULL;
    }

    if(!(result=(PyObject*)Pympz_new())) {
        return NULL;
    }

    mpz_sqrt(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self));
    return result;
}

static PyObject *
Pyxmpz_sqrt(PyObject *self, PyObject *args)
{
    if(mpz_sgn(Pympz_AS_MPZ(self)) < 0) {
        VALUE_ERROR("sqrt of negative number");
        return NULL;
    }

    mpz_sqrt(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self));
    Py_RETURN_NONE;
}

static char doc_sqrtremm[]="\
x.sqrtrem(): returns a 2-element tuple (s,t), such that\n\
s==x.sqrt() and x==s*s+t. x must be >= 0.\n\
";
static char doc_sqrtremg[]="\
sqrtrem(x): returns a 2-element tuple (s,t), such that\n\
s==sqrt(x) and x==s*s+t. x must be an mpz, or else gets\n\
coerced to one; further, x must be >= 0.\n\
";
static PyObject *
Pympz_sqrtrem(PyObject *self, PyObject *args)
{
    PyObject *root = 0, *rem = 0, *result = 0;

    PARSE_ONE_MPZ("sqrtrem() requires 'mpz' argument");

    if(mpz_sgn(Pympz_AS_MPZ(self)) < 0) {
        VALUE_ERROR("sqrt of negative number");
        Py_DECREF(self);
        return NULL;
    }

    if(Pyxmpz_Check(self)) {
        root = (PyObject*)Pyxmpz_new();
        rem = (PyObject*)Pyxmpz_new();
    } else {
        root = (PyObject*)Pympz_new();
        rem = (PyObject*)Pympz_new();
    }
    result = PyTuple_New(2);
    if(!root || !rem || !result) {
        Py_DECREF(self);
        Py_XDECREF(result);
        Py_XDECREF(root);
        Py_XDECREF(rem);
        return NULL;
    }

    mpz_sqrtrem(Pympz_AS_MPZ(root), Pympz_AS_MPZ(rem), Pympz_AS_MPZ(self));
    Py_DECREF(self);
    PyTuple_SET_ITEM(result, 0, (PyObject*)root);
    PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
    return result;
}

static char doc_removem[]="\
x.remove(f): returns a 2-element tuple (y,m) such that\n\
x==y*(f**m), and y%f==0; i.e., y is x with any factor f\n\
removed, and m (an ordinary Python int) is the multiplicity\n\
of the factor f in x (m=0, and y=x, unless x%f==0). f must\n\
be > 0.\n\
";
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
    unsigned long multiplicity;

    PARSE_TWO_MPZ(factor, "remove() requires 'mpz','mpz' arguments");

    if(mpz_sgn(Pympz_AS_MPZ(factor)) <= 0) {
        VALUE_ERROR("factor must be > 0");
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }

    if(!(result = Pympz_new())) {
        Py_DECREF(self);
        Py_DECREF(factor);
        return NULL;
    }
    multiplicity = mpz_remove(result->z, Pympz_AS_MPZ(self), Pympz_AS_MPZ(factor));
    Py_DECREF(self);
    Py_DECREF(factor);
    return Py_BuildValue("(Nk)", result, multiplicity);
}

static char doc_invertm[]="\
x.invert(m): returns the inverse of x modulo m, i.e., that y\n\
such that x*y==1 modulo m, or 0 if no such y exists.\n\
m must be an ordinary Python int, !=0.\n\
";
static char doc_invertg[]="\
invert(x,m): returns the inverse of x modulo m, i.e., that y\n\
such that x*y==1 modulo m, or 0 if no such y exists.\n\
m must be an ordinary Python int, !=0; x must be an mpz,\n\
or else gets converted to one.\n\
";
static PyObject *
Pygmpy_invert(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *tempx, *tempy;
    int success;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF(result);
            return NULL;
        }
        success = mpz_invert(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
        if(!success)
            mpz_set_ui(Pympz_AS_MPZ(result), 0);
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        success = mpz_invert(Pympz_AS_MPZ(result), tempx->z, tempy->z);
        if(!success)
            mpz_set_ui(Pympz_AS_MPZ(result), 0);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return result;
}

static PyObject *
Pympz_invert(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;
    int success;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("invert() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            success = mpz_invert(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            if(!success)
                mpz_set_ui(Pympz_AS_MPZ(self), 0);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            success = mpz_invert(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            if(!success)
                mpz_set_ui(Pympz_AS_MPZ(result), 0);
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("invert() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("invert() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            success = mpz_invert(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            if(!success)
                mpz_set_ui(Pympz_AS_MPZ(self), 0);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            success = mpz_invert(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            if(!success)
                mpz_set_ui(Pympz_AS_MPZ(result), 0);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
}

static char doc_hamdistm[]="\
x.hamdist(y): returns the Hamming distance (number of bit-positions\n\
where the bits differ) between x and y.  y must be an mpz, or else\n\
gets coerced to one.\n\
";
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

    result = PyIntOrLong_FromLong(
            mpz_hamdist(Pympz_AS_MPZ(self),Pympz_AS_MPZ(other)));
    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

static char doc_divexactm[]="\
x.divexact(y): returns the quotient of x divided by y. Faster than\n\
standard division but requires the remainder is zero!  y must be an\n\
mpz, or else gets coerced to one.\n\
";
static char doc_divexactg[]="\
divexact(x,y): returns the quotient of x divided by y. Faster than\n\
standard division but requires the remainder is zero!  x and y must\n\
be mpz, or else get coerced to mpz.\n\
";
static PyObject *
Pygmpy_divexact(PyObject *self, PyObject *args)
{
    PyObject *x, *y, *result;
    PympzObject *tempx, *tempy;

    if(PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
        return NULL;
    }

    x = PyTuple_GET_ITEM(args, 0);
    y = PyTuple_GET_ITEM(args, 1);
    CREATE1_ONE_MPZANY(x, result);

    if(CHECK_MPZANY(x) && CHECK_MPZANY(y)) {
        if(mpz_sgn(Pympz_AS_MPZ(y)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF(result);
            return NULL;
        }
        mpz_divexact(Pympz_AS_MPZ(result), Pympz_AS_MPZ(x), Pympz_AS_MPZ(y));
    } else {
        tempx = Pympz_From_Integer(x);
        tempy = Pympz_From_Integer(y);
        if(!tempx || !tempy) {
            TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
            Py_XDECREF((PyObject*)tempx);
            Py_XDECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempy)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)tempx);
            Py_DECREF((PyObject*)tempy);
            Py_DECREF(result);
            return NULL;
        }
        mpz_divexact(Pympz_AS_MPZ(result), tempx->z, tempy->z);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
    }
    return result;
}

static PyObject *
Pympz_divexact(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(CHECK_MPZANY(other)) {
        if(mpz_sgn(Pympz_AS_MPZ(other)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_divexact(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                return NULL;
            }
            mpz_divexact(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), Pympz_AS_MPZ(other));
            return result;
        }
    } else {
        if(!(tempx = Pympz_From_Integer(other))) {
            TYPE_ERROR("divexact() requires 'mpz','mpz' arguments");
            return NULL;
        }
        if(mpz_sgn(Pympz_AS_MPZ(tempx)) == 0) {
            ZERO_ERROR("divexact() division by 0");
            Py_DECREF((PyObject*)tempx);
            return NULL;
        }
        if(Pyxmpz_Check(self)) {
            mpz_divexact(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            Py_RETURN_NONE;
        } else {
            if(!(result = (PyObject*)Pympz_new())) {
                Py_DECREF((PyObject*)tempx);
                return NULL;
            }
            mpz_divexact(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self), tempx->z);
            Py_DECREF((PyObject*)tempx);
            return result;
        }
    }
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
    PympzObject* newob;

    if(self && (CHECK_MPZANY(self))) {
        res = mpz_perfect_square_p(Pympz_AS_MPZ(self));
    } else if(CHECK_MPZANY(other)) {
        res = mpz_perfect_square_p(Pympz_AS_MPZ(other));
    } else {
        newob = Pympz_From_Integer(other);
        if(newob) {
            res = mpz_perfect_square_p(Pympz_AS_MPZ(newob));
            Py_DECREF((PyObject*)newob);
        } else {
            TYPE_ERROR("is_square() requires 'mpz' argument");
            return NULL;
        }
    }
    if(res)
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
    PympzObject* newob;

    if(self && (CHECK_MPZANY(self))) {
        res = mpz_perfect_power_p(Pympz_AS_MPZ(self));
    } else if(CHECK_MPZANY(other)) {
        res = mpz_perfect_power_p(Pympz_AS_MPZ(other));
    } else {
        newob = Pympz_From_Integer(other);
        if(newob) {
            res = mpz_perfect_power_p(Pympz_AS_MPZ(newob));
            Py_DECREF((PyObject*)newob);
        } else {
            TYPE_ERROR("is_power() requires 'mpz' argument");
            return NULL;
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

    if(reps <= 0) {
        VALUE_ERROR("repetition count for is_prime must be positive");
        Py_DECREF(self);
        return NULL;
    }
    i = (long) mpz_probab_prime_p(Pympz_AS_MPZ(self), reps);
    Py_DECREF(self);
    if(i)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

static char doc_next_primem[]="\
x.next_prime(): returns the smallest prime number > x.  Note that\n\
GMP may use a probabilistic definition of 'prime', and also that\n\
if x<0 GMP considers x 'prime' iff -x is prime; gmpy reflects these\n\
GMP design choices.\n\
";
static char doc_next_primeg[]="\
next_prime(x): returns the smallest prime number > x.  Note that\n\
GMP may use a probabilistic definition of 'prime', and also that\n\
if x<0 GMP considers x 'prime' iff -x is prime; gmpy reflects these\n\
GMP design choices. x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_next_prime(PyObject *self, PyObject *other)
{
    PyObject *result;
    PympzObject *tempx;

    if(self && Pyxmpz_Check(self)) {
        mpz_nextprime(Pympz_AS_MPZ(self), Pympz_AS_MPZ(self));
        Py_RETURN_NONE;
    } else if(self && Pympz_Check(self)) {
        if(!(result = (PyObject*)Pympz_new()))
            return NULL;
        mpz_nextprime(Pympz_AS_MPZ(result), Pympz_AS_MPZ(self));
    } else if(CHECK_MPZANY(other)) {
        CREATE1_ONE_MPZANY(other, result);
        mpz_nextprime(Pympz_AS_MPZ(result), Pympz_AS_MPZ(other));
    } else {
        CREATE0_ONE_MPZANY(result);
        tempx = Pympz_From_Integer(other);
        if(tempx) {
            mpz_nextprime(Pympz_AS_MPZ(result), tempx->z);
            Py_DECREF((PyObject*)tempx);
        } else {
            TYPE_ERROR("next_prime() requires 'mpz' argument");
            return NULL;
        }
    }
    return result;
}

static char doc_jacobim[]="\
x.jacobi(y): returns the Jacobi symbol (x|y) (y should be odd\n\
and must be positive).\n\
";
static char doc_jacobig[]="\
jacobi(x,y): returns the Jacobi symbol (x|y) (y should be odd and\n\
must be positive); x and y must be mpz, or else get coerced to mpz.\n\
";
static PyObject *
Pympz_jacobi(PyObject *self, PyObject *args)
{
    PyObject* other;
    long i;

    PARSE_TWO_MPZ(other, "jacobi() requires 'mpz','mpz' arguments");

    if(mpz_sgn(Pympz_AS_MPZ(other))<=0) {
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

static char doc_legendrem[]="\
x.legendre(y): returns the Legendre symbol (x|y) (y should be odd\n\
and must be positive).\n\
";
static char doc_legendreg[]="\
legendre(x,y): returns the Legendre symbol (x|y) (y should be odd\n\
and must be positive); x must be an mpz, or else gets coerced to one.\n\
";
static PyObject *
Pympz_legendre(PyObject *self, PyObject *args)
{
    PyObject* other;
    long i;

    PARSE_TWO_MPZ(other, "legendre() requires 'mpz','mpz' arguments");

    if(mpz_sgn(Pympz_AS_MPZ(other))<=0) {
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

static char doc_kroneckerm[]="\
x.kronecker(y): returns the Kronecker-Jacobi symbol (x|y).\n\
(At least one of x and y must fit in a plain int).\n\
";
static char doc_kroneckerg[]="\
kronecker(x,y): returns the Kronecker-Jacobi symbol (x|y).\n\
x and y must be mpz, or else get coerced to mpz (at least\n\
one of x and y, however, must also fit in a plain int).\n\
";
static PyObject *
Pympz_kronecker(PyObject *self, PyObject *args)
{
    PyObject* other;
    int ires;

    PARSE_TWO_MPZ(other, "kronecker() requires 'mpz','mpz' arguments");

    if(mpz_fits_ulong_p(Pympz_AS_MPZ(self))) {
        ires = mpz_ui_kronecker(
                mpz_get_ui(Pympz_AS_MPZ(self)),Pympz_AS_MPZ(other));
    } else if(mpz_fits_ulong_p(Pympz_AS_MPZ(other))) {
        ires = mpz_kronecker_ui(
                Pympz_AS_MPZ(self),mpz_get_ui(Pympz_AS_MPZ(other)));
    } else if(mpz_fits_slong_p(Pympz_AS_MPZ(self))) {
        ires = mpz_si_kronecker(
                mpz_get_si(Pympz_AS_MPZ(self)),Pympz_AS_MPZ(other));
    } else if(mpz_fits_slong_p(Pympz_AS_MPZ(other))) {
        ires = mpz_kronecker_si(
                Pympz_AS_MPZ(self),mpz_get_si(Pympz_AS_MPZ(other)));
    } else {
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
    PympzObject* newob;

    if(self && (CHECK_MPZANY(self))) {
        res = mpz_even_p(Pympz_AS_MPZ(self));
    } else if(CHECK_MPZANY(other)) {
        res = mpz_even_p(Pympz_AS_MPZ(other));
    } else {
        newob = Pympz_From_Integer(other);
        if(newob) {
            res = mpz_even_p(Pympz_AS_MPZ(newob));
            Py_DECREF((PyObject*)newob);
        } else {
            TYPE_ERROR("is_even() requires 'mpz' argument");
            return NULL;
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
    PympzObject* newob;

    if(self && (CHECK_MPZANY(self))) {
        res = mpz_odd_p(Pympz_AS_MPZ(self));
    } else if(CHECK_MPZANY(other)) {
        res = mpz_odd_p(Pympz_AS_MPZ(other));
    } else {
        newob = Pympz_From_Integer(other);
        if(newob) {
            res = mpz_odd_p(Pympz_AS_MPZ(newob));
            Py_DECREF((PyObject*)newob);
        } else {
            TYPE_ERROR("is_odd() requires 'mpz' argument");
            return NULL;
        }
    }
    if(res)
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
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
    } else if (PySlice_Check(item)) {
        Py_ssize_t start, stop, step, slicelength, cur, i;
        PyObject* result;

        if (PySlice_GetIndicesEx((PySliceObject*)item, mpz_sizeinbase(self->z, 2),
                         &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }

        if ((step < 0 && start < stop) ||
            (step > 0 && start > stop))
            stop = start;

        if (!(result = (PyObject*)Pyxmpz_new())) return NULL;
        mpz_set_ui(Pympz_AS_MPZ(result), 0);
        if (slicelength > 0) {
            for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                if(mpz_tstbit(self->z, cur)) {
                    mpz_setbit(Pympz_AS_MPZ(result), i);
                }
            }
        }
        return result;
    } else {
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
        } else if (bit_value == 0) {
            mpz_clrbit(self->z, i);
            return 0;
        } else {
            VALUE_ERROR("bit value must be 0 or 1");
            return -1;
        }
    } else if (PySlice_Check(item)) {
        Py_ssize_t bit_value, start, stop, step, slicelength;

        if (PySlice_GetIndicesEx((PySliceObject*)item, mpz_sizeinbase(self->z, 2),
                         &start, &stop, &step, &slicelength) < 0) {
            return -1;
        }

        if ((step < 0 && start < stop) ||
            (step > 0 && start > stop))
            stop = start;

        if (value == NULL) {
            TYPE_ERROR("deleting bits not supported");
            return -1;
        } else {
            Py_ssize_t cur, i;

            bit_value = PyNumber_AsSsize_t(value, PyExc_ValueError);
            if (bit_value == -1 && PyErr_Occurred()) {
                VALUE_ERROR("bit value must be 0 or 1");
                return -1;
            }
            if (bit_value == 0) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_clrbit(self->z, cur);
                }
            } else if (bit_value == 1) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_setbit(self->z, cur);
                }
            } else {
                VALUE_ERROR("bit value must be 0 or 1");
                return -1;
            }
#if 0
        } else {
            Py_ssize_t cur, i;
            int bit;
            PympzObject *tempx;
            if(!(tempx=Pympz_From_Integer(value))) {
                VALUE_ERROR("must specify bits sequence as an mpz");
                return -1;
            }
            if (mpz_sgn(tempx->z) == 0) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_clrbit(self->z, cur);
                }
            } else if (!(mpz_cmp_si(tempx->z, -1))) {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    mpz_setbit(self->z, cur);
                }
            } else {
                for (cur = start, i = 0; i < slicelength; cur += step, i++) {
                    bit = mpz_tstbit(tempx->z, i);
                    if(bit)
                        mpz_setbit(self->z, cur);
                    else
                        mpz_clrbit(self->z, cur);
                }
            }
            Py_DECREF(tempx);
#endif
            return 0;
        }
    } else {
        TYPE_ERROR("bit positions must be integers");
        return -1;
    }
    return -1;
}
