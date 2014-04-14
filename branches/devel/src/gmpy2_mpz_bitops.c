/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz_bitops.c                                                      *
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

PyDoc_STRVAR(doc_bit_length_function,
"bit_length(x) -> int\n\n"
"Return the number of significant bits in the radix-2\n"
"representation of x. Note: bit_length(0) returns 0.");

PyDoc_STRVAR(doc_bit_length_method,
"x.bit_length() -> int\n\n"
"Return the number of significant bits in the radix-2\n"
"representation of x. Note: mpz(0).bit_length() returns 0.");

static PyObject *
GMPy_MPZ_bit_length(PyObject *self, PyObject *other)
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
        if (!(tempx = GMPy_MPZ_From_Integer(other, NULL))) {
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

PyDoc_STRVAR(doc_bit_mask,
"bit_mask(n) -> mpz\n\n"
"Return an 'mpz' exactly n bits in length with all bits set.\n");

static PyObject *
GMPy_MPZ_bit_mask(PyObject *self, PyObject *other)
{
    Py_ssize_t i = 0;
    MPZ_Object* result;
    CTXT_Object *context = NULL;

    i = ssize_t_From_Integer(other);

    if (i == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_mask() requires 'int' argument");
        return NULL;
    }

    if (i < 0) {
        VALUE_ERROR("mask length must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    mpz_set_ui(result->z, 1);
    mpz_mul_2exp(result->z, result->z, i);
    mpz_sub_ui(result->z, result->z, 1);

    return (PyObject*)result;
}

/* return scan0/scan1 for an mpz */
PyDoc_STRVAR(doc_bit_scan0_method,
"x.bit_scan0(n=0) -> int\n\n"
"Return the index of the first 0-bit of x with index >= n. n >= 0.\n"
"If there are no more 0-bits in x at or above index n (which can\n"
"only happen for x<0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

PyDoc_STRVAR(doc_bit_scan0_function,
"bit_scan0(x, n=0) -> int\n\n"
"Return the index of the first 0-bit of x with index >= n. n >= 0.\n"
"If there are no more 0-bits in x at or above index n (which can\n"
"only happen for x<0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

static PyObject *
GMPy_MPZ_bit_scan0(PyObject *self, PyObject *args)
{
    Py_ssize_t maxbit, starting_bit = 0;
    PyObject *result;
    CTXT_Object *context = NULL;

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

PyDoc_STRVAR(doc_bit_scan1_method,
"x.bit_scan1(n=0) -> int\n\n"
"Return the index of the first 1-bit of x with index >= n. n >= 0.\n"
"If there are no more 1-bits in x at or above index n (which can\n"
"only happen for x>=0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

PyDoc_STRVAR(doc_bit_scan1_function,
"bit_scan1(x, n=0) -> int\n\n"
"Return the index of the first 1-bit of x with index >= n. n >= 0.\n"
"If there are no more 1-bits in x at or above index n (which can\n"
"only happen for x>=0, assuming an infinitely long 2's complement\n"
"format), then None is returned.");

static PyObject *
GMPy_MPZ_bit_scan1(PyObject *self, PyObject *args)
{
    Py_ssize_t maxbit, starting_bit = 0;
    PyObject *result;
    CTXT_Object *context = NULL;

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

PyDoc_STRVAR(doc_popcount,
"popcount(x) -> int\n\n"
"Return the number of 1-bits set in x. If x<0, the number of\n"
"1-bits is infinite so -1 is returned in that case.");

static PyObject *
GMPy_MPZ_popcount(PyObject *self, PyObject *other)
{
    Py_ssize_t temp;
    MPZ_Object *tempx;
    CTXT_Object *context = NULL;

    if (self && (CHECK_MPZANY(self)))
        return PyIntOrLong_FromSsize_t(mpz_popcount(MPZ(self)));
    else if(CHECK_MPZANY(other))
        return PyIntOrLong_FromSsize_t(mpz_popcount(MPZ(other)));
    else {
        if ((tempx = GMPy_MPZ_From_Integer(other, context))) {
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
PyDoc_STRVAR(doc_bit_test_function,
"bit_test(x, n) -> bool\n\n"
"Return the value of the n-th bit of x.");

static PyObject *
GMPy_MPZ_bit_test_function(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    int temp;
    PyObject *x;
    MPZ_Object *tempx;
    CTXT_Object *context = NULL;

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
        if (!(tempx = GMPy_MPZ_From_Integer(x, context))) {
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

PyDoc_STRVAR(doc_bit_test_method,
"x.bit_test(n) -> bool\n\n"
"Return the value of the n-th bit of x.");

static PyObject *
GMPy_MPZ_bit_test_method(PyObject *self, PyObject *other)
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

PyDoc_STRVAR(doc_bit_clear_function,
"bit_clear(x, n) -> mpz\n\n"
"Return a copy of x with the n-th bit cleared.");

static PyObject *
GMPy_MPZ_bit_clear_function(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    PyObject *x;
    MPZ_Object *result;
    CTXT_Object *context = NULL;

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
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_set(result->z, MPZ(x));
        mpz_clrbit(result->z, bit_index);
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(x, context))) {
            TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_clrbit(result->z, bit_index);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_clear_method,
"x.bit_clear(n) -> mpz\n\n"
"Return a copy of x with the n-th bit cleared.");

static PyObject *
GMPy_MPZ_bit_clear_method(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_clear() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;
    mpz_set(result->z, MPZ(self));
    mpz_clrbit(result->z, bit_index);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_set_function,
"bit_set(x, n) -> mpz\n\n"
"Return a copy of x with the n-th bit set.");

static PyObject *
GMPy_MPZ_bit_set_function(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    PyObject *x;
    MPZ_Object *result;
    CTXT_Object *context = NULL;

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
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_set(result->z, MPZ(x));
        mpz_setbit(result->z, bit_index);
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(x, context))) {
            TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_setbit(result->z, bit_index);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_set_method,
"x.bit_set(n) -> mpz\n\n"
"Return a copy of x with the n-th bit set.");

static PyObject *
GMPy_MPZ_bit_set_method(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_set() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;
    mpz_set(result->z, MPZ(self));
    mpz_setbit(result->z, bit_index);
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_flip_function,
"bit_flip(x, n) -> mpz\n\n"
"Return a copy of x with the n-th bit inverted.");

static PyObject *
GMPy_MPZ_bit_flip_function(PyObject *self, PyObject *args)
{
    Py_ssize_t bit_index;
    PyObject *x;
    MPZ_Object *result;
    CTXT_Object *context = NULL;

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
        if (!(result = GMPy_MPZ_New(context)))
            return NULL;
        mpz_set(result->z, MPZ(x));
        mpz_combit(result->z, bit_index);
    }
    else {
        if (!(result = GMPy_MPZ_From_Integer(x, context))) {
            TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
            return NULL;
        }
        mpz_combit(result->z, bit_index);
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_bit_flip_method,
"x.bit_flip(n) -> mpz\n\n"
"Return a copy of x with the n-th bit inverted.");

static PyObject *
GMPy_MPZ_bit_flip_method(PyObject *self, PyObject *other)
{
    Py_ssize_t bit_index;
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    bit_index = ssize_t_From_Integer(other);
    if (bit_index == -1 && PyErr_Occurred()) {
        TYPE_ERROR("bit_flip() requires 'mpz','int' arguments");
        return NULL;
    }

    if (bit_index < 0) {
        VALUE_ERROR("bit_index must be >= 0");
        return NULL;
    }

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;
    mpz_set(result->z, MPZ(self));
    mpz_combit(result->z, bit_index);
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Invert_Slot(MPZ_Object *self)
{
    MPZ_Object *result;

    if ((result = GMPy_MPZ_New(NULL)))
        mpz_com(result->z, MPZ(self));

    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_And_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(self)) {
        if (CHECK_MPZANY(other)) {
            if (!(result = GMPy_MPZ_New(NULL)))
                return NULL;
            mpz_and(result->z, MPZ(self), MPZ(other));
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL)))
                return NULL;
            mpz_and(result->z, MPZ(self), result->z);
        }
    }
    else if (CHECK_MPZANY(other)) {
        if (!(result = GMPy_MPZ_From_Integer(self, NULL)))
            return NULL;
        mpz_and(result->z, result->z, MPZ(other));
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Ior_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(self)) {
        if (CHECK_MPZANY(other)) {
            if (!(result = GMPy_MPZ_New(NULL)))
                return NULL;
            mpz_ior(result->z, MPZ(self), MPZ(other));
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL)))
                return NULL;
            mpz_ior(result->z, MPZ(self), result->z);
        }
    }
    else if (CHECK_MPZANY(other)) {
        if (!(result = GMPy_MPZ_From_Integer(self, NULL)))
            return NULL;
        mpz_ior(result->z, result->z, MPZ(other));
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Xor_Slot(PyObject *self, PyObject *other)
{
    MPZ_Object *result;

    if (CHECK_MPZANY(self)) {
        if (CHECK_MPZANY(other)) {
            if (!(result = GMPy_MPZ_New(NULL)))
                return NULL;
            mpz_xor(result->z, MPZ(self), MPZ(other));
        }
        else {
            if (!(result = GMPy_MPZ_From_Integer(other, NULL)))
                return NULL;
            mpz_xor(result->z, MPZ(self), result->z);
        }
    }
    else if (CHECK_MPZANY(other)) {
        if (!(result = GMPy_MPZ_From_Integer(self, NULL)))
            return NULL;
        mpz_xor(result->z, result->z, MPZ(other));
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }
    return (PyObject*)result;
}

static PyObject *
GMPy_MPZ_Rshift_Slot(PyObject *self, PyObject *other)
{
    mpir_si count_si;
    int overflow;
    MPZ_Object *result, *tempa, *tempb;
    CTXT_Object *context = NULL;

    if (!(result = GMPy_MPZ_New(context)))
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

    tempa = GMPy_MPZ_From_Integer(self, context);
    tempb = GMPy_MPZ_From_Integer(other, context);
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
GMPy_MPZ_Lshift_Slot(PyObject *self, PyObject *other)
{
    mpir_si count_si;
    int overflow;
    MPZ_Object *result, *tempa, *tempb;
    CTXT_Object *context = NULL;

    if (!(result = GMPy_MPZ_New(context)))
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

    tempa = GMPy_MPZ_From_Integer(self, context);
    tempb = GMPy_MPZ_From_Integer(other, context);
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

PyDoc_STRVAR(doc_hamdist,
"hamdist(x, y) -> int\n\n"
"Return the Hamming distance (number of bit-positions where the\n"
"bits differ) between integers x and y.");

static PyObject *
GMPy_MPZ_hamdist(PyObject *self, PyObject *args)
{
    PyObject *result, *other;
    CTXT_Object *context = NULL;

    PARSE_TWO_MPZ(other, "hamdist() requires 'mpz','mpz' arguments");

    result = PyIntOrLong_FromSize_t(
            mpz_hamdist(MPZ(self),MPZ(other)));
    Py_DECREF(self);
    Py_DECREF(other);
    return (PyObject*)result;
}

