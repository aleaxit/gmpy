/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq.c                                                              *
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

PyDoc_STRVAR(doc_mpq,
"mpq() -> mpq(0,1)\n\n"
"     If no argument is given, return mpq(0,1).\n\n"
"mpq(n) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n. Decimal and\n"
"     Fraction values are converted exactly.\n\n"
"mpq(n,m) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n/m.\n\n"
"mpq(s[, base=10]) -> mpq\n\n"
"     Return an 'mpq' object from a string s made up of digits in\n"
"     the given base. s may be made up of two numbers in the same\n"
"     base separated by a '/' character.\n");

static PyObject *
Pygmpy_mpq(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPQ_Object *result = NULL, *temp;
    PyObject *n = NULL, *m = NULL;
    int base = 10;
    Py_ssize_t argc;
    static char *kwlist[] = {"s", "base", NULL };
    CTXT_Object *context = NULL;

    argc = PyTuple_Size(args);
    if (argc > 2) {
        TYPE_ERROR("mpq() requires 0, 1 or 2 arguments");
        return NULL;
    }

    /* Handle 0 arguments. */
    if (argc == 0) {
        if ((result = GMPy_MPQ_New(context))) {
            mpq_set_ui(result->q, 0, 1);
        }
        return (PyObject*)result;
    }

    /* Handle the case where the first argument is a string. */
    n = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(n)) {
        /* keyword base is legal */
        if (PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, &n, &base)) {
            if ((base != 0) && ((base < 2) || (base > 62))) {
                VALUE_ERROR("base for mpq() must be 0 or in the "
                            "interval 2 ... 62");
            }
            else {
                result = GMPy_MPQ_From_PyStr(n, base, context);
            }
        }
        return (PyObject*)result;
    }

    /* Handle 1 argument. It can be non-complex number. */
    if (argc == 1) {
        if (IS_REAL(n)) {
            return (PyObject*)GMPy_MPQ_From_Number_Temp(n, context);
        }
        else {
            goto type_error;
        }
    }

    /* Handle 2 arguments. Both arguments must be integer or rational. */
    m = PyTuple_GetItem(args, 1);

    if (!IS_RATIONAL(n) || !IS_RATIONAL(m)) {
        goto type_error;
    }

    if (!(result = GMPy_MPQ_From_Number_New(n, context))) {
        goto type_error;
    }

    if (!(temp = GMPy_MPQ_From_Number_Temp(m, context))) {
        Py_DECREF((PyObject*)result);
        goto type_error;
    }

    if (mpq_sgn(temp->q) == 0) {
        ZERO_ERROR("zero denominator in mpq()");
        Py_DECREF((PyObject*)result);
        Py_DECREF((PyObject*)temp);
        return NULL;
    }

    mpq_div(result->q, result->q, temp->q);
    Py_DECREF((PyObject*)temp);
    return (PyObject*)result;

  type_error:
    TYPE_ERROR("mpq() requires numeric or string argument");
    return NULL;
}

/* Functions that operate strictly on mpq. */

/* produce digits for an mpq in requested base, default 10 */
PyDoc_STRVAR(doc_qdigitsm,
"x.digits([base=10]) -> string\n\n"
"Return a Python string representing x in the given base (2 to 62,\n"
"default is 10). A leading '-' is present if x<0, but no leading '+'\n"
"is present if x>=0.\n");

static PyObject *
Pympq_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    PyObject *result;
    CTXT_Object *context = NULL;

    SELF_MPQ_ONE_ARG("|i", &base);
    result = GMPy_PyStr_From_MPQ((MPQ_Object*)self, base, 0, context);
    Py_DECREF(self);
    return result;
}

/* Since Pympq_sign() is called by Pympany_sign(), we know that 'other' is
 * a Rational type.
 */

static PyObject *
Pympq_sign(PyObject *self, PyObject *other)
{
    long res;
    MPQ_Object* tempx;
    CTXT_Object *context = NULL;

    if (MPQ_Check(other)) {
        res = mpq_sgn(MPQ(other));
    }
    else {
        if (!(tempx = GMPy_MPQ_From_Number_Temp(other, context))) {
            TYPE_ERROR("sign() requires 'mpq' argument");
            return NULL;
        }
        else {
            res = mpq_sgn(tempx->q);
            Py_DECREF((PyObject*)tempx);
        }
    }
    return PyIntOrLong_FromLong(res);
}

PyDoc_STRVAR(doc_numerg,
"numer(x) -> mpz\n\n"
"Return the numerator of x.");

static PyObject *
Pympq_numer(PyObject *self, PyObject *args)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    SELF_MPQ_NO_ARG;
    assert(MPQ_Check(self));
    mpz_set(result->z, mpq_numref(MPQ(self)));
    Py_DECREF(self);
    return (PyObject*)result;
}

static PyObject *
Pympq_getnumer(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, mpq_numref(MPQ(self)));
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_denomg,
"denom(x) -> mpz\n\n"
"Return the denominator of x.");

static PyObject *
Pympq_denom(PyObject *self, PyObject *args)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if (!(result = GMPy_MPZ_New(context)))
        return NULL;

    SELF_MPQ_NO_ARG;
    assert(MPQ_Check(self));
    mpz_set(result->z, mpq_denref(MPQ(self)));
    Py_DECREF(self);
    return (PyObject*)result;
}

static PyObject *
Pympq_getdenom(MPQ_Object *self, void *closure)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context)))
        mpz_set(result->z, mpq_denref(MPQ(self)));
    return (PyObject*)result;
}

/* TODO: add error checking for arguments. */

PyDoc_STRVAR(doc_qdivg,
"qdiv(x[, y=1]) -> number\n\n"
"Return x/y as 'mpz' if possible, or as 'mpq' if x is not exactly\n"
"divisible by y.");

static PyObject *
Pympq_qdiv(PyObject *self, PyObject *args)
{
    Py_ssize_t argc;
    int isOne = 0;
    PyObject *result, *x, *y;
    MPQ_Object *tempx = NULL, *tempy = NULL;
    CTXT_Object *context = NULL;

    argc = PyTuple_Size(args);
    if (argc == 0 || argc > 2)
        goto arg_error;

    /* Validate the first argument. */

    x = PyTuple_GET_ITEM(args, 0);
    if (!IS_RATIONAL(x))
        goto arg_error;

    /* Convert the second argument first and see if it is 1. If the second
     * argument is 1 (or didn't exist), then isOne is true and tempy does not
     * contain a valid reference.
     */

    if (argc == 2) {
        y = PyTuple_GET_ITEM(args, 1);
        if (!IS_RATIONAL(y)) {
            goto arg_error;
        }
        tempy = GMPy_MPQ_From_Rational_Temp(y, context);
        if (!tempy) {
            return NULL;
        }
        if (mpq_cmp_ui(tempy->q, 1, 1) == 0) {
            isOne = 1;
            Py_DECREF((PyObject*)tempy);
        }
    }
    else {
        isOne = 1;
    }

    if (isOne) {
        if (IS_INTEGER(x))
            return (PyObject*)GMPy_MPZ_From_Integer_Temp(x, context);

        if (IS_RATIONAL(x)) {
            tempx = GMPy_MPQ_From_Rational_Temp(x, context);
            if (!tempx)
                return NULL;
            if (mpz_cmp_ui(mpq_denref(tempx->q), 1) == 0) {
                if ((result = (PyObject*)GMPy_MPZ_New(context))) {
                    mpz_set(MPZ(result), mpq_numref(tempx->q));
                }
                Py_DECREF((PyObject*)tempx);
                return result;
            }
            else {
                return (PyObject*)tempx;
            }
        }
        goto arg_error;
    }

    if (mpq_sgn(tempy->q) == 0) {
        Py_DECREF((PyObject*)tempy);
        ZERO_ERROR("division by zero in qdiv()");
        return NULL;
    }

    if (!(tempx = GMPy_MPQ_From_Number_New(x, context))) {
        Py_DECREF((PyObject*)tempy);
        return NULL;
    }

    mpq_div(tempx->q, tempx->q, tempy->q);
    Py_DECREF((PyObject*)tempy);

    if (mpz_cmp_ui(mpq_denref(tempx->q), 1) == 0) {
        if ((result = (PyObject*)GMPy_MPZ_New(context))) {
            mpz_set(MPZ(result), mpq_numref(tempx->q));
        }
        Py_DECREF((PyObject*)tempx);
        return result;
    }

    return (PyObject*)tempx;

  arg_error:
    TYPE_ERROR("qdiv() takes one or two Integer or Rational argument(s)");
    return NULL;
}

static PyObject *
Pympq_neg(MPQ_Object *self)
{
    MPQ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPQ_New(context))) {
        mpq_neg(result->q, self->q);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpq_floor,
             "Return greatest integer less than or equal to an mpq.");

static PyObject *
Pympq_floor(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context))) {
        mpz_fdiv_q(result->z,
                   mpq_numref(MPQ(self)),
                   mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpq_ceil,
             "Return least integer greater than or equal to an mpq.");

static PyObject *
Pympq_ceil(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context))) {
        mpz_cdiv_q(result->z,
                   mpq_numref(MPQ(self)),
                   mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpq_trunc,
             "Return integer portion of an mpq.");

static PyObject *
Pympq_trunc(PyObject *self, PyObject *other)
{
    MPZ_Object *result;
    CTXT_Object *context = NULL;

    if ((result = GMPy_MPZ_New(context))) {
        mpz_tdiv_q(result->z,
                   mpq_numref(MPQ(self)),
                   mpq_denref(MPQ(self)));
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpq_round, "Round an mpq to power of 10.");

static PyObject *
Pympq_round(PyObject *self, PyObject *args)
{
    Py_ssize_t round_digits = 0;
    MPQ_Object *resultq;
    MPZ_Object *resultz;
    mpz_t temp, rem;
    CTXT_Object *context = NULL;

    /* If args is NULL or the size of args is 0, we just return an mpz. */

    if (!args || PyTuple_GET_SIZE(args) == 0) {
        if (!(resultz = GMPy_MPZ_New(context)))
            return NULL;

        mpz_inoc(rem);
        mpz_fdiv_qr(resultz->z, rem, mpq_numref(MPQ(self)),
                    mpq_denref(MPQ(self)));
        mpz_mul_2exp(rem, rem, 1);
        if (mpz_cmp(rem, mpq_denref(MPQ(self))) > 0) {
            mpz_add_ui(resultz->z, resultz->z, 1);
        }
        else if (mpz_cmp(rem, mpq_denref(MPQ(self))) == 0) {
            if (mpz_odd_p(resultz->z)) {
                mpz_add_ui(resultz->z, resultz->z, 1);
            }
        }
        mpz_cloc(rem);
        return (PyObject*)resultz;
    }

    if (PyTuple_GET_SIZE(args) > 1) {
        TYPE_ERROR("Too many arguments for __round__().");
        return NULL;
    }

    if (PyTuple_GET_SIZE(args) == 1) {
        round_digits = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 0));
        if (round_digits == -1 && PyErr_Occurred()) {
            TYPE_ERROR("__round__() requires 'int' argument");
            return NULL;
        }
    }

    if (!(resultq = GMPy_MPQ_New(context)))
        return NULL;

    mpz_inoc(temp);
    mpz_ui_pow_ui(temp, 10, round_digits > 0 ? round_digits : -round_digits);

    mpq_set(resultq->q, MPQ(self));
    if (round_digits > 0) {
        mpz_mul(mpq_numref(resultq->q), mpq_numref(resultq->q), temp);
        mpq_canonicalize(resultq->q);
        if (!(resultz = (MPZ_Object*)Pympq_round((PyObject*)resultq, NULL))) {
            mpz_cloc(temp);
            return NULL;
        }
        mpz_set(mpq_numref(resultq->q), resultz->z);
        Py_DECREF((PyObject*)resultz);
        mpz_set(mpq_denref(resultq->q), temp);
        mpz_cloc(temp);
        mpq_canonicalize(resultq->q);
    }
    else {
        mpz_mul(mpq_denref(resultq->q), mpq_denref(resultq->q), temp);
        mpq_canonicalize(resultq->q);
        if (!(resultz = (MPZ_Object*)Pympq_round((PyObject*)resultq, NULL))) {
            mpz_cloc(temp);
            return NULL;
        }
        mpq_set_ui(resultq->q, 0, 1);
        mpz_mul(mpq_numref(resultq->q), resultz->z, temp);
        Py_DECREF((PyObject*)resultz);
        mpz_cloc(temp);
        mpq_canonicalize(resultq->q);
    }
    return (PyObject*)resultq;
}

static PyObject *
Pympq_pos(MPQ_Object *self)
{
    Py_INCREF((PyObject*)self);
    return (PyObject*)self;
}

static PyObject *
Pympq_square(PyObject *self, PyObject *other)
{
    MPQ_Object *tempx, *result;
    CTXT_Object *context = NULL;

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (self && (MPQ_Check(self))) {
        mpq_mul(result->q, MPQ(self), MPQ(self));
    }
    else if (MPQ_Check(other)) {
        mpq_mul(result->q, MPQ(other), MPQ(other));
    }
    else {
        if (!(tempx = GMPy_MPQ_From_Rational_Temp(other, context))) {
            TYPE_ERROR("square() requires 'mpq' argument");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        else {
            mpq_mul(result->q, MPQ(tempx), MPQ(tempx));
            Py_DECREF((PyObject*)tempx);
        }
    }
    return (PyObject*)result;
}

static int
Pympq_nonzero(MPQ_Object *self)
{
    return mpq_sgn(self->q) != 0;
}

static Py_hash_t
Pympq_hash(MPQ_Object *self)
{
#ifdef _PyHASH_MODULUS
    Py_hash_t hash = 0;
    mpz_t temp, temp1, mask;

    if (self->hash_cache != -1)
        return self->hash_cache;

    mpz_inoc(temp);
    mpz_inoc(temp1);
    mpz_inoc(mask);
    mpz_set_si(mask, 1);
    mpz_mul_2exp(mask, mask, _PyHASH_BITS);
    mpz_sub_ui(mask, mask, 1);

    if (!mpz_invert(temp, mpq_denref(self->q), mask)) {
        mpz_cloc(temp);
        mpz_cloc(temp1);
        mpz_cloc(mask);
        hash = _PyHASH_INF;
        if (mpz_sgn(mpq_numref(self->q))<0)
            hash = -hash;
        self->hash_cache = hash;
        return hash;
    }
    mpz_set(temp1, mask);
    mpz_sub_ui(temp1, temp1, 2);
    mpz_powm(temp, mpq_denref(self->q), temp1, mask);

    mpz_tdiv_r(temp1, mpq_numref(self->q), mask);
    mpz_mul(temp, temp, temp1);
    hash = (Py_hash_t)mpn_mod_1(temp->_mp_d, mpz_size(temp), _PyHASH_MODULUS);

    if (mpz_sgn(mpq_numref(self->q))<0)
        hash = -hash;
    if (hash==-1) hash = -2;
    mpz_cloc(temp);
    mpz_cloc(temp1);
    mpz_cloc(mask);
    self->hash_cache = hash;
    return hash;
#else
    PyObject *temp;

    if (self->hash_cache != -1)
        return self->hash_cache;

    if (!(temp = Pympq_To_PyFloat(self))) {
        SYSTEM_ERROR("Could not convert 'mpq' to float.");
        return -1;
    }
    self->hash_cache = PyObject_Hash(temp);
    Py_DECREF(temp);
    return self->hash_cache;
#endif
}

/* Divide two Rational objects (see convert.c/IS_RATIONAL) and return the
 * remainder. Returns None and raises TypeError if both objects are not valid
 * rationals. */

static PyObject *
Pympq_Mod_Rational(PyObject *x, PyObject *y, CTXT_Object *context)
{
    mpz_t tempz;
    MPQ_Object *tempx, *tempy, *result;

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        tempx = GMPy_MPQ_From_Number_Temp(x, context);
        tempy = GMPy_MPQ_From_Number_Temp(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Rational to mpq.");
            goto error;
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }

        mpz_inoc(tempz);
        mpq_div(result->q, tempx->q, tempy->q);
        mpz_fdiv_q(tempz, mpq_numref(result->q), mpq_denref(result->q));
        /* Need to calculate x - tempz * y. */
        mpq_set_z(result->q, tempz);
        mpq_mul(result->q, result->q, tempy->q);
        mpq_sub(result->q, tempx->q, result->q);
        mpz_cloc(tempz);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        return (PyObject*)result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)result);
    return NULL;
}

static PyObject *
Pympq_mod_fast(PyObject *x, PyObject *y)
{
    CTXT_Object *context = NULL;

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return Pympq_Mod_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        return Pympfr_Mod_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return Pympc_Mod_Complex(x, y, context);

    Py_RETURN_NOTIMPLEMENTED;
}

/* Return the quotient and remainder from dividing two Rational objects (see
 * convert.c/IS_RATIONAL). Returns None and raises TypeError if both objects
 * are not valid rationals. */

static PyObject *
Pympq_DivMod_Rational(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPQ_Object *tempx, *tempy, *rem;
    MPZ_Object *quo;
    PyObject *result;

    result = PyTuple_New(2);
    rem = GMPy_MPQ_New(context);
    quo = GMPy_MPZ_New(context);
    if (!result || !rem || !quo) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)rem);
        Py_XDECREF((PyObject*)quo);
        return NULL;
    }

    if (IS_RATIONAL(x) && IS_RATIONAL(y)) {
        tempx = GMPy_MPQ_From_Number_Temp(x, context);
        tempy = GMPy_MPQ_From_Number_Temp(y, context);
        if (!tempx || !tempy) {
            SYSTEM_ERROR("Could not convert Rational to mpq.");
            goto error;
        }
        if (mpq_sgn(tempy->q) == 0) {
            ZERO_ERROR("division or modulo by zero");
            goto error;
        }

        mpq_div(rem->q, tempx->q, tempy->q);
        mpz_fdiv_q(quo->z, mpq_numref(rem->q), mpq_denref(rem->q));
        /* Need to calculate x - quo * y. */
        mpq_set_z(rem->q, quo->z);
        mpq_mul(rem->q, rem->q, tempy->q);
        mpq_sub(rem->q, tempx->q, rem->q);
        Py_DECREF((PyObject*)tempx);
        Py_DECREF((PyObject*)tempy);
        PyTuple_SET_ITEM(result, 0, (PyObject*)quo);
        PyTuple_SET_ITEM(result, 1, (PyObject*)rem);
        return result;
    }

    Py_DECREF((PyObject*)result);
    Py_RETURN_NOTIMPLEMENTED;

  error:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_DECREF((PyObject*)rem);
    Py_DECREF((PyObject*)quo);
    Py_DECREF(result);
    return NULL;
}

static PyObject *
Pympq_divmod_fast(PyObject *x, PyObject *y)
{
    CTXT_Object *context= NULL;

    if (IS_RATIONAL(x) && IS_RATIONAL(y))
        return Pympq_DivMod_Rational(x, y, context);
    else if (IS_REAL(x) && IS_REAL(y))
        return Pympfr_DivMod_Real(x, y, context);
    else if (IS_COMPLEX(x) && IS_COMPLEX(y))
        return Pympc_DivMod_Complex(x, y, context);

    Py_RETURN_NOTIMPLEMENTED;
}

PyDoc_STRVAR(doc_mpq_sizeof,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x. Note: deleted mpq objects\n"
"are reused and may or may not be resized when a new value is assigned.");

static PyObject *
Pympq_sizeof(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPQ_Object) + \
        (mpq_numref(MPQ(self))->_mp_alloc * sizeof(mp_limb_t)) + \
        (mpq_denref(MPQ(self))->_mp_alloc * sizeof(mp_limb_t)));
}

#ifdef PY3
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) GMPy_mpq_add_fast,         /* nb_add                  */
    (binaryfunc) GMPy_mpq_sub_fast,         /* nb_subtract             */
    (binaryfunc) GMPy_mpq_mul_fast,         /* nb_multiply             */
    (binaryfunc) Pympq_mod_fast,            /* nb_remainder            */
    (binaryfunc) Pympq_divmod_fast,         /* nb_divmod               */
    (ternaryfunc) GMPy_mpany_pow_fast,      /* nb_power                */
    (unaryfunc) Pympq_neg,                  /* nb_negative             */
    (unaryfunc) Pympq_pos,                  /* nb_positive             */
    (unaryfunc) GMPy_mpq_abs_fast,          /* nb_absolute             */
    (inquiry) Pympq_nonzero,                /* nb_bool                 */
        0,                                  /* nb_invert               */
        0,                                  /* nb_lshift               */
        0,                                  /* nb_rshift               */
        0,                                  /* nb_and                  */
        0,                                  /* nb_xor                  */
        0,                                  /* nb_or                   */
    (unaryfunc) GMPy_PyLong_From_MPQ,       /* nb_int                  */
        0,                                  /* nb_reserved             */
    (unaryfunc) GMPy_PyFloat_From_MPQ,      /* nb_float                */
        0,                                  /* nb_inplace_add          */
        0,                                  /* nb_inplace_subtract     */
        0,                                  /* nb_inplace_multiply     */
        0,                                  /* nb_inplace_remainder    */
        0,                                  /* nb_inplace_power        */
        0,                                  /* nb_inplace_lshift       */
        0,                                  /* nb_inplace_rshift       */
        0,                                  /* nb_inplace_and          */
        0,                                  /* nb_inplace_xor          */
        0,                                  /* nb_inplace_or           */
    (binaryfunc) GMPy_mpq_floordiv_fast,    /* nb_floor_divide         */
    (binaryfunc) GMPy_mpq_truediv_fast,     /* nb_true_divide          */
        0,                                  /* nb_inplace_floor_divide */
        0,                                  /* nb_inplace_true_divide  */
        0,                                  /* nb_index                */
};
#else
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) GMPy_mpq_add_fast,         /* nb_add                  */
    (binaryfunc) GMPy_mpq_sub_fast,         /* nb_subtract             */
    (binaryfunc) GMPy_mpq_mul_fast,         /* nb_multiply             */
    (binaryfunc) GMPy_mpq_truediv_fast,     /* nb_divide               */
    (binaryfunc) Pympq_mod_fast,            /* nb_remainder            */
    (binaryfunc) Pympq_divmod_fast,         /* nb_divmod               */
    (ternaryfunc) GMPy_mpany_pow_fast,      /* nb_power                */
    (unaryfunc) Pympq_neg,                  /* nb_negative             */
    (unaryfunc) Pympq_pos,                  /* nb_positive             */
    (unaryfunc) GMPy_mpq_abs_fast,          /* nb_absolute             */
    (inquiry) Pympq_nonzero,                /* nb_bool                 */
        0,                                  /* nb_invert               */
        0,                                  /* nb_lshift               */
        0,                                  /* nb_rshift               */
        0,                                  /* nb_and                  */
        0,                                  /* nb_xor                  */
        0,                                  /* nb_or                   */
        0,                                  /* nb_coerce               */
    (unaryfunc) GMPy_PyIntOrLong_From_MPQ,  /* nb_int                  */
    (unaryfunc) Pympq_To_PyLong,            /* nb_long                 */
    (unaryfunc) Pympq_To_PyFloat,           /* nb_float                */
        0,                                  /* nb_oct                  */
        0,                                  /* nb_hex                  */
        0,                                  /* nb_inplace_add;         */
        0,                                  /* nb_inplace_subtract     */
        0,                                  /* nb_inplace_multiply     */
        0,                                  /* nb_inplace_divide       */
        0,                                  /* nb_inplace_remainder    */
        0,                                  /* nb_inplace_power        */
        0,                                  /* nb_inplace_lshift       */
        0,                                  /* nb_inplace_rshift       */
        0,                                  /* nb_inplace_and          */
        0,                                  /* nb_inplace_xor          */
        0,                                  /* nb_inplace_or           */
    (binaryfunc) GMPy_mpq_floordiv_fast,    /* nb_floor_divide         */
    (binaryfunc) GMPy_mpq_truediv_fast,     /* nb_true_divide          */
        0,                                  /* nb_inplace_floor_divide */
        0,                                  /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympq_getseters[] =
{
    { "numerator", (getter)Pympq_getnumer, NULL, "numerator", NULL },
    { "denominator", (getter)Pympq_getdenom, NULL, "denominator", NULL },
    {NULL}
};

static PyMethodDef Pympq_methods [] =
{
    { "__ceil__", Pympq_ceil, METH_NOARGS, doc_mpq_ceil },
    { "__floor__", Pympq_floor, METH_NOARGS, doc_mpq_floor },
    { "__round__", Pympq_round, METH_VARARGS, doc_mpq_round },
    { "__sizeof__", Pympq_sizeof, METH_NOARGS, doc_mpq_sizeof },
    { "__trunc__", Pympq_trunc, METH_NOARGS, doc_mpq_trunc },
    { "digits", Pympq_digits, METH_VARARGS, doc_qdigitsm },
    { NULL, NULL, 1 }
};

static PyTypeObject MPQ_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpq",                                  /* tp_name          */
    sizeof(MPQ_Object),                     /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPQ_Dealloc,          /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPQ_Repr_Slot,          /* tp_repr          */
    &mpq_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympq_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPQ_Str_Slot,           /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE |
        Py_TPFLAGS_CHECKTYPES,              /* tp_flags         */
#endif
    "Multiple precision rational",          /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympq_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympq_getseters,                        /* tp_getset        */
};

