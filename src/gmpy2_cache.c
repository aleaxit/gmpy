/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_cache.c                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017 Case Van Horsen                              *
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


/* gmpy2 caches objects so they can be reused quickly without involving a new
 * memory allocation or object construction.
 *
 * The "py???cache" is used to cache Py??? objects. The cache is accessed
 * via Py???_new/Py???_dealloc. The functions set_py???cache and
 * set_py???cache are used to change the size of the array used to the store
 * the cached objects.
 */

/* Caching logic for Pympz. */

static void
set_gmpympzcache(void)
{
    if (global.in_gmpympzcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < global.in_gmpympzcache; ++i) {
            mpz_clear(global.gmpympzcache[i]->z);
            PyObject_Del(global.gmpympzcache[i]);
        }
        global.in_gmpympzcache = global.cache_size;
    }
    global.gmpympzcache = realloc(global.gmpympzcache, sizeof(MPZ_Object)*global.cache_size);
}

/* GMPy_MPZ_New returns a reference to a new MPZ_Object. Its value
 * is initialized to 0.
 */

static MPZ_Object *
GMPy_MPZ_New(CTXT_Object *context)
{
    MPZ_Object *result = NULL;

    if (global.in_gmpympzcache) {
        result = global.gmpympzcache[--(global.in_gmpympzcache)];
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)result);
        mpz_set_ui(result->z, 0);
        result->hash_cache = -1;
    }
    else {
        if ((result = PyObject_New(MPZ_Object, &MPZ_Type))) {
            mpz_init(result->z);
            result->hash_cache = -1;
        }
    }
    return result;
}

/* GMPy_MPZ_NewInit returns a reference to an initialized MPZ_Object. It is
 * used by mpz.__new__ to replace the old mpz() factory function.
 */

static PyObject *
GMPy_MPZ_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds)
{
    MPZ_Object *result = NULL;
    PyObject *n = NULL;
    PyObject *temp = NULL;
    int base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"s", "base", NULL };
    CTXT_Object *context = NULL;

    if (type != &MPZ_Type) {
        TYPE_ERROR("mpz.__new__() requires mpz type");
        return NULL;
    }

    /* Optimize the most common use cases first; either 0 or 1 argument */

    argc = PyTuple_GET_SIZE(args);

    if (argc == 0) {
        return (PyObject*)GMPy_MPZ_New(context);
    }

    if (argc == 1 && !keywds) {
        n = PyTuple_GET_ITEM(args, 0);

        if (MPZ_Check(n)) {
            Py_INCREF(n);
            return n;
        }

        if (PyIntOrLong_Check(n)) {
            return (PyObject*)GMPy_MPZ_From_PyIntOrLong(n, context);
        }

        if (MPQ_Check(n)) {
            return (PyObject*)GMPy_MPZ_From_MPQ((MPQ_Object*)n, context);
        }

        if (MPFR_Check(n)) {
            return (PyObject*)GMPy_MPZ_From_MPFR((MPFR_Object*)n, context);
        }

        if (PyFloat_Check(n)) {
            return (PyObject*)GMPy_MPZ_From_PyFloat(n, context);
        }

        if (XMPZ_Check(n)) {
            return (PyObject*)GMPy_MPZ_From_XMPZ((XMPZ_Object*)n, context);
        }

        if (IS_FRACTION(n)) {
            MPQ_Object *temp = GMPy_MPQ_From_Fraction(n, context);

            if (temp) {
                result = GMPy_MPZ_From_MPQ(temp, context);
                Py_DECREF((PyObject*)temp);
            }
            return (PyObject*)result;
        }

        if (PyStrOrUnicode_Check(n)) {
            return (PyObject*)GMPy_MPZ_From_PyStr(n, base, context);
        }

        if (PyObject_HasAttrString(n, "__mpz__")) {
             return (PyObject*) PyObject_CallMethod(n, "__mpz__", NULL);
        }

        /* Try converting to integer. */
        temp = PyNumber_Long(n);
        if (temp) {
            result = GMPy_MPZ_From_PyIntOrLong(temp, context);
            Py_DECREF(temp);
            return (PyObject*)result;
        }

        TYPE_ERROR("mpz() requires numeric or string argument");
        return NULL;
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, &n, &base)) {
        return NULL;
    }

    if ((base != 0) && ((base < 2)|| (base > 62))) {
        VALUE_ERROR("base for mpz() must be 0 or in the interval [2, 62]");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        return (PyObject*)GMPy_MPZ_From_PyStr(n, base, context);
    }

    if (IS_REAL(n)) {
        TYPE_ERROR("mpz() with number argument only takes 1 argument");
    }
    else {
        TYPE_ERROR("mpz() requires numeric or string (and optional base) arguments");
    }
    return NULL;
}

static void
GMPy_MPZ_Dealloc(MPZ_Object *self)
{
    if (global.in_gmpympzcache < global.cache_size &&
        self->z->_mp_alloc <= global.cache_obsize) {
        global.gmpympzcache[(global.in_gmpympzcache)++] = self;
    }
    else {
        mpz_clear(self->z);
        PyObject_Del(self);
    }
}

/* Caching logic for Pyxmpz. */

static void
set_gmpyxmpzcache(void)
{
    if (global.in_gmpyxmpzcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < global.in_gmpyxmpzcache; ++i) {
            mpz_clear(global.gmpyxmpzcache[i]->z);
            PyObject_Del(global.gmpyxmpzcache[i]);
        }
        global.in_gmpyxmpzcache = global.cache_size;
    }
    global.gmpyxmpzcache = realloc(global.gmpyxmpzcache, sizeof(XMPZ_Object)*global.cache_size);
}

static XMPZ_Object *
GMPy_XMPZ_New(CTXT_Object *context)
{
    XMPZ_Object *result = NULL;

    if (global.in_gmpyxmpzcache) {
        result = global.gmpyxmpzcache[--(global.in_gmpyxmpzcache)];
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)result);
        mpz_set_ui(result->z, 0);
    }
    else {
        if ((result = PyObject_New(XMPZ_Object, &XMPZ_Type))) {
            mpz_init(result->z);
        }
    }
    return result;
}

static PyObject *
GMPy_XMPZ_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds)
{
    XMPZ_Object *result = NULL;
    PyObject *n = NULL;
    PyObject *temp = NULL;
    int base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"s", "base", NULL };
    CTXT_Object *context = NULL;

    if (type != &XMPZ_Type) {
        TYPE_ERROR("xmpz.__new__() requires xmpz type");
        return NULL;
    }

    /* Optimize the most common use cases first; either 0 or 1 argument */

    argc = PyTuple_GET_SIZE(args);

    if (argc == 0) {
        return (PyObject*)GMPy_XMPZ_New(context);
    }

    if (argc == 1 && !keywds) {
        n = PyTuple_GET_ITEM(args, 0);

        if (XMPZ_Check(n)) {
            Py_INCREF(n);
            return n;
        }

        if (PyIntOrLong_Check(n)) {
            return (PyObject*)GMPy_XMPZ_From_PyIntOrLong(n, context);
        }

        if (MPQ_Check(n)) {
            return (PyObject*)GMPy_XMPZ_From_MPQ((MPQ_Object*)n, context);
        }

        if (MPFR_Check(n)) {
            return (PyObject*)GMPy_XMPZ_From_MPFR((MPFR_Object*)n, context);
        }

        if (PyFloat_Check(n)) {
            return (PyObject*)GMPy_XMPZ_From_PyFloat(n, context);
        }

        if (MPZ_Check(n)) {
            return (PyObject*)GMPy_XMPZ_From_MPZ((MPZ_Object*)n, context);
        }

        if (IS_FRACTION(n)) {
            MPQ_Object *temp = GMPy_MPQ_From_Fraction(n, context);

            if (temp) {
                result = GMPy_XMPZ_From_MPQ(temp, context);
                Py_DECREF((PyObject*)temp);
            }
            return (PyObject*)result;
        }

        if (PyStrOrUnicode_Check(n)) {
            return (PyObject*)GMPy_XMPZ_From_PyStr(n, base, context);
        }

        /* Try converting to integer. */
        temp = PyNumber_Long(n);
        if (temp) {
            result = GMPy_XMPZ_From_PyIntOrLong(temp, context);
            Py_DECREF(temp);
            return (PyObject*)result;
        }

        TYPE_ERROR("xmpz() requires numeric or string argument");
        return NULL;
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, &n, &base)) {
        return NULL;
    }

    if ((base != 0) && ((base < 2)|| (base > 62))) {
        VALUE_ERROR("base for xmpz() must be 0 or in the interval [2, 62]");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        return (PyObject*)GMPy_XMPZ_From_PyStr(n, base, context);
    }

    if (IS_REAL(n)) {
        TYPE_ERROR("xmpz() with number argument only takes 1 argument");
    }
    else {
        TYPE_ERROR("xmpz() requires numeric or string (and optional base) arguments");
    }
    return NULL;
}

static void
GMPy_XMPZ_Dealloc(XMPZ_Object *obj)
{
    if (global.in_gmpyxmpzcache < global.cache_size &&
        obj->z->_mp_alloc <= global.cache_obsize) {
        global.gmpyxmpzcache[(global.in_gmpyxmpzcache)++] = obj;
    }
    else {
        mpz_clear(obj->z);
        PyObject_Del((PyObject*)obj);
    }
}

/* Caching logic for Pympq. */

static void
set_gmpympqcache(void)
{
    if (global.in_gmpympqcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < global.in_gmpympqcache; ++i) {
            mpq_clear(global.gmpympqcache[i]->q);
            PyObject_Del(global.gmpympqcache[i]);
        }
        global.in_gmpympqcache = global.cache_size;
    }
    global.gmpympqcache = realloc(global.gmpympqcache, sizeof(MPQ_Object)*global.cache_size);
}

static MPQ_Object *
GMPy_MPQ_New(CTXT_Object *context)
{
    MPQ_Object *result = NULL;

    if (global.in_gmpympqcache) {
        result = global.gmpympqcache[--(global.in_gmpympqcache)];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)result);
    }
    else {
        if (!(result = PyObject_New(MPQ_Object, &MPQ_Type))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpq_init(result->q);
    }
    result->hash_cache = -1;
    return result;
}

static PyObject *
GMPy_MPQ_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds)
{
    MPQ_Object *result = NULL, *temp = NULL;
    PyObject *n = NULL, *m = NULL;
    int base = 10;
    Py_ssize_t argc, keywdc = 0;
    static char *kwlist[] = {"s", "base", NULL };
    CTXT_Object *context = NULL;

    argc = PyTuple_Size(args);
    if (keywds) {
        keywdc = PyDict_Size(keywds);
    }

    if (argc + keywdc > 2) {
        TYPE_ERROR("mpq() takes at most 2 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPQ_New(context))) {
            mpq_set_ui(result->q, 0, 1);
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpq() requires at least one non-keyword argument");
        return NULL;
    }

    n = PyTuple_GetItem(args, 0);

    /* Handle the case where the first argument is a string. */
    if (PyStrOrUnicode_Check(n)) {
        /* keyword base is legal */
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, &n, &base))) {
                return NULL;
            }
        }

        if ((base != 0) && ((base < 2) || (base > 62))) {
            VALUE_ERROR("base for mpq() must be 0 or in the interval [2, 62]");
            return NULL;
        }

        return (PyObject*)GMPy_MPQ_From_PyStr(n, base, context);
    }

    /* Handle 1 argument. It must be non-complex number or an object with a __mpq__ method. */
    if (argc == 1) {
        if (IS_REAL(n)) {
            return (PyObject*)GMPy_MPQ_From_Number(n, context);
        }

        if (PyObject_HasAttrString(n, "__mpq__")) {
             return (PyObject*) PyObject_CallMethod(n, "__mpq__", NULL);
        }

    }

    /* Handle 2 arguments. Both arguments must be integer or rational. */
    if (argc == 2) {
        m = PyTuple_GetItem(args, 1);

        if (IS_RATIONAL(n) && IS_RATIONAL(m)) {
           result = GMPy_MPQ_From_Rational(n, context);
           temp = GMPy_MPQ_From_Rational(m, context);
           if (!result || !temp) {
               Py_XDECREF((PyObject*)result);
               Py_XDECREF((PyObject*)temp);
               return NULL;
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
        }
    }

    TYPE_ERROR("mpq() requires numeric or string argument");
    return NULL;
}

static void
GMPy_MPQ_Dealloc(MPQ_Object *self)
{
    if (global.in_gmpympqcache<global.cache_size &&
        mpq_numref(self->q)->_mp_alloc <= global.cache_obsize &&
        mpq_denref(self->q)->_mp_alloc <= global.cache_obsize) {
        global.gmpympqcache[(global.in_gmpympqcache)++] = self;
    }
    else {
        mpq_clear(self->q);
        PyObject_Del(self);
    }
}

/* Caching logic for Pympfr. */

static void
set_gmpympfrcache(void)
{
    if (global.in_gmpympfrcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < global.in_gmpympfrcache; ++i) {
            mpfr_clear(global.gmpympfrcache[i]->f);
            PyObject_Del(global.gmpympfrcache[i]);
        }
        global.in_gmpympfrcache = global.cache_size;
    }
    global.gmpympfrcache = realloc(global.gmpympfrcache, sizeof(MPFR_Object)*global.cache_size);
}

static MPFR_Object *
GMPy_MPFR_New(mpfr_prec_t bits, CTXT_Object *context)
{
    MPFR_Object *result;

    if (bits < 2) {
        bits = GET_MPFR_PREC(context);
    }

    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }

    if (global.in_gmpympfrcache) {
        result = global.gmpympfrcache[--(global.in_gmpympfrcache)];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)result);
        mpfr_set_prec(result->f, bits);
    }
    else {
        if (!(result = PyObject_New(MPFR_Object, &MPFR_Type))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpfr_init2(result->f, bits);
    }
    result->hash_cache = -1;
    result->rc = 0;
    return result;
}

static void
GMPy_MPFR_Dealloc(MPFR_Object *self)
{
    size_t msize;

    /* Calculate the number of limbs in the mantissa. */
    msize = (self->f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    if (global.in_gmpympfrcache < global.cache_size &&
        msize <= (size_t)global.cache_obsize) {
        global.gmpympfrcache[(global.in_gmpympfrcache)++] = self;
    }
    else {
        mpfr_clear(self->f);
        PyObject_Del(self);
    }
}

static void
set_gmpympccache(void)
{
    if (global.in_gmpympccache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < global.in_gmpympccache; ++i) {
            mpc_clear(global.gmpympccache[i]->c);
            PyObject_Del(global.gmpympccache[i]);
        }
        global.in_gmpympccache = global.cache_size;
    }
    global.gmpympccache = realloc(global.gmpympccache, sizeof(MPC_Object)*global.cache_size);
}


static MPC_Object *
GMPy_MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context)
{
    MPC_Object *self;

    CHECK_CONTEXT(context);

    if (rprec < 2) {
        rprec = GET_REAL_PREC(context);
    }

    if (iprec < 2) {
        iprec = GET_IMAG_PREC(context);
    }

    if (rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
        iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (global.in_gmpympccache) {
        self = global.gmpympccache[--(global.in_gmpympccache)];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
        if (rprec == iprec) {
            mpc_set_prec(self->c, rprec);
        }
        else {
            mpc_clear(self->c);
            mpc_init3(self->c, rprec, iprec);
        }
    }
    else {
        if (!(self = PyObject_New(MPC_Object, &MPC_Type))) {
            /* LCOV_EXCL_START */
            return NULL;
            /* LCOV_EXCL_STOP */
        }
        mpc_init3(self->c, rprec, iprec);
    }
    self->hash_cache = -1;
    self->rc = 0;
    return self;
}

static void
GMPy_MPC_Dealloc(MPC_Object *self)
{
    size_t msize;

    /* Calculate the number of limbs in the mantissa. */
    msize = (mpc_realref(self->c)->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    msize += (mpc_imagref(self->c)->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    if (global.in_gmpympccache < global.cache_size &&
        msize <= (size_t)global.cache_obsize) {
        global.gmpympccache[(global.in_gmpympccache)++] = self;
    }
    else {
        mpc_clear(self->c);
        PyObject_Del(self);
    }
}

