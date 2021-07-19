/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_cache.c                                                           *
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
    PyObject *out = NULL;
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

        if (HAS_MPZ_CONVERSION(n)) {
            out = (PyObject *) PyObject_CallMethod(n, "__mpz__", NULL);

            if (out == NULL)
                return out;
            if (!MPZ_Check(out)) {
                PyErr_Format(PyExc_TypeError,
                             "object of type '%.200s' can not be interpreted as mpz",
                             out->ob_type->tp_name);
                Py_DECREF(out);
                return NULL;
            }
            return out;
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

    if (type != &MPQ_Type) {
        TYPE_ERROR("mpq.__new__() requires mpq type");
        return NULL;
    }

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
            return (PyObject *) GMPy_MPQ_From_Number(n, context);
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
        CHECK_CONTEXT(context);
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

static PyObject *
GMPy_MPFR_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = NULL;
    Py_ssize_t argc, keywdc = 0;
    PyObject *arg0 = NULL;
    PyObject *out = NULL;
    int base = 0;


    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t prec = 0;

    static char *kwlist_s[] = {"s", "precision", "base", "context", NULL};
    static char *kwlist_n[] = {"n", "precision", "context", NULL};

    if (type != &MPFR_Type) {
        TYPE_ERROR("mpfr.__new__() requires mpfr type");
        return NULL;
    }

    CHECK_CONTEXT(context);

    argc = PyTuple_Size(args);
    if (keywds) {
        keywdc = PyDict_Size(keywds);
    }

    if (argc + keywdc > 4) {
        TYPE_ERROR("mpfr() takes at most 4 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPFR_New(0, context))) {
            mpfr_set_ui(result->f, 0, MPFR_RNDN);
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpfr() requires at least one non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    /* A string can have precision, base, and context as additional arguments. */
    if (PyStrOrUnicode_Check(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|liO", kwlist_s,
                                              &arg0, &prec, &base, &context)))
                return NULL;
        }

        if (!CTXT_Check(context)) {
            TYPE_ERROR("context argument is not a valid context");
            return NULL;
        }

        if (prec < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
            return NULL;
        }

        if (base != 0 && (base < 2 || base > 62)) {
            VALUE_ERROR("base for mpfr() must be 0 or in the interval [2, 62]");
            return NULL;
        }

        return (PyObject*)GMPy_MPFR_From_PyStr(arg0, base, prec, context);
    }

    if (HAS_MPFR_CONVERSION(arg0)) {
        out = (PyObject *) PyObject_CallMethod(arg0, "__mpfr__", NULL);

        if(out == NULL)
            return out;
        if (!MPFR_Check(out)) {
            PyErr_Format(PyExc_TypeError,
                         "object of type '%.200s' can not be interpreted as mpfr",
                         out->ob_type->tp_name);
            Py_DECREF(out);
            return NULL;
        }
        return out;
    }

    /* A number can only have precision and context as additional arguments. */
    if (IS_REAL(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|lO", kwlist_n,
                                              &arg0, &prec, &context)))
                return NULL;
        }

        if (!CTXT_Check(context)) {
            TYPE_ERROR("context argument is not a valid context");
            return NULL;
        }

        if (prec < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
            return NULL;
        }

        return (PyObject*)GMPy_MPFR_From_Real(arg0, prec, context);
    }

    TYPE_ERROR("mpfr() requires numeric or string argument");
    return NULL;
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

    if (rprec < 2) {
        CHECK_CONTEXT(context);
        rprec = GET_REAL_PREC(context);
    }

    if (iprec < 2) {
        CHECK_CONTEXT(context);
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

static PyObject *
GMPy_MPC_NewInit(PyTypeObject *type, PyObject *args, PyObject *keywds)
{
    MPC_Object *result = NULL;
    MPFR_Object *tempreal = NULL, *tempimag = NULL;
    PyObject *arg0 = NULL, *arg1 = NULL, *prec = NULL, *out = NULL;
    int base = 10;
    Py_ssize_t argc = 0, keywdc = 0;
    CTXT_Object *context = NULL;

    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t rprec = 0, iprec = 0;

    static char *kwlist_c[] = {"c", "precision", "context", NULL};
    static char *kwlist_r[] = {"real", "imag", "precision", "context", NULL};
    static char *kwlist_s[] = {"s", "precision", "base", "context", NULL};

    if (type != &MPC_Type) {
        TYPE_ERROR("mpc.__new__() requires mpc type");
        return NULL;
    }

    CHECK_CONTEXT(context);

    argc = PyTuple_Size(args);
    if (keywds) {
        keywdc = PyDict_Size(keywds);
    }

    if (argc + keywdc > 4) {
        TYPE_ERROR("mpc() takes at most 4 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPC_New(0, 0, context))) {
            mpc_set_ui(result->c, 0, GET_MPC_ROUND(context));
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpc() requires at least one non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    /* A string can have precision, base, and context as additional arguments.
     */

    if (PyStrOrUnicode_Check(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|OiO", kwlist_s,
                                              &arg0, &prec, &base, &context)))
                return NULL;
        }

        if (!CTXT_Check(context)) {
            TYPE_ERROR("context argument is not a valid context");
            return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            else {
                TYPE_ERROR("precision for mpc() must be integer or tuple");
                return NULL;
            }

            if (rprec < 0 || iprec < 0) {
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc()");
                }
                else {
                    VALUE_ERROR("precision for mpc() must be >= 0");
                }
                return NULL;
            }
        }

        if (base < 2 || base > 36) {
            VALUE_ERROR("base for mpc() must be in the interval [2,36]");
            return NULL;
        }

        return (PyObject*)GMPy_MPC_From_PyStr(arg0, base, rprec, iprec, context);
    }

    if (HAS_MPC_CONVERSION(arg0)) {
        out = (PyObject*) PyObject_CallMethod(arg0, "__mpc__", NULL);
        if(out == NULL)
            return out;
        if (!MPC_Check(out)) {
            PyErr_Format(PyExc_TypeError,
                         "object of type '%.200s' can not be interpreted as mpc",
                         out->ob_type->tp_name);
            Py_DECREF(out);
            return NULL;
        }
        return out;
    }

    /* Should special case PyFLoat to avoid double rounding. */

    if (IS_REAL(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|OOO", kwlist_r,
                                            &arg0, &arg1, &prec, &context)))
                return NULL;
        }

        if (!CTXT_Check(context)) {
            TYPE_ERROR("context argument is not a valid context");
            return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            else {
                TYPE_ERROR("precision for mpc() must be integer or tuple");
                return NULL;
            }

            if (rprec < 0 || iprec < 0) {
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc()");
                }
                else {
                    VALUE_ERROR("precision for mpc() must be >= 0");
                }
                return NULL;
            }
        }

        if (arg1 && !IS_REAL(arg1)) {
            TYPE_ERROR("invalid type for imaginary component in mpc()");
            return NULL;
        }

        tempreal = GMPy_MPFR_From_Real(arg0, rprec, context);
        if (arg1) {
            tempimag = GMPy_MPFR_From_Real(arg1, iprec, context);
        }
        else {
            if ((tempimag = GMPy_MPFR_New(iprec, context))) {
                mpfr_set_ui(tempimag->f, 0, MPFR_RNDN);
            }
        }

        result = GMPy_MPC_New(rprec, iprec, context);
        if (!tempreal || !tempimag || !result) {
            Py_XDECREF(tempreal);
            Py_XDECREF(tempimag);
            Py_XDECREF(result);
            TYPE_ERROR("mpc() requires string or numeric argument.");
            return NULL;
        }

        mpc_set_fr_fr(result->c, tempreal->f, tempimag->f, GET_MPC_ROUND(context));
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
        return (PyObject*)result;
    }

    if (IS_COMPLEX_ONLY(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|O", kwlist_c,
                                          &arg0, &prec)))
            return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            else {
                TYPE_ERROR("precision for mpc() must be integer or tuple");
                return NULL;
            }

            if (rprec < 0 || iprec < 0) {
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc()");
                }
                else {
                    VALUE_ERROR("precision for mpc() must be >= 0");
                }
                return NULL;
            }
        }

        if (PyComplex_Check(arg0)) {
            result = GMPy_MPC_From_PyComplex(arg0, rprec, iprec, context);
        }
        else {
            result = GMPy_MPC_From_MPC((MPC_Object*)arg0, rprec, iprec, context);
        }
        return (PyObject*)result;
    }

    TYPE_ERROR("mpc() requires numeric or string argument");
    return NULL;
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

