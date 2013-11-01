/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_cache.c                                                            *
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


/* gmpy2 caches objects so they can be reused quickly without involving a new
 * memory allocation or object construction. There are two different types of
 * object caches used in gmpy2.
 *
 * "zcache" is used to cache mpz_t objects. The cache is accessed via the
 * functions mpz_inoc/mpz_cloc. The function set_zcache is used to change
 * the size of the array used to store the cached objects.
 *
 * The "py???cache" is used to cache Py??? objects. The cache is accessed
 * via Py???_new/Py???_dealloc. The functions set_py???cache and
 * set_py???cache are used to change the size of the array used to the store
 * the cached objects.
 */

static void
set_zcache(void)
{
    if (in_zcache > global.cache_size) {
        int i;
        for(i = global.cache_size; i < in_zcache; ++i)
            mpz_clear(zcache[i]);
        in_zcache = global.cache_size;
    }
    zcache = GMPY_REALLOC(zcache, sizeof(mpz_t) * global.cache_size);
}

static void
mpz_inoc(mpz_t newo)
{
    if (in_zcache) {
        newo[0] = (zcache[--in_zcache])[0];
    }
    else {
        mpz_init(newo);
    }
}

static void
mpz_cloc(mpz_t oldo)
{
    if (in_zcache<global.cache_size && oldo->_mp_alloc <= global.cache_obsize) {
        (zcache[in_zcache++])[0] = oldo[0];
    }
    else {
        mpz_clear(oldo);
    }
}

/* Caching logic for Pympz. */

static void
set_gmpympzcache(void)
{
    if (in_gmpympzcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_gmpympzcache; ++i) {
            mpz_cloc(gmpympzcache[i]->z);
            PyObject_Del(gmpympzcache[i]);
        }
        in_gmpympzcache = global.cache_size;
    }
    gmpympzcache = GMPY_REALLOC(gmpympzcache, sizeof(MPZ_Object)*global.cache_size);
}

static PyObject *
Pympz_new(void)
{
    MPZ_Object *self;

    if (in_gmpympzcache) {
        self = gmpympzcache[--in_gmpympzcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        if (!(self = PyObject_New(MPZ_Object, &MPZ_Type)))
            return NULL;
        mpz_inoc(self->z);
    }
    self->hash_cache = -1;
    return (PyObject*)self;
}

static void
Pympz_dealloc(MPZ_Object *self)
{
    if (in_gmpympzcache < global.cache_size &&
        self->z->_mp_alloc <= global.cache_obsize) {
        gmpympzcache[in_gmpympzcache++] = self;
    }
    else {
        mpz_cloc(self->z);
        PyObject_Del(self);
    }
}

/* Caching logic for Pyxmpz. */

static void
set_gmpyxmpzcache(void)
{
    if (in_gmpyxmpzcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_gmpyxmpzcache; ++i) {
            mpz_cloc(gmpyxmpzcache[i]->z);
            PyObject_Del(gmpyxmpzcache[i]);
        }
        in_gmpyxmpzcache = global.cache_size;
    }
    gmpyxmpzcache = GMPY_REALLOC(gmpyxmpzcache, sizeof(XMPZ_Object)*global.cache_size);
}

static PyObject *
Pyxmpz_new(void)
{
    XMPZ_Object *self;

    if (in_gmpyxmpzcache) {
        self = gmpyxmpzcache[--in_gmpyxmpzcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        if (!(self = PyObject_New(XMPZ_Object, &XMPZ_Type)))
            return NULL;
        mpz_inoc(self->z);
    }
    return (PyObject*)self;
}

static void
Pyxmpz_dealloc(XMPZ_Object *self)
{
    if (in_gmpyxmpzcache < global.cache_size &&
        self->z->_mp_alloc <= global.cache_obsize) {
        gmpyxmpzcache[in_gmpyxmpzcache++] = self;
    }
    else {
        mpz_cloc(self->z);
        PyObject_Del(self);
    }
}

/* Caching logic for Pympq. */

static void
set_gmpympqcache(void)
{
    if (in_gmpympqcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_gmpympqcache; ++i) {
            mpq_clear(gmpympqcache[i]->q);
            PyObject_Del(gmpympqcache[i]);
        }
        in_gmpympqcache = global.cache_size;
    }
    gmpympqcache = GMPY_REALLOC(gmpympqcache, sizeof(MPQ_Object)*global.cache_size);
}

static PyObject *
Pympq_new(void)
{
    MPQ_Object *self;

    if (in_gmpympqcache) {
        self = gmpympqcache[--in_gmpympqcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        if (!(self = PyObject_New(MPQ_Object, &MPQ_Type)))
            return NULL;
        mpq_init(self->q);
    }
    self->hash_cache = -1;
    return (PyObject*)self;
}

static void
Pympq_dealloc(MPQ_Object *self)
{
    if (in_gmpympqcache<global.cache_size &&
        mpq_numref(self->q)->_mp_alloc <= global.cache_obsize &&
        mpq_denref(self->q)->_mp_alloc <= global.cache_obsize) {
        gmpympqcache[in_gmpympqcache++] = self;
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
    if (in_gmpympfrcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_gmpympfrcache; ++i) {
            mpfr_clear(gmpympfrcache[i]->f);
            PyObject_Del(gmpympfrcache[i]);
        }
        in_gmpympfrcache = global.cache_size;
    }
    gmpympfrcache = GMPY_REALLOC(gmpympfrcache, sizeof(MPFR_Object)*global.cache_size);
}

static PyObject *
Pympfr_new(mpfr_prec_t bits)
{
    MPFR_Object *self;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (!bits)
        bits = GET_MPFR_PREC(context);

    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_gmpympfrcache) {
        self = gmpympfrcache[--in_gmpympfrcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
        mpfr_set_prec(self->f, bits);
    }
    else {
        if (!(self = PyObject_New(MPFR_Object, &MPFR_Type)))
            return NULL;
        mpfr_init2(self->f, bits);
    }
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = context->ctx.mpfr_round;
    return (PyObject*)self;
}

static PyObject *
Pympfr_new_bits_context(mpfr_prec_t bits, GMPyContextObject *context)
{
    MPFR_Object *self;

    if (!bits)
        bits = context->ctx.mpfr_prec;

    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_gmpympfrcache) {
        self = gmpympfrcache[--in_gmpympfrcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
        mpfr_set_prec(self->f, bits);
    }
    else {
        if (!(self = PyObject_New(MPFR_Object, &MPFR_Type)))
            return NULL;
        mpfr_init2(self->f, bits);
    }
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = context->ctx.mpfr_round;
    return (PyObject*)self;
}

static PyObject *
Pympfr_new_context(GMPyContextObject *context)
{
    mpfr_prec_t bits;
    MPFR_Object *result;

    bits = GET_MPFR_PREC(context);

    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_gmpympfrcache) {
        result = gmpympfrcache[--in_gmpympfrcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)result);
        mpfr_set_prec(result->f, bits);
    }
    else {
        if (!(result = PyObject_New(MPFR_Object, &MPFR_Type)))
            return NULL;
        mpfr_init2(result->f, bits);
    }
    result->hash_cache = -1;
    result->rc = 0;
    result->round_mode = GET_MPFR_ROUND(context);
    return (PyObject*)result;
}

static void
Pympfr_dealloc(MPFR_Object *self)
{
    size_t msize;

    /* Calculate the number of limbs in the mantissa. */
    msize = (self->f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    if (in_gmpympfrcache < global.cache_size &&
        msize <= (size_t)global.cache_obsize) {
        gmpympfrcache[in_gmpympfrcache++] = self;
    }
    else {
        mpfr_clear(self->f);
        PyObject_Del(self);
    }
}

static void
set_gmpympccache(void)
{
    if (in_gmpympccache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_gmpympccache; ++i) {
            mpc_clear(gmpympccache[i]->c);
            PyObject_Del(gmpympccache[i]);
        }
        in_gmpympccache = global.cache_size;
    }
    gmpympccache = GMPY_REALLOC(gmpympccache, sizeof(MPC_Object)*global.cache_size);
}


static PyObject *
Pympc_new(mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *self;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (!rprec) rprec = GET_REAL_PREC(context);
    if (!iprec) iprec = GET_IMAG_PREC(context);
    if (rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
        iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_gmpympccache) {
        self = gmpympccache[--in_gmpympccache];
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
        if (!(self = PyObject_New(MPC_Object, &MPC_Type)))
            return NULL;
        mpc_init3(self->c, rprec, iprec);
    }
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = GET_MPC_ROUND(context);
    return (PyObject*)self;
}

static PyObject *
Pympc_new_bits_context(mpfr_prec_t rprec, mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *self;

    if (!rprec) rprec = GET_REAL_PREC(context);
    if (!iprec) iprec = GET_IMAG_PREC(context);
    if (rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
        iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_gmpympccache) {
        self = gmpympccache[--in_gmpympccache];
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
        if (!(self = PyObject_New(MPC_Object, &MPC_Type)))
            return NULL;
        mpc_init3(self->c, rprec, iprec);
    }
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = GET_MPC_ROUND(context);
    return (PyObject*)self;
}

static PyObject *
Pympc_new_context(GMPyContextObject *context)
{
    mpfr_prec_t rprec, iprec;
    MPC_Object *self;

    rprec = GET_REAL_PREC(context);
    iprec = GET_IMAG_PREC(context);
    if (rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
        iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_gmpympccache) {
        self = gmpympccache[--in_gmpympccache];
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
        if (!(self = PyObject_New(MPC_Object, &MPC_Type)))
            return NULL;
        mpc_init3(self->c, rprec, iprec);
    }
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = GET_MPC_ROUND(context);
    return (PyObject*)self;
}

static void
Pympc_dealloc(MPC_Object *self)
{
    size_t msize;

    /* Calculate the number of limbs in the mantissa. */
    msize = (mpc_realref(self->c)->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    msize += (mpc_imagref(self->c)->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    if (in_gmpympccache < global.cache_size &&
        msize <= (size_t)global.cache_obsize) {
        gmpympccache[in_gmpympccache++] = self;
    }
    else {
        mpc_clear(self->c);
        PyObject_Del(self);
    }
}

