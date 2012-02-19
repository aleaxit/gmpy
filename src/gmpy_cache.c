/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_cache.c                                                            *
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
    TRACE("Entering set_zcache\n");
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
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "Getting %d from zcache\n", in_zcache);
#endif
        newo[0] = (zcache[--in_zcache])[0];
    }
    else {
        TRACE("Initing new not in zcache\n");
        mpz_init(newo);
    }
}

static void
mpz_cloc(mpz_t oldo)
{
    if (in_zcache<global.cache_size && oldo->_mp_alloc <= global.cache_obsize) {
        (zcache[in_zcache++])[0] = oldo[0];
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "Stashed %d to zcache\n", in_zcache);
#endif
    }
    else {
#ifdef DEBUG
        if (global.debug)
            fprintf(stderr, "Not placing in full zcache(%d/%d)\n",
                    in_zcache, global.cache_size);
#endif
        mpz_clear(oldo);
    }
}

/* Caching logic for Pympz. */

static void
set_pympzcache(void)
{
    TRACE("Entering set_pympzcache\n");
    if (in_pympzcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_pympzcache; ++i) {
            mpz_cloc(pympzcache[i]->z);
            PyObject_Del(pympzcache[i]);
        }
        in_pympzcache = global.cache_size;
    }
    pympzcache = GMPY_REALLOC(pympzcache, sizeof(PympzObject)*global.cache_size);
}

static PympzObject *
Pympz_new(void)
{
    PympzObject *self;

    TRACE("Entering Pympz_new\n");
    if (in_pympzcache) {
        TRACE("Pympz_new is reusing an old object\n");
        self = pympzcache[--in_pympzcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        TRACE("Pympz_new is creating a new object\n");
        if (!(self = PyObject_New(PympzObject, &Pympz_Type)))
            return NULL;
        mpz_inoc(self->z);
    }
    self->hash_cache = -1;
    return self;
}

static void
Pympz_dealloc(PympzObject *self)
{
    TRACE("Pympz_dealloc\n");
    if (in_pympzcache < global.cache_size &&
        self->z->_mp_alloc <= global.cache_obsize) {
        pympzcache[in_pympzcache++] = self;
    }
    else {
        mpz_cloc(self->z);
        PyObject_Del(self);
    }
}

/* Caching logic for Pyxmpz. */

static void
set_pyxmpzcache(void)
{
    TRACE("Entering set_pyxmpzcache\n");
    if (in_pyxmpzcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_pyxmpzcache; ++i) {
            mpz_cloc(pyxmpzcache[i]->z);
            PyObject_Del(pyxmpzcache[i]);
        }
        in_pyxmpzcache = global.cache_size;
    }
    pyxmpzcache = GMPY_REALLOC(pyxmpzcache, sizeof(PyxmpzObject)*global.cache_size);
}

static PyxmpzObject *
Pyxmpz_new(void)
{
    PyxmpzObject *self;

    TRACE("Entering Pyxmpz_new\n");
    if (in_pyxmpzcache) {
        TRACE("Pyxmpz_new is reusing an old object\n");
        self = pyxmpzcache[--in_pyxmpzcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
         * _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        TRACE("Pyxmpz_new is creating a new object\n");
        if (!(self = PyObject_New(PyxmpzObject, &Pyxmpz_Type)))
            return NULL;
        mpz_inoc(self->z);
    }
    return self;
}

static void
Pyxmpz_dealloc(PyxmpzObject *self)
{
    TRACE("Pyxmpz_dealloc\n");
    if (in_pyxmpzcache < global.cache_size &&
        self->z->_mp_alloc <= global.cache_obsize) {
        pyxmpzcache[in_pyxmpzcache++] = self;
    }
    else {
        mpz_cloc(self->z);
        PyObject_Del(self);
    }
}

/* Caching logic for Pympq. */

static void
set_pympqcache(void)
{
    TRACE("Entering set_pympqcache\n");
    if (in_pympqcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_pympqcache; ++i) {
            mpq_clear(pympqcache[i]->q);
            PyObject_Del(pympqcache[i]);
        }
        in_pympqcache = global.cache_size;
    }
    pympqcache = GMPY_REALLOC(pympqcache, sizeof(PympqObject)*global.cache_size);
}

static PympqObject *
Pympq_new(void)
{
    PympqObject *self;

    TRACE("Entering Pympq_new\n");
    if (in_pympqcache) {
        TRACE("Pympq_new is reusing an old object\n");
        self = pympqcache[--in_pympqcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
    }
    else {
        TRACE("Pympq_new is creating a new object\n");
        if (!(self = PyObject_New(PympqObject, &Pympq_Type)))
            return NULL;
        mpq_init(self->q);
    }
    self->hash_cache = -1;
    return self;
}

static void
Pympq_dealloc(PympqObject *self)
{
    TRACE("Pympq_dealloc\n");
    if (in_pympqcache<global.cache_size &&
        mpq_numref(self->q)->_mp_alloc <= global.cache_obsize &&
        mpq_denref(self->q)->_mp_alloc <= global.cache_obsize) {
        pympqcache[in_pympqcache++] = self;
    }
    else {
        mpq_clear(self->q);
        PyObject_Del(self);
    }
}

#ifdef WITHMPFR
/* Caching logic for Pympfr. */

static void
set_pympfrcache(void)
{
    TRACE("Entering set_pympfrcache\n");
    if (in_pympfrcache > global.cache_size) {
        int i;
        for (i = global.cache_size; i < in_pympfrcache; ++i) {
            mpfr_clear(pympfrcache[i]->f);
            PyObject_Del(pympfrcache[i]);
        }
        in_pympfrcache = global.cache_size;
    }
    pympfrcache = GMPY_REALLOC(pympfrcache, sizeof(PympfrObject)*global.cache_size);
}

static PympfrObject *
Pympfr_new(mpfr_prec_t bits)
{
    PympfrObject *self;

    TRACE("Entering Pympfr_new\n");
    if (bits == 0)
        bits = context->now.mpfr_prec;
    if (bits < MPFR_PREC_MIN || bits > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (in_pympfrcache) {
        TRACE("Pympfr_new is reusing an old object\n");
        self = pympfrcache[--in_pympfrcache];
        /* Py_INCREF does not set the debugging pointers, so need to use
           _Py_NewReference instead. */
        _Py_NewReference((PyObject*)self);
        mpfr_set_prec(self->f, bits);
    }
    else {
        TRACE("Pympfr_new is creating a new object\n");
        if (!(self = PyObject_New(PympfrObject, &Pympfr_Type)))
            return NULL;
        mpfr_init2(self->f, bits);
    }
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = context->now.mpfr_round;
    return self;
}

static void
Pympfr_dealloc(PympfrObject *self)
{
    size_t msize;

    TRACE("Pympfr_dealloc\n");
    /* Calculate the number of limbs in the mantissa. */
    msize = (self->f->_mpfr_prec + mp_bits_per_limb - 1) / mp_bits_per_limb;
    if (in_pympfrcache < global.cache_size &&
        msize <= (size_t)global.cache_obsize) {
        pympfrcache[in_pympfrcache++] = self;
    }
    else {
        mpfr_clear(self->f);
        PyObject_Del(self);
    }
}
#endif

#ifdef WITHMPC
/* Caching logic for Pympc.
   No caching is done for Pympc since support for setting real and imaginary
   precision independantly forces a new memory allocation.
*/

static PympcObject *
Pympc_new(mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *self;

    TRACE("Entering Pympc_new\n");

    if (rprec == 0)
        rprec = GET_MPC_RPREC(context);
    if (iprec == 0)
        iprec = GET_MPC_IPREC(context);
    if (rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
        iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    if (!(self = PyObject_New(PympcObject, &Pympc_Type)))
        return NULL;
    mpc_init3(self->c, rprec, iprec);
    self->hash_cache = -1;
    self->rc = 0;
    self->round_mode = GET_MPC_ROUND(context);
    return self;
}

static void
Pympc_dealloc(PympcObject *self)
{
    mpc_clear(self->c);
    PyObject_Del(self);
}
#endif

