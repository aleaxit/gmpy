/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_cache.h                                                            *
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

#ifndef GMPY_CACHE_H
#define GMPY_CACHE_H

#ifdef __cplusplus
extern "C" {
#endif

static void          set_zcache(void);
static void          mpz_inoc(mpz_t newo);
static void          mpz_cloc(mpz_t oldo);

static void          set_gmpympzcache(void);
static MPZ_Object *  GMPy_MPZ_New(CTXT_Object *context);
static void          GMPy_MPZ_Dealloc(MPZ_Object *self);

static void          set_gmpyxmpzcache(void);
static XMPZ_Object * GMPy_XMPZ_New(CTXT_Object *context);
static void          GMPy_XMPZ_Dealloc(XMPZ_Object *self);

static void          set_gmpympqcache(void);
static MPQ_Object *  GMPy_MPQ_New(CTXT_Object *context);
static void          GMPy_MPQ_Dealloc(MPQ_Object *self);

static void          set_gmpympfrcache(void);
static MPFR_Object * GMPy_MPFR_New(mpfr_prec_t bits, CTXT_Object *context);
static void          GMPy_MPFR_Dealloc(MPFR_Object *self);

static void          set_gmpympccache(void);
static MPC_Object *  GMPy_MPC_New(mpfr_prec_t rprec, mpfr_prec_t iprec, CTXT_Object *context);
static void          GMPy_MPC_Dealloc(MPC_Object *self);

#ifdef __cplusplus
}
#endif
#endif
