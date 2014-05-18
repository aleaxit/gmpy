/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpc.h                                                             *
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

#ifndef GMPY_MPC_H
#define GMPY_MPC_H

#ifdef __cplusplus
extern "C" {
#endif


/* gmpy_mpc C API extension header file.
 *
 * Provide interface to the MPC (Multiple Precision Complex) library.
 *
 * Version 2.00, April 2011 (created) casevh
 *
 * This file is expected to be included from gmpy.h
 */

#if defined(MS_WIN32) && defined(_MSC_VER)
#  pragma comment(lib,"mpc.lib")
#endif

typedef struct {
    PyObject_HEAD
    mpc_t c;
    Py_hash_t hash_cache;
    int rc;
    int round_mode;
} MPC_Object;

static PyTypeObject MPC_Type;
#define MPC(obj) (((MPC_Object *)(obj))->c)
#define MPC_Check(v) (((PyObject*)v)->ob_type == &MPC_Type)

/*
 * Define macros for comparing with zero, checking if either component is
 * 'nan' or 'inf', etc.
 */

#define MPC_IS_ZERO_P(x) \
    (mpfr_zero_p(mpc_realref(MPC(x))) && \
     mpfr_zero_p(mpc_imagref(MPC(x))))

#define MPC_IS_NAN_P(x) \
    (mpfr_nan_p(mpc_realref(MPC(x))) || \
     mpfr_nan_p(mpc_imagref(MPC(x))))

#define MPC_IS_INF_P(x) \
    (mpfr_inf_p(mpc_realref(MPC(x))) || \
     mpfr_inf_p(mpc_imagref(MPC(x))))

#define MPC_IS_FINITE_P(x) \
    (mpfr_number_p(mpc_realref(MPC(x))) && \
     mpfr_number_p(mpc_imagref(MPC(x))))

/* V is the value that is expected to be returned.
 * CTX is the context.
 * NAME is prepended to the error message.
 */

#define GMPY_MPC_CHECK_RANGE(V, CTX) \
    { \
        int rcr, rci; \
        rcr = MPC_INEX_RE(V->rc); \
        rci = MPC_INEX_IM(V->rc); \
        if (mpfr_regular_p(mpc_realref(V->c)) && \
            (!((mpc_realref(V->c)->_mpfr_exp >= CTX->ctx.emin) && \
               (mpc_realref(V->c)->_mpfr_exp <= CTX->ctx.emax)))) { \
            mpfr_exp_t _oldemin, _oldemax; \
            _oldemin = mpfr_get_emin(); \
            _oldemax = mpfr_get_emax(); \
            mpfr_set_emin(CTX->ctx.emin); \
            mpfr_set_emax(CTX->ctx.emax); \
            rcr = mpfr_check_range(mpc_realref(V->c), rcr, GET_REAL_ROUND(CTX)); \
            mpfr_set_emin(_oldemin); \
            mpfr_set_emax(_oldemax); \
        } \
        if (mpfr_regular_p(mpc_imagref(V->c)) && \
            (!((mpc_imagref(V->c)->_mpfr_exp >= CTX->ctx.emin) && \
               (mpc_imagref(V->c)->_mpfr_exp <= CTX->ctx.emax)))) { \
            mpfr_exp_t _oldemin, _oldemax; \
            _oldemin = mpfr_get_emin(); \
            _oldemax = mpfr_get_emax(); \
            mpfr_set_emin(CTX->ctx.emin); \
            mpfr_set_emax(CTX->ctx.emax); \
            rci = mpfr_check_range(mpc_imagref(V->c), rci, GET_IMAG_ROUND(CTX)); \
            mpfr_set_emin(_oldemin); \
            mpfr_set_emax(_oldemax); \
        } \
        V->rc = MPC_INEX(rcr, rci); \
    }

#define GMPY_MPC_SUBNORMALIZE(V, CTX) \
    { \
        int rcr, rci; \
        rcr = MPC_INEX_RE(V->rc); \
        rci = MPC_INEX_IM(V->rc); \
        if (CTX->ctx.subnormalize && \
            (!((mpc_realref(V->c)->_mpfr_exp >= CTX->ctx.emin) && \
               (mpc_realref(V->c)->_mpfr_exp <= CTX->ctx.emin + mpfr_get_prec(mpc_realref(V->c)) - 2)))) { \
            mpfr_exp_t _oldemin, _oldemax; \
            _oldemin = mpfr_get_emin(); \
            _oldemax = mpfr_get_emax(); \
            mpfr_set_emin(CTX->ctx.emin); \
            mpfr_set_emax(CTX->ctx.emax); \
            rcr = mpfr_subnormalize(mpc_realref(V->c), rcr, GET_REAL_ROUND(CTX)); \
            mpfr_set_emin(_oldemin); \
            mpfr_set_emax(_oldemax); \
        } \
        if (CTX->ctx.subnormalize && \
            (!((mpc_imagref(V->c)->_mpfr_exp >= CTX->ctx.emin) && \
               (mpc_imagref(V->c)->_mpfr_exp <= CTX->ctx.emin + mpfr_get_prec(mpc_imagref(V->c)) - 2)))) { \
            mpfr_exp_t _oldemin, _oldemax; \
            _oldemin = mpfr_get_emin(); \
            _oldemax = mpfr_get_emax(); \
            mpfr_set_emin(CTX->ctx.emin); \
            mpfr_set_emax(CTX->ctx.emax); \
            rci = mpfr_check_range(mpc_imagref(V->c), rci, GET_IMAG_ROUND(CTX)); \
            mpfr_set_emin(_oldemin); \
            mpfr_set_emax(_oldemax); \
        } \
        V->rc = MPC_INEX(rcr, rci); \
    }

#define GMPY_MPC_EXCEPTIONS(V, CTX, NAME) \
    do { \
        int _invalid = 0, _underflow = 0, _overflow = 0, _inexact = 0; \
        int rcr, rci; \
        rcr = MPC_INEX_RE(V->rc); \
        rci = MPC_INEX_IM(V->rc); \
        if (MPC_IS_NAN_P(V)) { \
            CTX->ctx.invalid = 1; \
            _invalid = 1; \
        } \
        if (V->rc) { \
            CTX->ctx.inexact = 1; \
            _inexact = 1; \
        } \
        if ((rcr && mpfr_zero_p(mpc_realref(V->c))) || (rci && mpfr_zero_p(mpc_imagref(V->c)))) { \
            CTX->ctx.underflow = 1; \
            _underflow = 1; \
        } \
        if ((rcr && mpfr_inf_p(mpc_realref(V->c))) || (rci && mpfr_inf_p(mpc_imagref(V->c)))) { \
            CTX->ctx.overflow = 1; \
            _overflow = 1; \
        } \
        if (CTX->ctx.traps) { \
            if ((CTX->ctx.traps & TRAP_UNDERFLOW) && _underflow) { \
                GMPY_UNDERFLOW(NAME" underflow"); \
                Py_XDECREF((PyObject*)V); \
                V = NULL; \
            } \
            if ((CTX->ctx.traps & TRAP_OVERFLOW) && _overflow) { \
                GMPY_OVERFLOW(NAME" overflow"); \
                Py_XDECREF((PyObject*)V); \
                V = NULL; \
            } \
            if ((CTX->ctx.traps & TRAP_INEXACT) && _inexact) { \
                GMPY_INEXACT(NAME" inexact result"); \
                Py_XDECREF((PyObject*)V); \
                V = NULL; \
            } \
            if ((CTX->ctx.traps & TRAP_INVALID) && _invalid) { \
                GMPY_INVALID(NAME" invalid operation"); \
                Py_XDECREF((PyObject*)V); \
                V = NULL; \
            } \
        } \
    } while(0); \

#define GMPY_MPC_CLEANUP(V, CTX, NAME) \
    GMPY_MPC_CHECK_RANGE(V, CTX); \
    GMPY_MPC_SUBNORMALIZE(V, CTX); \
    GMPY_MPC_EXCEPTIONS(V, CTX, NAME); \

static PyObject * GMPy_MPC_Factory(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args);
static PyObject * GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure);
static PyObject * GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure);
static int        GMPy_MPC_NonZero_Slot(MPC_Object *self);

#ifdef __cplusplus
}
#endif
#endif
