/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.h                                                              *
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

/* Verify that an object is an mpc and that both components have valid exp */
#define MPC_CheckAndExp(v) \
    (MPC_Check(v) && \
        (mpfr_zero_p(mpc_realref(MPC(v))) || \
            (mpfr_regular_p(mpc_realref(MPC(v))) && \
                ((mpc_realref(MPC(v))->_mpfr_exp >= context->ctx.emin)) && \
                ((mpc_realref(MPC(v))->_mpfr_exp <= context->ctx.emax)) \
            ) \
        ) && \
        (mpfr_zero_p(mpc_imagref(MPC(v))) || \
            (mpfr_regular_p(mpc_imagref(MPC(v))) && \
                ((mpc_imagref(MPC(v))->_mpfr_exp >= context->ctx.emin)) && \
                ((mpc_imagref(MPC(v))->_mpfr_exp <= context->ctx.emax)) \
            ) \
        ) \
    )

#define MPC_CHECK_UNDERFLOW(mpct, msg) \
    if (MPC_IS_ZERO_P(mpct) && mpct->rc) { \
        context->ctx.underflow = 1; \
        if (context->ctx.traps & TRAP_UNDERFLOW) { \
            GMPY_UNDERFLOW(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_OVERFLOW(mpct, msg) \
    if (MPC_IS_INF_P(mpct)) { \
        context->ctx.overflow = 1; \
        if (context->ctx.traps & TRAP_OVERFLOW) { \
            GMPY_OVERFLOW(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_INVALID(mpct, msg) \
    if (MPC_IS_NAN_P(mpct)) { \
        context->ctx.invalid = 1; \
        if (context->ctx.traps & TRAP_INVALID) { \
            GMPY_INVALID(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_INEXACT(mpct, msg) \
    if (mpct->rc) { \
        context->ctx.inexact = 1; \
        if (context->ctx.traps & TRAP_INEXACT) { \
            GMPY_INEXACT(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_FLAGS(mpct, NAME) \
    MPC_CHECK_INVALID(mpct, "'mpc' invalid operation in "NAME); \
    MPC_CHECK_UNDERFLOW(mpct, "'mpc' underflow in "NAME); \
    MPC_CHECK_OVERFLOW(mpct, "'mpc' overflow in "NAME); \
    MPC_CHECK_INEXACT(mpct, "'mpc' inexact result in "NAME);

#define MPC_CHECK_FLAGS_RESULT(NAME) \
    if (MPC_IS_NAN_P(result)) { \
        context->ctx.invalid = 1; \
        if (context->ctx.traps & TRAP_INVALID) { \
            GMPY_INVALID("'mpc' invalid operation "NAME); \
            Py_DECREF((PyObject*)result); \
            return NULL; \
        } \
    } \
    if (MPC_IS_ZERO_P(result) && result->rc) { \
        context->ctx.underflow = 1; \
        if (context->ctx.traps & TRAP_UNDERFLOW) { \
            GMPY_UNDERFLOW("'mpc' underflow in "NAME); \
            Py_DECREF((PyObject*)result); \
            return NULL; \
        } \
    } \
    if (MPC_IS_INF_P(result)) { \
        context->ctx.overflow = 1; \
        if (context->ctx.traps & TRAP_OVERFLOW) { \
            GMPY_OVERFLOW("'mpc' overflow in "NAME); \
            Py_DECREF((PyObject*)result); \
            return NULL; \
        } \
    } \
    if (result->rc) { \
        context->ctx.inexact = 1; \
        if (context->ctx.traps & TRAP_INEXACT) { \
            GMPY_INEXACT("'mpc' inexact result in "NAME); \
            Py_DECREF((PyObject*)result); \
            return NULL; \
        } \
    } \


#define MPC_SUBNORMALIZE(mpct) \
    if (context->ctx.subnormalize) { \
        int rcr, rci; \
        rcr = MPC_INEX_RE(mpct->rc); \
        rci = MPC_INEX_IM(mpct->rc); \
        rcr = mpfr_subnormalize(mpc_realref(mpct->c), rcr, GET_REAL_ROUND(context)); \
        rci = mpfr_subnormalize(mpc_imagref(mpct->c), rci, GET_IMAG_ROUND(context)); \
        mpct->rc = MPC_INEX(rcr, rci); \
    } \

#define MPC_CLEANUP(mpct, NAME) \
    MPC_SUBNORMALIZE(mpct); \
    MPC_CHECK_FLAGS(mpct, NAME); \
  done:\
    if (PyErr_Occurred()) { \
        Py_DECREF((PyObject*)mpct); \
        mpct = NULL; \
    } \
    return (PyObject*)mpct;

#define MPC_CLEANUP_RESULT(NAME) \
    MPC_SUBNORMALIZE(result); \
    MPC_CHECK_FLAGS_RESULT(NAME);

#define MPC_SUBNORMALIZE_2(V, CTX) \
    if (CTX->ctx.subnormalize) { \
        int rcr, rci; \
        rcr = MPC_INEX_RE(V->rc); \
        rci = MPC_INEX_IM(V->rc); \
        rcr = mpfr_subnormalize(mpc_realref(V->c), rcr, GET_REAL_ROUND(CTX)); \
        rci = mpfr_subnormalize(mpc_imagref(V->c), rci, GET_IMAG_ROUND(CTX)); \
        V->rc = MPC_INEX(rcr, rci); \
    } \

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

#define MPC_CLEANUP_2(V, CTX, NAME) \
    MPC_SUBNORMALIZE_2(V, CTX); \
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
        if ((rcr && mpfr_zero_p(mpc_realref(MPC(V)))) || (rci && mpfr_zero_p(mpc_imagref(MPC(V))))) { \
            CTX->ctx.underflow = 1; \
            _underflow = 1; \
        } \
        if ((rcr && mpfr_inf_p(mpc_realref(MPC(V)))) || (rci && mpfr_inf_p(mpc_imagref(MPC(V))))) { \
            CTX->ctx.overflow = 1; \
            _overflow = 1; \
        } \
        if (_invalid && (CTX->ctx.traps & TRAP_INVALID)) { \
            GMPY_INVALID(NAME" invalid operation"); \
            Py_DECREF((PyObject*)V); \
            return NULL; \
        } \
        if (_inexact && (CTX->ctx.traps & TRAP_INEXACT)) { \
            GMPY_INEXACT(NAME" inexact result"); \
            Py_DECREF((PyObject*)V); \
            return NULL; \
        } \
        if (_underflow && (CTX->ctx.traps & TRAP_UNDERFLOW)) { \
            GMPY_UNDERFLOW(NAME" underflow"); \
            Py_DECREF((PyObject*)V); \
            return NULL; \
        } \
        if (_overflow && (CTX->ctx.traps & TRAP_OVERFLOW)) { \
            GMPY_OVERFLOW(NAME" overflow"); \
            Py_DECREF((PyObject*)V); \
            return NULL; \
        } \
    } while(0); \

/*
 * Parses one, and only one, argument into "self" and converts it to an
 * mpc. Is faster, but not as generic, as using PyArg_ParseTuple. It
 * supports either gmpy.fname(z) or z.fname(). "self" must be decref'ed.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments.
 */

#define PARSE_ONE_MPC_ARGS(msg) \
    if(self && MPC_Check(self)) { \
        if (PyTuple_GET_SIZE(args) != 0) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        if (MPC_CheckAndExp(self)) { \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)GMPy_MPC_From_Complex(self, 0, 0, context))) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
        } \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) != 1) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        self = PyTuple_GET_ITEM(args, 0);\
        if (MPC_CheckAndExp(self)) { \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)GMPy_MPC_From_Complex(self, 0, 0, context))) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
        } \
    }

/*
 * Parses one, and only one, argument into "self" and converts it to an
 * mpc. Is faster, but not as generic, as using PyArg_ParseTuple. It
 * supports either gmpy.fname(z) or z.fname(). "self" must be decref'ed.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. It assumes the functions is declared
 * as either METH_O or METH_NOARGS. It is faster than PARSE_ONE_MPFR and
 * passing a tuple as args.
 */

#define PARSE_ONE_MPC_OTHER(msg) \
    if(self && MPC_Check(self)) { \
        if (MPC_CheckAndExp(self)) { \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)GMPy_MPC_From_Complex(self, 0, 0, context))) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
        } \
    } \
    else { \
        if (MPC_CheckAndExp(other)) { \
            self = other; \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)GMPy_MPC_From_Complex(other, 0, 0, context))) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
        } \
    }

/*
 * Parses two, and only two, arguments into "self" and "var" and converts
 * them both to mpC. Is faster, but not as generic, as using PyArg_ParseTuple.
 * It supports either gmpy.fname(f,f) or f.fname(f). "self" & "var" must be
 * decref'ed after use. "msg" should be an error message that includes the
 * function name and describes the required arguments. Replaces
 * SELF_MPF_ONE_ARG_CONVERTED(var).
 */

#define PARSE_TWO_MPC_ARGS(var, msg) \
    if(self && MPC_Check(self)) { \
        if (PyTuple_GET_SIZE(args) != 1) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        self = (PyObject*)GMPy_MPC_From_Complex_Temp(self, 0, 0, context); \
        var = (PyObject*)GMPy_MPC_From_Complex_Temp(PyTuple_GET_ITEM(args, 0), 0, 0, context); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) != 2) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        self = (PyObject*)GMPy_MPC_From_Complex_Temp(PyTuple_GET_ITEM(args, 0), 0, 0, context); \
        var = (PyObject*)GMPy_MPC_From_Complex_Temp(PyTuple_GET_ITEM(args, 1), 0, 0, context); \
    } \
    if (!self || !var) { \
        TYPE_ERROR(msg); \
        Py_XDECREF((PyObject*)var); \
        Py_XDECREF((PyObject*)self); \
        return NULL; \
    }

static PyObject * GMPy_MPC_Factory(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * Pympc_conjugate(PyObject *self, PyObject *args);
static PyObject * Pympc_getprec_attrib(MPC_Object *self, void *closure);
static PyObject * Pympc_getrc_attrib(MPC_Object *self, void *closure);
static PyObject * Pympc_getimag_attrib(MPC_Object *self, void *closure);
static PyObject * Pympc_getreal_attrib(MPC_Object *self, void *closure);
static int        Pympc_nonzero(MPC_Object *self);
static PyObject * Pympc_phase(PyObject *self, PyObject *other);
static PyObject * Pympc_norm(PyObject *self, PyObject *other);
static PyObject * Pympc_polar(PyObject *self, PyObject *other);
static PyObject * Pympc_rect(PyObject *self, PyObject *args);
static PyObject * Pympc_proj(PyObject *self, PyObject *other);

#ifdef __cplusplus
}
#endif
#endif
