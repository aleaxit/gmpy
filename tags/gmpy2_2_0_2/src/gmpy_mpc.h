/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.h                                                              *
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
} PympcObject;

static PyTypeObject Pympc_Type;
#define Pympc_AS_MPC(obj) (((PympcObject *)(obj))->c)
#define Pympc_Check(v) (((PyObject*)v)->ob_type == &Pympc_Type)

/*
 * Define macros for comparing with zero, checking if either component is
 * 'nan' or 'inf', etc.
 */

#define MPC_IS_ZERO_P(x) \
    (mpfr_zero_p(mpc_realref(Pympc_AS_MPC(x))) && \
     mpfr_zero_p(mpc_imagref(Pympc_AS_MPC(x))))

#define MPC_IS_NAN_P(x) \
    (mpfr_nan_p(mpc_realref(Pympc_AS_MPC(x))) || \
     mpfr_nan_p(mpc_imagref(Pympc_AS_MPC(x))))

#define MPC_IS_INF_P(x) \
    (mpfr_inf_p(mpc_realref(Pympc_AS_MPC(x))) || \
     mpfr_inf_p(mpc_imagref(Pympc_AS_MPC(x))))

#define MPC_IS_FINITE_P(x) \
    (mpfr_number_p(mpc_realref(Pympc_AS_MPC(x))) && \
     mpfr_number_p(mpc_imagref(Pympc_AS_MPC(x))))

/* Verify that an object is an mpc and that both components have valid exp */
#define Pympc_CheckAndExp(v) \
    (Pympc_Check(v) && \
        (mpfr_zero_p(mpc_realref(Pympc_AS_MPC(v))) || \
            (mpfr_regular_p(mpc_realref(Pympc_AS_MPC(v))) && \
                ((mpc_realref(Pympc_AS_MPC(v))->_mpfr_exp >= context->ctx.emin)) && \
                ((mpc_realref(Pympc_AS_MPC(v))->_mpfr_exp <= context->ctx.emax)) \
            ) \
        ) && \
        (mpfr_zero_p(mpc_imagref(Pympc_AS_MPC(v))) || \
            (mpfr_regular_p(mpc_imagref(Pympc_AS_MPC(v))) && \
                ((mpc_imagref(Pympc_AS_MPC(v))->_mpfr_exp >= context->ctx.emin)) && \
                ((mpc_imagref(Pympc_AS_MPC(v))->_mpfr_exp <= context->ctx.emax)) \
            ) \
        ) \
    )

#define MPC_CHECK_UNDERFLOW(mpct, msg) \
    if (MPC_IS_ZERO_P(mpct) && mpct->rc) { \
        context->ctx.underflow = 1; \
        if (context->ctx.trap_underflow) { \
            GMPY_UNDERFLOW(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_OVERFLOW(mpct, msg) \
    if (MPC_IS_INF_P(mpct)) { \
        context->ctx.overflow = 1; \
        if (context->ctx.trap_overflow) { \
            GMPY_OVERFLOW(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_INVALID(mpct, msg) \
    if (MPC_IS_NAN_P(mpct)) { \
        context->ctx.invalid = 1; \
        if (context->ctx.trap_invalid) { \
            GMPY_INVALID(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_INEXACT(mpct, msg) \
    if (mpct->rc) { \
        context->ctx.inexact = 1; \
        if (context->ctx.trap_inexact) { \
            GMPY_INEXACT(msg); \
            goto done; \
        } \
    }

#define MPC_CHECK_FLAGS(mpct, NAME) \
    MPC_CHECK_INVALID(mpct, "'mpc' invalid operation in "NAME); \
    MPC_CHECK_UNDERFLOW(mpct, "'mpc' underflow in "NAME); \
    MPC_CHECK_OVERFLOW(mpct, "'mpc' overflow in "NAME); \
    MPC_CHECK_INEXACT(mpct, "'mpc' inexact result in "NAME);

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

/*
 * Parses one, and only one, argument into "self" and converts it to an
 * mpc. Is faster, but not as generic, as using PyArg_ParseTuple. It
 * supports either gmpy.fname(z) or z.fname(). "self" must be decref'ed.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments.
 */

#define PARSE_ONE_MPC_ARGS(msg) \
    if(self && Pympc_Check(self)) { \
        if (PyTuple_GET_SIZE(args) != 0) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        if (Pympc_CheckAndExp(self)) { \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)Pympc_From_Complex(self, 0, 0))) { \
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
        if (Pympc_CheckAndExp(self)) { \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)Pympc_From_Complex(self, 0, 0))) { \
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
    if(self && Pympc_Check(self)) { \
        if (Pympc_CheckAndExp(self)) { \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)Pympc_From_Complex(self, 0, 0))) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
        } \
    } \
    else { \
        if (Pympc_CheckAndExp(other)) { \
            self = other; \
            Py_INCREF(self); \
        } \
        else { \
            if (!(self = (PyObject*)Pympc_From_Complex(other, 0, 0))) { \
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
    if(self && Pympc_Check(self)) { \
        if (PyTuple_GET_SIZE(args) != 1) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        self = (PyObject*)Pympc_From_Complex(self, 0, 0); \
        var = (PyObject*)Pympc_From_Complex(PyTuple_GET_ITEM(args, 0), 0, 0); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) != 2) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        self = (PyObject*)Pympc_From_Complex(PyTuple_GET_ITEM(args, 0), 0, 0); \
        var = (PyObject*)Pympc_From_Complex(PyTuple_GET_ITEM(args, 1), 0, 0); \
    } \
    if (!self || !var) { \
        TYPE_ERROR(msg); \
        Py_XDECREF((PyObject*)var); \
        Py_XDECREF((PyObject*)self); \
        return NULL; \
    }

static PyObject * Pympc_digits(PyObject *self, PyObject *args);
static PyObject * Pygmpy_mpc(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject * Pympc_format(PyObject *self, PyObject *args);
static PyObject * Pympc_abs(PyObject *self);
static PyObject * Pympc_neg(PympcObject *self);
static PyObject * Pympc_pos(PympcObject *self);
static PyObject * Pympc_sqr(PyObject* self, PyObject *other);
static PyObject * Pympc_pow(PyObject *base, PyObject *exp, PyObject *m);
static PyObject * Pympc_conjugate(PyObject *self, PyObject *args);
static PyObject * Pympc_getprec_attrib(PympcObject *self, void *closure);
static PyObject * Pympc_getrc_attrib(PympcObject *self, void *closure);
static PyObject * Pympc_getimag_attrib(PympcObject *self, void *closure);
static PyObject * Pympc_getreal_attrib(PympcObject *self, void *closure);
static int Pympc_nonzero(PympcObject *self);
static PyObject * Pympc_is_NAN(PyObject *self, PyObject *other);
static PyObject * Pympc_is_ZERO(PyObject *self, PyObject *other);
static PyObject * Pympc_is_INF(PyObject *self, PyObject *other);
static PyObject * Pympc_is_FINITE(PyObject *self, PyObject *other);
static PyObject * Pympc_phase(PyObject *self, PyObject *other);
static PyObject * Pympc_norm(PyObject *self, PyObject *other);
static PyObject * Pympc_polar(PyObject *self, PyObject *other);
static PyObject * Pympc_rect(PyObject *self, PyObject *args);
static PyObject * Pympc_proj(PyObject *self, PyObject *other);
static PyObject * Pympc_log(PyObject *self, PyObject *other);
static PyObject * Pympc_log10(PyObject *self, PyObject *other);
static PyObject * Pympc_exp(PyObject *self, PyObject *other);
static PyObject * Pympc_sin(PyObject *self, PyObject *other);
static PyObject * Pympc_cos(PyObject *self, PyObject *other);
static PyObject * Pympc_tan(PyObject *self, PyObject *other);
static PyObject * Pympc_sinh(PyObject *self, PyObject *other);
static PyObject * Pympc_cosh(PyObject *self, PyObject *other);
static PyObject * Pympc_tanh(PyObject *self, PyObject *other);
static PyObject * Pympc_asin(PyObject *self, PyObject *other);
static PyObject * Pympc_acos(PyObject *self, PyObject *other);
static PyObject * Pympc_atan(PyObject *self, PyObject *other);
static PyObject * Pympc_asinh(PyObject *self, PyObject *other);
static PyObject * Pympc_acosh(PyObject *self, PyObject *other);
static PyObject * Pympc_atanh(PyObject *self, PyObject *other);
static PyObject * Pympc_sqrt(PyObject *self, PyObject *other);
static PyObject * Pympc_sin_cos(PyObject *self, PyObject *other);
static PyObject * Pympc_fma(PyObject *self, PyObject *args);
static PyObject * Pympc_fms(PyObject *self, PyObject *args);
static PyObject * Pympc_div_2exp(PyObject *self, PyObject *args);
static PyObject * Pympc_mul_2exp(PyObject *self, PyObject *args);
static Py_hash_t Pympc_hash(PympcObject *self);
static PyObject * Pympc_add_fast(PyObject *x, PyObject *y);
static PyObject * Pympc_add(PyObject *self, PyObject *args);
static PyObject * Pympc_sub_fast(PyObject *x, PyObject *y);
static PyObject * Pympc_sub(PyObject *self, PyObject *args);
static PyObject * Pympc_mul_fast(PyObject *x, PyObject *y);
static PyObject * Pympc_mul(PyObject *self, PyObject *args);
static PyObject * Pympc_truediv_fast(PyObject *x, PyObject *y);
#ifdef PY2
static PyObject * Pympc_div2_fast(PyObject *x, PyObject *y);
#endif
static PyObject * Pympc_div(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
