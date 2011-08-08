/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.h                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008 - 2011 Case Van Horsen                              *
 *                                                                         *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2.1 of the License, or  *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This library is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * Lesser General Public License for more details.                         *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with this library; if not, write to the Free Software     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA           *
 * 02110-1301  USA                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


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
    mpob ob;
    mpc_t c;
    Py_hash_t hash_cache;
    int rc;
    int round_mode;
} PympcObject;

static PyTypeObject Pympc_Type;
#define Pympc_AS_MPC(obj) (((PympcObject *)(obj))->c)
#define Pympc_Check(v) (((PyObject*)v)->ob_type == &Pympc_Type)
/* Verify that an object is an mpc and that both components have valid exp */
#define Pympc_CheckAndExp(v) \
    (Pympc_Check(v) && \
        (mpfr_zero_p(mpc_realref(Pympc_AS_MPC(v))) || \
            (mpfr_regular_p(mpc_realref(Pympc_AS_MPC(v))) && \
                ((mpc_realref(Pympc_AS_MPC(v))->_mpfr_exp >= context->now.emin)) && \
                ((mpc_realref(Pympc_AS_MPC(v))->_mpfr_exp <= context->now.emax)) \
            ) \
        ) && \
        (mpfr_zero_p(mpc_imagref(Pympc_AS_MPC(v))) || \
            (mpfr_regular_p(mpc_imagref(Pympc_AS_MPC(v))) && \
                ((mpc_imagref(Pympc_AS_MPC(v))->_mpfr_exp >= context->now.emin)) && \
                ((mpc_imagref(Pympc_AS_MPC(v))->_mpfr_exp <= context->now.emax)) \
            ) \
        ) \
    )

#define MPC_CHECK_OVERFLOW(mpct, msg) \
    if ((mpfr_inf_p(mpc_realref(mpct->c)) || mpfr_inf_p(mpc_imagref(mpct->c))) \
        && context->now.trap_overflow) { \
        GMPY_OVERFLOW(msg); \
        goto done; \
    }

#define MPC_CHECK_INVALID(mpct, msg) \
    if ((mpfr_nan_p(mpc_realref(mpct->c)) || mpfr_nan_p(mpc_imagref(mpct->c))) \
        && context->now.trap_invalid) { \
        GMPY_INVALID(msg); \
        goto done; \
    }

#define MPC_CHECK_FLAGS(mpct, NAME) \
    MPC_CHECK_INVALID(mpct, "invalid operation in 'mpc' "NAME); \
    MPC_CHECK_OVERFLOW("overflow in 'mpc' "NAME); \

#define MPC_SUBNORMALIZE(NAME) \
    if (context->now.subnormalize) { \
        int rcr, rci; \
        rcr = MPC_INEX_RE(NAME->rc); \
        rci = MPC_INEX_IM(NAME->rc); \
        rcr = mpfr_subnormalize(mpc_realref(NAME->c), rcr, GET_MPC_RROUND(context)); \
        rci = mpfr_subnormalize(mpc_imagref(NAME->c), rci, GET_MPC_IROUND(context)); \
        NAME->rc = MPC_INEX(rcr, rci); \
    } \

#define MPC_CLEANUP_SELF(NAME) \
    MPC_SUBNORMALIZE(result); \
    MPC_CHECK_INVALID("invalid operation in 'mpc' "NAME); \
    MPC_CHECK_OVERFLOW("overflow in 'mpc' "NAME); \
  done: \
    Py_DECREF(self); \
    if (PyErr_Occurred()) { \
        Py_XDECREF((PyObject*)result); \
        result = NULL; \
    } \
    return (PyObject*)result;

#define MPC_CLEANUP_SELF_OTHER(NAME) \
    MPC_SUBNORMALIZE(result); \
    MPC_CHECK_INVALID("invalid operation in 'mpc' "NAME); \
    MPC_CHECK_OVERFLOW("overflow in 'mpc' "NAME); \
  done: \
    Py_DECREF(self); \
    Py_DECREF(other); \
    if (PyErr_Occurred()) { \
        Py_XDECREF((PyObject*)result); \
        result = NULL; \
    } \
    return (PyObject*)result;

#define MPC_CLEANUP_RC(NAME) \
    MPC_SUBNORMALIZE(rc); \
    if ((mpfr_inf_p(mpc_realref(rc->c)) || mpfr_inf_p(mpc_imagref(rc->c))) \
        && context->now.trap_overflow) { \
        GMPY_OVERFLOW("overflow in 'mpc' " #NAME); \
        return NULL; \
    } \
    return (PyObject*)rc;

#define MPC_CLEANUP_RESULT(NAME) \
    MPC_SUBNORMALIZE(result); \
    if ((mpfr_nan_p(mpc_realref(result->c)) || \
         mpfr_nan_p(mpc_imagref(result->c))) && \
        context->now.trap_invalid) { \
        GMPY_INVALID("invalid operation in 'mpc' "NAME); \
        Py_DECREF((PyObject*)result); \
        result = NULL; \
        goto done; \
    } \
    if ((mpfr_inf_p(mpc_realref(result->c)) || \
         mpfr_inf_p(mpc_imagref(result->c))) && \
        context->now.trap_overflow) { \
        GMPY_OVERFLOW("overflow in 'mpc' "NAME); \
        Py_DECREF((PyObject*)result); \
        result = NULL; \
        goto done; \
    }
