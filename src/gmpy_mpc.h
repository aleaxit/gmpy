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
