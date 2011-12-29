/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy.h                                                                  *
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

typedef struct {
    mpfr_prec_t mpfr_prec;   /* current precision in bits, for MPFR */
#ifdef WITHMPC
    mpfr_prec_t mpc_rprec;   /* current precision in bits, for Re(MPC) */
    mpfr_prec_t mpc_iprec;   /* current precision in bits, for Im(MPC) */
#endif
    mpfr_rnd_t mpfr_round;   /* current rounding mode for float (MPFR) */
#ifdef WITHMPC
    mpfr_rnd_t mpc_rround;   /* current rounding mode for Re(MPC) */
    mpfr_rnd_t mpc_iround;   /* current rounding mode for Im(MPC) */
#endif
    mpfr_exp_t emax;         /* maximum exponent */
    mpfr_exp_t emin;         /* minimum exponent */
    int subnormalize;        /* if 1, subnormalization is performed */
    int underflow;           /* did an underflow occur? */
    int overflow;            /* did an overflow occur? */
    int inexact;             /* was the result inexact? */
    int invalid;             /* invalid operation (i.e. NaN)? */
    int erange;              /* did a range error occur? */
    int divzero;             /* divided by zero? */
    int trap_underflow;      /* if 1, raise exception for underflow */
    int trap_overflow;       /* if 1, raise exception for overflow */
    int trap_inexact;        /* if 1, raise exception for inexact */
    int trap_invalid;        /* if 1, raise exception for invalid (NaN) */
    int trap_erange;         /* if 1, raise exception for range error */
    int trap_divzero;        /* if 1, raise exception for divide by zero */
    int trap_expbound;       /* if 1, raise exception if mpfr/mpc exponents */
                             /*       are out of bounds */
#ifdef WITHMPC
    int allow_complex;       /* if 1, allow mpfr functions to return an mpc */
#endif
} gmpy_context;

typedef struct {
    PyObject_HEAD
    gmpy_context now;        /* The "new" values, used by __enter__ */
    PyObject *orig;          /* Original context, restored by __exit__*/
} GMPyContextObject;

static PyTypeObject GMPyContext_Type;
#define GMPyContext_Check(v) (((PyObject*)v)->ob_type == &GMPyContext_Type)
