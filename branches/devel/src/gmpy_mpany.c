/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpany.c                                                            *
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


/* Generic module-level methods for gmpy types.
 *
 * These methods are designed to accept any number type as input and call
 * the appropriate type-specific method. For example, gmpy2.digits(n) will
 * call gmpy2.mpz(n).digits() if n is an integer, gmpy2.mpq(n).digits() if
 * n is a rational, or gmpy2.mpf(n).digits() is n is a float.
 */

/* gmpy_sign is only intended to be used at the module level!
 * gmpy_sign uses the METH_O/METH_NOARGS calling convention!
 * gmpy_sign assumes mpX_sign also use the METH_O/METH_NOARGS convention!
 */

PyDoc_STRVAR(doc_g_mpany_sign,
"sign(x) -> number\n\n"
"Return -1 if x < 0, 0 if x == 0, or +1 if x >0.");

static PyObject *
Pympany_sign(PyObject *self, PyObject *other)
{
    if (IS_INTEGER(other))
        return Pympz_sign(self, other);
    else if (IS_RATIONAL(other))
        return Pympq_sign(self, other);
    else if (IS_REAL(other))
        return Pympfr_sign(self, other);

    TYPE_ERROR("sign() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_printf,
"_printf(fmt, x) -> string\n\n"
"Return a Python string by formatting 'x' using the format string\n"
"'fmt'.\n\n"
"WARNING: Invalid format strings will cause a crash. Please see the\n"
"         GMP and MPFR manuals for details on the format code. 'mpc'\n"
"         objects are not supported.");

static PyObject *
Pympany_printf(PyObject *self, PyObject *args)
{
    PyObject *result = 0, *x = 0;
    char *buffer = 0, *fmtcode = 0;
    void *generic;
    int buflen;

    if (!PyArg_ParseTuple(args, "sO", &fmtcode, &x))
        return NULL;

    if (CHECK_MPZANY(x) || MPQ_Check(x)) {
        if (CHECK_MPZANY(x))
            generic = MPZ(x);
        else
            generic = MPQ(x);
        buflen = gmp_asprintf(&buffer, fmtcode, generic);
        if (buflen < 0) {
            VALUE_ERROR("printf() could not format the 'mpz' or 'mpq' object");
        }
        else {
            result = Py_BuildValue("s", buffer);
            GMPY_FREE(buffer);
        }
        return result;
    }
    else if(MPFR_Check(x)) {
        generic = MPFR(x);
        buflen = mpfr_asprintf(&buffer, fmtcode, generic);
        if (buflen < 0) {
            VALUE_ERROR("printf() could not format the 'mpfr' object");
        }
        else {
            result = Py_BuildValue("s", buffer);
            GMPY_FREE(buffer);
        }
        return result;
    }
    else if(MPC_Check(x)) {
        TYPE_ERROR("printf() does not support 'mpc'");
        return NULL;
    }
    else {
        TYPE_ERROR("printf() argument type not supported");
        return NULL;
    }
}

#define MPANY_MPFR_MPC(NAME) \
static PyObject * \
Pympany_##NAME(PyObject *self, PyObject *other) \
{ \
    if (IS_REAL(other)) \
        return Pympfr_##NAME(self, other); \
    else if (IS_COMPLEX(other)) \
        return Pympc_##NAME(self, other); \
    TYPE_ERROR(#NAME"() argument type not supported"); \
    return NULL; \
}

PyDoc_STRVAR(doc_mpany_is_nan,
"is_nan(x) -> boolean\n\n"
"Return True if x is NaN (Not-A-Number).");

static PyObject *
Pympany_is_nan(PyObject *self, PyObject *other)
{
    if (IS_REAL(other))
        return Pympfr_is_nan(self, other);
    else if (IS_COMPLEX(other))
        return Pympc_is_NAN(self, other);
    TYPE_ERROR("is_nan() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_is_inf,
"is_inf(x) -> boolean\n\n"
"Return True if x is +Infinity or -Infinity.\n"
"Note: is_inf() is deprecated; please use is_infinite().");

static PyObject *
Pympany_is_inf(PyObject *self, PyObject *other)
{
    if (IS_REAL(other))
        return Pympfr_is_inf(self, other);
    else if (IS_COMPLEX(other))
        return Pympc_is_INF(self, other);
    TYPE_ERROR("is_inf() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_is_infinite,
"is_infinite(x) -> boolean\n\n"
"Return True if x is +Infinity or -Infinity. If x is an mpc, return True\n"
"if either x.real or x.imag is infinite.");

static PyObject *
Pympany_is_infinite(PyObject *self, PyObject *other)
{
    if (IS_REAL(other))
        return Pympfr_is_inf(self, other);
    else if (IS_COMPLEX(other))
        return Pympc_is_INF(self, other);
    TYPE_ERROR("is_infinite() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_is_finite,
"is_finite(x) -> boolean\n\n"
"Return True if x is an actual number (i.e. not NaN or Infinity). If x is\n"
"an mpc, return True if both x.real and x.imag are finite.");

static PyObject *
Pympany_is_finite(PyObject *self, PyObject *other)
{
    if (IS_REAL(other))
        return Pympfr_is_number(self, other);
    else if (IS_COMPLEX(other))
        return Pympc_is_FINITE(self, other);
    TYPE_ERROR("is_finite() argument type not supported");
    return NULL;
}

PyDoc_STRVAR(doc_mpany_is_zero,
"is_zero(x) -> boolean\n\n"
"Return True if x is zero.");

static PyObject *
Pympany_is_zero(PyObject *self, PyObject *other)
{
    if (IS_REAL(other))
        return Pympfr_is_zero(self, other);
    else if (IS_COMPLEX(other))
        return Pympc_is_ZERO(self, other);
    TYPE_ERROR("is_zero() argument type not supported");
    return NULL;
}

/* COMPARING */

static PyObject *_cmp_to_object(int c, int op)
{
    PyObject *result;
    switch (op) {
    case Py_LT: c = c <  0; break;
    case Py_LE: c = c <= 0; break;
    case Py_EQ: c = c == 0; break;
    case Py_NE: c = c != 0; break;
    case Py_GT: c = c >  0; break;
    case Py_GE: c = c >= 0; break;
    }
    result = c ? Py_True : Py_False;
    Py_INCREF(result);
    return result;
}
static PyObject *
mpany_richcompare(PyObject *a, PyObject *b, int op)
{
    int c, overflow;
    mpir_si temp_si;
    mpz_t tempz;
    PyObject *tempa = 0, *tempb = 0;
    PyObject *result = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            temp_si = PyLong_AsSIAndOverflow(b, &overflow);
            if (overflow) {
                mpz_inoc(tempz);
                mpz_set_PyIntOrLong(tempz, b);
                c = mpz_cmp(MPZ(a), tempz);
                mpz_cloc(tempz);
            }
            else {
                c = mpz_cmp_si(MPZ(a), temp_si);
            }
            return _cmp_to_object(c, op);
        }
        if (CHECK_MPZANY(b)) {
            return _cmp_to_object(mpz_cmp(MPZ(a), MPZ(b)), op);
        }
        if (IS_INTEGER(b)) {
            tempb = (PyObject*)GMPy_MPZ_From_Integer(b, context);
            if (!tempb)
                return NULL;
            c = mpz_cmp(MPZ(a), MPZ(tempb));
            Py_DECREF(tempb);
            return _cmp_to_object(c, op);
        }
        if (IS_RATIONAL(b)) {
            tempa = (PyObject*)GMPy_MPQ_From_Rational(a, context);
            tempb = (PyObject*)GMPy_MPQ_From_Rational(b, context);
            if (!tempa || !tempb) {
                Py_XDECREF(a);
                Py_XDECREF(b);
                return NULL;
            }
            c = mpq_cmp(MPQ(tempa), MPQ(tempb));
            Py_DECREF(tempa);
            Py_DECREF(tempb);
            return _cmp_to_object(c, op);
        }
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            if (Py_IS_NAN(d)) {
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else if (Py_IS_INFINITY(d)) {
                if (d < 0.0)
                    return _cmp_to_object(1, op);
                else
                    return _cmp_to_object(-1, op);
            }
            else {
                return _cmp_to_object(mpz_cmp_d(MPZ(a), d), op);
            }
        }
        if (IS_DECIMAL(b)) {
            tempa = (PyObject*)GMPy_MPQ_From_Rational(a, context);
            tempb = (PyObject*)GMPy_MPQ_From_Decimal(b, context);
            if (!tempa || !tempb) {
                Py_XDECREF(a);
                Py_XDECREF(b);
                return NULL;
            }
            if (!mpz_cmp_si(mpq_denref(MPQ(tempb)), 0)) {
                if (!mpz_cmp_si(mpq_numref(MPQ(tempb)), 0)) {
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_DECREF(tempa);
                    Py_DECREF(tempb);
                    Py_INCREF(result);
                    return result;
                }
                else if (mpz_cmp_si(mpq_numref(MPQ(tempb)), 0) < 0) {
                    Py_DECREF(tempa);
                    Py_DECREF(tempb);
                    return _cmp_to_object(1, op);
                }
                else {
                    Py_DECREF(tempa);
                    Py_DECREF(tempb);
                    return _cmp_to_object(-1, op);
                }
            }
            else {
                c = mpq_cmp(MPQ(tempa), MPQ(tempb));
                Py_DECREF(tempa);
                Py_DECREF(tempb);
                return _cmp_to_object(c, op);
            }
        }
    }
    if (MPQ_Check(a)) {
        if (MPQ_Check(b)) {
            return _cmp_to_object(mpq_cmp(MPQ(a), MPQ(b)), op);
        }
        if (IS_RATIONAL(b)) {
            tempb = (PyObject*)GMPy_MPQ_From_Rational(b, context);
            c = mpq_cmp(MPQ(a), MPQ(tempb));
            Py_DECREF(tempb);
            return _cmp_to_object(c, op);
        }
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            if (Py_IS_NAN(d)) {
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else if (Py_IS_INFINITY(d)) {
                if (d < 0.0)
                    return _cmp_to_object(1, op);
                else
                    return _cmp_to_object(-1, op);
            }
            else {
                tempb = (PyObject*)GMPy_MPQ_New(context);
                if (!tempb)
                    return NULL;
                mpq_set_d(MPQ(tempb), d);
                c = mpq_cmp(MPQ(a), MPQ(tempb));
                Py_DECREF(tempb);
                return _cmp_to_object(c, op);
            }
        }
        if (IS_DECIMAL(b)) {
            if (!(tempb = (PyObject*)GMPy_MPQ_From_Decimal(b, context)))
                return NULL;
            if (!mpz_cmp_si(mpq_denref(MPQ(tempb)), 0)) {
                if (!mpz_cmp_si(mpq_numref(MPQ(tempb)), 0)) {
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_DECREF(tempb);
                    Py_INCREF(result);
                    return result;
                }
                else if (mpz_cmp_si(mpq_numref(MPQ(tempb)), 0) < 0) {
                    Py_DECREF(tempb);
                    return _cmp_to_object(1, op);
                }
                else {
                    Py_DECREF(tempb);
                    return _cmp_to_object(-1, op);
                }
            }
            else {
                c = mpq_cmp(MPQ(a), MPQ(tempb));
                Py_DECREF(tempb);
                return _cmp_to_object(c, op);
            }
        }
    }

    if (MPFR_Check(a)) {
        if (MPFR_Check(b)) {
            mpfr_clear_flags();
            c = mpfr_cmp(MPFR(a), MPFR(b));
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (PyFloat_Check(b)) {
            double d = PyFloat_AS_DOUBLE(b);
            mpfr_clear_flags();
            c = mpfr_cmp_d(MPFR(a), d);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (IS_INTEGER(b)) {
            tempb = (PyObject*)GMPy_MPZ_From_Integer(b, context);
            if (!tempb)
                return NULL;
            mpfr_clear_flags();
            c = mpfr_cmp_z(MPFR(a), MPZ(tempb));
            Py_DECREF(tempb);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (IS_RATIONAL(b)) {
            tempb = (PyObject*)GMPy_MPQ_From_Rational(b, context);
            if (!tempb)
                return NULL;
            mpfr_clear_flags();
            c = mpfr_cmp_q(MPFR(a), MPQ(tempb));
            Py_DECREF(tempb);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (IS_DECIMAL(b)) {
            tempb = (PyObject*)GMPy_MPQ_From_Decimal(b, context);
            if (!tempb)
                return NULL;
            if (!mpz_cmp_si(mpq_denref(MPQ(tempb)), 0)) {
                if (!mpz_cmp_si(mpq_numref(MPQ(tempb)), 0)) {
                    context->ctx.erange = 1;
                    if (context->ctx.traps & TRAP_ERANGE) {
                        GMPY_ERANGE("comparison with NaN");
                        return NULL;
                    }
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_DECREF(tempb);
                    Py_INCREF(result);
                    return result;
                }
                else if (mpz_cmp_si(mpq_numref(MPQ(tempb)), 0) < 0) {
                    Py_DECREF(tempb);
                    return _cmp_to_object(1, op);
                }
                else {
                    Py_DECREF(tempb);
                    return _cmp_to_object(-1, op);
                }
            }
            else {
                mpfr_clear_flags();
                c = mpfr_cmp_q(MPFR(a), MPQ(tempb));
                Py_DECREF(tempb);
                if (mpfr_erangeflag_p()) {
                    /* Set erange and check if an exception should be raised. */
                    context->ctx.erange = 1;
                    if (context->ctx.traps & TRAP_ERANGE) {
                        GMPY_ERANGE("comparison with NaN");
                        return NULL;
                    }
                    result = (op == Py_NE) ? Py_True : Py_False;
                    Py_INCREF(result);
                    return result;
                }
                else {
                    return _cmp_to_object(c, op);
                }
            }
        }
        if (IS_REAL(b)) {
            tempb = (PyObject*)GMPy_MPFR_From_Real(b, 1, context);
            if (!tempb)
                return NULL;
            mpfr_clear_flags();
            c = mpfr_cmp(MPFR(a), MPFR(tempb));
            Py_DECREF(tempb);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
    }

    if (MPC_Check(a)) {
        if (!(op == Py_EQ || op == Py_NE)) {
            TYPE_ERROR("no ordering relation is defined for complex numbers");
            return NULL;
        }
        if (MPC_Check(b)) {
            mpfr_clear_flags();
            c = mpc_cmp(MPC(a), MPC(b));
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        if (PyComplex_Check(b)) {
            MPC_Object *tempmpc;

            if (!(tempmpc = GMPy_MPC_From_PyComplex(b, 53, 53, context)))
                return NULL;
            mpfr_clear_flags();
            c = mpc_cmp(MPC(a), MPC(tempmpc));
            Py_DECREF((PyObject*)tempmpc);
            if (mpfr_erangeflag_p()) {
                /* Set erange and check if an exception should be raised. */
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
                result = (op == Py_NE) ? Py_True : Py_False;
                Py_INCREF(result);
                return result;
            }
            else {
                return _cmp_to_object(c, op);
            }
        }
        /* a.imag must be 0 or else all further comparisons will be NE */
        if (!mpfr_zero_p(mpc_imagref(MPC(a)))) {
            /* if a.real is NaN, possibly raise exception */
            if (mpfr_nan_p(mpc_realref(MPC(a)))) {
                context->ctx.erange = 1;
                if (context->ctx.traps & TRAP_ERANGE) {
                    GMPY_ERANGE("comparison with NaN");
                    return NULL;
                }
            }
            result = (op == Py_NE) ? Py_True : Py_False;
            Py_INCREF(result);
            return result;
        }
        else {
            MPFR_Object *tempmpfr;

            tempmpfr = GMPy_MPFR_New(mpfr_get_prec(mpc_realref(MPC(a))), context);
            if (!tempmpfr)
                return NULL;
            mpc_real(tempmpfr->f, MPC(a), GET_MPFR_ROUND(context));
            result = mpany_richcompare((PyObject*)tempmpfr, b, op);
            Py_DECREF((PyObject*)tempmpfr);
            return result;
        }
    }

    Py_RETURN_NOTIMPLEMENTED;
}

