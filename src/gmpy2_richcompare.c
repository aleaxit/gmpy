/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_richcompare.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
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
GMPy_RichCompare_Slot(PyObject *a, PyObject *b, int op)
{
    int c;
    PyObject *tempa = NULL, *tempb = NULL, *result = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (CHECK_MPZANY(a)) {
        if (PyIntOrLong_Check(b)) {
            int error;
            long temp = PyLong_AsLongAndOverflow(b, &error);

            if (!error) {
                c = mpz_cmp_si(MPZ(a), temp);
            }
            else {
                mpz_set_PyIntOrLong(global.tempz, b);
                c = mpz_cmp(MPZ(a), global.tempz);
            }
            return _cmp_to_object(c, op);
        }

        if (CHECK_MPZANY(b)) {
            return _cmp_to_object(mpz_cmp(MPZ(a), MPZ(b)), op);
        }

        if (IS_INTEGER(b)) {
            if (!(tempb = (PyObject*)GMPy_MPZ_From_Integer(b, context))) {
                return NULL;
            }
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
    }

    if (MPQ_Check(a)) {
        if (MPQ_Check(b)) {
            return _cmp_to_object(mpq_cmp(MPQ(a), MPQ(b)), op);
        }

        if (IS_RATIONAL(b)) {
            if (!(tempb = (PyObject*)GMPy_MPQ_From_Rational(b, context))) {
                return NULL;
            }
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
                if (!(tempb = (PyObject*)GMPy_MPQ_New(context))) {
                    return NULL;
                }
                mpq_set_d(MPQ(tempb), d);
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
            if (!(tempb = (PyObject*)GMPy_MPZ_From_Integer(b, context)))  {
                return NULL;
            }
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
            if (!(tempb = (PyObject*)GMPy_MPQ_From_Rational(b, context))) {
                return NULL;
            }
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

        if (IS_REAL(b)) {
            if (!(tempb = (PyObject*)GMPy_MPFR_From_Real(b, 1, context))) {
                return NULL;
            }
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
            if (!(tempb = (PyObject*)GMPy_MPC_From_PyComplex(b, 1, 1, context))) {
                return NULL;
            }
            mpfr_clear_flags();
            c = mpc_cmp(MPC(a), MPC(tempb));
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
            if (!(tempb = (PyObject*)GMPy_MPFR_New(mpfr_get_prec(mpc_realref(MPC(a))), context))) {
                return NULL;
            }
            mpc_real(MPFR(tempb), MPC(a), GET_MPFR_ROUND(context));
            result = GMPy_RichCompare_Slot(tempb, b, op);
            Py_DECREF(tempb);
            return result;
        }
    }

    Py_RETURN_NOTIMPLEMENTED;
}

