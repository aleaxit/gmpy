/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert_mpfr.c                                                     *
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

/* This file contains all the conversion functions for MPFR data types.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. The basic operations (+, -, *, /) are optimized to
 * directly work with the basic types such as C longs or doubles.
 *
 * gmpy2 supports two different strategies for creating new references to
 * an mpfr instance. If bits (or prec) is set to 0, the precison of the
 * result exactly matches the precision of the context. This is the default
 * behavior of the mpfr() function. If bits (or prec) is set to 1, the
 * precision of the result depends on the type of the source.
 *
 *   1) if the source number is already a radix-2 floating point number,
 *      the precision is not changed. In practical terms, this only applies
 *      to sources operands that either an mpfr or Python double.
 *
 *   2) otherwise, the source is converted to an mpfr using a precision of
 *      (context.precison + context.guard_bits). guard_bits is
 *      set to 0 by default but can be changed to provide additional
 *      precision for creating temporary mpfr instances that are used
 *      are arguments for MPFR functions.
 *
 * Even though the precision of an mpfr or PyDouble can remain unchanged,
 * the exponent of the result is always guaranteed to be valid.
 *
 * Support for the Decimal type is a challenge. For the basic operations, it
 * is most accurate to convert a Decimal instance into an mpq and then use
 * MPFR's functions to accurately operate on an mpfr and mpq. This approach is
 * challenging because (1) a large exponent can create a very large mpq and
 * (2) the changes made to C-coded version of Decimal in Python 3.3. See
 * gmpy2_convert_gmp.c for the code to convert a Decimal to an mpq exactly.
 */

/* If prec == 0:
 *   Return an additional reference to an existing mpfr if the exponent is
 *   valid in the current context and the precision of the existing mpfr
 *   equals the current context precision.
 *
 *   If the exponent is not valid or if the precision does not equal the
 *   current context precision, then a reference to a new, valid instance is
 *   returned. The precision of the current context will be used.
 *
 * If prec == 1:
 *   Return an additional reference to an existing mpfr if the exponent is
 *   valid in the current context. The precision is not changed.
 *
 *   If the exponent is not valid, then a reference to a new, valid instance
 *   is returned. The precision is not changed.
 *
 * If prec >= 2:
 *   Return a reference to a valid instance with the specified precision.
 *
 * All mpfr arguments to functions in the MPFR library should go through this
 * function to guarantee that the exponent is valid. prec=1 should be used for
 * validating all arguments. To maintain maximum accuracy, the precision is
 * not changed. References returned when prec=1 should not be returned to the
 * user (although they must still be decremented like normal references).
 *
 * prec=0 should be used for creating references that may be returned to the
 * user. Note: do not mutate the returned reference; they are not guaranteed
 * to be new objects.
 */

static MPFR_Object *
GMPy_MPFR_From_MPFR(MPFR_Object *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL;

    assert(MPFR_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0)
        prec = GET_MPFR_PREC(context);
    else if (prec == 1)
        prec = mpfr_get_prec(obj->f);

    if (MPFR_CheckAndExp(obj) && !mpfr_nan_p(MPFR(obj))) {
        if (prec == mpfr_get_prec(obj->f)) {
            Py_INCREF((PyObject*)obj);
            return obj;
        }
        else {
            if ((result = GMPy_MPFR_New(prec, context))) {
                mpfr_clear_flags();
                result->rc = mpfr_set(result->f, obj->f, GET_MPFR_ROUND(context));
                MPFR_CLEANUP_2(result, context, "mpfr()");
            }
            return result;
        }
    }

    if (context->ctx.traps & TRAP_EXPBOUND) {
        GMPY_EXPBOUND("exponent of existing mpfr incompatible with current context");
        return NULL;
    }

    if ((result = GMPy_MPFR_New(mpfr_get_prec(obj->f), context))) {
        /* First make the exponent valid. */
        mpfr_set(result->f, obj->f, GET_MPFR_ROUND(context));
        mpfr_clear_flags();
        result->rc = mpfr_check_range(result->f, obj->rc, obj->round_mode);
        /* Then round to the desired precision. */
        result->rc = mpfr_prec_round(result->f, prec, GET_MPFR_ROUND(context));
        /* Removing any exception checks. Need to think about what should happen.
         * MPFR_CLEANUP_2(result, context, "mpfr()");
         */
    }

    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of the current context plus guard bits is
 * used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_PyIntOrLong(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPZ_Object *tempz;

    assert(PyIntOrLong_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0 || prec == 1)
        prec = GET_MPFR_PREC(context) + prec * GET_GUARD_BITS(context);

    if ((tempz = GMPy_MPZ_From_PyIntOrLong(obj, context))) {
        result = GMPy_MPFR_From_MPZ(tempz, prec, context);
        Py_DECREF((PyObject*)tempz);
    }
    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of Python float type is used (typically 53).
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_PyFloat(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;

    assert(PyFloat_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0)
        prec = GET_MPFR_PREC(context);
    else if (prec == 1)
        prec = DBL_MANT_DIG;

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set_d(result->f, PyFloat_AS_DOUBLE(obj),
                                GET_MPFR_ROUND(context));
        MPFR_CLEANUP_2(result, context, "mpfr()");
    }

    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of the current context plus guard bits is
 * used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_MPZ(MPZ_Object *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;

    assert(MPZ_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0 || prec == 1)
        prec = GET_MPFR_PREC(context) + prec * GET_GUARD_BITS(context);

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set_z(result->f, obj->z, GET_MPFR_ROUND(context));
        MPFR_CLEANUP_2(result, context, "mpfr()");
    }

    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of the current context plus guard bits is
 * used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_MPQ(MPQ_Object *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;

    assert(MPQ_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0 || prec == 1)
        prec = GET_MPFR_PREC(context) + prec * GET_GUARD_BITS(context);

    if ((result = GMPy_MPFR_New(prec, context))) {
        mpfr_clear_flags();
        result->rc = mpfr_set_q(result->f, obj->q, GET_MPFR_ROUND(context));
        MPFR_CLEANUP_2(result, context, "mpfr()");
    }

    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of the current context plus guard bits is
 * used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_Fraction(PyObject *obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result = NULL;
    MPQ_Object *tempq;

    assert(IS_RATIONAL(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0 || prec == 1)
        prec = GET_MPFR_PREC(context) + prec * GET_GUARD_BITS(context);

    if ((tempq = GMPy_MPQ_From_Fraction(obj, context))) {
        result = GMPy_MPFR_From_MPQ(tempq, prec, context);
        Py_DECREF((PyObject*)tempq);
    }

    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of the current context plus guard bits is
 * used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_Decimal(PyObject* obj, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;
    MPQ_Object *temp;

    assert(IS_DECIMAL(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0 || prec == 1)
        prec = GET_MPFR_PREC(context) + prec * GET_GUARD_BITS(context);

    result = GMPy_MPFR_New(prec, context);
    temp = GMPy_MPQ_From_DecimalRaw(obj, context);

    if (!temp || !result) {
        Py_XDECREF((PyObject*)temp);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    if (!mpz_cmp_si(mpq_numref(temp->q), 0)) {
        if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
            mpfr_set_nan(result->f);
        }
        else {
            mpfr_set_zero(result->f, mpz_sgn(mpq_denref(temp->q)));
        }
    }
    else if (!mpz_cmp_si(mpq_denref(temp->q), 0)) {
        if (mpz_cmp_si(mpq_numref(temp->q), 0) < 0) {
            mpfr_set_inf(result->f, -1);
        }
        else {
            mpfr_set_inf(result->f, 1);
        }
    }
    else {
        Py_DECREF((PyObject*)result);
        result = GMPy_MPFR_From_MPQ(temp, prec, context);
    }
    Py_DECREF((PyObject*)temp);
    return result;
}

/* If prec==0, then the precision of the current context is used.
 *
 * If prec==1, then the precision of the current context plus guard bits is
 * used.
 *
 * If prec>=2, then the specified precision is used.
 */

static MPFR_Object *
GMPy_MPFR_From_PyStr(PyObject *s, int base, mpfr_prec_t prec, CTXT_Object *context)
{
    MPFR_Object *result;
    char *cp, *endptr;
    Py_ssize_t len;
    PyObject *ascii_str = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (prec == 0 || prec == 1)
        prec = GET_MPFR_PREC(context) + prec * GET_GUARD_BITS(context);

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = PyBytes_AsString(s);
    }
    else if (PyUnicode_Check(s)) {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = PyBytes_AsString(ascii_str);
    }
    else {
        TYPE_ERROR("object is not string or Unicode");
        return NULL;
    }

    if (!(result = GMPy_MPFR_New(prec, context))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* delegate the rest to MPFR */
    mpfr_clear_flags();
    result->rc = mpfr_strtofr(result->f, cp, &endptr, base, GET_MPFR_ROUND(context));
    Py_XDECREF(ascii_str);

    MPFR_CLEANUP_2(result, context, "mpfr()");

    if (len != (Py_ssize_t)(endptr - cp)) {
        VALUE_ERROR("invalid digits");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    return result;
}

/* GMPy_MPFR_From_Real() converts a real number (see IS_REAL()) to an mpfr.
 *
 * If prec==0, then the result has the precision of the current context.
 *
 * If prec==1 and the value can be converted exactly (i.e. the input value is
 * a floating-point number using radix-2 representation), then the conversion
 * is done with the maximum possible precision. If the input value can't be
 * converted exactly, then the context precision plus guard bits is used.
 *
 * If prec >=2, then the specified precision is used.
 *
 * The return value is guaranteed to have a valid exponent.
 */

static MPFR_Object *
GMPy_MPFR_From_Real(PyObject *obj, mp_prec_t prec, CTXT_Object *context)
{
    CHECK_CONTEXT_SET_EXPONENT(context);

    if (MPFR_Check(obj))
        return GMPy_MPFR_From_MPFR((MPFR_Object*)obj, prec, context);

    if (PyFloat_Check(obj))
        return GMPy_MPFR_From_PyFloat(obj, prec, context);

    if (MPQ_Check(obj))
        return GMPy_MPFR_From_MPQ((MPQ_Object*)obj, prec, context);

    if (MPZ_Check(obj) || XMPZ_Check(obj))
        return GMPy_MPFR_From_MPZ((MPZ_Object*)obj, prec, context);

    if (PyIntOrLong_Check(obj))
        return GMPy_MPFR_From_PyIntOrLong(obj, prec, context);

    if (IS_DECIMAL(obj))
        return GMPy_MPFR_From_Decimal(obj, prec, context);

    if (IS_FRACTION(obj))
        return GMPy_MPFR_From_Fraction(obj, prec, context);

    TYPE_ERROR("object could not be converted to 'mpfr'");
    return NULL;
}

static MPZ_Object *
GMPy_MPZ_From_MPFR(MPFR_Object *obj, CTXT_Object *context)
{
    MPZ_Object *result;

    assert(MPFR_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPZ_New(context))) {
        if (mpfr_nan_p(MPFR(obj))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'mpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(MPFR(obj))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'mpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, MPFR(obj), GET_MPFR_ROUND(context));
    }

    return result;
}

static XMPZ_Object *
GMPy_XMPZ_From_MPFR(MPFR_Object *self, CTXT_Object *context)
{
    XMPZ_Object *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_XMPZ_New(context))) {
        if (mpfr_nan_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            VALUE_ERROR("'xmpz' does not support NaN");
            return NULL;
        }
        if (mpfr_inf_p(MPFR(self))) {
            Py_DECREF((PyObject*)result);
            OVERFLOW_ERROR("'xmpz' does not support Infinity");
            return NULL;
        }
        /* return code is ignored */
        mpfr_get_z(result->z, MPFR(self), GET_MPFR_ROUND(context));
    }

    return result;
}

/* Return the simpliest rational number that approximates 'self' to the
 * requested precision 'err'. If 'err' is negative, then the requested
 * precision is -2**abs(int(err)). If 'err' is NULL, then the requested
 * precision is -2**prec. If 'prec' is 0, then the requested precision is
 * the precision of 'self'.
 */

static MPQ_Object *
stern_brocot(MPFR_Object* self, MPFR_Object *err, mpfr_prec_t prec, int mayz, CTXT_Object *context)
{
    MPQ_Object *result = 0;
    int i, negative, errsign;
    mpfr_t f, al, a, r1[3], r2[3], minerr, curerr, newerr, temp;

    CHECK_CONTEXT_SET_EXPONENT(context);

#define F2Q_PREC 20

    if (mpfr_nan_p(self->f)) {
        VALUE_ERROR("Cannot convert NaN to a number.");
        return NULL;
    }

    if (mpfr_inf_p(self->f)) {
        OVERFLOW_ERROR("Cannot convert Infinity to a number.");
        return NULL;
    }

    if (prec == 0)
        prec = mpfr_get_prec(self->f);

    errsign = err ? mpfr_sgn(err->f) : 0;
    if (errsign < 0)
        prec = (mpfr_prec_t)(-mpfr_get_si(err->f, GET_MPFR_ROUND(context)));

    if (errsign <= 0 && (prec < 2 || prec > mpfr_get_prec(self->f))) {
        VALUE_ERROR("Requested precision out-of-bounds.");
        return NULL;
    }

    if (!(result = GMPy_MPQ_New(context)))
        return NULL;

    mpfr_init2(minerr, F2Q_PREC);
    if (errsign <= 0) {
        mpfr_set_ui(minerr, 1, MPFR_RNDN);
        mpfr_div_2si(minerr, minerr, prec, GET_MPFR_ROUND(context));
    }
    else {
        mpfr_set(minerr, err->f, GET_MPFR_ROUND(context));
    }

    mpfr_init2(f, prec);
    if (mpfr_sgn(self->f) < 0) {
        negative = 1;
        mpfr_abs(f, self->f, GET_MPFR_ROUND(context));
    }
    else {
        negative = 0;
        mpfr_set(f, self->f, GET_MPFR_ROUND(context));
    }

    mpfr_init2(al, prec);
    mpfr_set(al, f, GET_MPFR_ROUND(context));
    mpfr_init2(a, prec);
    mpfr_floor(a, al);
    mpfr_init2(temp, prec);
    for (i=0; i<3; ++i) {
        mpfr_init2(r1[i], prec);
        mpfr_init2(r2[i], prec);
    }
    mpfr_set_si(r1[0], 0, MPFR_RNDN);
    mpfr_set_si(r1[1], 0, MPFR_RNDN);
    mpfr_set_si(r1[2], 1, MPFR_RNDN);
    mpfr_set_si(r2[0], 0, MPFR_RNDN);
    mpfr_set_si(r2[1], 1, MPFR_RNDN);
    mpfr_set(r2[2], a, GET_MPFR_ROUND(context));
    mpfr_init2(curerr, F2Q_PREC);
    mpfr_init2(newerr, F2Q_PREC);
    mpfr_reldiff(curerr, f, a, GET_MPFR_ROUND(context));
    while (mpfr_cmp(curerr, minerr) > 0) {
        mpfr_sub(temp, al, a, GET_MPFR_ROUND(context));
        mpfr_ui_div(al, 1, temp, GET_MPFR_ROUND(context));
        mpfr_floor(a, al);
        mpfr_swap(r1[0], r1[1]);
        mpfr_swap(r1[1], r1[2]);
        mpfr_mul(r1[2], r1[1], a, GET_MPFR_ROUND(context));
        mpfr_add(r1[2], r1[2], r1[0], GET_MPFR_ROUND(context));
        mpfr_swap(r2[0], r2[1]);
        mpfr_swap(r2[1], r2[2]);
        mpfr_mul(r2[2], r2[1], a, GET_MPFR_ROUND(context));
        mpfr_add(r2[2], r2[2], r2[0], GET_MPFR_ROUND(context));
        mpfr_div(temp, r2[2], r1[2], GET_MPFR_ROUND(context));
        mpfr_reldiff(newerr, f, temp, GET_MPFR_ROUND(context));
        if (mpfr_cmp(curerr, newerr) <= 0) {
            mpfr_swap(r1[1],r1[2]);
            mpfr_swap(r2[1],r2[2]);
            break;
        }
        mpfr_swap(curerr, newerr);
    }

    if (mayz && (mpfr_cmp_ui(r1[2],1) == 0)) {
        Py_DECREF((PyObject*)result);
        result = (MPQ_Object*)GMPy_MPZ_New(context);
        mpfr_get_z(MPZ(result), r2[2], GET_MPFR_ROUND(context));
        if (negative)
            mpz_neg(MPZ(result), MPZ(result));
    }
    else {
        mpfr_get_z(mpq_numref(result->q), r2[2], GET_MPFR_ROUND(context));
        mpfr_get_z(mpq_denref(result->q), r1[2], GET_MPFR_ROUND(context));
        if (negative)
            mpz_neg(mpq_numref(result->q), mpq_numref(result->q));
    }

    mpfr_clear(minerr);
    mpfr_clear(al);
    mpfr_clear(a);
    mpfr_clear(f);
    for (i=0; i<3; ++i) {
        mpfr_clear(r1[i]);
        mpfr_clear(r2[i]);
    }
    mpfr_clear(curerr);
    mpfr_clear(newerr);
    mpfr_clear(temp);
    return result;
}

static MPQ_Object *
GMPy_MPQ_From_MPFR(MPFR_Object *self, CTXT_Object *context)
{
    return stern_brocot((MPFR_Object*)self, 0, 0, 0, context);
}

static PyObject *
GMPy_PyIntOrLong_From_MPFR(MPFR_Object *obj, CTXT_Object *context)
{
    PyObject *result;
    MPZ_Object *tempz;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(tempz = GMPy_MPZ_From_MPFR(obj, context)))
        return NULL;

    result = GMPy_PyIntOrLong_From_MPZ(tempz, context);
    Py_DECREF((PyObject*)tempz);

    return result;
}

static PyObject *
GMPy_MPFR_Int_Slot(MPFR_Object *self) {
    return GMPy_PyIntOrLong_From_MPFR(self, NULL);
}

static PyObject *
GMPy_PyLong_From_MPFR(MPFR_Object *obj, CTXT_Object *context)
{
    PyObject *result;
    MPZ_Object *tempz;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(tempz = GMPy_MPZ_From_MPFR(obj, context)))
        return NULL;

    result = GMPy_PyLong_From_MPZ(tempz, context);
    Py_DECREF((PyObject*)tempz);

    return result;
}

static PyObject *
GMPy_MPFR_Long_Slot(MPFR_Object *self) {
    return GMPy_PyLong_From_MPFR(self, NULL);
}

static PyObject *
GMPy_PyFloat_From_MPFR(MPFR_Object *self, CTXT_Object *context)
{
    double res;

    CHECK_CONTEXT_SET_EXPONENT(context);

    res = mpfr_get_d(self->f, GET_MPFR_ROUND(context));

    return PyFloat_FromDouble(res);
}

static PyObject *
GMPy_MPFR_Float_Slot(MPFR_Object *self) {
    return GMPy_PyFloat_From_MPFR(self, NULL);
}

static PyObject*
GMPy_PyStr_From_MPFR(MPFR_Object *self, int base, int digits, CTXT_Object *context)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;

    CHECK_CONTEXT_SET_EXPONENT(context);

    /* check arguments are valid */
    assert(MPFR_Check((PyObject*)self));
    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    /* Process special cases first */
    if (!(mpfr_regular_p(self->f))) {
        if (mpfr_nan_p(self->f)) {
            return Py_BuildValue("(sii)", "nan", 0, 0);
        }
        else if (mpfr_inf_p(self->f) && !mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "inf", 0, 0);
        }
        else if (mpfr_inf_p(self->f) && mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "-inf", 0, 0);
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self->f)) {
            return Py_BuildValue("(sii)", "-0", 0, mpfr_get_prec(self->f));
        }
        else {
            return Py_BuildValue("(sii)", "0", 0, mpfr_get_prec(self->f));
        }
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self->f, GET_MPFR_ROUND(context));
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in Pympfr_To_PyStr");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self->f));
    mpfr_free_str(buffer);
    return result;
}

/*
 * coerce any number to a mpf
 */

int
GMPy_MPFR_convert_arg(PyObject *arg, PyObject **ptr)
{
    MPFR_Object* newob = GMPy_MPFR_From_Real(arg, 1, NULL);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("argument can not be converted to 'mpfr'");
        return 0;
    }
}

/* str and repr implementations for mpfr */
static PyObject *
GMPy_MPFR_Str_Slot(MPFR_Object *self)
{
    PyObject *result, *temp;
    long precision;
    char fmtstr[30];

    precision = (long)(log10(2) * (double)mpfr_get_prec(MPFR(self))) + 2;

    sprintf(fmtstr, "{0:.%ldg}", precision);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
GMPy_MPFR_Repr_Slot(MPFR_Object *self)
{
    PyObject *result, *temp;
    long precision, bits;
    char fmtstr[30];

    bits = mpfr_get_prec(MPFR(self));
    precision = (long)(log10(2) * (double)bits) + 2;

    if (mpfr_number_p(MPFR(self)) && bits != DBL_MANT_DIG)
        sprintf(fmtstr, "mpfr('{0:.%ldg}',%ld)", precision, bits);
    else
        sprintf(fmtstr, "mpfr('{0:.%ldg}')", precision);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
mpfr_ascii(mpfr_t self, int base, int digits, int round)
{
    PyObject *result;
    char *buffer;
    mpfr_exp_t the_exp;

    /* Process special cases first */
    if (!(mpfr_regular_p(self))) {
        if (mpfr_nan_p(self)) {
            return Py_BuildValue("(sii)", "nan", 0, 0);
        }
        else if (mpfr_inf_p(self) && !mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "inf", 0, 0);
        }
        else if (mpfr_inf_p(self) && mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "-inf", 0, 0);
        }
        /* 0 is not considered a 'regular" number */
        else if (mpfr_signbit(self)) {
            return Py_BuildValue("(sii)", "-0", 0, mpfr_get_prec(self));
        }
        else {
            return Py_BuildValue("(sii)", "0", 0, mpfr_get_prec(self));
        }
    }

    /* obtain digits-string and exponent */
    buffer = mpfr_get_str(0, &the_exp, base, digits, self, round);
    if (!*buffer) {
        SYSTEM_ERROR("Internal error in mpfr_ascii");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self));
    mpfr_free_str(buffer);
    return result;
}
