/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_mpc.c                                                     *
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

/* Return a copy of an mpc. If the value for rprec and iprec is 0, then the
 * context's precision is used. If the values for rprec and iprec are >= 2,
 * then their value will be used. This function will always return a new
 * instance.
 */

static MPC_Object *
GMPy_MPC_From_MPC_New(MPC_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;
    mpfr_prec_t tempr = 0, tempi = 0;
    int rr, ri, dr, di;

    assert(MPC_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!rprec)
        rprec = GET_REAL_PREC(context);

    if (!iprec)
        iprec = GET_IMAG_PREC(context);

    if (MPC_CheckAndExp(obj)) {
        /* The exponents are valid in the current context. */
        if ((result = GMPy_MPC_New(rprec, iprec, context)))
            result->rc = mpc_set(result->c, obj->c, GET_MPC_ROUND(context));
        return result;
    }
    else {
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing mpc incompatible with current context");
            return NULL;
        }

        /* Get the real & imaginary precisions of the source. */
        mpc_get_prec2(&tempr, &tempi, obj->c);
        /* Get the real & imaginary ternary result codes of the source. */
        rr = MPC_INEX_RE(obj->rc);
        ri = MPC_INEX_IM(obj->rc);
        /* Get the real & imaginary rounding modes of the source. */
        dr = MPC_RND_RE(obj->round_mode);
        di = MPC_RND_IM(obj->round_mode);

        if ((result = GMPy_MPC_New(tempr, tempi, context))) {
            /* First make the exponent valid. */
            mpc_set(result->c, obj->c, GET_MPC_ROUND(context));
            rr = mpfr_check_range(mpc_realref(result->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(result->c), ri, di);
            /* Round to the desired precision. */
            rr = mpfr_prec_round(mpc_realref(result->c), rprec, GET_REAL_ROUND(context));
            ri = mpfr_prec_round(mpc_imagref(result->c), iprec, GET_IMAG_ROUND(context));
            result->rc = MPC_INEX(rr, ri);
        }
        return result;
    }
}

/* Return a new reference to an existing mpc if the exponents are valid in the
 * current context. If the exponents are not valid, a reference to a new, valid
 * instance is returned.
 *
 * Note: the precision will not be changed.
 *
 * All mpc arguments to functions in the MPC library should go through this
 * function to guarantee that the exponents are valid. References returned by
 * function should not be returned to the user (although they must still be
 * decremented like normal references).
 */

static MPC_Object *
GMPy_MPC_From_MPC_Temp(MPC_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;
    mpfr_prec_t tempr = 0, tempi = 0;
    int rr, ri, dr, di;

    assert(MPC_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (MPC_CheckAndExp(obj)) {
        /* The exponents are valid in the current context. */
        Py_INCREF((PyObject*)obj);
        return result;
    }
    else {
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing mpc incompatible with current context");
            return NULL;
        }

        mpc_get_prec2(&tempr, &tempi, obj->c);
        rr = MPC_INEX_RE(obj->rc);
        ri = MPC_INEX_IM(obj->rc);
        dr = MPC_RND_RE(obj->round_mode);
        di = MPC_RND_IM(obj->round_mode);

        if ((result = GMPy_MPC_New(tempr, tempi, context))) {
            /* First make the exponent valid. */
            mpc_set(result->c, obj->c, GET_MPC_ROUND(context));
            rr = mpfr_check_range(mpc_realref(result->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(result->c), ri, di);
            /* Round to the desired precision. */
            result->rc = MPC_INEX(rr, ri);
        }
        return result;
    }
}

static MPC_Object *
GMPy_MPC_From_PyComplex(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                        CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPC_New(rprec, iprec, context)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(obj),
                    PyComplex_ImagAsDouble(obj), GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
GMPy_MPC_From_MPFR(MPFR_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                   CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!rprec)
        rprec = mpfr_get_prec(obj->f);

    if ((result = GMPy_MPC_New(rprec, iprec, context)))
        result->rc = mpc_set_fr(result->c, obj->f,GET_MPC_ROUND(context));

    return result;
}

static MPC_Object *
GMPy_MPC_From_PyFloat(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                   CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!rprec)
        rprec = DBL_MANT_DIG;

    if ((result = GMPy_MPC_New(rprec, iprec, context)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(obj),
                               GET_MPC_ROUND(context));

    return result;
}

static MPC_Object *
GMPy_MPC_From_MPZ(MPZ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;

    assert(MPZ_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPC_New(rprec, iprec, context)))
        result->rc = mpc_set_z(result->c, obj->z, GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
GMPy_MPC_From_MPQ(MPQ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;

    assert(MPQ_Object(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPC_New(rprec, iprec, context)))
        result->rc = mpc_set_q(result->c, obj->q, GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
GMPy_MPC_From_Fraction(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                       CTXT_Object *context)
{
    MPC_Object *result = NULL;
    MPQ_Object *tempq;

    assert(IS_RATIONALt(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(tempq = GMPy_MPQ_From_Fraction(obj, context))) {
        result = GMPy_MPC_From_MPQ(tempq, rprec, iprec, context);
        Py_DECREF((PyObject*)tempq);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_Decimal(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                      CTXT_Object *context)
{
    MPC_Object *result = NULL;
    MPFR_Object *tempf;

    assert(IS_DECIMAL(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(tempf = GMPy_MPFR_From_Decimal(obj, rprec, context))) {
        result = GMPy_MPC_From_MPFR(tempf, rprec, iprec, context);
        Py_DECREF((PyObject*)tempf);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_PyIntOrLong(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                          CTXT_Object *context)
{
    MPC_Object *result = NULL;
    MPZ_Object *tempz;

    assert(PyIntOrLong_Check(obj));

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!(tempz = GMPy_MPZ_From_PyIntOrLong(obj, context))) {
        result = GMPy_MPC_From_MPZ(tempz, rprec, iprec, context);
        Py_DECREF((PyObject*)tempz);
    }
    return result;
}

/* Python's string representation of a complex number differs from the format
 * used by MPC. Both MPC and Python surround the complex number with '(' and
 * ')' but Python adds a 'j' after the imaginary component and MPC requires a
 * space between the real and imaginery components. PyStr2Pympc tries to work
 * around the differences as follows reading two MPFR-compatible numbers from
 * the string and storing into the real and imaginary components respectively.
 */

static MPC_Object *
GMPy_MPC_From_PyStr(PyObject *s, int base, mpfr_prec_t rprec, mpfr_prec_t iprec,
                    CTXT_Object *context)
{
    MPC_Object *result;
    PyObject *ascii_str = NULL;
    Py_ssize_t len;
    char *cp, *unwind, *tempchar, *lastchar;
    int firstp = 0, lastp = 0, real_rc = 0, imag_rc = 0;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (PyBytes_Check(s)) {
        len = PyBytes_Size(s);
        cp = (char*)PyBytes_AsString(s);
    }
    else if (PyUnicode_Check(s)) {
        ascii_str = PyUnicode_AsASCIIString(s);
        if (!ascii_str) {
            VALUE_ERROR("string contains non-ASCII characters");
            return NULL;
        }
        len = PyBytes_Size(ascii_str);
        cp = (char*)PyBytes_AsString(ascii_str);
    }
    else {
        TYPE_ERROR("string required");
        return NULL;
    }

    /* Don't allow NULL characters */
    if (strlen(cp) != len) {
        VALUE_ERROR("string without NULL characters expected");
        Py_XDECREF(ascii_str);
        return NULL;
    }

    if (!(result = GMPy_MPC_New(rprec, iprec, context))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Get a pointer to the last valid character (ignoring trailing
     * whitespace.) */
    lastchar = cp + len - 1;
    while (isspace(*lastchar))
        lastchar--;

    /* Skip trailing ). */
    if (*lastchar == ')') {
        lastp = 1;
        lastchar--;
    }

    /* Skip trailing j. */
    if (*lastchar == 'j')
        lastchar--;

    /* Skip leading whitespace. */
    while (isspace(*cp))
        cp++;

    /* Skip a leading (. */
    if (*cp == '(') {
        firstp = 1;
        cp++;
    }

    if (firstp != lastp) goto invalid_string;

    /* Read the real component first. */
    unwind = cp;
    real_rc = mpfr_strtofr(mpc_realref(result->c), cp, &tempchar, base,
                           GET_REAL_ROUND(context));

    /* Verify that at least one valid character was read. */
    if (cp == tempchar) goto invalid_string;

    /* If the next character is a j, then the real component is 0 and
     * we just read the imaginary componenet.
     */
    if (*tempchar == 'j') {
        mpfr_set_zero(mpc_realref(result->c), +1);
        cp = unwind;
    }
    else {
        /* Read the imaginary component next. */
        cp = tempchar;
    }

    imag_rc = mpfr_strtofr(mpc_imagref(result->c), cp, &tempchar, base,
                           GET_IMAG_ROUND(context));

    if (cp == tempchar && tempchar > lastchar)
        goto valid_string;

    if (*tempchar != 'j' && *cp != ' ')
        goto invalid_string;

    if (tempchar <= lastchar)
        goto invalid_string;

  valid_string:
    Py_XDECREF(ascii_str);
    result->rc = MPC_INEX(real_rc, imag_rc);
    return result;

  invalid_string:
    VALUE_ERROR("invalid string in mpc()");
    Py_DECREF((PyObject*)result);
    Py_XDECREF(ascii_str);
    return NULL;
}

/* See the comments for GMPy_MPFR_From_Real_Temp. */

static MPC_Object *
GMPy_MPC_From_Complex_Temp(PyObject* obj, mp_prec_t rprec, mp_prec_t iprec,
                           CTXT_Object *context)
{
    CHECK_CONTEXT_SET_EXPONENT(context);

    if (MPC_Check(obj))
        return GMPy_MPC_From_MPC_Temp((MPC_Object*)obj, rprec, iprec, context);

    if (MPFR_Check(obj))
        return GMPy_MPC_From_MPFR((MPFR_Object*)obj,
                                  mpfr_get_prec(MPFR(obj)),
                                  mpfr_get_prec(MPFR(obj)),
                                  context
                                 );

    if (PyFloat_Check(obj))
        return GMPy_MPC_From_PyFloat(obj, 53, 53, context);

    if (PyComplex_Check(obj))
        return GMPy_MPC_From_PyComplex(obj, 53, 53, context);

    if (MPQ_Check(obj))
        return GMPy_MPC_From_MPQ((MPQ_Object*)obj, rprec, iprec, context);

    if (MPZ_Check(obj) || XMPZ_Check(obj))
        return GMPy_MPC_From_MPZ((MPZ_Object*)obj, rprec, iprec, context);

    if (PyIntOrLong_Check(obj))
        return GMPy_MPC_From_PyIntOrLong(obj, rprec, iprec, context);

    if (IS_DECIMAL(obj))
        return GMPy_MPC_From_Decimal(obj, rprec, iprec, context);

    if (IS_FRACTION(obj))
        return GMPy_MPC_From_Fraction(obj, rprec, iprec, context);

    TYPE_ERROR("object could not be converted to 'mpc'");
    return NULL;
}

static MPC_Object *
GMPy_MPC_From_Complex_New(PyObject* obj, mp_prec_t rprec, mp_prec_t iprec,
                          CTXT_Object *context)
{
    CHECK_CONTEXT_SET_EXPONENT(context);

    if (MPC_Check(obj))
        return GMPy_MPC_From_MPC_New((MPC_Object*)obj, rprec, rprec, context);

    if (MPFR_Check(obj))
        return GMPy_MPC_From_MPFR((MPFR_Object*)obj,
                                  mpfr_get_prec(MPFR(obj)),
                                  mpfr_get_prec(MPFR(obj)),
                                  context
                                 );

    if (PyFloat_Check(obj))
        return GMPy_MPC_From_PyFloat(obj, rprec, iprec, context);

    if (PyComplex_Check(obj))
        return GMPy_MPC_From_PyComplex(obj, rprec, iprec, context);

    if (MPQ_Check(obj))
        return GMPy_MPC_From_MPQ((MPQ_Object*)obj, rprec, iprec, context);

    if (MPZ_Check(obj) || XMPZ_Check(obj))
        return GMPy_MPC_From_MPZ((MPZ_Object*)obj, rprec, iprec, context);

    if (PyIntOrLong_Check(obj))
        return GMPy_MPC_From_PyIntOrLong(obj, rprec, iprec, context);

    if (IS_DECIMAL(obj))
        return GMPy_MPC_From_Decimal(obj, rprec, iprec, context);

    if (IS_FRACTION(obj))
        return GMPy_MPC_From_Fraction(obj, rprec, iprec, context);

    TYPE_ERROR("object could not be converted to 'mpc'");
    return NULL;
}

static PyObject *
GMPy_PyStr_From_MPC(MPC_Object *self, int base, int digits, CTXT_Object *context)
{
    PyObject *tempreal = 0, *tempimag = 0, *result;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    tempreal = mpfr_ascii(mpc_realref(self->c), base, digits,
                          MPC_RND_RE(GET_MPC_ROUND(context)));
    tempimag = mpfr_ascii(mpc_imagref(self->c), base, digits,
                          MPC_RND_IM(GET_MPC_ROUND(context)));

    if (!tempreal || !tempimag) {
        Py_XDECREF(tempreal);
        Py_XDECREF(tempimag);
        return NULL;
    }

    result = Py_BuildValue("(NN)", tempreal, tempimag);
    if (!result) {
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
    }
    return result;
}

static PyObject *
Pympc_To_PyFloat(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'float'");
    return NULL;
}

PyDoc_STRVAR(doc_mpc_complex, "Convert 'mpc' to 'complex'.");

static PyObject *
Pympc_To_PyComplex(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    double real = mpfr_get_d(mpc_realref(MPC(self)), GET_REAL_ROUND(context));
    double imag = mpfr_get_d(mpc_imagref(MPC(self)), GET_IMAG_ROUND(context));

    return PyComplex_FromDoubles(real, imag);
}

static PyObject *
Pympc_To_PyComplex_context(PyObject *self, PyObject *other,
                           CTXT_Object *context)
{
    CHECK_CONTEXT_SET_EXPONENT(context);

    double real = mpfr_get_d(mpc_realref(MPC(self)), GET_REAL_ROUND(context));
    double imag = mpfr_get_d(mpc_imagref(MPC(self)), GET_IMAG_ROUND(context));

    return PyComplex_FromDoubles(real, imag);
}

static PyObject *
Pympc_To_PyLong(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'long'");
    return NULL;
}

/*
 * coerce any number to a mpc
 */

int
GMPy_MPC_convert_arg(PyObject *arg, PyObject **ptr)
{
    MPC_Object *newob = GMPy_MPC_From_Complex_Temp(arg, 0, 0, NULL);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("can't convert argument to 'mpc'");
        return 0;
    }
}

/* str and repr implementations for mpc */
static PyObject *
Pympc_To_Str(MPC_Object *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[30];

    mpc_get_prec2(&rbits, &ibits, MPC(self));
    rprec = (long)(log10(2) * (double)rbits) + 2;
    iprec = (long)(log10(2) * (double)ibits) + 2;

    sprintf(fmtstr, "{0:.%ld.%ldg}", rprec, iprec);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympc_To_Repr(MPC_Object *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[30];

    mpc_get_prec2(&rbits, &ibits, MPC(self));
    rprec = (long)(log10(2) * (double)rbits) + 2;
    iprec = (long)(log10(2) * (double)ibits) + 2;

    if (rbits != DBL_MANT_DIG || ibits !=DBL_MANT_DIG)
        sprintf(fmtstr, "mpc('{0:.%ld.%ldg}',(%ld,%ld))",
                rprec, iprec, rbits, ibits);
    else
        sprintf(fmtstr, "mpc('{0:.%ld.%ldg}')", rprec, iprec);

    temp = Py_BuildValue("s", fmtstr);
    if (!temp)
        return NULL;
    result = PyObject_CallMethod(temp, "format", "O", self);
    Py_DECREF(temp);
    return result;
}



