/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_convert_mpc.c                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018 Case Van Horsen                        *
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

/* Please see the documentation for GMPy_MPFR_From_MPFR for detailed
 * description of this function.
 */


static MPC_Object *
GMPy_MPC_From_MPC(MPC_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;

    assert(MPC_Check(obj));

    /* Optimize the critical case when prec==1 or obj is NaN or Inf. */

    if ((rprec == 1 && iprec == 1) ||
        (!mpfr_number_p(mpc_realref(obj->c)) && !mpfr_number_p(mpc_imagref(obj->c)))) {
        Py_INCREF((PyObject*)obj);
        return obj;
    }

    CHECK_CONTEXT(context);

    if (rprec == 0)
        rprec = GET_REAL_PREC(context);
    else if (rprec == 1)
        rprec = mpfr_get_prec(mpc_realref(obj->c));

    if (iprec == 0)
        iprec = GET_IMAG_PREC(context);
    else if (iprec == 1)
        iprec = mpfr_get_prec(mpc_imagref(obj->c));

    /* Try to identify when an additional reference to existing instance can
     * be returned. It is possible when (1) the precision matches, (2) the
     * exponent is valid and not in the range that might require subnormal-
     * ization, and (3) subnormalize is not enabled.
     */

    if ((rprec == mpfr_get_prec(mpc_realref(obj->c))) &&
        (iprec == mpfr_get_prec(mpc_imagref(obj->c))) &&
        (!context->ctx.subnormalize) &&
        (mpc_realref(obj->c)->_mpfr_exp >= (context->ctx.emin + mpfr_get_prec(mpc_realref(obj->c)) - 1)) &&
        (mpc_realref(obj->c)->_mpfr_exp <= context->ctx.emax) &&
        (mpc_imagref(obj->c)->_mpfr_exp >= (context->ctx.emin + mpfr_get_prec(mpc_imagref(obj->c)) - 1)) &&
        (mpc_imagref(obj->c)->_mpfr_exp <= context->ctx.emax)
        ) {

        Py_INCREF((PyObject*)obj);
        return obj;
    }

    if ((result = GMPy_MPC_New(rprec, iprec, context))) {
        result->rc = mpc_set(result->c, obj->c, GET_MPC_ROUND(context));
        _GMPy_MPC_Cleanup(&result, context);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_PyComplex(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                        CTXT_Object *context)
{
    MPC_Object *result;

    assert(PyComplex_Check(obj));

    CHECK_CONTEXT(context);

    if (rprec == 0)
        rprec = GET_REAL_PREC(context);
    else if (rprec == 1)
        rprec = DBL_MANT_DIG;

    if (iprec == 0)
        iprec = GET_IMAG_PREC(context);
    else if (iprec == 1)
        rprec = DBL_MANT_DIG;

    if ((result = GMPy_MPC_New(rprec, iprec, context))) {
        result->rc = mpc_set_d_d(result->c, PyComplex_RealAsDouble(obj),
                                 PyComplex_ImagAsDouble(obj), GET_MPC_ROUND(context));
        if (rprec != 1 || iprec != 1) {
            GMPY_MPC_CHECK_RANGE(result, context);
        }
        GMPY_MPC_SUBNORMALIZE(result, context);
        GMPY_MPC_EXCEPTIONS(result, context);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_MPFR(MPFR_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                   CTXT_Object *context)
{
    MPC_Object *result;

    assert(MPFR_Check(obj));

    CHECK_CONTEXT(context);

    if (rprec == 0)
        rprec = GET_REAL_PREC(context);
    else if (rprec == 1)
        rprec = mpfr_get_prec(obj->f);

    if (iprec == 0)
        iprec = GET_IMAG_PREC(context);
    else if (iprec == 1)
        rprec = mpfr_get_prec(obj->f);

    if ((result = GMPy_MPC_New(rprec, iprec, context))) {
        result->rc = mpc_set_fr(result->c, obj->f, GET_MPC_ROUND(context));
        if (rprec != 1) {
            GMPY_MPC_CHECK_RANGE(result, context);
        }
        GMPY_MPC_SUBNORMALIZE(result, context);
        GMPY_MPC_EXCEPTIONS(result, context);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_PyFloat(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                   CTXT_Object *context)
{
    MPC_Object *result;

    CHECK_CONTEXT(context);

    assert(PyFloat_Check(obj));

    if (rprec == 0)
        rprec = GET_REAL_PREC(context);
    else if (rprec == 1)
        rprec = DBL_MANT_DIG;

    if (iprec == 0)
        iprec = GET_IMAG_PREC(context);
    else if (iprec == 1)
        rprec = DBL_MANT_DIG;

    if ((result = GMPy_MPC_New(rprec, iprec, context))) {
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(obj),
                               GET_MPC_ROUND(context));
        if (rprec != 1) {
            GMPY_MPC_CHECK_RANGE(result, context);
        }
        GMPY_MPC_SUBNORMALIZE(result, context);
        GMPY_MPC_EXCEPTIONS(result, context);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_MPZ(MPZ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;

    assert(CHECK_MPZANY(obj));

    CHECK_CONTEXT(context);

    if (rprec < 2) {
        rprec = GET_REAL_PREC(context);
    }

    if (iprec < 2) {
        iprec = GET_IMAG_PREC(context);
    }

    if ((result = GMPy_MPC_New(rprec, iprec, context))) {
        result->rc = mpc_set_z(result->c, obj->z, GET_MPC_ROUND(context));
        if (rprec != 1) {
            GMPY_MPC_CHECK_RANGE(result, context);
        }
        GMPY_MPC_SUBNORMALIZE(result, context);
        GMPY_MPC_EXCEPTIONS(result, context);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_MPQ(MPQ_Object *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                  CTXT_Object *context)
{
    MPC_Object *result = NULL;

    assert(MPQ_Check(obj));

    CHECK_CONTEXT(context);

    if (rprec < 2) {
        rprec = GET_REAL_PREC(context);
    }

    if (iprec < 2) {
        iprec = GET_IMAG_PREC(context);
    }

    if ((result = GMPy_MPC_New(rprec, iprec, context))) {
        result->rc = mpc_set_q(result->c, obj->q, GET_MPC_ROUND(context));
        if (rprec != 1) {
            GMPY_MPC_CHECK_RANGE(result, context);
        }
        GMPY_MPC_SUBNORMALIZE(result, context);
        GMPY_MPC_EXCEPTIONS(result, context);
    }
    return result;
}

static MPC_Object *
GMPy_MPC_From_Fraction(PyObject *obj, mpfr_prec_t rprec, mpfr_prec_t iprec,
                       CTXT_Object *context)
{
    MPC_Object *result = NULL;
    MPQ_Object *tempq;

    assert(IS_RATIONAL(obj));

    CHECK_CONTEXT(context);

    if ((tempq = GMPy_MPQ_From_Fraction(obj, context))) {
        result = GMPy_MPC_From_MPQ(tempq, rprec, iprec, context);
        Py_DECREF((PyObject*)tempq);
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

    CHECK_CONTEXT(context);

    if ((tempz = GMPy_MPZ_From_PyIntOrLong(obj, context))) {
        result = GMPy_MPC_From_MPZ(tempz, rprec, iprec, context);
        Py_DECREF((PyObject*)tempz);
    }

    return result;
}

/* Python's string representation of a complex number differs from the format
 * used by MPC. Both MPC and Python surround the complex number with '(' and
 * ')' but Python adds a 'j' after the imaginary component and MPC requires a
 * space between the real and imaginery components. GMPy_MPC_From_PyStr tries
 * to work around the differences as by reading two MPFR-compatible numbers
 * from the string and storing into the real and imaginary components
 * respectively.
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

    CHECK_CONTEXT(context);

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
    if ((Py_ssize_t)strlen(cp) != len) {
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
        mpfr_set_zero(mpc_realref(result->c), MPFR_RNDN);
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

    if (rprec != 1 || iprec != 1) {
        GMPY_MPC_CHECK_RANGE(result, context);
    }
    GMPY_MPC_SUBNORMALIZE(result, context);
    GMPY_MPC_EXCEPTIONS(result, context);

    return result;

  invalid_string:
    VALUE_ERROR("invalid string in mpc()");
    Py_DECREF((PyObject*)result);
    Py_XDECREF(ascii_str);
    return NULL;
}

/* See the comments for GMPy_MPFR_From_Real_Temp. */

static MPC_Object *
GMPy_MPC_From_Complex(PyObject* obj, mp_prec_t rprec, mp_prec_t iprec,
                           CTXT_Object *context)
{
    CHECK_CONTEXT(context);

    if (MPC_Check(obj))
        return GMPy_MPC_From_MPC((MPC_Object*)obj, rprec, iprec, context);

    if (MPFR_Check(obj))
        return GMPy_MPC_From_MPFR((MPFR_Object*)obj, rprec, iprec, context);

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

    if (IS_FRACTION(obj))
        return GMPy_MPC_From_Fraction(obj, rprec, iprec, context);

    TYPE_ERROR("object could not be converted to 'mpc'");
    return NULL;
}

static PyObject *
GMPy_PyStr_From_MPC(MPC_Object *self, int base, int digits, CTXT_Object *context)
{
    PyObject *tempreal = 0, *tempimag = 0, *result;

    CHECK_CONTEXT(context);

    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval [2,62]");
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
GMPy_MPC_Float_Slot(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'float'");
    return NULL;
}

PyDoc_STRVAR(GMPy_doc_mpc_complex, "Convert 'mpc' to 'complex'.");

static PyObject *
GMPy_PyComplex_From_MPC(PyObject *self, PyObject *other)
{
    CTXT_Object *context = NULL;
    double real, imag;

    CHECK_CONTEXT(context);

    real = mpfr_get_d(mpc_realref(MPC(self)), GET_REAL_ROUND(context));
    imag = mpfr_get_d(mpc_imagref(MPC(self)), GET_IMAG_ROUND(context));

    return PyComplex_FromDoubles(real, imag);
}

#ifdef PY2
static PyObject *
GMPy_MPC_Long_Slot(PyObject *self)
{
    TYPE_ERROR("can't covert mpc to long");
    return NULL;
}
#endif

static PyObject *
GMPy_MPC_Int_Slot(PyObject *self)
{
    TYPE_ERROR("can't covert mpc to int");
    return NULL;
}

#ifdef SHARED
/*
 * coerce any number to a mpc
 */

int
GMPy_MPC_ConvertArg(PyObject *arg, PyObject **ptr)
{
    MPC_Object *newob = GMPy_MPC_From_Complex(arg, 0, 0, NULL);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("can't convert argument to 'mpc'");
        return 0;
    }
}
#endif

/* str and repr implementations for mpc */
static PyObject *
GMPy_MPC_Str_Slot(MPC_Object *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[60];

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
GMPy_MPC_Repr_Slot(MPC_Object *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[60];

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



