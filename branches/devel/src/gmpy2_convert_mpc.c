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

/* This file contains all the conversion functions for gmpy2.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. gmpy2 uses a LBYL (Look Before You Leap) approach and
 * identifies the numeric type before conversion before conversion to a gmpy2
 * type. The basic operations (+, -, *, /) are optimized to directly work with
 * some basic types such as C longs or doubles.
 *
 * Support for the Decimal type is a challenge. For the basic operations, it
 * is most accurate to convert a Decimal instance into an mpq and then use
 * MPFR's functions to accurately operate on an mpfr and mpq. This approach is
 * challenging because (1) a large exponent can create a very large mpq and
 * (2) the changes made to C-coded version of Decimal in Python 3.3.
 *
 */

static MPC_Object *
Pympc_From_Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, MPC(self));
    if ((result = (MPC_Object*)Pympc_new(rprec, iprec)))
        mpc_set(result->c, MPC(self), GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympc_bits_context(PyObject *self, mpfr_prec_t rprec,
                              mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *result;

    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, MPC(self));
    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        mpc_set(result->c, MPC(self), GET_MPC_ROUND(context));
    return result;
}

/* Return an mpc instance based on the context. If the precision of self is
 * the same as the context's precision, then a new reference is created. If
 * the precisions are different, then a new object is created. */

static MPC_Object *
Pympc_From_Pympc_context(PyObject *self, GMPyContextObject *context)
{
    mpfr_prec_t rprec, iprec;
    MPC_Object *result;

    mpc_get_prec2(&rprec, &iprec, MPC(self));
    if ((rprec == GET_REAL_PREC(context)) && (iprec == GET_IMAG_PREC(context))) {
        Py_INCREF(self);
        return (MPC_Object*)self;
    }
    else {
        if ((result = (MPC_Object*)Pympc_new_context(context)))
            mpc_set(result->c, MPC(self), GET_MPC_ROUND(context));
        return result;
    }
}

static MPC_Object *
Pympc_From_PyComplex(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPC_Object*)Pympc_new(rprec, iprec)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(self),
                    PyComplex_ImagAsDouble(self), GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyComplex_bits_context(PyObject *self, mpfr_prec_t rprec,
                                  mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(self),
                    PyComplex_ImagAsDouble(self), GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyComplex_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_context(context)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(self),
                    PyComplex_ImagAsDouble(self), GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympfr(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (!rprec)
        rprec = mpfr_get_prec(MPFR(self));
    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_fr(result->c, MPFR(self),
                                GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympfr_bits_context(PyObject *self, mpfr_prec_t rprec,
                               mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *result;

    if (!rprec)
        rprec = mpfr_get_prec(MPFR(self));
    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_fr(result->c, MPFR(self),
                                GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympfr_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_context(context)))
        result->rc = mpc_set_fr(result->c, MPFR(self),
                                GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyFloat(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (!rprec)
        rprec = DBL_MANT_DIG;
    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(self),
                               GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyFloat_bits_context(PyObject *self, mpfr_prec_t rprec,
                                mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *result;

    if (!rprec)
        rprec = DBL_MANT_DIG;
    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(self),
                               GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyFloat_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_context(context)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(self),
                               GET_MPC_ROUND(context));
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
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    double real = mpfr_get_d(mpc_realref(MPC(self)),
                             GET_REAL_ROUND(context));
    double imag = mpfr_get_d(mpc_imagref(MPC(self)),
                             GET_IMAG_ROUND(context));

    return PyComplex_FromDoubles(real, imag);
}

static PyObject *
Pympc_To_PyComplex_context(PyObject *self, PyObject *other,
                           GMPyContextObject *context)
{
    double real = mpfr_get_d(mpc_realref(MPC(self)),
                             GET_REAL_ROUND(context));
    double imag = mpfr_get_d(mpc_imagref(MPC(self)),
                             GET_IMAG_ROUND(context));

    return PyComplex_FromDoubles(real, imag);
}

static MPC_Object *
Pympc_From_Pympz(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPC_Object*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_z(result->c, MPZ(self),
                                GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympz_bits_context(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec,
                         GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_z(result->c, MPZ(self),
                                GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympz_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_context(context)))
        result->rc = mpc_set_z(result->c, MPZ(self),
                                GET_MPC_ROUND(context));
    return result;
}

#define Pympc_From_Pyxmpz Pympc_From_Pympz
#define Pympc_From_Pyxmpz_bits_context Pympc_From_Pympz_bits_context
#define Pympc_From_Pyxmpz_context Pympc_From_Pympz_context

static MPC_Object *
Pympc_From_Pympq(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_q(result->c, MPQ(self),
                               GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympq_bits_context(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec,
                         GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context)))
        result->rc = mpc_set_q(result->c, MPQ(self),
                               GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_Pympq_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_context(context)))
        result->rc = mpc_set_q(result->c, MPQ(self),
                               GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyLong(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object *result;
    PyObject *temp = (PyObject*)GMPy_MPZ_From_PyLong(self);

    if (!temp)
        return NULL;
    result = Pympc_From_Pympz(temp, rprec, iprec);
    Py_DECREF(temp);
    return result;
}

static MPC_Object *
Pympc_From_PyLong_bits_context(PyObject *self, mpfr_prec_t rprec,
                               mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *result;
    PyObject *temp = (PyObject*)GMPy_MPZ_From_PyLong(self);

    if (!temp)
        return NULL;
    result = Pympc_From_Pympz_bits_context(temp, rprec, iprec, context);
    Py_DECREF(temp);
    return result;
}

static MPC_Object *
Pympc_From_PyLong_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;
    PyObject *temp = (PyObject*)GMPy_MPZ_From_PyLong(self);

    if (!temp)
        return NULL;
    result = Pympc_From_Pympz_context(temp, context);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympc_To_PyLong(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'long'");
    return NULL;
}

#ifdef PY2
static MPC_Object *
Pympc_From_PyInt_bits_context(PyObject *self, mpfr_prec_t rprec,
                              mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new(rprec, iprec)))
        result->rc = mpc_set_si(result->c, PyInt_AsLong(self),
                                GET_MPC_ROUND(context));
    return result;
}

static MPC_Object *
Pympc_From_PyInt_context(PyObject *self, GMPyContextObject *context)
{
    MPC_Object *result;

    if ((result = (MPC_Object*)Pympc_new_context(context)))
        result->rc = mpc_set_si(result->c, PyInt_AsLong(self),
                                GET_MPC_ROUND(context));
    return result;
}

static PyObject *
Pympc_To_PyIntOrLong(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'int'");
    return NULL;
}
#endif

/* Conversion to/from MPC
 * Python's string representation of a complex number differs from the format
 * used by MPC. Both MPC and Python surround the complex number with '(' and
 * ')' but Python adds a 'j' after the imaginary component and MPC requires a
 * space between the real and imaginery components. PyStr2Pympc tries to work
 * around the differences as follows reading two MPFR-compatible numbers from
 * the string and storing into the real and imaginary components respectively.
 */

static MPC_Object *
Pympc_From_PyStr(PyObject *s, int base, mpfr_prec_t rbits, mpfr_prec_t ibits)
{
    MPC_Object *newob;
    PyObject *ascii_str = NULL;
    Py_ssize_t len;
    char *cp, *unwind, *tempchar, *lastchar;
    int firstp = 0, lastp = 0, real_rc = 0, imag_rc = 0;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

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
        TYPE_ERROR("string required for PyStr2Pympc");
        return NULL;
    }

    if (!(newob = (MPC_Object*)Pympc_new(rbits, ibits))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Don't allow NULL characters */
    if (strlen(cp) != len) {
        VALUE_ERROR("string without NULL characters expected");
        Py_DECREF((PyObject*)newob);
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
    real_rc = mpfr_strtofr(mpc_realref(newob->c), cp, &tempchar, base,
                           GET_REAL_ROUND(context));
    /* Verify that at least one valid character was read. */
    if (cp == tempchar) goto invalid_string;
    /* If the next character is a j, then the real component is 0 and
     * we just read the imaginary componenet.
     */
    if (*tempchar == 'j') {
        mpfr_set_zero(mpc_realref(newob->c), +1);
        cp = unwind;
    }
    else {
        /* Read the imaginary component next. */
        cp = tempchar;
    }
    imag_rc = mpfr_strtofr(mpc_imagref(newob->c), cp, &tempchar, base,
                           GET_IMAG_ROUND(context));

    if (cp == tempchar && tempchar > lastchar)
        goto valid_string;

    if (*tempchar != 'j' && *cp != ' ')
        goto invalid_string;

    if (tempchar <= lastchar)
        goto invalid_string;

  valid_string:
    Py_XDECREF(ascii_str);
    newob->rc = MPC_INEX(real_rc, imag_rc);
    return newob;

  invalid_string:
    VALUE_ERROR("invalid string in mpc()");
    Py_DECREF((PyObject*)newob);
    Py_XDECREF(ascii_str);
    return NULL;
}

static MPC_Object *
Pympc_From_PyStr_bits_context(PyObject *s, int base, mpfr_prec_t rbits,
                              mpfr_prec_t ibits, GMPyContextObject *context)
{
    MPC_Object *newob;
    PyObject *ascii_str = NULL;
    Py_ssize_t len;
    char *cp, *unwind, *tempchar, *lastchar;
    int firstp = 0, lastp = 0, real_rc = 0, imag_rc = 0;

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
        TYPE_ERROR("string required for PyStr2Pympc");
        return NULL;
    }

    if (!(newob = (MPC_Object*)Pympc_new_bits_context(rbits, ibits, context))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Don't allow NULL characters */
    if (strlen(cp) != len) {
        VALUE_ERROR("string without NULL characters expected");
        Py_DECREF((PyObject*)newob);
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
    real_rc = mpfr_strtofr(mpc_realref(newob->c), cp, &tempchar, base,
                           GET_REAL_ROUND(context));
    /* Verify that at least one valid character was read. */
    if (cp == tempchar) goto invalid_string;
    /* If the next character is a j, then the real component is 0 and
     * we just read the imaginary componenet.
     */
    if (*tempchar == 'j') {
        mpfr_set_zero(mpc_realref(newob->c), +1);
        cp = unwind;
    }
    else {
        /* Read the imaginary component next. */
        cp = tempchar;
    }
    imag_rc = mpfr_strtofr(mpc_imagref(newob->c), cp, &tempchar, base,
                           GET_IMAG_ROUND(context));

    if (cp == tempchar && tempchar > lastchar)
        goto valid_string;

    if (*tempchar != 'j' && *cp != ' ')
        goto invalid_string;

    if (tempchar <= lastchar)
        goto invalid_string;

  valid_string:
    Py_XDECREF(ascii_str);
    newob->rc = MPC_INEX(real_rc, imag_rc);
    return newob;

  invalid_string:
    VALUE_ERROR("invalid string in mpc()");
    Py_DECREF((PyObject*)newob);
    Py_XDECREF(ascii_str);
    return NULL;
}

static MPC_Object *
Pympc_From_PyStr_context(PyObject *s, int base, GMPyContextObject *context)
{
    MPC_Object *newob;
    PyObject *ascii_str = NULL;
    Py_ssize_t len;
    char *cp, *unwind, *tempchar, *lastchar;
    int firstp = 0, lastp = 0, real_rc = 0, imag_rc = 0;

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
        TYPE_ERROR("string required for PyStr2Pympc");
        return NULL;
    }

    if (!(newob = (MPC_Object*)Pympc_new_context(context))) {
        Py_XDECREF(ascii_str);
        return NULL;
    }

    /* Don't allow NULL characters */
    if (strlen(cp) != len) {
        VALUE_ERROR("string without NULL characters expected");
        Py_DECREF((PyObject*)newob);
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
    real_rc = mpfr_strtofr(mpc_realref(newob->c), cp, &tempchar, base,
                           GET_REAL_ROUND(context));
    /* Verify that at least one valid character was read. */
    if (cp == tempchar) goto invalid_string;
    /* If the next character is a j, then the real component is 0 and
     * we just read the imaginary componenet.
     */
    if (*tempchar == 'j') {
        mpfr_set_zero(mpc_realref(newob->c), +1);
        cp = unwind;
    }
    else {
        /* Read the imaginary component next. */
        cp = tempchar;
    }
    imag_rc = mpfr_strtofr(mpc_imagref(newob->c), cp, &tempchar, base,
                           GET_IMAG_ROUND(context));

    if (cp == tempchar && tempchar > lastchar)
        goto valid_string;

    if (*tempchar != 'j' && *cp != ' ')
        goto invalid_string;

    if (tempchar <= lastchar)
        goto invalid_string;

  valid_string:
    Py_XDECREF(ascii_str);
    newob->rc = MPC_INEX(real_rc, imag_rc);
    return newob;

  invalid_string:
    VALUE_ERROR("invalid string in mpc()");
    Py_DECREF((PyObject*)newob);
    Py_XDECREF(ascii_str);
    return NULL;
}

static PyObject *
Pympc_To_PyStr(MPC_Object *self, int base, int digits)
{
    PyObject *tempreal = 0, *tempimag = 0, *result;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (!((base >= 2) && (base <= 62))) {
        VALUE_ERROR("base must be in the interval 2 ... 62");
        return NULL;
    }
    if ((digits < 0) || (digits == 1)) {
        VALUE_ERROR("digits must be 0 or >= 2");
        return NULL;
    }

    tempreal = raw_mpfr_ascii(mpc_realref(self->c), base, digits,
                            MPC_RND_RE(GET_MPC_ROUND(context)));
    tempimag = raw_mpfr_ascii(mpc_imagref(self->c), base, digits,
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

/*
 * If obj is a Pympc and rprec/iprec are 0/0 or the same as the precision of
 * obj, then a new reference is created.
 *
 * For all other numerical types with bits = 0, the conversion is rounded
 * according to the context.
 */

static MPC_Object *
Pympc_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    MPC_Object* newob = 0;
    MPQ_Object* temp = 0;
    mpfr_prec_t pr = 0, pi = 0;
    int rr, ri, dr, di;
    GMPyContextObject *context;

    CURRENT_CONTEXT(context);

    if (MPC_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpc is still
         * valid in the current context. */
        if (!rprec && !iprec) {
            Py_INCREF(obj);
            newob = (MPC_Object*)obj;
        }
        else {
            mpc_get_prec2(&pr, &pi, MPC(obj));
            if (rprec == pr && iprec == pi) {
                Py_INCREF(obj);
                newob = (MPC_Object*)obj;
            }
            else {
                newob = Pympc_From_Pympc((PyObject*)obj, rprec, iprec);
            }
        }
    }
    else if (MPC_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer
         * valid and mpfr_check_range needs to be called. */
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpc' incompatible with current context");
            return NULL;
        }
        /* Get the real and imaginary precisions. */
        mpc_get_prec2(&pr, &pi, MPC(obj));

        /* Get the real and imaginary inexact codes. */
        rr = MPC_INEX_RE( ((MPC_Object*)obj)->rc );
        ri = MPC_INEX_IM( ((MPC_Object*)obj)->rc );

        /* Get the real and imaginary rounding modes. */
        dr = MPC_RND_RE( ((MPC_Object*)obj)->round_mode );
        di = MPC_RND_IM( ((MPC_Object*)obj)->round_mode );

        if ((newob = (MPC_Object*)Pympc_new(pr, pi))) {
            mpc_set(newob->c, MPC(obj), GET_MPC_ROUND(context));
            newob->round_mode = ((MPC_Object*)obj)->round_mode;
            rr = mpfr_check_range(mpc_realref(newob->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(newob->c), ri, di);
            newob->rc = MPC_INEX(rr, ri);
        }
    }
    else if (MPFR_Check(obj)) {
            newob = Pympc_From_Pympfr((PyObject*)obj, rprec, iprec);
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympc_From_PyFloat(obj, rprec, iprec);
    }
    else if (PyComplex_Check(obj)) {
            newob = Pympc_From_PyComplex(obj, rprec, iprec);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympc_From_PyInt_bits_context(obj, rprec, iprec, context);
#endif
    }
    else if (MPQ_Check(obj)) {
        newob = Pympc_From_Pympq(obj, rprec, iprec);
    }
    else if (MPZ_Check(obj)) {
        newob = Pympc_From_Pympz(obj, rprec, iprec);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympc_From_PyLong(obj, rprec, iprec);
    }
    else if (XMPZ_Check(obj)) {
        newob = Pympc_From_Pyxmpz(obj, rprec, iprec);
    }
    else if (IS_DECIMAL(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = Pympc_From_PyStr(s, 10, rprec, iprec);
            if (!newob) {
                Py_DECREF(s);
                return NULL;
            }
            Py_DECREF(s);
        }
    }
    else if (IS_FRACTION(obj)) {
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympc_From_Pympq((PyObject *)temp, rprec, iprec);
            Py_DECREF((PyObject*)temp);
        }
    }
    return newob;
}

/* See the comments for GMPy_MPFR_From_Real_Temp. */

static MPC_Object *
GMPy_MPC_From_Complex_Temp(PyObject* obj, GMPyContextObject *context)
{
    MPC_Object* result = NULL;

    /* Check if obj is an mpc and the exponents are valid. */

    if (MPC_CheckAndExp(obj)) {
        /* Return a new reference with the precision of the input. */
        result = (MPC_Object*)obj;
        Py_INCREF(obj);
        return result;
    }

    /* The exponent is not valid. */

    if (MPC_Check(obj)) {
        mpfr_prec_t pr = 0, pi = 0;
        int rr, ri, dr, di;

        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpc' incompatible with current context");
            return NULL;
        }

        /* Get the real and imaginary precisions. */
        mpc_get_prec2(&pr, &pi, MPC(obj));

        /* Get the real and imaginary inexact codes. */
        rr = MPC_INEX_RE( ((MPC_Object*)obj)->rc );
        ri = MPC_INEX_IM( ((MPC_Object*)obj)->rc );

        /* Get the real and imaginary rounding modes. */
        dr = MPC_RND_RE( ((MPC_Object*)obj)->round_mode );
        di = MPC_RND_IM( ((MPC_Object*)obj)->round_mode );

        if ((result = (MPC_Object*)Pympc_new(pr, pi))) {
            mpc_set(result->c, MPC(obj), GET_MPC_ROUND(context));
            result->round_mode = ((MPC_Object*)obj)->round_mode;
            rr = mpfr_check_range(mpc_realref(result->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(result->c), ri, di);
            result->rc = MPC_INEX(rr, ri);
        }
        return result;
    }

    if (MPFR_Check(obj))
        return Pympc_From_Pympfr((PyObject*)obj, mpfr_get_prec(MPFR(obj)),
                                 mpfr_get_prec(MPFR(obj)));

    if (PyFloat_Check(obj))
        return Pympc_From_PyFloat(obj, 53, 53);

    if (PyComplex_Check(obj))
        return Pympc_From_PyComplex(obj, 53, 53);

#ifdef PY2
    if (PyInt_Check(obj))
        return Pympc_From_PyInt_bits_context(obj, 0, 0, context);
#endif

    if (MPQ_Check(obj))
        return Pympc_From_Pympq(obj, 0, 0);

    if (MPZ_Check(obj))
        return Pympc_From_Pympz(obj, 0, 0);

    if (PyLong_Check(obj))
        return Pympc_From_PyLong(obj, 0, 0);

    if (XMPZ_Check(obj))
        return Pympc_From_Pyxmpz(obj, 0, 0);

    if (IS_DECIMAL(obj)) {
        PyObject *temps = PyObject_Str(obj);

        if (temps) {
            result = Pympc_From_PyStr(temps, 10, 0, 0);
            Py_DECREF(temps);
        }
        return result;
    }

    if (IS_FRACTION(obj)) {
        MPQ_Object *tempq = Pympq_From_Fraction(obj);

        if (tempq) {
            result = Pympc_From_Pympq((PyObject *)tempq, 0, 0);
            Py_DECREF((PyObject*)tempq);
        }
        return result;
    }

    TYPE_ERROR("object could not be converted to 'mpc'");
    return NULL;
}

static MPC_Object *
Pympc_From_Complex_bits_context(PyObject* obj, mpfr_prec_t rprec,
                                mpfr_prec_t iprec, GMPyContextObject *context)
{
    MPC_Object* newob = 0;
    MPQ_Object* temp = 0;
    mpfr_prec_t pr = 0, pi = 0;
    int rr, ri, dr, di;

    if (MPC_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpc is still
         * valid in the current context. */
        if (!rprec && !iprec) {
            Py_INCREF(obj);
            newob = (MPC_Object*)obj;
        }
        else {
            mpc_get_prec2(&pr, &pi, MPC(obj));
            if (rprec == pr && iprec == pi) {
                Py_INCREF(obj);
                newob = (MPC_Object*)obj;
            }
            else {
                newob = Pympc_From_Pympc_bits_context((PyObject*)obj, rprec,
                                                      iprec, context);
            }
        }
    }
    else if (MPC_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer
         * valid and mpfr_check_range needs to be called. */
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpc' incompatible with current context");
            return NULL;
        }
        /* Get the real and imaginary precisions. */
        mpc_get_prec2(&pr, &pi, MPC(obj));

        /* Get the real and imaginary inexact codes. */
        rr = MPC_INEX_RE( ((MPC_Object*)obj)->rc );
        ri = MPC_INEX_IM( ((MPC_Object*)obj)->rc );

        /* Get the real and imaginary rounding modes. */
        dr = MPC_RND_RE( ((MPC_Object*)obj)->round_mode );
        di = MPC_RND_IM( ((MPC_Object*)obj)->round_mode );

        if ((newob = (MPC_Object*)Pympc_new_bits_context(pr, pi, context))) {
            mpc_set(newob->c, MPC(obj), GET_MPC_ROUND(context));
            newob->round_mode = ((MPC_Object*)obj)->round_mode;
            rr = mpfr_check_range(mpc_realref(newob->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(newob->c), ri, di);
            newob->rc = MPC_INEX(rr, ri);
        }
    }
    else if (MPFR_Check(obj)) {
            newob = Pympc_From_Pympfr_bits_context((PyObject*)obj, rprec, iprec,
                                                   context);
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympc_From_PyFloat_bits_context(obj, rprec, iprec, context);
    }
    else if (PyComplex_Check(obj)) {
            newob = Pympc_From_PyComplex_bits_context(obj, rprec, iprec, context);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympc_From_PyInt_bits_context(obj, rprec, iprec, context);
#endif
    }
    else if (MPQ_Check(obj)) {
        newob = Pympc_From_Pympq_bits_context(obj, rprec, iprec, context);
    }
    else if (MPZ_Check(obj)) {
        newob = Pympc_From_Pympz_bits_context(obj, rprec, iprec, context);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympc_From_PyLong_bits_context(obj, rprec, iprec, context);
    }
    else if (XMPZ_Check(obj)) {
        newob = Pympc_From_Pyxmpz_bits_context(obj, rprec, iprec, context);
    }
    else if (IS_DECIMAL(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = Pympc_From_PyStr_bits_context(s, 10, rprec, iprec, context);
            if (!newob) {
                Py_DECREF(s);
                return NULL;
            }
            Py_DECREF(s);
        }
    }
    else if (IS_FRACTION(obj)) {
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympc_From_Pympq_bits_context((PyObject *)temp, rprec,
                                                  iprec, context);
            Py_DECREF((PyObject*)temp);
        }
    }
    return newob;
}

/* Return an mpc instance from an integer/rational/real/complex number. If
 * obj is an mpc, and its precision matches the precision in context, then
 * just return a new reference to obj. Otherwise, create a new mpc. */

static MPC_Object *
Pympc_From_Complex_context(PyObject* obj, GMPyContextObject *context)
{
    MPC_Object *newob = NULL;
    MPQ_Object *temp;
    mpfr_prec_t rprec, iprec;
    int rr, ri, dr, di;

    /* Handle the likely case where the exponent of the mpc is still
     * valid in the current context. */
    if (MPC_CheckAndExp(obj)) {
        mpc_get_prec2(&rprec, &iprec, MPC(obj));
        if ((rprec == GET_REAL_PREC(context)) &&
            (iprec == GET_IMAG_PREC(context))) {
                Py_INCREF(obj);
                newob = (MPC_Object*)obj;
        }
        else {
            newob = Pympc_From_Pympc_context((PyObject*)obj, context);
        }
    }
    else if (MPC_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer
         * valid and mpfr_check_range needs to be called. */
        if (context->ctx.traps & TRAP_EXPBOUND) {
            GMPY_EXPBOUND("exponent of existing 'mpc' incompatible with current context");
            return NULL;
        }
        /* Get the real and imaginary precisions. */
        mpc_get_prec2(&rprec, &iprec, MPC(obj));

        /* Get the real and imaginary inexact codes. */
        rr = MPC_INEX_RE( ((MPC_Object*)obj)->rc );
        ri = MPC_INEX_IM( ((MPC_Object*)obj)->rc );

        /* Get the real and imaginary rounding modes. */
        dr = MPC_RND_RE( ((MPC_Object*)obj)->round_mode );
        di = MPC_RND_IM( ((MPC_Object*)obj)->round_mode );

        if ((newob = (MPC_Object*)Pympc_new_bits_context(rprec, iprec, context))) {
            mpc_set(newob->c, MPC(obj), GET_MPC_ROUND(context));
            newob->round_mode = ((MPC_Object*)obj)->round_mode;
            rr = mpfr_check_range(mpc_realref(newob->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(newob->c), ri, di);
            newob->rc = MPC_INEX(rr, ri);
        }
    }
    else if (MPFR_Check(obj)) {
            newob = Pympc_From_Pympfr_context((PyObject*)obj, context);
    }
    else if (PyFloat_Check(obj)) {
        newob = Pympc_From_PyFloat_context(obj, context);
    }
    else if (PyComplex_Check(obj)) {
            newob = Pympc_From_PyComplex_context(obj, context);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = Pympc_From_PyInt_context(obj, context);
#endif
    }
    else if (MPQ_Check(obj)) {
        newob = Pympc_From_Pympq_context(obj, context);
    }
    else if (MPZ_Check(obj)) {
        newob = Pympc_From_Pympz_context(obj, context);
    }
    else if (PyLong_Check(obj)) {
        newob = Pympc_From_PyLong_context(obj, context);
    }
    else if (XMPZ_Check(obj)) {
        newob = Pympc_From_Pyxmpz_context(obj, context);
    }
    else if (IS_DECIMAL(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = Pympc_From_PyStr_context(s, 10, context);
            Py_DECREF(s);
        }
    }
    else if (IS_FRACTION(obj)) {
        temp = Pympq_From_Fraction(obj);
        if (temp) {
            newob = Pympc_From_Pympq_context((PyObject*)temp, context);
            Py_DECREF((PyObject*)temp);
        }
    }
    return newob;
}

/*
 * coerce any number to a mpc
 */

int
Pympc_convert_arg(PyObject *arg, PyObject **ptr)
{
    MPC_Object* newob = Pympc_From_Complex(arg, 0, 0);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("can't convert argument 'mpc'");
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



