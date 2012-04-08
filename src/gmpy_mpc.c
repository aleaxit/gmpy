/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen                  *
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

/* Classify an object as a type of number. If an object is recognized as a
 * number, it must be properly converted by the routines below.
 */

static int isComplex(PyObject* obj)
{
    if (Pympz_Check(obj))       return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (Pympq_Check(obj))       return 1;
    if (Pympfr_Check(obj))      return 1;
    if (Pyxmpz_Check(obj))      return 1;
    if (Pympc_Check(obj))       return 1;
    if (PyFloat_Check(obj))     return 1;
    if (PyComplex_Check(obj))   return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

/* Verify that a valid rounding mode is specified for complex arithmetic.
 * Returns 0 (false) if the rounding mode is not valid else returns 1 (true).
 */

static PympcObject *
Pympc2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, Pympc_AS_MPC(self));
    if ((result = Pympc_new(rprec, iprec)))
        mpc_set(result->c, Pympc_AS_MPC(self), GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
PyComplex2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(self),
                    PyComplex_ImagAsDouble(self), GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
Pympfr2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (!rprec)
        rprec = mpfr_get_prec(Pympfr_AS_MPFR(self));
    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_fr(result->c, Pympfr_AS_MPFR(self),
                                GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
PyFloat2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (!rprec)
        rprec = DBL_MANT_DIG;
    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(self),
                               GET_MPC_ROUND(context));
    return result;
}

static PyObject *
Pympc2PyFloat(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'float'");
    return NULL;
}

static PympcObject *
Pympz2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_z(result->c, Pympz_AS_MPZ(self),
                                GET_MPC_ROUND(context));
    return result;
}

#define Pyxmpz2Pympc Pympz2Pympc

static PympcObject *
Pympq2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_q(result->c, Pympq_AS_MPQ(self),
                               GET_MPC_ROUND(context));
    return result;
}

static PympcObject *
PyLong2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;
    PyObject *temp = (PyObject*)PyLong2Pympz(self);

    if (!temp)
        return NULL;
    result = Pympz2Pympc(temp, rprec, iprec);
    Py_DECREF(temp);
    return result;
}

static PyObject *
Pympc2PyLong(PyObject *self)
{
    TYPE_ERROR("can't covert 'mpc' to 'long'");
    return NULL;
}

#ifdef PY2
static PympcObject *
PyInt2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_si(result->c, PyInt_AsLong(self),
                                GET_MPC_ROUND(context));
    return result;
}
static PyObject *
Pympc2PyInt(PyObject *self)
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

static PympcObject *
PyStr2Pympc(PyObject *s, long base, mpfr_prec_t rbits, mpfr_prec_t ibits)
{
    PympcObject *newob;
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

    if (!(newob = Pympc_new(rbits, ibits))) {
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
raw_mpfr_ascii(mpfr_t self, int base, int digits, int round)
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
        SYSTEM_ERROR("Internal error in raw_mpfr_ascii");
        return NULL;
    }

    result = Py_BuildValue("(sii)", buffer, the_exp, mpfr_get_prec(self));
    mpfr_free_str(buffer);
    return result;
}

static PyObject *
Pympc_ascii(PympcObject *self, int base, int digits)
{
    PyObject *tempreal = 0, *tempimag = 0;

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

    return Py_BuildValue("(NN)", tempreal, tempimag);
}

/*
 * If obj is a Pympc and rprec/iprec are 0/0 or the same as the precision of
 * obj, then a new reference is created.
 *
 * For all other numerical types with bits = 0, the conversion is rounded
 * according to the context.
 */

static PympcObject *
Pympc_From_Complex(PyObject* obj, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject* newob = 0;
    PympqObject* temp = 0;
    mpfr_prec_t pr = 0, pi = 0;
    int rr, ri, dr, di;

    if (Pympc_CheckAndExp(obj)) {
        /* Handle the likely case where the exponent of the mpc is still
         * valid in the current context. */
        if (!rprec && !iprec) {
            Py_INCREF(obj);
            newob = (PympcObject*)obj;
        }
        else {
            mpc_get_prec2(&pr, &pi, Pympc_AS_MPC(obj));
            if (rprec == pr && iprec == pi) {
                Py_INCREF(obj);
                newob = (PympcObject*)obj;
            }
            else {
                newob = Pympc2Pympc((PyObject*)obj, rprec, iprec);
            }
        }
    }
    else if (Pympc_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer
         * valid and mpfr_check_range needs to be called. */
        if (context->ctx.trap_expbound) {
            GMPY_EXPBOUND("exponent of existing 'mpc' incompatible with current context");
            return NULL;
        }
        /* Get the real and imaginary precisions. */
        mpc_get_prec2(&pr, &pi, Pympc_AS_MPC(obj));

        /* Get the real and imaginary inexact codes. */
        rr = MPC_INEX_RE( ((PympcObject*)obj)->rc );
        ri = MPC_INEX_IM( ((PympcObject*)obj)->rc );

        /* Get the real and imaginary rounding modes. */
        dr = MPC_RND_RE( ((PympcObject*)obj)->round_mode );
        di = MPC_RND_IM( ((PympcObject*)obj)->round_mode );

        if ((newob = Pympc_new(pr, pi))) {
            mpc_set(newob->c, Pympc_AS_MPC(obj), GET_MPC_ROUND(context));
            newob->round_mode = ((PympcObject*)obj)->round_mode;
            rr = mpfr_check_range(mpc_realref(newob->c), rr, dr);
            ri = mpfr_check_range(mpc_imagref(newob->c), ri, di);
            newob->rc = MPC_INEX(rr, ri);
        }
    }
    else if (Pympfr_Check(obj)) {
            newob = Pympfr2Pympc((PyObject*)obj, rprec, iprec);
    }
    else if (PyFloat_Check(obj)) {
        newob = PyFloat2Pympc(obj, rprec, iprec);
    }
    else if (PyComplex_Check(obj)) {
            newob = PyComplex2Pympc(obj, rprec, iprec);
#ifdef PY2
    }
    else if (PyInt_Check(obj)) {
        newob = PyInt2Pympc(obj, rprec, iprec);
#endif
    }
    else if (Pympq_Check(obj)) {
        newob = Pympq2Pympc(obj, rprec, iprec);
    }
    else if (Pympz_Check(obj)) {
        newob = Pympz2Pympc(obj, rprec, iprec);
    }
    else if (PyLong_Check(obj)) {
        newob = PyLong2Pympc(obj, rprec, iprec);
    }
    else if (Pyxmpz_Check(obj)) {
        newob = Pyxmpz2Pympc(obj, rprec, iprec);
    }
    else if (isDecimal(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            newob = PyStr2Pympc(s, 10, rprec, iprec);
            if (!newob) {
                Py_DECREF(s);
                return NULL;
            }
            Py_DECREF(s);
        }
    }
    else if (isFraction(obj)) {
        PyObject *s = PyObject_Str(obj);
        if (s) {
            temp = PyStr2Pympq(s, 10);
            newob = Pympq2Pympc((PyObject *)temp, rprec, iprec);
            Py_DECREF(s);
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
    PympcObject* newob = Pympc_From_Complex(arg, 0, 0);

    if (newob) {
        *ptr = (PyObject*)newob;
        return 1;
    }
    else {
        TYPE_ERROR("can't convert argument 'mpc'");
        return 0;
    }
}

PyDoc_STRVAR(doc_mpc_digits,
"c.digits(base=10, prec=0) -> ((mant, exp, prec), (mant, exp, prec))\n\n"
"Returns up to 'prec' digits in the given base. If 'prec' is 0, as many\n"
"digits that are available given c's precision are returned. 'base' must\n"
"be between 2 and 62. The result consists of 2 three-element tuples that\n"
"contain the mantissa, exponent, and number of bits of precision of the\n"
"real and imaginary components.");

/* TODO: support keyword arguments. */

static PyObject *
Pympc_digits(PyObject *self, PyObject *args)
{
    int base = 10;
    int prec = 0;
    PyObject *result;

    if (self && Pympc_Check(self)) {
        if (!PyArg_ParseTuple(args, "|ii", &base, &prec))
            return NULL;
        Py_INCREF(self);
    }
    else {
        if(!PyArg_ParseTuple(args, "O&|ii", Pympc_convert_arg, &self,
                            &base, &prec))
        return NULL;
    }
    result = Pympc_ascii((PympcObject*)self, base, prec);
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_g_mpc,
"mpc(c[, precision=0]) -> mpc\n\n"
"      Return a new 'mpc' object from an existing complex number\n"
"      (either a Python complex object or another 'mpc' object). If\n"
"      the precision is not specified, then the precision is taken\n"
"      from the current context. The rounding mode is always taken\n"
"      from the current context.\n\n"
"mpc(r[, i=0[, precision=0]]) -> mpc\n\n"
"      Return a new 'mpc' object by converting two non-complex numbers\n"
"      into the real and imaginary components of an 'mpc' object. If\n"
"      the precision is not specified, then the precision is taken from\n"
"      the current context. The rounding mode is always taken from the\n"
"      current context.\n\n"
"mpc(s[, [precision=0[, base=10]]) -> mpc\n\n"
"      Return a new 'mpc' object by converting a string s into a complex\n"
"      number. If base is omitted, then a base-10 representation is\n"
"      assumed otherwise a base between 2 and 36 can be specified. If\n"
"      the precision is not specified, then the precision is taken from\n"
"      the current context. The rounding mode is always taken from the\n"
"      current context.\n\n"
"Note: The precision can be specified either a single number that\n"
"      is used for both the real and imaginary components, or as a\n"
"      tuple that can specify different precisions for the real\n"
"      and imaginary components.");
static PyObject *
Pygmpy_mpc(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PympcObject *result = NULL;
    PympfrObject *tempreal = NULL, *tempimag = NULL;
    PyObject *arg0 = NULL, *arg1 = NULL, *prec = NULL;
    long base = 10;
    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t rbits = 0, ibits = 0;
    Py_ssize_t argc;
    static char *kwlist_c[] = {"c", "precision", NULL};
    static char *kwlist_r[] = {"r", "i", "precision", NULL};
    static char *kwlist_s[] = {"s", "precision", "base", NULL};

    argc = PyTuple_Size(args);
    if (argc < 1) {
        TYPE_ERROR("mpc() requires at least 1 non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GetItem(args, 0);
    if (PyStrOrUnicode_Check(arg0)) {
        /* First argument is a string */
        if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|Oi", kwlist_s,
                                          &arg0, &prec, &base)))
            return NULL;

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rbits = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                ibits = rbits;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rbits = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GetItem(prec, 0));
                ibits = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GetItem(prec, 1));
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in gmpy2.mpc().");
                    return NULL;
                }
            }
        }

        if (base < 2 || base > 36) {
            VALUE_ERROR("base for mpc() must be in the interval 2 ... 36.");
            return NULL;
        }

        result = PyStr2Pympc(arg0, base, rbits, ibits);
    }
    else if (PyComplex_Check(arg0) || Pympc_Check(arg0)) {
        /* First argument is a complex number */
        if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist_c,
                                          &arg0, &prec)))
            return NULL;

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rbits = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                ibits = rbits;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rbits = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GetItem(prec, 0));
                ibits = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GetItem(prec, 1));
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc().");
                    return NULL;
                }
            }
        }

        if (PyComplex_Check(arg0)) {
            result = PyComplex2Pympc(arg0, rbits, ibits);
        }
        else {
            result = Pympc2Pympc(arg0, rbits, ibits);
        }
    }
    else if (isReal(arg0)) {
        /* First argument is a real number */
        if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|OO", kwlist_r,
                                          &arg0, &arg1, &prec)))
            return NULL;

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rbits = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                ibits = rbits;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rbits = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GetItem(prec, 0));
                ibits = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GetItem(prec, 1));
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc().");
                    return NULL;
                }
            }
        }

        if (arg0) tempreal = Pympfr_From_Real(arg0, rbits);
        if (arg1) tempimag = Pympfr_From_Real(arg1, ibits);

        if (!tempreal) {
            if ((tempreal = Pympfr_new(rbits)))
                mpfr_set_ui(Pympfr_AS_MPFR(tempreal), 0, context->ctx.mpfr_round);
        }

        if (!tempimag) {
            if ((tempimag = Pympfr_new(ibits)))
                mpfr_set_ui(Pympfr_AS_MPFR(tempimag), 0, context->ctx.mpfr_round);
        }

        result = Pympc_new(rbits, ibits);
        if (!tempreal || !tempimag || !result) {
            Py_XDECREF(tempreal);
            Py_XDECREF(tempimag);
            Py_XDECREF(result);
            TYPE_ERROR("mpc() require string or numeric argument.");
            return NULL;
        }

        mpc_set_fr_fr(Pympc_AS_MPC(result), Pympfr_AS_MPFR(tempreal),
                      Pympfr_AS_MPFR(tempimag), GET_MPC_ROUND(context));
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
    }
    else {
        TYPE_ERROR("mpc() requires numeric or string argument");
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpc_format,
"x.__format__(fmt) -> string\n\n"
"Return a Python string by formatting 'x' using the format string\n"
"'fmt'. A valid format string consists of:\n"
"     optional alignment code:\n"
"        '<' -> left shifted in field\n"
"        '>' -> right shifted in field\n"
"        '^' -> centered in field\n"
"     optional leading sign code\n"
"        '+' -> always display leading sign\n"
"        '-' -> only display minus for negative values\n"
"        ' ' -> minus for negative values, space for positive values\n"
"     optional width.real_precision.imag_precision\n"
"     optional rounding mode:\n"
"        'U' -> round toward plus infinity\n"
"        'D' -> round toward minus infinity\n"
"        'Z' -> round toward zero\n"
"        'N' -> round to nearest\n"
"     optional output style:\n"
"        'P' -> Python style, 1+2j, (default)\n"
"        'M' -> MPC style, (1 2)\n"
"     optional conversion code:\n"
"        'a','A' -> hex format\n"
"        'b'     -> binary format\n"
"        'e','E' -> scientific format\n"
"        'f','F' -> fixed point format\n"
"        'g','G' -> fixed or scientific format\n\n"
"The default format is 'f'.");

static PyObject *
Pympc_format(PyObject *self, PyObject *args)
{
    PyObject *result = 0, *tempstr = 0;
    char *realbuf = 0, *imagbuf = 0, *tempbuf = 0, *fmtcode = 0;
    char *p, *rfmtptr, *ifmtptr, *fmtptr;
    char rfmt[100], ifmt[100], fmt[30];
    int rbuflen, ibuflen;
    int seensign = 0, seenalign = 0, seendecimal = 0, seendigits = 0;
    int seenround = 0, seenconv = 0, seenstyle = 0, mpcstyle = 0;

    if (!Pympc_Check(self)) {
        TYPE_ERROR("requires 'mpc' object");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "s", &fmtcode))
        return NULL;

    rfmtptr = rfmt;
    ifmtptr = ifmt;
    fmtptr = fmt;
    *(rfmtptr++) = '%';
    *(ifmtptr++) = '%';

    for (p = fmtcode; *p != '\00'; p++) {
        if (*p == '<' || *p == '>' || *p == '^') {
            if (seenalign || seensign || seendecimal || seendigits ||
                seenround || seenstyle) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(fmtptr++) = *p;
                seenalign = 1;
                continue;
            }
        }
        if (*p == '+' || *p == ' ' || *p == '-') {
            if (seensign || seendecimal || seendigits || seenround ||
                seenstyle) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(rfmtptr++) = *p;
                *(ifmtptr++) = *p;
                seensign = 1;
                continue;
            }
        }
        if (!seensign) {
            *(rfmtptr++) = '-';
            *(ifmtptr++) = '-';
            seensign = 1;
        }
        if (*p == '.') {
            if (seendecimal == 2 || seendigits || seenround || seenstyle) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                if (!seendecimal) {
                    *(rfmtptr++) = *p;
                    *(ifmtptr++) = *p;
                }
                seendecimal++;
                if (seendecimal == 2) {
                    while (isdigit(*(ifmtptr-1)))
                        ifmtptr--;
                }
                continue;
            }
        }
        if (isdigit(*p)) {
            if (seendigits || seenround || seenstyle) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else if (seendecimal == 1) {
                *(rfmtptr++) = *p;
                *(ifmtptr++) = *p;
                continue;
            }
            else if (seendecimal == 2) {
                *(ifmtptr++) = *p;
                continue;
            }
            else {
                if (fmtptr == fmt) {
                    *(fmtptr++) = '>';
                    seenalign = 1;
                }
                *(fmtptr++) = *p;
                continue;
            }
        }
        if (!seendigits) {
            seendigits = 1;
            *(rfmtptr++) = 'R';
            *(ifmtptr++) = 'R';
        }
        if (*p == 'U' || *p == 'D' || *p == 'Y' || *p == 'Z' ||
            *p == 'N' ) {
            if (seenround || seenstyle) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(rfmtptr++) = *p;
                *(ifmtptr++) = *p;
                seenround = 1;
                continue;
            }
        }
        if (*p == 'P' || *p == 'M') {
            if (seenstyle) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                if (*p == 'M')
                    mpcstyle = 1;
                seenstyle = 1;
                continue;
            }
        }
        if (*p == 'a' || *p == 'A' || *p == 'b' || *p == 'e' ||
            *p == 'E' || *p == 'f' || *p == 'F' || *p == 'g' ||
            *p == 'G' ) {
            *(rfmtptr++) = *p;
            *(ifmtptr++) = *p;
            seenconv = 1;
            break;
        }
        VALUE_ERROR("Invalid conversion specification");
        return NULL;
    }

    if (!seensign) {
        *(rfmtptr++) = '-';
        *(ifmtptr++) = '-';
    }
    if (!seendigits) {
        *(rfmtptr++) = 'R';
        *(ifmtptr++) = 'R';
    }
    if (!seenconv) {
        *(rfmtptr++) = 'f';
        *(ifmtptr++) = 'f';
    }

    *(rfmtptr) = '\00';
    *(ifmtptr) = '\00';
    *(fmtptr) = '\00';

    /* Format the real part.... */

    rbuflen = mpfr_asprintf(&realbuf, rfmt,
                           mpc_realref(Pympc_AS_MPC(self)));

    if (rbuflen < 0) {
        mpfr_free_str(realbuf);
        SYSTEM_ERROR("Internal error in mpfr_asprintf");
        return NULL;
    }

    /* Format the imaginary part. If Python style is wanted, convert the '-'
     * or ' ' sign indicator to '+'. */

    if (!mpcstyle) {
        if (ifmt[1] == ' ' || ifmt[1] == '-' || ifmt[1] == '+') {
            ifmt[1] = '+';
        }
        else {
            mpfr_free_str(realbuf);
            VALUE_ERROR("Invalid conversion specification for imag");
            return NULL;
        }
    }

    ibuflen = mpfr_asprintf(&imagbuf, ifmt,
                           mpc_imagref(Pympc_AS_MPC(self)));

    if (ibuflen < 0) {
        mpfr_free_str(realbuf);
        mpfr_free_str(imagbuf);
        SYSTEM_ERROR("Internal error in mpfr_asprintf");
        return NULL;
    }

    /* Combine the real and imaginary components into a single buffer.
     * Include space for '(', ' ', and 'j)' and possibly appending '.0' twice.
     */

    tempbuf = GMPY_MALLOC(rbuflen + ibuflen + 10);
    if (!tempbuf) {
        mpfr_free_str(realbuf);
        mpfr_free_str(imagbuf);
        return PyErr_NoMemory();
    }
    tempbuf[0] = '\00';
    if (mpcstyle)
        strcat(tempbuf, "(");
    strcat(tempbuf, realbuf);

    /* If there isn't a decimal point in the output and the output
     * is short and only consists of digits, then append .0 */
    if (strlen(realbuf) < 50 &&
        strlen(realbuf) == strspn(realbuf, "+- 0123456789")) {
        strcat(tempbuf, ".0");
    }

    if (mpcstyle)
        strcat(tempbuf, " ");
    else {
        /* Need to insert + if imag is nan or +inf. */
        if (mpfr_nan_p(mpc_imagref(Pympc_AS_MPC(self))) ||
            (mpfr_inf_p(mpc_imagref(Pympc_AS_MPC(self))) &&
             mpfr_sgn(mpc_imagref(Pympc_AS_MPC(self))) > 0)) {
            strcat(tempbuf, "+");
        }
    }
    strcat(tempbuf, imagbuf);
    if (strlen(imagbuf) < 50 &&
        strlen(imagbuf) == strspn(imagbuf, "+- 0123456789")) {
        strcat(tempbuf, ".0");
    }

    if (mpcstyle)
        strcat(tempbuf, ")");
    else
        strcat(tempbuf, "j");

    mpfr_free_str(realbuf);
    mpfr_free_str(imagbuf);

    tempstr = Py_BuildValue("s", tempbuf);
    if (!tempstr) {
        GMPY_FREE(tempbuf);
        return NULL;
    }

    result = PyObject_CallMethod(tempstr, "__format__", "(s)", fmt);

    Py_DECREF(tempstr);
    return result;
}

/* str and repr implementations for mpc */
static PyObject *
Pympc2str(PympcObject *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[30];

    mpc_get_prec2(&rbits, &ibits, Pympc_AS_MPC(self));
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
Pympc2repr(PympcObject *self)
{
    PyObject *result, *temp;
    mpfr_prec_t rbits, ibits;
    long rprec, iprec;
    char fmtstr[30];

    mpc_get_prec2(&rbits, &ibits, Pympc_AS_MPC(self));
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

static PyObject *
Pympc_abs(PyObject *self)
{
    PympfrObject *result = 0;
    PympcObject *tempx = 0;

    result = Pympfr_new(0);
    tempx = Pympc_From_Complex(self, 0, 0);
    if (!tempx || !result) {
        SYSTEM_ERROR("Can't convert argument to 'mpc'.");
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    result->rc = mpc_abs(result->f, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    MPFR_SUBNORMALIZE(result);
    MPFR_CHECK_INVALID(result, "invalid operation in 'mpc' __abs__");
    MPFR_CHECK_UNDERFLOW(result, "underflow in 'mpc' __abs__");
    MPFR_CHECK_OVERFLOW(result, "overflow in 'mpc' __abs__");
    MPFR_CHECK_INEXACT(result, "inexact result in 'mpc' __abs__");
  done:
    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

static PyObject *
Pympc_neg(PympcObject *self)
{
    PympcObject *result = 0;

    if (!(result = Pympc_new(0, 0)))
        return NULL;

    if (!(self = Pympc_From_Complex((PyObject*)self, 0, 0))) {
        SYSTEM_ERROR("__neg__() requires 'mpc' argument");
        Py_DECREF(result);
        return NULL;
    }

    result->rc = mpc_neg(result->c, self->c, GET_MPC_ROUND(context));

    MPC_CLEANUP(result, "__neg__");
}

static PyObject *
Pympc_pos(PympcObject *self)
{
    PympcObject *result = 0;

    if (!(result = Pympc_From_Complex((PyObject*)self, 0, 0))) {
        SYSTEM_ERROR("__pos__ requires 'mpc' argument");
        return NULL;
    }

    MPC_CLEANUP(result, "__pos__");
}

/* Support Pympany_square */

static PyObject *
Pympc_sqr(PyObject* self, PyObject *other)
{
    PympcObject *result;

    PARSE_ONE_MPC_OTHER("square() requires 'mpc' argument");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_sqr(result->c, Pympc_AS_MPC(self),
                         GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "square()");
}

static PyObject *
Pympc_pow(PyObject *base, PyObject *exp, PyObject *m)
{
    PympcObject *tempb, *tempe, *result;

    if (m != Py_None) {
        TYPE_ERROR("pow() 3rd argument not allowed unless all arguments are integers");
        return NULL;
    }

    tempb = Pympc_From_Complex(base, 0, 0);
    tempe = Pympc_From_Complex(exp, 0, 0);

    if (!tempe || !tempb) {
        Py_XDECREF((PyObject*)tempe);
        Py_XDECREF((PyObject*)tempb);
        Py_RETURN_NOTIMPLEMENTED;
    }

    result = Pympc_new(0, 0);

    if (!result) {
        Py_DECREF((PyObject*)tempe);
        Py_DECREF((PyObject*)tempb);
        return NULL;
    }

    if (MPC_IS_ZERO_P(tempb) && MPC_IS_ZERO_P(tempe)) {
        mpc_set_ui(result->c, 1, GET_MPC_ROUND(context));
        Py_DECREF((PyObject*)tempe);
        Py_DECREF((PyObject*)tempb);
        return (PyObject*)result;
    }

    if (MPC_IS_ZERO_P(tempb) &&
        (!mpfr_zero_p(mpc_imagref(tempe->c)) ||
         mpfr_sgn(mpc_realref(tempe->c)) < 0)) {

        context->ctx.divzero = 1;
        if (context->ctx.trap_divzero) {
            GMPY_DIVZERO("zero cannot be raised to a negative or complex power");
            Py_DECREF((PyObject*)tempe);
            Py_DECREF((PyObject*)tempb);
            Py_DECREF((PyObject*)result);
            return NULL;
        }
    }

    result->rc = mpc_pow(result->c, tempb->c,
                         tempe->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempe);
    Py_DECREF((PyObject*)tempb);

    MPC_CLEANUP(result, "pow()");
}

/* Implement the conjugate() method. */

PyDoc_STRVAR(doc_mpc_conjugate,
"x.conjugate() -> mpc\n\n"
"Returns the conjugate of x.");

static PyObject *
Pympc_conjugate(PyObject *self, PyObject *args)
{
    PympcObject *result;

    PARSE_ONE_MPC_ARGS("conjugate() requires 'mpc' argument");

    if (!(result = Pympc_new(0,0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_conj(result->c, Pympc_AS_MPC(self),
                          GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "conjugate()");
}

/* Implement the .precision attribute of an mpfr. */

static PyObject *
Pympc_getprec_attrib(PympcObject *self, void *closure)
{
    mpfr_prec_t rprec = 0, iprec = 0;

    mpc_get_prec2(&rprec, &iprec, self->c);
    return Py_BuildValue("(nn)", rprec, iprec);
}

/* Implement the .rc attribute of an mpfr. */

static PyObject *
Pympc_getrc_attrib(PympcObject *self, void *closure)
{
    return Py_BuildValue("(ii)", MPC_INEX_RE(self->rc), MPC_INEX_IM(self->rc));
}

/* Implement the .imag attribute of an mpfr. */

static PyObject *
Pympc_getimag_attrib(PympcObject *self, void *closure)
{
    PympfrObject *result;

    if ((result = Pympfr_new(0)))
        mpc_imag(result->f, self->c, context->ctx.mpfr_round);
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpfr. */

static PyObject *
Pympc_getreal_attrib(PympcObject *self, void *closure)
{
    PympfrObject *result;

    if ((result = Pympfr_new(0)))
        mpc_real(result->f, self->c, context->ctx.mpfr_round);
    return (PyObject*)result;
}

/* Implement the nb_bool slot. */

static int
Pympc_nonzero(PympcObject *self)
{
    return !MPC_IS_ZERO_P(self->c);
}

/* To work with the MPC_IS_ macros, NAN, INF, and ZERO are all upper-case. */

#define MPC_TEST_OTHER(NAME, msg) \
static PyObject * \
Pympc_is_##NAME(PyObject *self, PyObject *other)\
{\
    int res;\
    if(self && Pympc_Check(self)) {\
        Py_INCREF(self);\
    }\
    else if(Pympc_Check(other)) {\
        self = other;\
        Py_INCREF((PyObject*)self);\
    }\
    else if (!(self = (PyObject*)Pympc_From_Complex(other, 0, 0))) {\
        PyErr_SetString(PyExc_TypeError, msg);\
        return NULL;\
    }\
    res = MPC_IS_##NAME##_P(Pympc_AS_MPC(self));\
    Py_DECREF(self);\
    if (res)\
        Py_RETURN_TRUE;\
    else\
        Py_RETURN_FALSE;\
}

MPC_TEST_OTHER(NAN, "is_nan() requires 'mpc' argument");

MPC_TEST_OTHER(INF, "is_inf() requires 'mpc' argument");

MPC_TEST_OTHER(ZERO, "is_zero() requires 'mpc' argument");

PyDoc_STRVAR(doc_mpc_phase,
"phase(x) -> mpfr\n\n"
"Return the phase angle, also known as argument, of a complex x.");

static PyObject *
Pympc_phase(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPC_OTHER("phase() requires 'mpc' argument");

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_arg(result->f, Pympc_AS_MPC(self),
                         context->ctx.mpfr_round);
    Py_DECREF((PyObject*)self);

    MPFR_SUBNORMALIZE(result);
    MPFR_CHECK_OVERFLOW(result, "overflow in 'mpc' phase()");
    MPFR_CHECK_INVALID(result, "invalid operation 'mpc' phase()");
    MPFR_CHECK_UNDERFLOW(result, "underflow in 'mpc' phase()");
    MPFR_CHECK_INEXACT(result, "inexact operation in 'mpc' phase()");
  done:
    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpc_norm,
"norm(x) -> mpfr\n\n"
"Return the norm of a complex x. The norm(x) is defined as\n"
"x.real**2 + x.imag**2. abs(x) is the square root of norm(x).\n");

static PyObject *
Pympc_norm(PyObject *self, PyObject *other)
{
    PympfrObject *result;

    PARSE_ONE_MPC_OTHER("norm() requires 'mpc' argument");

    if (!(result = Pympfr_new(0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_norm(result->f, Pympc_AS_MPC(self),
                          context->ctx.mpfr_round);
    Py_DECREF((PyObject*)self);

    MPFR_SUBNORMALIZE(result);
    MPFR_CHECK_OVERFLOW(result, "overflow in 'mpc' norm()");
    MPFR_CHECK_INVALID(result, "invalid operation 'mpc' norm()");
    MPFR_CHECK_UNDERFLOW(result, "underflow in 'mpc' norm()");
    MPFR_CHECK_INEXACT(result, "inexact operation in 'mpc' norm()");
  done:
    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpc_polar,
"polar(x) -> (abs(x), phase(x))\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

static PyObject *
Pympc_polar(PyObject *self, PyObject *other)
{
    PyObject *abs, *phase;

    PARSE_ONE_MPC_OTHER("norm() requires 'mpc' argument");

    if (!(abs = Pympc_abs(self))) {
        Py_DECREF(self);
        return NULL;
    }
    if (!(phase = Pympc_phase(self, other))) {
        Py_DECREF(abs);
        Py_DECREF(self);
        return NULL;
    }

    return Py_BuildValue("(NN)", abs, phase);
}

PyDoc_STRVAR(doc_mpc_rect,
"rect(x) -> mpc\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

/* Note: does not properly check for inexact or underflow */

static PyObject *
Pympc_rect(PyObject *self, PyObject *args)
{
    PyObject *other;
    PympcObject *result;

    PARSE_TWO_MPFR_ARGS(other, "rect() requires 'mpfr','mpfr' arguments");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_cos(mpc_realref(result->c), Pympfr_AS_MPFR(other),
             GET_REAL_ROUND(context));
    mpfr_mul(mpc_realref(result->c), mpc_realref(result->c),
             Pympfr_AS_MPFR(self), GET_REAL_ROUND(context));
    mpfr_sin(mpc_imagref(result->c), Pympfr_AS_MPFR(other),
             GET_IMAG_ROUND(context));
    mpfr_mul(mpc_imagref(result->c), mpc_imagref(result->c),
             Pympfr_AS_MPFR(self), GET_IMAG_ROUND(context));
    Py_DECREF(self);
    Py_DECREF(other);

    MPC_CLEANUP(result, "rect()");
}

PyDoc_STRVAR(doc_mpc_proj,
"proj(x) -> mpc\n\n"
"Returns the projection of a complex x on to the Riemann sphere.");

static PyObject *
Pympc_proj(PyObject *self, PyObject *other)
{
    PympcObject *result;

    PARSE_ONE_MPC_OTHER("proj() requires 'mpc' argument");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_proj(result->c, Pympc_AS_MPC(self),
                          GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "proj()");
}

#define MPC_UNIOP(NAME) \
static PyObject * \
Pympc_##NAME(PyObject* self, PyObject *other) \
{ \
    PympcObject *result; \
    PARSE_ONE_MPC_OTHER(#NAME "() requires 'mpc' argument"); \
    if (!(result = Pympc_new(0, 0))) { \
        Py_DECREF(self); \
        return NULL; \
    } \
    result->rc = mpc_##NAME(result->c, Pympc_AS_MPC(self), GET_MPC_ROUND(context)); \
    Py_DECREF(self); \
    MPC_CLEANUP(result, #NAME"()"); \
}

MPC_UNIOP(log)

MPC_UNIOP(exp)

MPC_UNIOP(sin)

MPC_UNIOP(cos)

MPC_UNIOP(tan)

MPC_UNIOP(sinh)

MPC_UNIOP(cosh)

MPC_UNIOP(tanh)

MPC_UNIOP(asin)

MPC_UNIOP(acos)

MPC_UNIOP(atan)

MPC_UNIOP(asinh)

MPC_UNIOP(acosh)

MPC_UNIOP(atanh)

MPC_UNIOP(sqrt)

static PyObject *
Pympc_sin_cos(PyObject *self, PyObject *other)
{
    PympcObject *s, *c;
    PyObject *result;
    int code;

    PARSE_ONE_MPC_OTHER("sin_cos() requires 'mpc' argument");

    s = Pympc_new(0, 0);
    c = Pympc_new(0, 0);
    result = PyTuple_New(2);
    if (!s || !c || !result) {
        Py_DECREF(self);
        return NULL;
    }

    code = mpc_sin_cos(s->c, c->c, Pympc_AS_MPC(self),
                       GET_MPC_ROUND(context), GET_MPC_ROUND(context));
    s->rc = MPC_INEX1(code);
    c->rc = MPC_INEX2(code);
    MPC_SUBNORMALIZE(s);
    MPC_SUBNORMALIZE(c);
    MPC_CHECK_FLAGS(s, "sin_cos()");
    MPC_CHECK_FLAGS(c, "sin_cos()");

  done:
    Py_DECREF(self);
    if (PyErr_Occurred()) {
        Py_XDECREF((PyObject*)s);
        Py_XDECREF((PyObject*)c);
        Py_XDECREF(result);
        result = NULL;
    }
    else {
        PyTuple_SET_ITEM(result, 0, (PyObject*)s);
        PyTuple_SET_ITEM(result, 1, (PyObject*)c);
    }
    return result;
}

static PyObject *
Pympc_fma(PyObject *self, PyObject *args)
{
    PympcObject *result, *x, *y, *z;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fma() requires 'mpc','mpc','mpc' arguments.");
        return NULL;
    }

    result = Pympc_new(0, 0);
    x = Pympc_From_Complex(PyTuple_GET_ITEM(args, 0), 0, 0);
    y = Pympc_From_Complex(PyTuple_GET_ITEM(args, 1), 0, 0);
    z = Pympc_From_Complex(PyTuple_GET_ITEM(args, 2), 0, 0);
    if (!result || !x || !y || !z) {
        TYPE_ERROR("fma() requires 'mpc','mpc','mpc' arguments.");
        goto done;
    }

    result->rc = mpc_fma(result->c, x->c, y->c, z->c,
                         context->ctx.mpfr_round);
    MPC_SUBNORMALIZE(result);
    MPC_CHECK_FLAGS(result, "fma()");

  done:
    Py_XDECREF((PyObject*)x);
    Py_XDECREF((PyObject*)y);
    Py_XDECREF((PyObject*)z);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        result = NULL;
    }
    return (PyObject*)result;
}

static PyObject *
Pympc_fms(PyObject *self, PyObject *args)
{
    PympcObject *result, *x, *y, *z;

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fms() requires 'mpc','mpc','mpc' arguments.");
        return NULL;
    }

    result = Pympc_new(0, 0);
    x = Pympc_From_Complex(PyTuple_GET_ITEM(args, 0), 0, 0);
    y = Pympc_From_Complex(PyTuple_GET_ITEM(args, 1), 0, 0);
    z = Pympc_From_Complex(PyTuple_GET_ITEM(args, 2), 0, 0);
    if (!result || !x || !y || !z) {
        TYPE_ERROR("fms() requires 'mpc','mpc','mpc' arguments.");
        goto done;
    }

    mpc_neg(z->c, z->c, GET_MPC_ROUND(context));
    result->rc = mpc_fma(result->c, x->c, y->c, z->c,
                         context->ctx.mpfr_round);
    MPC_SUBNORMALIZE(result);
    MPC_CHECK_FLAGS(result, "fms()");

  done:
    Py_XDECREF((PyObject*)x);
    Py_XDECREF((PyObject*)y);
    Py_XDECREF((PyObject*)z);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        result = NULL;
    }
    return (PyObject*)result;
}

static PyObject *
Pympc_div_2exp(PyObject *self, PyObject *args)
{
    PympcObject *result = 0;
    mp_bitcnt_t exp = 0;

    if (!PyArg_ParseTuple(args, "O&k", Pympc_convert_arg, &self, &exp)) {
        TYPE_ERROR("div_2exp() requires 'mpc', 'int' arguments");
        return NULL;
    }

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_div_2exp(Pympc_AS_MPC(result), Pympc_AS_MPC(self),
                              exp, GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "div_2exp()");
}

static PyObject *
Pympc_mul_2exp(PyObject *self, PyObject *args)
{
    PympcObject *result = 0;
    mp_bitcnt_t exp = 0;

    if (!PyArg_ParseTuple(args, "O&k", Pympc_convert_arg, &self, &exp)) {
        TYPE_ERROR("mul_2exp() requires 'mpc', 'int' arguments");
        return NULL;
    }

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_mul_2exp(Pympc_AS_MPC(result), Pympc_AS_MPC(self),
                             exp, GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "mul_2exp()");
}

static Py_hash_t
Pympc_hash(PympcObject *self)
{
    Py_uhash_t hashreal, hashimag, combined;

    if (self->hash_cache != -1)
        return self->hash_cache;

    hashreal = (Py_uhash_t)_mpfr_hash(mpc_realref(self->c));
    if (hashreal == (Py_uhash_t)-1)
        return -1;
    hashimag = (Py_uhash_t)_mpfr_hash(mpc_imagref(self->c));
    if (hashimag == (Py_uhash_t)-1)
        return -1;
    combined = hashreal + _PyHASH_IMAG * hashimag;
    if (combined == (Py_uhash_t)-1)
        combined = (Py_uhash_t)-2;
    self->hash_cache = combined;
    return (Py_hash_t)combined;
}

static PyObject *
Pympc_add(PyObject *self, PyObject *args)
{
    PympcObject *result;
    PyObject *other;

    PARSE_TWO_MPC_ARGS(other, "add() requires 'mpc','mpc' arguments");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    result->rc = mpc_add(result->c, Pympc_AS_MPC(self),
                         Pympc_AS_MPC(other), GET_MPC_ROUND(context));
    Py_DECREF(self);
    Py_DECREF(other);
    MPC_CLEANUP(result, "add()");
}

static PyObject *
Pympc_sub(PyObject *self, PyObject *args)
{
    PympcObject *result;
    PyObject *other;

    PARSE_TWO_MPC_ARGS(other, "sub() requires 'mpc','mpc' arguments");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    result->rc = mpc_sub(result->c, Pympc_AS_MPC(self),
                         Pympc_AS_MPC(other), GET_MPC_ROUND(context));
    Py_DECREF(self);
    Py_DECREF(other);
    MPC_CLEANUP(result, "sub()");
}

static PyObject *
Pympc_mul(PyObject *self, PyObject *args)
{
    PympcObject *result;
    PyObject *other;

    PARSE_TWO_MPC_ARGS(other, "mul() requires 'mpc','mpc' arguments");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    result->rc = mpc_mul(result->c, Pympc_AS_MPC(self),
                         Pympc_AS_MPC(other), GET_MPC_ROUND(context));
    Py_DECREF(self);
    Py_DECREF(other);
    MPC_CLEANUP(result, "mul()");
}

static PyObject *
Pympc_div(PyObject *self, PyObject *args)
{
    PympcObject *result;
    PyObject *other;

    PARSE_TWO_MPC_ARGS(other, "div() requires 'mpc','mpc' arguments");

    if (!(result = Pympc_new(0, 0))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    if (MPC_IS_ZERO_P(Pympc_AS_MPC(other))) {
        context->ctx.divzero = 1;
        if (context->ctx.trap_divzero) {
            GMPY_DIVZERO("'mpc' division by zero");
            Py_DECREF(self);
            Py_DECREF(other);
            return NULL;
        }
    }

    result->rc = mpc_div(result->c, Pympc_AS_MPC(self),
                         Pympc_AS_MPC(other), GET_MPC_ROUND(context));
    Py_DECREF(self);
    Py_DECREF(other);
    MPC_CLEANUP(result, "div()");
}

static PyMethodDef Pympc_methods[] =
{
    { "__format__", Pympc_format, METH_VARARGS, doc_mpc_format },
    { "conjugate", Pympc_conjugate, METH_NOARGS, doc_mpc_conjugate },
    { "digits", Pympc_digits, METH_VARARGS, doc_mpc_digits },
    { NULL, NULL, 1 }
};


#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympc_neg,               /* nb_negative             */
    (unaryfunc) Pympc_pos,               /* nb_positive             */
    (unaryfunc) Pympc_abs,               /* nb_absolute             */
    (inquiry) Pympc_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) Pympc2PyLong,            /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) Pympc2PyFloat,           /* nb_float                */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) Pybasic_add,            /* nb_add                  */
    (binaryfunc) Pybasic_sub,            /* nb_subtract             */
    (binaryfunc) Pybasic_mul,            /* nb_multiply             */
    (binaryfunc) Pybasic_div2,           /* nb_divide               */
    (binaryfunc) Pybasic_rem,            /* nb_remainder            */
    (binaryfunc) Pybasic_divmod,         /* nb_divmod               */
    (ternaryfunc) Pympany_pow,           /* nb_power                */
    (unaryfunc) Pympc_neg,               /* nb_negative             */
    (unaryfunc) Pympc_pos,               /* nb_positive             */
    (unaryfunc) Pympc_abs,               /* nb_absolute             */
    (inquiry) Pympc_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) Pympc2PyInt,             /* nb_int                  */
    (unaryfunc) Pympc2PyLong,            /* nb_long                 */
    (unaryfunc) Pympc2PyFloat,           /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    (binaryfunc) Pybasic_floordiv,       /* nb_floor_divide         */
    (binaryfunc) Pybasic_truediv,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympc_getseters[] =
{
    {"precision", (getter)Pympc_getprec_attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)Pympc_getrc_attrib, NULL, "return code", NULL},
    {"imag", (getter)Pympc_getimag_attrib, NULL, "imaginary component", NULL},
    {"real", (getter)Pympc_getreal_attrib, NULL, "real component", NULL},
    {NULL}
};

static PyTypeObject Pympc_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpc",                                  /* tp_name          */
    sizeof(PympcObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympc_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pympc2repr,                  /* tp_repr          */
    &mpc_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympc_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) Pympc2str,                   /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "MPC-based complex number",             /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&mpany_richcompare,        /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympc_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympc_getseters,                        /* tp_getset        */
};


