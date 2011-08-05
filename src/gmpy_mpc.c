/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.c                                                              *
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

/* Verify that a valid rounding mode is specified for complex arithmetic.
 * Returns 0 (false) if the rounding mode is not valid else returns 1 (true).
 */

static int
Pymisc_verify_mpc_round(int rmode)
{
    if ( rmode == MPC_RNDNN || rmode == MPC_RNDNZ ||
         rmode == MPC_RNDNU || rmode == MPC_RNDND ||
         rmode == MPC_RNDZN || rmode == MPC_RNDZZ ||
         rmode == MPC_RNDZU || rmode == MPC_RNDZD ||
         rmode == MPC_RNDUN || rmode == MPC_RNDUZ ||
         rmode == MPC_RNDUU || rmode == MPC_RNDUD ||
         rmode == MPC_RNDDN || rmode == MPC_RNDDZ ||
         rmode == MPC_RNDDU || rmode == MPC_RNDDD )
        return 1;
    else
        return 0;
}

/* Verify that valid precisions are requested for complex arithmetic.
 * Returns 0 if the precisions are not valid else returns 1.
 */

static int
Pymisc_verify_mpc_precision(Py_ssize_t rprec, Py_ssize_t iprec)
{
    if ( rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
         iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX )
        return 0;
    else
        return 1;
}

static PympcObject *
Pympc2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, Pympc_AS_MPC(self));
    if ((result = Pympc_new(rprec, iprec)))
        mpc_set(result->c, Pympc_AS_MPC(self), context->now.mpc_round);
    return result;
}

static PympcObject *
PyComplex2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        mpc_set_d_d(result->c, PyComplex_RealAsDouble(self),
                    PyComplex_ImagAsDouble(self), context->now.mpc_round);
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
                                context->now.mpc_round);
    return result;
}

static PympfrObject *
Pympc2Pympfr(PyObject *self)
{
    TYPE_ERROR("can not covert 'mpc' to 'mpfr'");
    return NULL;
}

static PympcObject *
PyFloat2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if (!rprec)
        rprec = DBL_MANT_DIG;
    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_d(result->c, PyFloat_AS_DOUBLE(self),
                               context->now.mpc_round);
    return result;
}

static PyObject *
Pympc2PyFloat(PyObject *self)
{
    TYPE_ERROR("can not covert 'mpc' to 'float'");
    return NULL;
}

static PympcObject *
Pympz2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_z(result->c, Pympz_AS_MPZ(self),
                                context->now.mpc_round);
    return result;
}

static PyObject *
Pympc2Pympz(PyObject *self)
{
    TYPE_ERROR("can not covert 'mpc' to 'mpz'");
    return NULL;
}

#define Pyxmpz2Pympc Pympz2Pympc

static PyObject *
Pympc2Pyxmpz(PyObject *self)
{
    TYPE_ERROR("can not covert 'mpc' to 'xmpz'");
    return NULL;
}

static PympcObject *
Pympq2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_q(result->c, Pympq_AS_MPQ(self),
                               context->now.mpc_round);
    return result;
}

static PyObject *
Pympc2Pympq(PyObject *self)
{
    TYPE_ERROR("can not covert 'mpc' to 'mpq'");
    return NULL;
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
    TYPE_ERROR("can not covert 'mpc' to 'long'");
    return NULL;
}

#ifdef PY2
static PympcObject *
PyInt2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *result;

    if ((result = Pympc_new(rprec, iprec)))
        result->rc = mpc_set_si(result->c, PyInt_AsLong(self),
                                context->now.mpc_round);
    return result;
}
static PyObject *
Pympc2PyInt(PyObject *self)
{
    TYPE_ERROR("can not covert 'mpc' to 'int'");
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
    int firstp = 0, lastp = 0;

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
    mpfr_strtofr(mpc_realref(newob->c), cp, &tempchar, base,
                 GET_MPC_RROUND(context));
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
    mpfr_strtofr(mpc_imagref(newob->c), cp, &tempchar, base,
                 GET_MPC_IROUND(context));

    if (cp == tempchar && tempchar > lastchar)
        goto valid_string;

    if (*tempchar != 'j' && *cp != ' ')
        goto invalid_string;

    if (tempchar <= lastchar)
        goto invalid_string;

  valid_string:
    Py_XDECREF(ascii_str);
    return newob;

  invalid_string:
    VALUE_ERROR("invalid string in gmpy2.mpc()");
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

    return Py_BuildValue("(OO)", tempreal, tempimag);
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
        mpc_get_prec2(&pr, &pi, Pympc_AS_MPC(obj));
        if ((!rprec && !iprec) || (rprec == pr && iprec == pi)) {
            newob = (PympcObject*) obj;
            Py_INCREF(obj);
        }
        else {
            newob = Pympc2Pympc((PyObject*)obj, rprec, iprec);
        }
    }
    else if (Pympc_Check(obj)) {
        /* Handle the unlikely case where the exponent is no longer
         * valid and mpfr_check_range needs to be called. */

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
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) {
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
    else if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) {
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
        TYPE_ERROR("argument can not be converted to mpc");
        return 0;
    }
}

PyDoc_STRVAR(doc_mpc_digits,
"c.digits(base=10, prec=0) -> ((mant, exp, prec), (mant, exp, prec))\n\n"
"Returns up to 'prec' digits in the given base. If 'prec' is 0, as many\n"
"digits that are available are returned. No more digits than available\n"
"given c's precision are returned. 'base' must be between 2 and 62,\n"
"inclusive. The result consists of 2 three-element tuples containing the\n"
"mantissa, exponent, and number of bits of precision of the real and\n"
"imaginary components.");

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
"mpc(c, [precision=0]) -> mpc object\n"
"    Return a new mpc object from an existing complex number (either\n"
"    a Python complex object or another mpc object). If the precision\n"
"    is not specified, then the precision is taken from the current\n"
"    context. The rounding mode is always taken from the current\n"
"    context.\n\n"
"mpc(r, [i=0], [precision=0]) -> mpc object\n"
"    Return a new mpc object by converting two non-complex numbers into\n"
"    the real and imaginary components of an mpc object. If the prec-\n"
"    ision is not specified, then the precision is taken from the cur-\n"
"    rent context. The rounding mode is always taken from the current\n"
"    context.\n\n"
"mpc(s, [precision=0], [base=10]) -> mpc object\n"
"    Return a new mpc object by converting a string 's' into a complex\n"
"    number. If 'base' is omitted, then a base 10 representation is\n"
"    assumed otherwise a base between 2 and 36 can be specified. The\n"
"    precision and rounding modes are taken from the current context.\n\n"
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
        TYPE_ERROR("gmpy2.mpc() requires at least 1 non-keyword argument");
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
            VALUE_ERROR("base for gmpy2.mpc() must be in the interval 2..36.");
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
                    VALUE_ERROR("invalid value for precision in gmpy2.mpc().");
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
                    VALUE_ERROR("invalid value for precision in gmpy2.mpc().");
                    return NULL;
                }
            }
        }

        if (arg0) tempreal = Pympfr_From_Real(arg0, rbits);
        if (arg1) tempimag = Pympfr_From_Real(arg1, ibits);

        if (!tempreal) {
            if ((tempreal = Pympfr_new(rbits)))
                mpfr_set_ui(Pympfr_AS_MPFR(tempreal), 0, context->now.mpfr_round);
        }

        if (!tempimag) {
            if ((tempimag = Pympfr_new(ibits)))
                mpfr_set_ui(Pympfr_AS_MPFR(tempimag), 0, context->now.mpfr_round);
        }

        result = Pympc_new(rbits, ibits);
        if (!tempreal || !tempimag || !result) {
            Py_XDECREF(tempreal);
            Py_XDECREF(tempimag);
            Py_XDECREF(result);
            TYPE_ERROR("gmpy2.mpc() require string or numeric argument.");
            return NULL;
        }

        mpc_set_fr_fr(Pympc_AS_MPC(result), Pympfr_AS_MPFR(tempreal),
                      Pympfr_AS_MPFR(tempimag), GET_MPC_ROUND(context));
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
    }
    else {
        TYPE_ERROR("gmpy2.mpc() requires numeric or string argument");
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
"     optional width.precision\n"
"     optional rounding mode:\n"
"        'U' -> round toward plus infinity\n"
"        'D' -> round toward minus infinity\n"
"        'Y' -> round away from zero\n"
"        'Z' -> round toward zero\n"
"        'N' -> round to nearest\n"
"     optional output style:\n"
"        'P' -> Python style, 1+2j\n"
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
        TYPE_ERROR("requires mpc type");
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
        if (ifmt[1] == ' ' || ifmt[1] == '-') {
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

    tempbuf = PyMem_Malloc(rbuflen + ibuflen + 10);
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
        PyMem_Free(tempbuf);
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


