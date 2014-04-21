/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.c                                                              *
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

PyDoc_STRVAR(doc_g_mpc,
"mpc() -> mpc(0.0+0.0j)\n\n"
"      If no argument is given, return mpc(0.0+0.0j).\n\n"
"mpc(c [, precision=0]) -> mpc\n\n"
"      Return a new 'mpc' object from an existing complex number (either\n"
"      a Python complex object or another 'mpc' object).\n\n"
"mpc(real [,imag=0 [, precision=0]]) -> mpc\n\n"
"      Return a new 'mpc' object by converting two non-complex numbers\n"
"      into the real and imaginary components of an 'mpc' object.\n\n"
"mpc(s [, precision=0 [, base=10]]) -> mpc\n\n"
"      Return a new 'mpc' object by converting a string s into a complex\n"
"      number. If base is omitted, then a base-10 representation is\n"
"      assumed otherwise the base must be in the interval [2,36].\n\n"
"Note: The precision can be specified either a single number that\n"
"      is used for both the real and imaginary components, or as a\n"
"      tuple that can specify different precisions for the real\n"
"      and imaginary components.\n\n"
"      If a precision greater than or equal to 2 is specified, then it\n"
"      is used.\n\n"
"      A precision of 0 (the default) implies the precision of the\n"
"      current context is used.\n\n"
"      A precision of 1 minimizes the loss of precision by following\n"
"      these rules:\n"
"        1) If n is a radix-2 floating point number, then the full\n"
"           precision of n is retained.\n"
"        2) For all other n, the precision of the result is the context\n"
"           precision + guard_bits.\n" );

static PyObject *
Pygmpy_mpc(PyObject *self, PyObject *args, PyObject *kwargs)
{
    MPC_Object *result = NULL;
    MPFR_Object *tempreal = NULL, *tempimag = NULL;
    PyObject *arg0 = NULL, *arg1 = NULL, *prec = NULL;
    int base = 10;
    Py_ssize_t argc = 0, keywdc = 0;
    CTXT_Object *context = NULL;

    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t rprec = 0, iprec = 0;

    static char *kwlist_c[] = {"c", "precision", NULL};
    static char *kwlist_r[] = {"real", "imag", "precision", NULL};
    static char *kwlist_s[] = {"s", "precision", "base", NULL};

    CHECK_CONTEXT_SET_EXPONENT(context);

    argc = PyTuple_Size(args);
    if (kwargs)
        keywdc = PyDict_Size(kwargs);

    if (argc + keywdc > 3) {
        TYPE_ERROR("mpc() takes at most 3 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPC_New(0, 0, context))) {
            mpc_set_ui(result->c, 0, GET_MPC_ROUND(context));
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpc() requires at least one non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    /* If building an mpc from a string, there can be upto two additional
     * arguments. Note that precision can be either a single integer or
     * a tuple.
     */
    if (PyStrOrUnicode_Check(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|Oi", kwlist_s,
                                            &arg0, &prec, &base)))
                return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }

            if (rprec < 0 || iprec < 0) {
                VALUE_ERROR("precision for mpc() must be >= 0");
                return NULL;
            }

            if (PyErr_Occurred()) {
                VALUE_ERROR("invalid value for precision in mpc()");
                return NULL;
            }
        }

        if (base < 2 || base > 36) {
            VALUE_ERROR("base for mpc() must be in the interval [2,36]");
            return NULL;
        }

        return (PyObject*)GMPy_MPC_From_PyStr(arg0, base, rprec, iprec, context);
    }

    if (IS_REAL(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|OO", kwlist_r,
                                            &arg0, &arg1, &prec)))
                return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }

            if (rprec < 0 || iprec < 0) {
                VALUE_ERROR("precision for mpc() must be >= 0");
                return NULL;
            }

            if (PyErr_Occurred()) {
                VALUE_ERROR("invalid value for precision in mpc()");
                return NULL;
            }
        }

        if (arg1 && !IS_REAL(arg1)) {
            TYPE_ERROR("invalid type for imaginary component in mpc()");
            return NULL;
        }

        tempreal = GMPy_MPFR_From_Real(arg0, rprec, context);
        if (arg1) {
            tempimag = GMPy_MPFR_From_Real(arg1, iprec, context);
        }
        else {
            if ((tempimag = GMPy_MPFR_New(iprec, context))) {
                mpfr_set_ui(tempimag->f, 0, MPFR_RNDN);
            }
        }

        result = GMPy_MPC_New(rprec, iprec, context);
        if (!tempreal || !tempimag || !result) {
            Py_XDECREF(tempreal);
            Py_XDECREF(tempimag);
            Py_XDECREF(result);
            TYPE_ERROR("mpc() requires string or numeric argument.");
            return NULL;
        }

        mpc_set_fr_fr(result->c, tempreal->f, tempimag->f, GET_MPC_ROUND(context));
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
        return (PyObject*)result;
    }

    if (IS_COMPLEX_ONLY(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist_c,
                                          &arg0, &prec)))
            return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            if (PyErr_Occurred()) {
                VALUE_ERROR("invalid value for precision in mpc()");
                return NULL;
            }
        }

        if (PyComplex_Check(arg0)) {
            result = GMPy_MPC_From_PyComplex(arg0, rprec, iprec, context);
        }
        else {
            result = GMPy_MPC_From_MPC((MPC_Object*)arg0, rprec, iprec, context);
        }
        return (PyObject*)result;
    }

    TYPE_ERROR("mpc() requires numeric or string argument");
    return NULL;
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

    if (!MPC_Check(self)) {
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
                           mpc_realref(MPC(self)));

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
                           mpc_imagref(MPC(self)));

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
        if (mpfr_nan_p(mpc_imagref(MPC(self))) ||
            (mpfr_inf_p(mpc_imagref(MPC(self))) &&
             mpfr_sgn(mpc_imagref(MPC(self))) > 0)) {
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

    if (self && MPC_Check(self)) {
        if (!PyArg_ParseTuple(args, "|ii", &base, &prec))
            return NULL;
        Py_INCREF(self);
    }
    else {
        if(!PyArg_ParseTuple(args, "O&|ii", GMPy_MPC_convert_arg, &self,
                            &base, &prec))
        return NULL;
    }
    result = GMPy_PyStr_From_MPC((MPC_Object*)self, base, prec, NULL);
    Py_DECREF(self);
    return result;
}

/* Support Pympany_square */

static PyObject *
Pympc_sqr(PyObject* self, PyObject *other)
{
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("square() requires 'mpc' argument");

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_sqr(result->c, MPC(self),
                         GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "square()");
}

/* Implement the conjugate() method. */

PyDoc_STRVAR(doc_mpc_conjugate,
"x.conjugate() -> mpc\n\n"
"Returns the conjugate of x.");

static PyObject *
Pympc_conjugate(PyObject *self, PyObject *args)
{
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_ARGS("conjugate() requires 'mpc' argument");

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_conj(result->c, MPC(self),
                          GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "conjugate()");
}

/* Implement the .precision attribute of an mpfr. */

static PyObject *
Pympc_getprec_attrib(MPC_Object *self, void *closure)
{
    mpfr_prec_t rprec = 0, iprec = 0;

    mpc_get_prec2(&rprec, &iprec, self->c);
    return Py_BuildValue("(nn)", rprec, iprec);
}

/* Implement the .rc attribute of an mpfr. */

static PyObject *
Pympc_getrc_attrib(MPC_Object *self, void *closure)
{
    return Py_BuildValue("(ii)", MPC_INEX_RE(self->rc), MPC_INEX_IM(self->rc));
}

/* Implement the .imag attribute of an mpfr. */

static PyObject *
Pympc_getimag_attrib(MPC_Object *self, void *closure)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPFR_New(0, context)))
        mpc_imag(result->f, self->c, context->ctx.mpfr_round);
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpfr. */

static PyObject *
Pympc_getreal_attrib(MPC_Object *self, void *closure)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if ((result = GMPy_MPFR_New(0, context)))
        mpc_real(result->f, self->c, context->ctx.mpfr_round);
    return (PyObject*)result;
}

/* Implement the nb_bool slot. */

static int
Pympc_nonzero(MPC_Object *self)
{
    return !MPC_IS_ZERO_P(self->c);
}

/* To work with the MPC_IS_ macros, NAN, INF, FINITE, and ZERO are
 * all upper-case.
 */

#define MPC_TEST_OTHER(NAME, msg) \
static PyObject * \
Pympc_is_##NAME(PyObject *self, PyObject *other)\
{\
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT_SET_EXPONENT(context); \
    if(self && MPC_Check(self)) {\
        Py_INCREF(self);\
    }\
    else if(MPC_Check(other)) {\
        self = other;\
        Py_INCREF((PyObject*)self);\
    }\
    else if (!(self = (PyObject*)GMPy_MPC_From_Complex(other, 0, 0, context))) {\
        PyErr_SetString(PyExc_TypeError, msg);\
        return NULL;\
    }\
    if (MPC_IS_##NAME##_P(self)) {\
        Py_DECREF(self);\
        Py_RETURN_TRUE;\
    }\
    else {\
        Py_DECREF(self);\
        Py_RETURN_FALSE;\
    }\
}

MPC_TEST_OTHER(NAN, "is_nan() requires 'mpc' argument");

MPC_TEST_OTHER(INF, "is_infinite() requires 'mpc' argument");

MPC_TEST_OTHER(FINITE, "is_finite() requires 'mpc' argument");

MPC_TEST_OTHER(ZERO, "is_zero() requires 'mpc' argument");

PyDoc_STRVAR(doc_mpc_phase,
"phase(x) -> mpfr\n\n"
"Return the phase angle, also known as argument, of a complex x.");

static PyObject *
Pympc_phase(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("phase() requires 'mpc' argument");

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_arg(result->f, MPC(self),
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
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("norm() requires 'mpc' argument");

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_norm(result->f, MPC(self),
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
    PyObject *abs, *phase, *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("norm() requires 'mpc' argument");

    if (!(abs = GMPy_Complex_Abs(self, context))) {
        Py_DECREF(self);
        return NULL;
    }
    if (!(phase = Pympc_phase(self, other))) {
        Py_DECREF(abs);
        Py_DECREF(self);
        return NULL;
    }

    result = Py_BuildValue("(NN)", abs, phase);
    if (!result) {
        Py_DECREF(abs);
        Py_DECREF(phase);
    }
    Py_DECREF(self);
    return result;
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
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "rect() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_cos(mpc_realref(result->c), MPFR(other),
             GET_REAL_ROUND(context));
    mpfr_mul(mpc_realref(result->c), mpc_realref(result->c),
             MPFR(self), GET_REAL_ROUND(context));
    mpfr_sin(mpc_imagref(result->c), MPFR(other),
             GET_IMAG_ROUND(context));
    mpfr_mul(mpc_imagref(result->c), mpc_imagref(result->c),
             MPFR(self), GET_IMAG_ROUND(context));
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
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("proj() requires 'mpc' argument");

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_proj(result->c, MPC(self),
                          GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "proj()");
}

#define MPC_UNIOP(NAME) \
static PyObject * \
Pympc_##NAME(PyObject* self, PyObject *other) \
{ \
    MPC_Object *result; \
    CTXT_Object *context = NULL; \
    CHECK_CONTEXT_SET_EXPONENT(context); \
    PARSE_ONE_MPC_OTHER(#NAME "() requires 'mpc' argument"); \
    if (!(result = GMPy_MPC_New(0, 0, context))) { \
        Py_DECREF(self); \
        return NULL; \
    } \
    result->rc = mpc_##NAME(result->c, MPC(self), GET_MPC_ROUND(context)); \
    Py_DECREF(self); \
    MPC_CLEANUP(result, #NAME"()"); \
}

MPC_UNIOP(log)

MPC_UNIOP(log10)

MPC_UNIOP(exp)

MPC_UNIOP(sqrt)

static PyObject *
Pympc_sin_cos(PyObject *self, PyObject *other)
{
    MPC_Object *s, *c;
    PyObject *result;
    int code;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("sin_cos() requires 'mpc' argument");

    s = GMPy_MPC_New(0, 0, context);
    c = GMPy_MPC_New(0, 0, context);
    result = PyTuple_New(2);
    if (!s || !c || !result) {
        Py_DECREF(self);
        return NULL;
    }

    code = mpc_sin_cos(s->c, c->c, MPC(self),
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
    MPC_Object *result, *tempx, *tempy, *tempz;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fma() requires 'mpc','mpc','mpc' arguments.");
        return NULL;
    }

    result = GMPy_MPC_New(0, 0, context);
    tempx = GMPy_MPC_From_Complex(PyTuple_GET_ITEM(args, 0), 0, 0, context);
    tempy = GMPy_MPC_From_Complex(PyTuple_GET_ITEM(args, 1), 0, 0, context);
    tempz = GMPy_MPC_From_Complex(PyTuple_GET_ITEM(args, 2), 0, 0, context);
    if (!result || !tempx || !tempy || !tempz) {
        TYPE_ERROR("fma() requires 'mpc','mpc','mpc' arguments.");
        goto done;
    }

    result->rc = mpc_fma(result->c, tempx->c, tempy->c, tempz->c,
                         GET_MPC_ROUND(context));
    MPC_SUBNORMALIZE(result);
    MPC_CHECK_FLAGS(result, "fma()");

  done:
    Py_XDECREF((PyObject*)tempx);
    Py_XDECREF((PyObject*)tempy);
    Py_XDECREF((PyObject*)tempz);
    if (PyErr_Occurred()) {
        Py_XDECREF(result);
        result = NULL;
    }
    return (PyObject*)result;
}

static PyObject *
Pympc_fms(PyObject *self, PyObject *args)
{
    MPC_Object *result, *x, *y, *z;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (PyTuple_GET_SIZE(args) != 3) {
        TYPE_ERROR("fms() requires 'mpc','mpc','mpc' arguments.");
        return NULL;
    }

    result = GMPy_MPC_New(0, 0, context);
    x = GMPy_MPC_From_Complex(PyTuple_GET_ITEM(args, 0), 0, 0, context);
    y = GMPy_MPC_From_Complex(PyTuple_GET_ITEM(args, 1), 0, 0, context);
    z = GMPy_MPC_From_Complex(PyTuple_GET_ITEM(args, 2), 0, 0, context);
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
    MPC_Object *result = 0;
    unsigned long exp = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!PyArg_ParseTuple(args, "O&k", GMPy_MPC_convert_arg, &self, &exp)) {
        TYPE_ERROR("div_2exp() requires 'mpc', 'int' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_div_2ui(MPC(result), MPC(self),
                              exp, GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "div_2exp()");
}

static PyObject *
Pympc_mul_2exp(PyObject *self, PyObject *args)
{
    MPC_Object *result = 0;
    unsigned long exp = 0;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    if (!PyArg_ParseTuple(args, "O&k", GMPy_MPC_convert_arg, &self, &exp)) {
        TYPE_ERROR("mul_2exp() requires 'mpc', 'int' arguments");
        return NULL;
    }

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_mul_2ui(MPC(result), MPC(self),
                             exp, GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "mul_2exp()");
}

static Py_hash_t
Pympc_hash(MPC_Object *self)
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

PyDoc_STRVAR(doc_mpc_sizeof,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x.");

static PyObject *
Pympc_sizeof(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPC_Object) + \
        (((mpc_realref(MPC(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t)) + \
        (((mpc_imagref(MPC(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t)));
}

static PyMethodDef Pympc_methods[] =
{
    { "__complex__", GMPy_PyComplex_From_MPC, METH_O, doc_mpc_complex },
    { "__format__", Pympc_format, METH_VARARGS, doc_mpc_format },
    { "__sizeof__", Pympc_sizeof, METH_NOARGS, doc_mpc_sizeof },
    { "conjugate", Pympc_conjugate, METH_NOARGS, doc_mpc_conjugate },
    { "digits", Pympc_digits, METH_VARARGS, doc_mpc_digits },
    { NULL, NULL, 1 }
};


#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) GMPy_MPC_Add_Slot,      /* nb_add                  */
    (binaryfunc) GMPy_MPC_Sub_Slot,      /* nb_subtract             */
    (binaryfunc) GMPy_MPC_Mul_Slot,      /* nb_multiply             */
    (binaryfunc) GMPy_MPC_Mod_Slot,      /* nb_remainder            */
    (binaryfunc) GMPy_MPC_DivMod_Slot,   /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,   /* nb_power                */
    (unaryfunc) GMPy_MPC_Minus_Slot,     /* nb_negative             */
    (unaryfunc) GMPy_MPC_Plus_Slot,      /* nb_positive             */
    (unaryfunc) GMPy_MPC_Abs_Slot,       /* nb_absolute             */
    (inquiry) Pympc_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) GMPy_MPC_Int_Slot,       /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) GMPy_MPC_Float_Slot,     /* nb_float                */
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
    (binaryfunc) GMPy_MPC_FloorDiv_Slot, /* nb_floor_divide         */
    (binaryfunc) GMPy_MPC_TrueDiv_Slot,  /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) GMPy_MPC_Add_Slot,      /* nb_add                  */
    (binaryfunc) GMPy_MPC_Sub_Slot,      /* nb_subtract             */
    (binaryfunc) GMPy_MPC_Mul_Slot,      /* nb_multiply             */
    (binaryfunc) GMPy_MPC_TrueDiv_Slot,  /* nb_divide               */
    (binaryfunc) GMPy_MPC_Mod_Slot,      /* nb_remainder            */
    (binaryfunc) GMPy_MPC_DivMod_Slot,   /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,   /* nb_power                */
    (unaryfunc) GMPy_MPC_Minus_Slot,     /* nb_negative             */
    (unaryfunc) GMPy_MPC_Plus_Slot,      /* nb_positive             */
    (unaryfunc) GMPy_MPC_Abs_Slot,       /* nb_absolute             */
    (inquiry) Pympc_nonzero,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) GMPy_MPC_Int_Slot,       /* nb_int                  */
    (unaryfunc) GMPy_MPC_Long_Slot,      /* nb_long                 */
    (unaryfunc) GMPy_MPC_Float_Slot,     /* nb_float                */
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
    (binaryfunc) GMPy_MPC_FloorDiv_Slot, /* nb_floor_divide         */
    (binaryfunc) GMPy_MPC_TrueDiv_Slot,  /* nb_true_divide          */
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

static PyTypeObject MPC_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpc",                                  /* tp_name          */
    sizeof(MPC_Object),                     /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPC_Dealloc,          /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPC_Repr_Slot,          /* tp_repr          */
    &mpc_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) Pympc_hash,                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPC_Str_Slot,           /* tp_str           */
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


