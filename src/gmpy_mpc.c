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
    PympcObject *newob;

    assert(Pympc_Check(self));
    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, Pympc_AS_MPC(self));
    if ((newob = Pympc_new(rprec, iprec)))
        mpc_set(newob->c, Pympc_AS_MPC(self), context->now.mpc_round);
    return newob;
}

/* Conversion to/from MPC
 * Python's string representation of a complex number differs from the format
 * used by MPC. Python concatenates the real and imaginary components and adds
 * a 'j' character to the end of string. MPC surrounds the complex number with
 * parentheses, requires a space between the real and imaginary components,
 * and does not tolerate the trailing 'j'. PyStr2Pympc tries to work around
 * the differences as follows:
 *
 *  1) If the first non-whitespace character is a '(', the MPC format is
 *     assumed.
 *  2) Otherwise, two MPFR-compatible numbers are read from the string and
 *     stored into the real and imaginary components respectively.
 */

static PympcObject *
PyStr2Pympc(PyObject *s, long base)
{
    PympcObject *newob;
    PyObject *ascii_str = NULL;
    Py_ssize_t len;
    char *cp, *tempchar, *lastchar;

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

    if (!(newob = Pympc_new(GET_MPC_RPREC(context),
                            GET_MPC_IPREC(context)))) {
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

    /* Skip leading whitespace. */
    while (isspace(*cp))
        cp++;

    if (*cp == '(') {
        newob->rc = mpc_set_str(newob->c, cp, base,
                                GET_MPC_ROUND(context));
        if (newob->rc == -1) goto invalid_string;
    }
    else {
        /* Read the read component first. */
        mpfr_strtofr(mpc_realref(newob->c), cp, &tempchar, base,
                     GET_MPC_RROUND(context));
        /* Verify that at least one valid character was read. */
        if (cp == tempchar) goto invalid_string;
        /* Read the imaginary component next. */
        cp = tempchar;
        mpfr_strtofr(mpc_imagref(newob->c), cp, &tempchar, base,
                     GET_MPC_IROUND(context));
        /* For a valid string, either nothing must be read and we consumed
         * the string or there must be a trailing j optionally followed by
         * whitespace.
         */
        if (!((cp == tempchar && *tempchar == '\00') ||
              (*tempchar == 'j' && tempchar == lastchar))) {
            goto invalid_string;
        }
    }
    Py_XDECREF(ascii_str);
    return newob;

  invalid_string:
    VALUE_ERROR("invalid string in gmpy2.mpc()");
    Py_DECREF((PyObject*)newob);
    Py_XDECREF(ascii_str);
    return NULL;
}

PyDoc_STRVAR(doc_g_mpc,
"mpc(c) -> mpc object\n"
"    Return a new mpc object from an existing complex number (either\n"
"    a Python complex object or another mpc object). The precision and\n"
"    rounding modes are taken from the current context.\n\n"
"mpc(r, [i=0]) -> mpc object\n"
"    Return a new mpc object by converting two non-complex numbers into\n"
"    the real and imaginary components of an mpc object. The precision\n"
"    and rounding modes are taken from the current context.\n\n"
"mpc(s, [base=10]) -> mpc object\n"
"    Return a new mpc object by converting a string 's' into a complex\n"
"    number. If 'base' is omitted, then a base 10 representation is\n"
"    assumed otherwise a base between 2 and 36 can be specified. The\n"
"    precision and rounding modes are taken from the current context.");
static PyObject *
Pygmpy_mpc(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PympcObject *result = NULL;
    PyObject *arg0 = NULL;
    int base = 10;
    static char *kwlist[] = {"n", "base", NULL};

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|i", kwlist,
                                      &arg0, &base)))
        return NULL;

    if (base < 2 || base > 36) {
        VALUE_ERROR("base for gmpy2.mpc() must be in the interval 2..36.");
        return NULL;
    }

    if (arg0 == NULL)
        TYPE_ERROR("gmpy2.mpc() requires 1 or 2 arguments");
    else if (isComplex(arg0))
        TYPE_ERROR("gmpy2.mpc() can't convert a number yet");
    else if (PyStrOrUnicode_Check(arg0))
        result = PyStr2Pympc(arg0, base);
    else
        TYPE_ERROR("gmpy2.mpc() requires numeric or string argument");

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

    if (rbits != DBL_MANT_DIG && ibits !=DBL_MANT_DIG)
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

#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
        0,                               /* nb_add                  */
        0,                               /* nb_subtract             */
        0,                               /* nb_multiply             */
        0,                               /* nb_remaider             */
        0,                               /* nb_divmod               */
        0,                               /* nb_power                */
        0,                               /* nb_negative             */
        0,                               /* nb_positive             */
        0,                               /* nb_absolute             */
        0,                               /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_int                  */
        0,                               /* nb_reserved             */
        0,                               /* nb_float                */
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
        0,                               /* nb_floor_divide         */
        0,                               /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
        0,                               /* nb_add                  */
        0,                               /* nb_subtract             */
        0,                               /* nb_multiply             */
        0,                               /* nb_divide               */
        0,                               /* nb_remaider             */
        0,                               /* nb_divmod               */
        0,                               /* nb_power                */
        0,                               /* nb_negative             */
        0,                               /* nb_positive             */
        0,                               /* nb_absolute             */
        0,                               /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
        0,                               /* nb_int                  */
        0,                               /* nb_long                 */
        0,                               /* nb_float                */
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
        0,                               /* nb_floor_divide         */
        0,                               /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyMethodDef Pympc_methods[] =
{
    { "__format__", Pympc_format, METH_VARARGS, doc_mpc_format },
    { NULL, NULL, 1 }
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
        0,                                  /* tp_hash          */
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
        0,                                  /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympc_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
        0,                                  /* tp_getset        */
};


