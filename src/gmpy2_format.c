/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_format.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015, 2016, 2017, 2018, 2019 Case Van Horsen                  *
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

PyDoc_STRVAR(GMPy_doc_mpz_format,
"x.__format__(fmt) -> string\n\n"
"Return a Python string by formatting mpz 'x' using the format string\n"
"'fmt'. A valid format string consists of:\n"
"     optional alignment code:\n"
"        '<' -> left shifted in field\n"
"        '>' -> right shifted in field\n"
"        '^' -> centered in field\n"
"     optional leading sign code:\n"
"        '+' -> always display leading sign\n"
"        '-' -> only display minus sign\n"
"        ' ' -> minus for negative values, space for positive values\n"
"     optional base indicator\n"
"        '#' -> precede binary, octal, or hex with 0b, 0o or 0x\n"
"     optional width\n"
"     optional conversion code:\n"
"        'd' -> decimal format\n"
"        'b' -> binary format\n"
"        'o' -> octal format\n"
"        'x' -> hex format\n"
"        'X' -> upper-case hex format\n"
"The default format is 'd'.");

/* Formatting occurs in two phases. Pympz_ascii() is used to create a string
 * with the appropriate binary/octal/decimal/hex formatting, including the
 * leading sign character (+ , -, or space) and base encoding (0b, 0o, or 0x).
 * Left/right/centering using the specified width is done by creating a
 * format string and calling the __format__() method of the string object
 * returned by Pympz_ascii().
 */

static PyObject *
GMPy_MPZ_Format(PyObject *self, PyObject *args)
{
    PyObject *result = NULL, *mpzstr = NULL;
    char *fmtcode = 0, *p1, *p2;
    char fmt[30];
    int base = 10, option = 16;
    int seensign = 0, seenindicator = 0, seenalign = 0, seendigits = 0;

    if (!CHECK_MPZANY(self)) {
        TYPE_ERROR("requires mpz type");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "s", &fmtcode))
        return NULL;

    p2 = fmt;
    for (p1 = fmtcode; *p1 != '\00'; p1++) {
        if (*p1 == '<' || *p1 == '>' || *p1 == '^') {
            if (seenalign || seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seenalign = 1;
                continue;
            }
        }
        if (*p1 == '+') {
            if (seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                option |= 2;
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '-') {
            if (seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                seensign = 1;
                continue;
            }
        }
        if (*p1 == ' ') {
            if (seensign || seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                option |= 4;
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '#') {
            if (seenindicator || seendigits) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                option |= 8;
                seenindicator = 1;
                continue;
            }
        }
        if (isdigit(*p1)) {
            if (!seenalign) {
                *(p2++) = '>';
                seenalign = 1;
            }
            *(p2++) = *p1;
            seendigits = 1;
            continue;
        }
        if (*p1 == 'b') {
            base = 2;
            break;
        }
        if (*p1 == 'o') {
            base = 8;
            break;
        }
        if (*p1 == 'x') {
            base = 16;
            break;
        }
        if (*p1 == 'd') {
            base = 10;
            break;
        }
        if (*p1 == 'X') {
            base = -16;
            break;
        }
        VALUE_ERROR("Invalid conversion specification");
        return NULL;
    }
    *(p2++) = '\00';

    if (!(mpzstr = mpz_ascii(MPZ(self), base, option, 0)))
        return NULL;

    result = PyObject_CallMethod(mpzstr, "__format__", "(s)", fmt);
    Py_DECREF(mpzstr);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpfr_format,
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
"        'U' -> round toward plus Infinity\n"
"        'D' -> round toward minus Infinity\n"
"        'Y' -> round away from zero\n"
"        'Z' -> round toward zero\n"
"        'N' -> round to nearest\n"
"     optional conversion code:\n"
"        'a','A' -> hex format\n"
"        'b'     -> binary format\n"
"        'e','E' -> scientific format\n"
"        'f','F' -> fixed point format\n"
"        'g','G' -> fixed or float format\n\n"
"The default format is '.6f'.");

static PyObject *
GMPy_MPFR_Format(PyObject *self, PyObject *args)
{
    PyObject *result = NULL, *mpfrstr = NULL;
    char *buffer = 0, *newbuf = 0, *fmtcode = 0, *p1, *p2, *p3;
    char mpfrfmt[100], fmt[30];
    int buflen;
    int seensign = 0, seenalign = 0, seendecimal = 0, seendigits = 0;
    int seenround = 0, seenconv = 0;

    if (!MPFR_Check(self)) {
        TYPE_ERROR("requires mpfr type");
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "s", &fmtcode))
        return NULL;

    p2 = mpfrfmt;
    p3 = fmt;
    *(p2++) = '%';

    for (p1 = fmtcode; *p1 != '\00'; p1++) {
        if (*p1 == '<' || *p1 == '>' || *p1 == '^') {
            if (seenalign || seensign || seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p3++) = *p1;
                seenalign = 1;
                continue;
            }
        }
        if (*p1 == '+' || *p1 == ' ') {
            if (seensign || seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '-') {
            if (seensign || seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                seensign = 1;
                continue;
            }
        }
        if (*p1 == '.') {
            if (seendecimal || seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seendecimal = 1;
                continue;
            }
        }
        if (isdigit(*p1)) {
            if (seendigits || seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else if (seendecimal) {
                *(p2++) = *p1;
                continue;
            }
            else {
                if (p3 == fmt) {
                    *(p3++) = '>';
                    seenalign = 1;
                }
                *(p3++) = *p1;
                continue;
            }
        }
        if (!seendigits) {
            seendigits = 1;
            *(p2++) = 'R';
        }
        if (*p1 == 'U' || *p1 == 'D' || *p1 == 'Y' || *p1 == 'Z' ||
            *p1 == 'N' ) {
            if (seenround) {
                VALUE_ERROR("Invalid conversion specification");
                return NULL;
            }
            else {
                *(p2++) = *p1;
                seenround = 1;
                continue;
            }
        }
        if (*p1 == 'a' || *p1 == 'A' || *p1 == 'b' || *p1 == 'e' ||
            *p1 == 'E' || *p1 == 'f' || *p1 == 'F' || *p1 == 'g' ||
            *p1 == 'G' ) {
            *(p2++) = *p1;
            seenconv = 1;
            break;
        }
        VALUE_ERROR("Invalid conversion specification");
        return NULL;
    }

    if (!seendigits)
        *(p2++) = 'R';
    if (!seenconv)
        *(p2++) = 'f';

    *(p2) = '\00';
    *(p3) = '\00';

    buflen = mpfr_asprintf(&buffer, mpfrfmt, MPFR(self));

    /* If there isn't a decimal point in the output and the output
     * only consists of digits, then append .0 */
    if (strlen(buffer) == strspn(buffer, "+- 0123456789")) {
        newbuf = malloc(buflen + 3);
        if (!newbuf) {
            mpfr_free_str(buffer);
            return PyErr_NoMemory();
        }
        *newbuf = '\0';
        strcat(newbuf, buffer);
        strcat(newbuf, ".0");
        mpfr_free_str(buffer);
        mpfrstr = Py_BuildValue("s", newbuf);
        free(newbuf);
    }
    else {
        mpfrstr = Py_BuildValue("s", buffer);
        mpfr_free_str(buffer);
    }
    if (!mpfrstr) {
        return NULL;
    }

    result = PyObject_CallMethod(mpfrstr, "__format__", "(s)", fmt);
    Py_DECREF(mpfrstr);
    return result;
}

PyDoc_STRVAR(GMPy_doc_mpc_format,
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
GMPy_MPC_Format(PyObject *self, PyObject *args)
{
    PyObject *result = NULL, *tempstr = NULL;
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

    if (!PyArg_ParseTuple(args, "s", &fmtcode)) {
        return NULL;
    }

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

    tempbuf = malloc(rbuflen + ibuflen + 10);
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
        free(tempbuf);
        return NULL;
    }

    result = PyObject_CallMethod(tempstr, "__format__", "(s)", fmt);

    Py_DECREF(tempstr);
    return result;
}

/* produce digits for an mpz in requested base, default 10 */
PyDoc_STRVAR(GMPy_doc_mpz_digits_method,
"x.digits([base=10]) -> string\n\n"
"Return Python string representing x in the given base. Values for\n"
"base can range between 2 to 62. A leading '-' is present if x<0\n"
"but no leading '+' is present if x>=0.");

static PyObject *
GMPy_MPZ_Digits_Method(PyObject *self, PyObject *args)
{
    int base = 10;

    if (PyTuple_GET_SIZE(args) && !PyArg_ParseTuple(args, "|i", &base)) {
        return NULL;
    }

    return GMPy_PyStr_From_MPZ((MPZ_Object*)self, base, 16, NULL);
}

static PyObject *
GMPy_XMPZ_Digits_Method(PyObject *self, PyObject *args)
{
    int base = 10;

    if (PyTuple_GET_SIZE(args) && !PyArg_ParseTuple(args, "|i", &base)) {
        return NULL;
    }

    return  GMPy_PyStr_From_XMPZ((XMPZ_Object*)self, base, 0, NULL);
}

PyDoc_STRVAR(GMPy_doc_mpq_digits_method,
"x.digits([base=10]) -> string\n\n"
"Return a Python string representing x in the given base (2 to 62,\n"
"default is 10). A leading '-' is present if x<0, but no leading '+'\n"
"is present if x>=0.\n");

static PyObject *
GMPy_MPQ_Digits_Method(PyObject *self, PyObject *args)
{
    int base = 10;

    if (PyTuple_GET_SIZE(args) && !PyArg_ParseTuple(args, "|i", &base)) {
        return NULL;
    }

    return GMPy_PyStr_From_MPQ((MPQ_Object*)self, base, 0, NULL);
}


PyDoc_STRVAR(GMPy_doc_mpfr_digits_method,
"x.digits([base=10[, prec=0]]) -> (mantissa, exponent, bits)\n\n"
"Returns up to 'prec' digits in the given base. If 'prec' is 0, as many\n"
"digits that are available are returned. No more digits than available\n"
"given x's precision are returned. 'base' must be between 2 and 62,\n"
"inclusive. The result is a three element tuple containing the mantissa,\n"
"the exponent, and the number of bits of precision.");

static PyObject *
GMPy_MPFR_Digits_Method(PyObject *self, PyObject *args)
{
    int base = 10, prec = 0;

    if (PyTuple_GET_SIZE(args) && !PyArg_ParseTuple(args, "|ii", &base, &prec)) {
        return NULL;
    }

    return GMPy_PyStr_From_MPFR((MPFR_Object*)self, base, prec, NULL);
}

PyDoc_STRVAR(GMPy_doc_mpc_digits_method,
"c.digits(base=10, prec=0) -> ((mant, exp, prec), (mant, exp, prec))\n\n"
"Returns up to 'prec' digits in the given base. If 'prec' is 0, as many\n"
"digits that are available given c's precision are returned. 'base' must\n"
"be between 2 and 62. The result consists of 2 three-element tuples that\n"
"contain the mantissa, exponent, and number of bits of precision of the\n"
"real and imaginary components.");

static PyObject *
GMPy_MPC_Digits_Method(PyObject *self, PyObject *args)
{
    int base = 10, prec = 0;

    if (PyTuple_GET_SIZE(args) && !PyArg_ParseTuple(args, "|ii", &base, &prec)) {
        return NULL;
    }

    return GMPy_PyStr_From_MPC((MPC_Object*)self, base, prec, NULL);
}

PyDoc_STRVAR(GMPy_doc_context_digits,
"digits(x[, base[, prec]]) -> string\n\n"
"Return string representing x. Calls mpz.digits, mpq.digits,\n"
"mpfr.digits, or mpc.digits as appropriate.");

static PyObject *
GMPy_Context_Digits(PyObject *self, PyObject *args)
{
    PyObject *arg0, *tuple, *temp, *result;
    Py_ssize_t argc;

    argc = PyTuple_GET_SIZE(args);
    if (argc == 0) {
        TYPE_ERROR("digits() requires at least one argument");
        return NULL;
    }

    if (argc > 3) {
        TYPE_ERROR("digits() accepts at most three arguments");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);
    if (!(tuple = PyTuple_GetSlice(args, 1, argc))) {
        return NULL;
    }

    if (IS_INTEGER(arg0)) {
        temp = (PyObject*)GMPy_MPZ_From_Integer(arg0, NULL);
        result = GMPy_MPZ_Digits_Method(temp, tuple);
        Py_DECREF(temp);
        Py_DECREF(tuple);
        return result;
    }
    if (IS_RATIONAL(arg0)) {
        temp = (PyObject*)GMPy_MPQ_From_Rational(arg0, NULL);
        result = GMPy_MPQ_Digits_Method(temp, tuple);
        Py_DECREF(temp);
        Py_DECREF(tuple);
        return result;
    }
    if (IS_REAL(arg0)) {
        temp = (PyObject*)GMPy_MPFR_From_Real(arg0, 1, NULL);
        result = GMPy_MPFR_Digits_Method(temp, tuple);
        Py_DECREF(temp);
        Py_DECREF(tuple);
        return result;
    }
    if (IS_COMPLEX(arg0)) {
        temp = (PyObject*)GMPy_MPC_From_Complex(arg0, 1, 1, NULL);
        result = GMPy_MPC_Digits_Method(temp, tuple);
        Py_DECREF(temp);
        Py_DECREF(tuple);
        return result;
    }

    Py_DECREF(tuple);
    TYPE_ERROR("digits() argument type not supported");
    return NULL;
}


