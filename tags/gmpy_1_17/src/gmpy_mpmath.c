/* gmpy_mpmath.c
 *
 * Internal helper function for mpmath.
 *
 * This file should be considered part of gmpy.c.
 */

static PyObject *
mpmath_build_mpf(long sign, PympzObject *man, PyObject *exp, long bc)
{
    PyObject *tup, *tsign, *tbc;
    if(!(tup = PyTuple_New(4))){
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        return NULL;
    }
    if(!(tsign=Py2or3Int_FromLong(sign))){
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        Py_DECREF(tup);
        return NULL;
    }
    if(!(tbc=Py2or3Int_FromLong(bc))){
        Py_DECREF((PyObject*)man);
        Py_DECREF(exp);
        Py_DECREF(tup);
        Py_DECREF(tsign);
        return NULL;
    }
    PyTuple_SET_ITEM(tup, 0, tsign);
    PyTuple_SET_ITEM(tup, 1, (PyObject*)man);
    PyTuple_SET_ITEM(tup, 2, (exp)?exp:Py2or3Int_FromLong(0));
    PyTuple_SET_ITEM(tup, 3, tbc);
    return tup;
}

static char doc_mpmath_normalizeg[]="\
_mpmath_normalize(...): helper function for mpmath.\n\
";
static PyObject *
Pympz_mpmath_normalize(PyObject *self, PyObject *args)
{
    long sign = 0, bc = 0, prec = 0, shift, zbits, carry = 0;
    PyObject *exp = 0, *newexp = 0, *newexp2 = 0, *tmp = 0;
    PympzObject *man = 0;
    mpz_t upper, lower;
    char rnd = 0;

    if(PyTuple_GET_SIZE(args) == 6){
        /* Need better error-checking here. Under Python 3.0, overflow into
           C-long is possible. */
        sign = clong_From_Integer(PyTuple_GET_ITEM(args, 0));
        man = (PympzObject *)PyTuple_GET_ITEM(args, 1);
        exp = PyTuple_GET_ITEM(args, 2);
        bc = clong_From_Integer(PyTuple_GET_ITEM(args, 3));
        prec = clong_From_Integer(PyTuple_GET_ITEM(args, 4));
        rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 5))[0];
        if(PyErr_Occurred()){
            PyErr_SetString(PyExc_TypeError, "arguments long, PympzObject*,"
                "PyObject*, long, long, char needed");
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "6 arguments required");
        return NULL;
    }
    if(!Pympz_Check(man)){
        PyErr_SetString(PyExc_TypeError, "argument is not an mpz");
        return NULL;
    }

    /* If the mantissa is 0, return the normalized representation. */
    if(!mpz_sgn(man->z)) {
        Py_INCREF((PyObject*)man);
        return mpmath_build_mpf(0, man, 0, 0);
    }

    /* if bc <= prec and the number is odd return it */
    if ((bc <= prec) && mpz_odd_p(man->z)) {
        Py_INCREF((PyObject*)man);
        Py_INCREF((PyObject*)exp);
        return mpmath_build_mpf(sign, man, exp, bc);
    }

    mpz_inoc(upper);
    mpz_inoc(lower);

    shift = bc - prec;
    if(shift>0) {
        switch(rnd) {
            case 'f':
                if(sign) {
                    mpz_cdiv_q_2exp(upper, man->z, shift);
                } else {
                    mpz_fdiv_q_2exp(upper, man->z, shift);
                }
                break;
            case 'c':
                if(sign) {
                    mpz_fdiv_q_2exp(upper, man->z, shift);
                } else {
                    mpz_cdiv_q_2exp(upper, man->z, shift);
                }
                break;
            case 'd':
                mpz_fdiv_q_2exp(upper, man->z, shift);
                break;
            case 'u':
                mpz_cdiv_q_2exp(upper, man->z, shift);
                break;
            case 'n':
            default:
                mpz_tdiv_r_2exp(lower, man->z, shift);
                mpz_tdiv_q_2exp(upper, man->z, shift);
                if(mpz_sgn(lower)) {
                    /* lower is not 0 so it must have at least 1 bit set */
                    if(mpz_sizeinbase(lower, 2)==shift) {
                        /* lower is >= 1/2 */
                        if(mpz_scan1(lower, 0)==shift-1) {
                            /* lower is exactly 1/2 */
                            if(mpz_odd_p(upper))
                                carry = 1;
                        } else {
                            carry = 1;
                        }
                    }
                }
                if(carry)
                    mpz_add_ui(upper, upper, 1);
        }
        if (!(tmp = Py2or3Int_FromLong(shift))) {
            mpz_cloc(upper);
            mpz_cloc(lower);
            return NULL;
        }
        if (!(newexp = PyNumber_Add(exp, tmp))) {
            mpz_cloc(upper);
            mpz_cloc(lower);
            Py_DECREF(tmp);
            return NULL;
        }
        Py_DECREF(tmp);
        bc = prec;
    } else {
        mpz_set(upper, man->z);
        newexp = exp;
        Py_INCREF(newexp);
    }

    /* Strip trailing 0 bits. */
    if((zbits = mpz_scan1(upper, 0)))
        mpz_tdiv_q_2exp(upper, upper, zbits);

    if (!(tmp = Py2or3Int_FromLong(zbits))) {
        mpz_cloc(upper);
        mpz_cloc(lower);
        Py_DECREF(newexp);
        return NULL;
    }
    if (!(newexp2 = PyNumber_Add(newexp, tmp))) {
        mpz_cloc(upper);
        mpz_cloc(lower);
        Py_DECREF(tmp);
        Py_DECREF(newexp);
        return NULL;
    }
    Py_DECREF(newexp);
    Py_DECREF(tmp);

    bc -= zbits;
    /* Check if one less than a power of 2 was rounded up. */
    if(!mpz_cmp_ui(upper, 1))
        bc = 1;

    mpz_cloc(lower);
    return mpmath_build_mpf(sign, Pympz_FROM_MPZ(upper), newexp2, bc);
}

static char doc_mpmath_createg[]="\
_mpmath_create(...): helper function for mpmath.\n\
";
static PyObject *
Pympz_mpmath_create(PyObject *self, PyObject *args)
{
    long sign, bc, shift, zbits, carry = 0;
    PyObject *exp = 0, *newexp = 0, *newexp2 = 0, *tmp = 0;
    PympzObject *man = 0, *upper = 0, *lower = 0;

    const char *rnd = "f";
    long prec = 0;

    if(PyTuple_GET_SIZE(args) < 2) {
        PyErr_SetString(PyExc_TypeError,
                "mpmath_create() expects 'mpz','int'[,'int','str'] arguments");
        return NULL;
    }

    switch(PyTuple_GET_SIZE(args)) {
        case 4:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 3));
        case 3:
            prec = clong_From_Integer(PyTuple_GET_ITEM(args, 2));
            if(prec == -1 && PyErr_Occurred())
                return NULL;
            prec = abs(prec);
        case 2:
            exp = PyTuple_GET_ITEM(args, 1);
        case 1:
            man = Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
            if(!man) {
                PyErr_SetString(PyExc_TypeError,
                        "mpmath_create() expects 'mpz','int'[,'int','str'] arguments");
                return NULL;
            }
    }

    /* If the mantissa is 0, return the normalized representation. */
    if(!mpz_sgn(man->z)) {
        return mpmath_build_mpf(0, man, 0, 0);
    }

    upper = Pympz_new();
    lower = Pympz_new();
    if(!upper||!lower) {
        Py_DECREF((PyObject*)man);
        Py_XDECREF((PyObject*)upper);
        Py_XDECREF((PyObject*)lower);
        return NULL;
    }

    /* Extract sign, make man positive, and set bit count */
    sign = (mpz_sgn(man->z) == -1);
    mpz_abs(upper->z, man->z);
    bc = mpz_sizeinbase(upper->z, 2);

    if(!prec) prec = bc;

    shift = bc - prec;
    if(shift>0) {
        switch(rnd[0]) {
            case 'f':
                if(sign) {
                    mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                } else {
                    mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                }
                break;
            case 'c':
                if(sign) {
                    mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                } else {
                    mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                }
                break;
            case 'd':
                mpz_fdiv_q_2exp(upper->z, upper->z, shift);
                break;
            case 'u':
                mpz_cdiv_q_2exp(upper->z, upper->z, shift);
                break;
            case 'n':
            default:
                mpz_tdiv_r_2exp(lower->z, upper->z, shift);
                mpz_tdiv_q_2exp(upper->z, upper->z, shift);
                if(mpz_sgn(lower->z)) {
                    /* lower is not 0 so it must have at least 1 bit set */
                    if(mpz_sizeinbase(lower->z, 2)==shift) {
                        /* lower is >= 1/2 */
                        if(mpz_scan1(lower->z, 0)==shift-1) {
                            /* lower is exactly 1/2 */
                            if(mpz_odd_p(upper->z))
                                carry = 1;
                        } else {
                            carry = 1;
                        }
                    }
                }
                if(carry)
                    mpz_add_ui(upper->z, upper->z, 1);
        }
        if (!(tmp = Py2or3Int_FromLong(shift))) {
            Py_DECREF((PyObject*)upper);
            Py_DECREF((PyObject*)lower);
            return NULL;
        }
        if (!(newexp = PyNumber_Add(exp, tmp))) {
            Py_DECREF((PyObject*)man);
            Py_DECREF((PyObject*)upper);
            Py_DECREF((PyObject*)lower);
            Py_DECREF(tmp);
            return NULL;
        }
        Py_DECREF(tmp);
        bc = prec;
    } else {
        newexp = exp;
        Py_INCREF(newexp);
    }

    /* Strip trailing 0 bits. */
    if((zbits = mpz_scan1(upper->z, 0)))
        mpz_tdiv_q_2exp(upper->z, upper->z, zbits);

    if (!(tmp = Py2or3Int_FromLong(zbits))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF((PyObject*)upper);
        Py_DECREF((PyObject*)lower);
        Py_DECREF(newexp);
        return NULL;
    }
    if (!(newexp2 = PyNumber_Add(newexp, tmp))) {
        Py_DECREF((PyObject*)man);
        Py_DECREF((PyObject*)upper);
        Py_DECREF((PyObject*)lower);
        Py_DECREF(tmp);
        Py_DECREF(newexp);
        return NULL;
    }
    Py_DECREF(newexp);
    Py_DECREF(tmp);

    bc -= zbits;
    /* Check if one less than a power of 2 was rounded up. */
    if(!mpz_cmp_ui(upper->z, 1))
        bc = 1;

    Py_DECREF((PyObject*)lower);
    Py_DECREF((PyObject*)man);
    return mpmath_build_mpf(sign, upper, newexp2, bc);
}

/* Second version of helper functions for mpmath. See Issue 33 for details.
 */

static PyObject *
do_mpmath_trim(mpz_t xman, mpz_t xexp, long prec, char rnd) {
    PyObject *result = 0;
    PympzObject *rman = 0, *rexp = 0;
    long bc, shift, zbits, carry = 0;
    mpz_t lower;

    result = PyTuple_New(2);
    rman = Pympz_new();
    rexp = Pympz_new();
    if(!result || !rman || !rexp) {
        Py_XDECREF(result);
        Py_XDECREF((PyObject*)rman);
        Py_XDECREF((PyObject*)rexp);
        return NULL;
    }
    mpz_set(rman->z, xman);
    mpz_set(rexp->z, xexp);

    /* If the mantissa is 0, just return the canonical representation of 0. */
    if(!mpz_sgn(rman->z)) {
        mpz_set_ui(rexp->z, 0);
        goto return_result;
    }

    /* Remove trailing 0 bits and adjust exponenet. */
    if((zbits = mpz_scan1(rman->z, 0))) {
        mpz_tdiv_q_2exp(rman->z, rman->z, zbits);
        mpz_add_ui(rexp->z, rexp->z, zbits);
    }

    /* If prec is 0, return with trailing 0 bits removed. */
    if(prec == 0) goto return_result;

    bc = mpz_sizeinbase(rman->z, 2);

    /* If bc <= prec, just return. */
    if(bc <= prec) goto return_result;

    /* We need to round the mantissa. */
    shift = bc - prec;
    switch(rnd) {
        case 'f':
            mpz_fdiv_q_2exp(rman->z, rman->z, shift);
            break;
        case 'c':
            mpz_cdiv_q_2exp(rman->z, rman->z, shift);
            break;
        case 'd':
            if(mpz_sgn(rman->z) > 0) {
                mpz_fdiv_q_2exp(rman->z, rman->z, shift);
            } else {
                mpz_cdiv_q_2exp(rman->z, rman->z, shift);
            }
            break;
        case 'u':
            if(mpz_sgn(rman->z) > 0) {
                mpz_cdiv_q_2exp(rman->z, rman->z, shift);
            } else {
                mpz_fdiv_q_2exp(rman->z, rman->z, shift);
            }
            break;
        case 'n':
        default:
            mpz_inoc(lower);
            mpz_tdiv_r_2exp(lower, rman->z, shift);
            mpz_tdiv_q_2exp(rman->z, rman->z, shift);
            /* lower is not 0 so it must have at least 1 bit set */
            if(mpz_sizeinbase(lower, 2) == shift) {
                /* lower is >= 1/2 */
                if(mpz_scan1(lower, 0) == shift-1) {
                    /* lower is exactly 1/2 */
                    if(mpz_odd_p(rman->z))
                        carry = 1;
                } else {
                    carry = 1;
                }
            }
            mpz_cloc(lower);
            /* Add the carry bit. */
            if(carry) {
                if(mpz_sgn(rman->z) < 0) {
                    mpz_sub_ui(rman->z, rman->z, 1);
                } else {
                    mpz_add_ui(rman->z, rman->z, 1);
                }
            }
    }
    if((zbits = mpz_scan1(rman->z, 0))) {
        mpz_tdiv_q_2exp(rman->z, rman->z, zbits);
        mpz_add_ui(rexp->z, rexp->z, zbits);
    }
    mpz_add_ui(rexp->z, rexp->z, shift);

return_result:
    PyTuple_SET_ITEM(result, 0, (PyObject*)rman);
    PyTuple_SET_ITEM(result, 1, Pympz_To_Integer(rexp));
    Py_DECREF((PyObject*)rexp);
    return result;
}

static char doc_mpmath_trimg[]="\
_mpmath_trim(xman, xexp, prec, rounding):\n\
    Return (man, exp) by rounding xman*(2**xexp) to prec bits using the\n\
    specified rounding mode.\n\
";
static PyObject *
Pympz_mpmath_trim(PyObject *self, PyObject *args)
{
    PyObject *arg0 = 0, *arg1 = 0, *result;
    long prec = 0;
    const char *rnd = "d";

    switch(PyTuple_GET_SIZE(args)) {
        case 4:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 3));
        case 3:
            prec = clong_From_Integer(PyTuple_GET_ITEM(args, 2));
        case 2:
            arg1 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
        case 1:
            arg0 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    }
    if(!arg0 || !arg1 || (prec < 0) || PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "arguments mpz, mpz, long(>=0), char needed");
        Py_XDECREF(arg0);
        Py_XDECREF(arg1);
        return NULL;
    } else {
        result = do_mpmath_trim(Pympz_AS_MPZ(arg0), Pympz_AS_MPZ(arg1), prec, rnd[0]);
        Py_DECREF(arg0);
        Py_DECREF(arg1);
        return result;
    }
}

static char doc_mpmath_addg[]="\
_mpmath_add(xman, xexp, yman, yexp, prec, rounding):\n\
    Return (man, exp) by rounding xman*2**xexp + yman*2**yexp to prec\n\
    bits using the specified rounding mode.\n\
";
static PyObject *
Pympz_mpmath_add(PyObject *self, PyObject *args)
{
    PyObject *arg0 = 0, *arg1 = 0, *arg2 = 0, * arg3 = 0, *result, *temp;
    mpz_t man, exp, xbc_z, ybc_z, prec_z, offset_z, temp_z;
    long prec = 0, offset, zbits;
    const char *rnd = "d";

    switch(PyTuple_GET_SIZE(args)) {
        case 6:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 5));
        case 5:
            prec = clong_From_Integer(PyTuple_GET_ITEM(args, 4));
        case 4:
            arg3 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 3));
        case 3:
            arg2 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 2));
        case 2:
            arg1 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
        case 1:
            arg0 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    }
    if(!arg0 || !arg1 || !arg2 || !arg3 || (prec < 0) || PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "arguments mpz, mpz, mpz, mpz, long(>=0), char needed");
        Py_XDECREF(arg0);
        Py_XDECREF(arg1);
        Py_XDECREF(arg2);
        Py_XDECREF(arg3);
        return NULL;
    }
    /* Check if either argument is zero. */
    if(mpz_sgn(Pympz_AS_MPZ(arg0)) == 0) {
        result = do_mpmath_trim(Pympz_AS_MPZ(arg2), Pympz_AS_MPZ(arg3), prec, rnd[0]);
        goto return_result;
    }
    if(mpz_sgn(Pympz_AS_MPZ(arg2)) == 0) {
        result = do_mpmath_trim(Pympz_AS_MPZ(arg0), Pympz_AS_MPZ(arg1), prec, rnd[0]);
        goto return_result;
    }

    /* Remove trailing 0 bits. */
    if((zbits = mpz_scan1(Pympz_AS_MPZ(arg0), 0))) {
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(arg0), Pympz_AS_MPZ(arg0), zbits);
        mpz_add_ui(Pympz_AS_MPZ(arg1), Pympz_AS_MPZ(arg1), zbits);
    }
    if((zbits = mpz_scan1(Pympz_AS_MPZ(arg2), 0))) {
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(arg2), Pympz_AS_MPZ(arg2), zbits);
        mpz_add_ui(Pympz_AS_MPZ(arg3), Pympz_AS_MPZ(arg3), zbits);
    }

    /* Swap arguments to ensure arg1 >= arg3. Note: this does NOT imply that (arg0,arg1)
     * represents a number with a larger (in absolute terms) than (arg2,arg3).
     */
    if(mpz_cmp(Pympz_AS_MPZ(arg1), Pympz_AS_MPZ(arg3)) < 0) {
        temp = arg0;
        arg0 = arg2;
        arg2 = temp;
        temp = arg1;
        arg1 = arg3;
        arg3 = temp;
    }

    /* Get the bit lengths of the mantissas. */
    mpz_inoc(xbc_z);
    mpz_set_ui(xbc_z, mpz_sizeinbase(Pympz_AS_MPZ(arg0), 2));
    mpz_inoc(ybc_z);
    mpz_set_ui(ybc_z, mpz_sizeinbase(Pympz_AS_MPZ(arg2), 2));

    /* Calculate the amount arg0 must be shifted to line up with arg2. */
    mpz_inoc(offset_z);
    mpz_set(offset_z, Pympz_AS_MPZ(arg1));
    mpz_sub(offset_z, offset_z, Pympz_AS_MPZ(arg3));

    /* xbc_z now has the effective bitlength. It assumes the mantissa is
     * shifted.
     */
    mpz_add(xbc_z, xbc_z, offset_z);

    /* ybc_z is incremented by 2. If offset_z is greater than ybc_z, then
     * we only need to perturb the result.
     */
    mpz_add_ui(ybc_z, ybc_z, 2);

    mpz_inoc(prec_z);
    mpz_set_ui(prec_z, prec);
    mpz_add_ui(prec_z, prec_z, 3);

    mpz_inoc(temp_z);
    mpz_sub(temp_z, offset_z, ybc_z);

    mpz_inoc(man);
    mpz_inoc(exp);

    if(prec && mpz_cmp(temp_z, prec_z) > 0) {
        /* only need to perturb the result */
        if(!mpz_fits_slong_p(offset_z)) {
            PyErr_SetString(PyExc_ValueError, "offset too large");
            result = NULL;
            goto return_result;
        } else {
            offset = mpz_get_si(offset_z);
        }
        mpz_set(man, Pympz_AS_MPZ(arg0));
        mpz_mul_2exp(man, man, offset + 3);
        if(mpz_sgn(Pympz_AS_MPZ(arg2)) > 0) {
            mpz_add_ui(man, man, 1);
        } else {
            mpz_sub_ui(man, man, 1);
        }
        mpz_set(exp, Pympz_AS_MPZ(arg1));
        mpz_sub_ui(exp, exp, offset + 3);
        result = do_mpmath_trim(man, exp, prec, rnd[0]);
    } else {
        /* do a full addition */
        if(!mpz_fits_slong_p(offset_z)) {
            PyErr_SetString(PyExc_ValueError, "offset too large");
            result = NULL;
            goto return_result;
        } else {
            offset = mpz_get_si(offset_z);
        }
        mpz_set(man, Pympz_AS_MPZ(arg0));
        if(offset)
            mpz_mul_2exp(man, man, offset);
        mpz_add(man, man, Pympz_AS_MPZ(arg2));
        result = do_mpmath_trim(man, Pympz_AS_MPZ(arg3), prec, rnd[0]);
    }
    mpz_cloc(exp);
    mpz_cloc(man);
    mpz_cloc(offset_z);
    mpz_cloc(temp_z);
    mpz_cloc(prec_z);
    mpz_cloc(xbc_z);
    mpz_cloc(ybc_z);
return_result:
    Py_DECREF(arg0);
    Py_DECREF(arg1);
    Py_DECREF(arg2);
    Py_DECREF(arg3);
    return result;
}

static char doc_mpmath_multg[]="\
_mpmath_mult(xman, xexp, yman, yexp, prec, rounding):\n\
    Return (man, exp) by rounding xman*2**xexp * yman*2**yexp to prec\n\
    bits using the specified rounding mode.\n\
";
static PyObject *
Pympz_mpmath_mult(PyObject *self, PyObject *args)
{
    PyObject *arg0 = 0, *arg1 = 0, *arg2 = 0, * arg3 = 0, *result;
    mpz_t man, exp;
    long prec = 0;
    const char *rnd = "d";

    switch(PyTuple_GET_SIZE(args)) {
        case 6:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 5));
        case 5:
            prec = clong_From_Integer(PyTuple_GET_ITEM(args, 4));
        case 4:
            arg3 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 3));
        case 3:
            arg2 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 2));
        case 2:
            arg1 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
        case 1:
            arg0 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    }
    if(!arg0 || !arg1 || !arg2 || !arg3 || (prec < 0) || PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "arguments mpz, mpz, mpz, mpz, long(>=0), char needed");
        Py_XDECREF(arg0);
        Py_XDECREF(arg1);
        Py_XDECREF(arg2);
        Py_XDECREF(arg3);
        return NULL;
    } else {
        mpz_inoc(man);
        mpz_inoc(exp);
        mpz_mul(man, Pympz_AS_MPZ(arg0), Pympz_AS_MPZ(arg2));
        mpz_add(exp, Pympz_AS_MPZ(arg1), Pympz_AS_MPZ(arg3));
        result = do_mpmath_trim(man, exp, prec, rnd[0]);
        mpz_cloc(man);
        mpz_cloc(exp);
        Py_DECREF(arg0);
        Py_DECREF(arg1);
        Py_DECREF(arg2);
        Py_DECREF(arg3);
        return result;
    }
}

static char doc_mpmath_divg[]="\
_mpmath_div(xman, xexp, yman, yexp, prec, rounding):\n\
    Return (man, exp) by rounding xman*2**xexp / yman*2**yexp to prec\n\
    bits using the specified rounding mode.\n\
";
static PyObject *
Pympz_mpmath_div(PyObject *self, PyObject *args)
{
    PyObject *arg0 = 0, *arg1 = 0, *arg2 = 0, * arg3 = 0, *result;
    mpz_t quot, rem, exp, delta_z;
    long prec = 0, delta, zbits;
    const char *rnd = "d";

    switch(PyTuple_GET_SIZE(args)) {
        case 6:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 5));
        case 5:
            prec = clong_From_Integer(PyTuple_GET_ITEM(args, 4));
        case 4:
            arg3 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 3));
        case 3:
            arg2 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 2));
        case 2:
            arg1 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
        case 1:
            arg0 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    }
    if(!arg0 || !arg1 || !arg2 || !arg3 || (prec < 1) || PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "arguments mpz, mpz, mpz, mpz, long(>=1), char needed");
        Py_XDECREF(arg0);
        Py_XDECREF(arg1);
        Py_XDECREF(arg2);
        Py_XDECREF(arg3);
        return NULL;
    }
    /* Check if either argument is zero. */
    if(mpz_sgn(Pympz_AS_MPZ(arg2)) == 0) {
        PyErr_SetString(PyExc_ZeroDivisionError, "mpmath division by 0");
        result = NULL;
        goto return_result;
    }

    if(mpz_sgn(Pympz_AS_MPZ(arg0)) == 0) {
        result = do_mpmath_trim(Pympz_AS_MPZ(arg0), Pympz_AS_MPZ(arg1), prec, rnd[0]);
        goto return_result;
    }

    /* Remove trailing 0 bits. */
    if((zbits = mpz_scan1(Pympz_AS_MPZ(arg0), 0))) {
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(arg0), Pympz_AS_MPZ(arg0), zbits);
        mpz_add_ui(Pympz_AS_MPZ(arg1), Pympz_AS_MPZ(arg1), zbits);
    }
    if((zbits = mpz_scan1(Pympz_AS_MPZ(arg2), 0))) {
        mpz_tdiv_q_2exp(Pympz_AS_MPZ(arg2), Pympz_AS_MPZ(arg2), zbits);
        mpz_add_ui(Pympz_AS_MPZ(arg3), Pympz_AS_MPZ(arg3), zbits);
    }

    mpz_inoc(delta_z);
    mpz_set_ui(delta_z, prec);
    mpz_sub_ui(delta_z, delta_z, mpz_sizeinbase(Pympz_AS_MPZ(arg0), 2));
    mpz_add_ui(delta_z, delta_z, mpz_sizeinbase(Pympz_AS_MPZ(arg2), 2));
    mpz_add_ui(delta_z, delta_z, 5);
    if(mpz_cmp_ui(delta_z, 5) < 0) {
        mpz_set_ui(delta_z, 5);
    }

    mpz_inoc(quot);
    mpz_inoc(rem);
    mpz_inoc(exp);
    if(!mpz_fits_slong_p(delta_z)) {
        PyErr_SetString(PyExc_ValueError, "delta too large");
        result = NULL;
        goto return_result;
    } else {
        delta = mpz_get_si(delta_z);
    }
    mpz_set(quot, Pympz_AS_MPZ(arg0));
    mpz_mul_2exp(quot, quot, delta);
    mpz_tdiv_qr(quot, rem, quot, Pympz_AS_MPZ(arg2));

    if(mpz_sgn(rem)) {
        mpz_mul_2exp(quot, quot, 1);
        if(mpz_sgn(quot) < 0) {
            mpz_sub_ui(quot, quot, 1);
        } else {
            mpz_add_ui(quot, quot, 1);
        }
        mpz_add_ui(delta_z, delta_z, 1);
    }
    mpz_set(exp, Pympz_AS_MPZ(arg1));
    mpz_sub(exp, exp, Pympz_AS_MPZ(arg3));
    mpz_sub(exp, exp, delta_z);
    result = do_mpmath_trim(quot, exp, prec, rnd[0]);

    mpz_cloc(quot);
    mpz_cloc(rem);
    mpz_cloc(exp);
    mpz_cloc(delta_z);
return_result:
    Py_DECREF(arg0);
    Py_DECREF(arg1);
    Py_DECREF(arg2);
    Py_DECREF(arg3);
    return result;
}

static char doc_mpmath_sqrtg[]="\
_mpmath_sqrt(man, exp, prec, rounding):\n\
    Return (man, exp) by rounding square_root(xman*2**xexp)) to prec\n\
    bits using the specified rounding mode.\n\
";
static PyObject *
Pympz_mpmath_sqrt(PyObject *self, PyObject *args)
{
    PyObject *arg0 = 0, *arg1 = 0, *result;
    mpz_t man, exp, rem;
    long prec = 0, zbits;
    unsigned long shift, temp;
    const char *rnd = "d";

    switch(PyTuple_GET_SIZE(args)) {
        case 4:
            rnd = Py2or3String_AsString(PyTuple_GET_ITEM(args, 3));
        case 3:
            prec = clong_From_Integer(PyTuple_GET_ITEM(args, 2));
        case 2:
            arg1 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 1));
        case 1:
            arg0 = (PyObject*)Pympz_From_Integer(PyTuple_GET_ITEM(args, 0));
    }
    if(!arg0 || !arg1 || (prec < 1) || PyErr_Occurred()) {
        PyErr_SetString(PyExc_TypeError, "arguments mpz, mpz, long(>=1), char needed");
        Py_XDECREF(arg0);
        Py_XDECREF(arg1);
        return NULL;
    }

    mpz_inoc(man);
    mpz_inoc(exp);
    mpz_inoc(rem);
    mpz_set(man, Pympz_AS_MPZ(arg0));
    mpz_set(exp, Pympz_AS_MPZ(arg1));

    if(mpz_sgn(man) < 0) {
        PyErr_SetString(PyExc_ValueError, "square root of a negative number");
        result = NULL;
        goto return_result;
    }
    if(mpz_sgn(man) == 0) {
        result = do_mpmath_trim(man, exp, prec, rnd[0]);
        goto return_result;
    }
    if((zbits = mpz_scan1(man, 0))) {
        mpz_tdiv_q_2exp(man, man, zbits);
        mpz_add_ui(exp, exp, zbits);
    }
    if(mpz_odd_p(exp)) {
        mpz_sub_ui(exp, exp, 1);
        mpz_mul_2exp(man, man, 1);
    } else if(!mpz_cmp_ui(man, 1)) {
        /* Handle even powers of 2. */
        mpz_tdiv_q_2exp(exp, exp, 1);
        result = do_mpmath_trim(man, exp, prec, rnd[0]);
        goto return_result;
    }

    shift = (2 * prec) + 4;
    temp = mpz_sizeinbase(man, 2);
    if(temp >= shift) {
        shift = 4;
    } else {
        shift -= temp;
    }
    if(shift < 4) shift = 4;
    shift += shift & 1;
    mpz_mul_2exp(man, man, shift);
    if((rnd[0] == 'f') || (rnd[0] == 'd')) {
        mpz_sqrt(man, man);
    } else {
        mpz_sqrtrem(man, rem, man);
        if(mpz_sgn(rem)) {
            mpz_mul_2exp(man, man, 1);
            mpz_add_ui(man, man, 1);
            shift += 2;
        }
    }

    mpz_sub_ui(exp, exp, shift);
    mpz_tdiv_q_2exp(exp, exp, 1);
    result = do_mpmath_trim(man, exp, prec, rnd[0]);

return_result:
    mpz_cloc(man);
    mpz_cloc(exp);
    mpz_cloc(rem);
    Py_DECREF(arg0);
    Py_DECREF(arg1);
    return result;
}

