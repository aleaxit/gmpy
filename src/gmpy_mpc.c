PyDoc_STRVAR(doc_g_mpc_get_mpc_round,
"get_mpc_round() -> (integer, integer)\n\n"
"Return the rounding mode for 'mpc' arithmetic. The first value is"
"the rounding mode for the real portion. The second value is the"
"rounding mode for the imaginary portion. Rounding mode can be"
"RoundToNearest, RoundToZero, RoundUp, or RoundDown.");

static PyObject *
Pympc_get_mpc_round(PyObject *self, PyObject *args)
{
    return Py_BuildValue("(ii)",
                         MPC_RND_RE(global.mpc_round),
                         MPC_RND_IM(global.mpc_round));
}

PyDoc_STRVAR(doc_g_mpc_set_mpc_round,
"set_mpc_round((n,n))\n\n"
"Set the rounding mode for complex arithmetic. The first value is"
"the rounding more for the real portion. The second value is the"
"rounding mode for the imaginary portion. Valid rounding modes are"
"RoundToNearest, RoundToZero, RoundUp, or RoundDown.");

static PyObject *
Pympc_set_mpc_round(PyObject *self, PyObject *args)
{
    int rmode, imode;

    if(!PyArg_ParseTuple(args, "(ii)", &rmode, &imode))
        return NULL;

    if (rmode == MPFR_RNDN && imode == MPFR_RNDN)
        global.mpc_round = MPC_RNDNN;
    else if (rmode == MPFR_RNDN && imode == MPFR_RNDZ)
        global.mpc_round = MPC_RNDNZ;
    else if (rmode == MPFR_RNDN && imode == MPFR_RNDU)
        global.mpc_round = MPC_RNDNU;
    else if (rmode == MPFR_RNDN && imode == MPFR_RNDD)
        global.mpc_round = MPC_RNDND;

    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDN)
        global.mpc_round = MPC_RNDZN;
    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDZ)
        global.mpc_round = MPC_RNDZZ;
    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDU)
        global.mpc_round = MPC_RNDZU;
    else if (rmode == MPFR_RNDZ && imode == MPFR_RNDD)
        global.mpc_round = MPC_RNDZD;

    else if (rmode == MPFR_RNDU && imode == MPFR_RNDN)
        global.mpc_round = MPC_RNDUN;
    else if (rmode == MPFR_RNDU && imode == MPFR_RNDZ)
        global.mpc_round = MPC_RNDUZ;
    else if (rmode == MPFR_RNDU && imode == MPFR_RNDU)
        global.mpc_round = MPC_RNDUU;
    else if (rmode == MPFR_RNDU && imode == MPFR_RNDD)
        global.mpc_round = MPC_RNDUD;

    else if (rmode == MPFR_RNDD && imode == MPFR_RNDN)
        global.mpc_round = MPC_RNDDN;
    else if (rmode == MPFR_RNDD && imode == MPFR_RNDZ)
        global.mpc_round = MPC_RNDDZ;
    else if (rmode == MPFR_RNDD && imode == MPFR_RNDU)
        global.mpc_round = MPC_RNDDU;
    else if (rmode == MPFR_RNDD && imode == MPFR_RNDD)
        global.mpc_round = MPC_RNDDD;

    else {
        VALUE_ERROR("invalid rounding mode");
        return NULL;
    }
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g_mpc_get_mpc_precision,
"get_mpc_precision() -> (integer, integer)\n\n"
"Return the precision for 'mpc' arithmetic. The first value is the"
"precision for the real portion. The second value is the precision"
"for the imaginary portion.");

static PyObject *
Pympc_get_mpc_precision(PyObject *self, PyObject *args)
{
    return Py_BuildValue("(nn)",
                         (Py_ssize_t)(global.mpc_rprec),
                         (Py_ssize_t)(global.mpc_iprec));
}

PyDoc_STRVAR(doc_g_mpc_set_mpc_precision,
"set_mpc_precision((n,n))\n\n"
"Set the precision for complex arithmetic. The first value is the"
"precision for the real portion. The second value is the precision"
"for the imaginary portion.");

static PyObject *
Pympc_set_mpc_precision(PyObject *self, PyObject *args)
{
    Py_ssize_t rprec, iprec;

    if(!PyArg_ParseTuple(args, "(nn)", &rprec, &iprec))
        return NULL;
    if (rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
        iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }
    global.mpc_rprec = rprec;
    global.mpc_iprec = iprec;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_g_mpc_get_mpc_status,
"get_mpc_status() -> (integer, integer)\n\n"
"Return the ternary result code from the most recent MPC operation.\n"
"The values are for the real and imaginary components, respectively.\n"
"If the ternary value is 0, the result of the operation is exact.\n"
"If the ternary value is > 0, the result of the operation is greater\n"
"than the exact result. If the ternary value < 0, then the result\n"
"of the operation is less than the exact result.");

static PyObject *
Pympc_get_mpc_status(PyObject *self, PyObject *args)
{
    return Py_BuildValue("(ii)",
                         MPC_INEX_RE(global.mpc_rc),
                         MPC_INEX_IM(global.mpc_rc));
}
