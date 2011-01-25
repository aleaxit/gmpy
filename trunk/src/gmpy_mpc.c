PyDoc_STRVAR(doc_g_mpc_get_mpc_round,
"get_mpc_round() -> integer\n\n"
"Return the rounding mode for complex arithmetic. The value will\n"
"be one of the MPC_RNDxy constants. The values for 'x' and 'y'\n"
"indicate the rounding mode for the real and imaginary components\n"
"respectively. The values are interpreted as:\n"
"   N -> RoundToNearest\n"
"   Z -> RoundToZero\n"
"   U -> RoundUp\n"
"   D -> RoundDown\n");

static PyObject *
Pympc_get_mpc_round(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", context->now.mpc_round);
}

PyDoc_STRVAR(doc_g_mpc_set_mpc_round,
"set_mpc_round((n,n))\n\n"
"Set the rounding mode for complex arithmetic. The value will\n"
"be one of the MPC_RNDxy constants. The values for 'x' and 'y'\n"
"indicate the rounding mode for the real and imaginary components\n"
"respectively. The values are interpreted as:\n"
"   N -> RoundToNearest\n"
"   Z -> RoundToZero\n"
"   U -> RoundUp\n"
"   D -> RoundDown\n");

static PyObject *
Pympc_set_mpc_round(PyObject *self, PyObject *args)
{
    int temp = context->now.mpc_round;

    if (!PyArg_ParseTuple(args, "i", &temp))
        return NULL;

    if (!Pymisc_verify_mpc_round(temp)) {
        VALUE_ERROR("invalid rounding mode for complex arithmetic.");
        return NULL;
    }

    context->now.mpc_round = temp;
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
                         (Py_ssize_t)(context->now.mpc_rprec),
                         (Py_ssize_t)(context->now.mpc_iprec));
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

    if (!PyArg_ParseTuple(args, "(nn)", &rprec, &iprec))
        return NULL;

    if (!Pymisc_verify_mpc_precision(rprec, iprec)) {
        VALUE_ERROR("invalid precision for complex arithmetic");
        return NULL;
    }

    context->now.mpc_rprec = rprec;
    context->now.mpc_iprec = iprec;
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
