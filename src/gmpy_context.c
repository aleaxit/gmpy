/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_context.c                                                          *
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

/* A GMPyContextObject contains an instance of the C struct gmpy_context
 * and a PyObject* used to reference the enclosing instance when used as a
 * context manager in Python.
 *
 * gmpy2 uses a global pointer "context" to refer to the active
 * GMPyContextObject.
 *
 * WARNING: The context manager is not thread-safe. This may be fixed in a
 *          future version.
 */

/* Create and delete Context objects. */

static PyObject *
GMPyContext_new(void)
{
    GMPyContextObject *result;

    if ((result = PyObject_New(GMPyContextObject, &GMPyContext_Type))) {
        result->ctx.mpfr_prec = DBL_MANT_DIG;
        result->ctx.mpfr_round = MPFR_RNDN;
        result->ctx.emax = MPFR_EMAX_DEFAULT;
        result->ctx.emin = MPFR_EMIN_DEFAULT;
        result->ctx.subnormalize = 0;
        result->ctx.underflow = 0;
        result->ctx.overflow = 0;
        result->ctx.inexact = 0;
        result->ctx.invalid = 0;
        result->ctx.erange = 0;
        result->ctx.divzero = 0;
        result->ctx.trap_underflow = 0;
        result->ctx.trap_overflow = 0;
        result->ctx.trap_inexact = 0;
        result->ctx.trap_invalid = 0;
        result->ctx.trap_erange = 0;
        result->ctx.trap_divzero = 0;
        result->ctx.trap_expbound = 0;

#ifdef WITHMPC
        result->ctx.real_prec = -1;
        result->ctx.imag_prec = -1;
        result->ctx.real_round = -1;
        result->ctx.imag_round = -1;
        result->ctx.allow_complex = 0;
#endif
    }
    return (PyObject*)result;
};

static void
GMPyContext_dealloc(GMPyContextObject *self)
{
    PyObject_Del(self);
};

PyDoc_STRVAR(doc_context_ieee,
"ieee(bitwidth) -> context\n\n"
"Return a new context corresponding to a standard IEEE floating point\n"
"format. The currently supported precisions are 32, 64, and 128 bits.");

static PyObject *
GMPyContext_ieee(PyObject *self, PyObject *other)
{
    long bitwidth;
    GMPyContextObject *result;

    bitwidth = PyIntOrLong_AsLong(other);
    if (bitwidth == -1 && PyErr_Occurred()) {
        TYPE_ERROR("ieee() requires 'int' argument");
        return NULL;
    }

    if (bitwidth == 32) {
        result = (GMPyContextObject*)GMPyContext_new();
        if (result) {
            result->ctx.subnormalize = 1;
            result->ctx.mpfr_prec = 24;
            result->ctx.emax = 128;
            result->ctx.emin = -148;
        }
        return (PyObject*)result;
    }
    else if (bitwidth == 64) {
        result = (GMPyContextObject*)GMPyContext_new();
        if (result) {
            result->ctx.subnormalize = 1;
            result->ctx.mpfr_prec = 53;
            result->ctx.emax = 1024;
            result->ctx.emin = -1073;
        }
        return (PyObject*)result;
    }
    else if (bitwidth == 128) {
        result = (GMPyContextObject*)GMPyContext_new();
        if (result) {
            result->ctx.subnormalize = 1;
            result->ctx.mpfr_prec = 113;
            result->ctx.emax = 16384;
            result->ctx.emin = -16493;
        }
        return (PyObject*)result;
    }
    else {
        VALUE_ERROR("bitwidth must be 32, 64, or 128");
        return NULL;
    }
}

/* Create and delete ContextManager objects. */

static PyObject *
GMPyContextManager_new(void)
{
    return (PyObject*)PyObject_New(GMPyContextManagerObject,
                                   &GMPyContextManager_Type);
};

static void
GMPyContextManager_dealloc(GMPyContextManagerObject *self)
{
    PyObject_Del(self);
};

/* Helper function to convert to convert a rounding mode to a string. */

static PyObject *
_round_to_name(int val)
{
    if (val == MPFR_RNDN) return Py2or3String_FromString("RoundToNearest");
    if (val == MPFR_RNDZ) return Py2or3String_FromString("RoundToZero");
    if (val == MPFR_RNDU) return Py2or3String_FromString("RoundUp");
    if (val == MPFR_RNDD) return Py2or3String_FromString("RoundDown");
    if (val == MPFR_RNDA) return Py2or3String_FromString("RoundAwayZero");
    if (val == GMPY_DEFAULT) return Py2or3String_FromString("Default");
    return NULL;
};

static PyObject *
GMPyContext_repr(GMPyContextObject *self)
{
    PyObject *format;
    PyObject *tuple;
    PyObject *result = NULL;
    int i = 0;

#ifdef WITHMPC
    tuple = PyTuple_New(23);
#else
    tuple = PyTuple_New(18);
#endif
    if (!tuple) return NULL;

#ifdef WITHMPC
    format = Py2or3String_FromString(
            "context(precision=%s, real_prec=%s, imag_prec=%s,\n"
            "        round=%s, real_round=%s, imag_round=%s,\n"
            "        emax=%s, emin=%s,\n"
            "        subnormalize=%s,\n"
            "        trap_underflow=%s, underflow=%s,\n"
            "        trap_overflow=%s, overflow=%s,\n"
            "        trap_inexact=%s, inexact=%s,\n"
            "        trap_invalid=%s, invalid=%s,\n"
            "        trap_erange=%s, erange=%s,\n"
            "        trap_divzero=%s, divzero=%s,\n"
            "        trap_expbound=%s,\n"
            "        allow_complex=%s)"
            );
#else
    format = Py2or3String_FromString(
            "context(precision=%s,\n"
            "        round=%s,\n"
            "        emax=%s, emin=%s,\n"
            "        subnormalize=%s,\n"
            "        trap_underflow=%s, underflow=%s,\n"
            "        trap_overflow=%s, overflow=%s,\n"
            "        trap_inexact=%s, inexact=%s,\n"
            "        trap_invalid=%s, invalid=%s,\n"
            "        trap_erange=%s, erange=%s,\n"
            "        trap_divzero=%s, divzero=%s,\n"
            "        trap_expbound=%s)"
            );
#endif
    if (!format) {
        Py_DECREF(tuple);
        return NULL;
    }

    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.mpfr_prec));
#ifdef WITHMPC
    if (self->ctx.real_prec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, i++, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.real_prec));
    if (self->ctx.imag_prec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, i++, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.imag_prec));
#endif
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->ctx.mpfr_round));
#ifdef WITHMPC
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->ctx.real_round));
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->ctx.imag_round));
#endif
    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.emax));
    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.emin));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.subnormalize));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_underflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.underflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_overflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.overflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_inexact));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.inexact));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_invalid));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.invalid));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_erange));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.erange));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_divzero));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.divzero));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.trap_expbound));
#ifdef WITHMPC
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.allow_complex));
#endif

    if (!PyErr_Occurred())
        result = Py2or3String_Format(format, tuple);
    else
        SYSTEM_ERROR("internal error in GMPyContext_repr");

    Py_DECREF(format);
    Py_DECREF(tuple);
    return result;
};

static PyObject *
GMPyContextManager_repr(GMPyContextManagerObject *self)
{
    return Py_BuildValue("s", "<gmpy2.ContextManagerObject>");
}

PyDoc_STRVAR(doc_get_context,
"get_context([keywords]) -> gmpy2 context\n\n"
"Return a reference to the current context. Keyword arguments will\n"
"modify the current context. Changing attributes of the returned object\n"
"will modify the current context unless another context is activated by\n"
"set_context(). When used in 'with get_context(...):', the context as\n"
"modified by any keyword arguments will be restored at the end of the\n"
"with block. This differs from local_context() which will restore the\n"
"context before any changes by keyword arguments.");

/* Should parse keyword arguments. */

static PyObject *
GMPyContext_get_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    gmpy_context old;

#ifdef WITHMPC
    static char *kwlist[] = {
        "precision", "real_prec", "imag_prec", "round",
        "real_round", "imag_round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero",
        "trap_expbound", "allow_complex", NULL };
#else
    static char *kwlist[] = {
        "precision", "round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero", "trap_expbound",
        NULL };
#endif

    if (PyTuple_GET_SIZE(args)) {
        VALUE_ERROR("get_context() only supports keyword arguments");
        return NULL;
    }

    /* Temporarily save the options of the global context manager in case
     * there is an error while processing the arguments.
     */

    old = context->ctx;

#ifdef WITHMPC
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|llliiilliiiiiiiii", kwlist,
            &context->ctx.mpfr_prec,
            &context->ctx.real_prec,
            &context->ctx.imag_prec,
            &context->ctx.mpfr_round,
            &context->ctx.real_round,
            &context->ctx.imag_round,
            &context->ctx.emax,
            &context->ctx.emin,
            &context->ctx.subnormalize,
            &context->ctx.trap_underflow,
            &context->ctx.trap_overflow,
            &context->ctx.trap_inexact,
            &context->ctx.trap_invalid,
            &context->ctx.trap_erange,
            &context->ctx.trap_divzero,
            &context->ctx.trap_expbound,
            &context->ctx.allow_complex))) {
#else
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|lilliiiiiiii", kwlist,
            &context->ctx.mpfr_prec,
            &context->ctx.mpfr_round,
            &context->ctx.emax,
            &context->ctx.emin,
            &context->ctx.subnormalize,
            &context->ctx.trap_underflow,
            &context->ctx.trap_overflow,
            &context->ctx.trap_inexact,
            &context->ctx.trap_invalid,
            &context->ctx.trap_erange,
            &context->ctx.trap_divzero,
            &context->ctx.trap_expbound))) {
#endif
        VALUE_ERROR("invalid keyword arguments in get_context()");
        return NULL;
    }

    /* Sanity check for values. */
    if (context->ctx.mpfr_prec < MPFR_PREC_MIN ||
        context->ctx.mpfr_prec > MPFR_PREC_MAX) {
        context->ctx = old;
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }

#ifdef WITHMPC
    if (!(context->ctx.real_prec == GMPY_DEFAULT ||
        (context->ctx.real_prec >= MPFR_PREC_MIN &&
        context->ctx.real_prec <= MPFR_PREC_MAX))) {
        context->ctx = old;
        VALUE_ERROR("invalid value for real_prec");
        return NULL;
    }
    if (!(context->ctx.imag_prec == GMPY_DEFAULT ||
        (context->ctx.imag_prec >= MPFR_PREC_MIN &&
        context->ctx.imag_prec <= MPFR_PREC_MAX))) {
        context->ctx = old;
        VALUE_ERROR("invalid value for imag_prec");
        return NULL;
    }
#endif

    if (!(context->ctx.mpfr_round == MPFR_RNDN ||
        context->ctx.mpfr_round == MPFR_RNDZ ||
        context->ctx.mpfr_round == MPFR_RNDU ||
        context->ctx.mpfr_round == MPFR_RNDD ||
        context->ctx.mpfr_round == MPFR_RNDA)) {
        context->ctx = old;
        VALUE_ERROR("invalid value for round");
        return NULL;
    }

#ifdef WITHMPC
    if (context->ctx.mpfr_round == MPFR_RNDA) {
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
         * to MPFR_RNDN.
         */
        context->ctx.real_round = MPFR_RNDN;
        context->ctx.imag_round = MPFR_RNDN;
    }
    if (!(context->ctx.real_round == MPFR_RNDN ||
        context->ctx.real_round == MPFR_RNDZ ||
        context->ctx.real_round == MPFR_RNDU ||
        context->ctx.real_round == MPFR_RNDD ||
        context->ctx.real_round == GMPY_DEFAULT)) {
        context->ctx = old;
        VALUE_ERROR("invalid value for real_round");
        return NULL;
    }
    if (!(context->ctx.imag_round == MPFR_RNDN ||
        context->ctx.imag_round == MPFR_RNDZ ||
        context->ctx.imag_round == MPFR_RNDU ||
        context->ctx.imag_round == MPFR_RNDD ||
        context->ctx.imag_round == GMPY_DEFAULT)) {
        context->ctx = old;
        VALUE_ERROR("invalid value for imag_round");
        return NULL;
    }
#endif

    if (!(context->ctx.emin < 0 && context->ctx.emax > 0)) {
        VALUE_ERROR("invalid values for emin and/or emax");
        context->ctx = old;
        return NULL;
    }

    if (mpfr_set_emin(context->ctx.emin)) {
        VALUE_ERROR("invalid value for emin");
        context->ctx = old;
        return NULL;
    }
    if (mpfr_set_emax(context->ctx.emax)) {
        VALUE_ERROR("invalid value for emax");
        context->ctx = old;
        return NULL;
    }

    /* The values in the global context manager have been changed. Also return
     * another reference to that context manager.
     */

    Py_INCREF((PyObject*)context);
    return (PyObject*)context;
}

PyDoc_STRVAR(doc_local_context,
"local_context([context[,keywords]]) -> context manager\n\n"
"Create a context manager object that will restore either the context\n"
"that is passed as an argument, or the context that is active when\n"
"called when the 'with ...' block terminates. Keyword arguments are\n"
"supported and will modify the new context.");

static PyObject *
GMPyContext_local_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPyContextManagerObject *result;
    PyObject *local_args = args;
    int arg_context = 0;

#ifdef WITHMPC
    static char *kwlist[] = {
        "precision", "real_prec", "imag_prec", "round",
        "real_round", "imag_round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero",
        "trap_expbound", "allow_complex", NULL };
#else
    static char *kwlist[] = {
        "precision", "round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero", "trap_expbound",
        NULL };
#endif

    if (PyTuple_GET_SIZE(args) == 1 && GMPyContext_Check(PyTuple_GET_ITEM(args, 0))) {
        arg_context = 1;
        if (!(local_args = PyTuple_New(0)))
            return NULL;
    }
    else if (PyTuple_GET_SIZE(args)) {
        VALUE_ERROR("local_context() only supports [context[,keyword]] arguments");
        return NULL;
    }

    if (!(result = (GMPyContextManagerObject*)GMPyContextManager_new()))
        return NULL;

    if (arg_context) {
        result->new_ctx = ((GMPyContextObject*)PyTuple_GET_ITEM(args, 0))->ctx;
    }
    else {
        result->new_ctx = context->ctx;
    }
    result->old_ctx = context->ctx;

#ifdef WITHMPC
    if (!(PyArg_ParseTupleAndKeywords(local_args, kwargs,
            "|llliiilliiiiiiiii", kwlist,
            &result->new_ctx.mpfr_prec,
            &result->new_ctx.real_prec,
            &result->new_ctx.imag_prec,
            &result->new_ctx.mpfr_round,
            &result->new_ctx.real_round,
            &result->new_ctx.imag_round,
            &result->new_ctx.emax,
            &result->new_ctx.emin,
            &result->new_ctx.subnormalize,
            &result->new_ctx.trap_underflow,
            &result->new_ctx.trap_overflow,
            &result->new_ctx.trap_inexact,
            &result->new_ctx.trap_invalid,
            &result->new_ctx.trap_erange,
            &result->new_ctx.trap_divzero,
            &result->new_ctx.trap_expbound,
            &result->new_ctx.allow_complex))) {
#else
    if (!(PyArg_ParseTupleAndKeywords(local_args, kwargs,
            "|lilliiiiiiii", kwlist,
            &result->new_ctx.mpfr_prec,
            &result->new_ctx.mpfr_round,
            &result->new_ctx.emax,
            &result->new_ctx.emin,
            &result->new_ctx.subnormalize,
            &result->new_ctx.trap_underflow,
            &result->new_ctx.trap_overflow,
            &result->new_ctx.trap_inexact,
            &result->new_ctx.trap_invalid,
            &result->new_ctx.trap_erange,
            &result->new_ctx.trap_divzero,
            &result->new_ctx.trap_expbound))) {
#endif
        VALUE_ERROR("invalid keyword arguments in local_context()");
        goto error;
    }

    /* Sanity check for values. */
    if (result->new_ctx.mpfr_prec < MPFR_PREC_MIN ||
        result->new_ctx.mpfr_prec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        goto error;
    }

#ifdef WITHMPC
    if (!(result->new_ctx.real_prec == GMPY_DEFAULT ||
        (result->new_ctx.real_prec >= MPFR_PREC_MIN &&
        result->new_ctx.real_prec <= MPFR_PREC_MAX))) {
        VALUE_ERROR("invalid value for real_prec");
        goto error;
    }
    if (!(result->new_ctx.imag_prec == GMPY_DEFAULT ||
        (result->new_ctx.imag_prec >= MPFR_PREC_MIN &&
        result->new_ctx.imag_prec <= MPFR_PREC_MAX))) {
        VALUE_ERROR("invalid value for imag_prec");
        goto error;
    }
#endif

    if (!(result->new_ctx.mpfr_round == MPFR_RNDN ||
        result->new_ctx.mpfr_round == MPFR_RNDZ ||
        result->new_ctx.mpfr_round == MPFR_RNDU ||
        result->new_ctx.mpfr_round == MPFR_RNDD ||
        result->new_ctx.mpfr_round == MPFR_RNDA)) {
        VALUE_ERROR("invalid value for round");
        goto error;
    }

#ifdef WITHMPC
    if (result->new_ctx.mpfr_round == MPFR_RNDA) {
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
         * to MPFR_RNDN.
         */
        result->new_ctx.real_round = MPFR_RNDN;
        result->new_ctx.imag_round = MPFR_RNDN;
    }
    if (!(result->new_ctx.real_round == MPFR_RNDN ||
        result->new_ctx.real_round == MPFR_RNDZ ||
        result->new_ctx.real_round == MPFR_RNDU ||
        result->new_ctx.real_round == MPFR_RNDD ||
        result->new_ctx.real_round == GMPY_DEFAULT)) {
        VALUE_ERROR("invalid value for real_round");
        goto error;
    }
    if (!(result->new_ctx.imag_round == MPFR_RNDN ||
        result->new_ctx.imag_round == MPFR_RNDZ ||
        result->new_ctx.imag_round == MPFR_RNDU ||
        result->new_ctx.imag_round == MPFR_RNDD ||
        result->new_ctx.imag_round == GMPY_DEFAULT)) {
        VALUE_ERROR("invalid value for imag_round");
        goto error;
    }
#endif

    if (!(result->new_ctx.emin < 0 && result->new_ctx.emax > 0)) {
        VALUE_ERROR("invalid values for emin and/or emax");
        goto error;
    }

    if (mpfr_set_emin(result->new_ctx.emin)) {
        VALUE_ERROR("invalid value for emin");
        goto error;
    }
    if (mpfr_set_emax(result->new_ctx.emax)) {
        VALUE_ERROR("invalid value for emax");
        goto error;
    }

    return (PyObject*)result;

  error:
    if (arg_context) {
        Py_DECREF(local_args);
    }
    Py_DECREF((PyObject*)result);
    return NULL;
}

#ifdef WITHMPC

PyDoc_STRVAR(doc_context,
"context() -> context manager\n\n"
"Return a new context for controlling MPFR and MPC arithmetic. To load\n"
"the new context, use set_context(). Options can only be specified as\n"
"keyword arguments. \n\n"
"    precision:      precision, in bits, of an MPFR result\n"
"    real_prec:      precision, in bits, of Re(MPC)\n"
"                      -1 implies use mpfr_prec\n"
"    imag_prec:      precision, in bits, of Im(MPC)\n"
"                      -1 implies use real_prec\n"
"    round:          rounding mode for MPFR\n"
"    real_round:     rounding mode for Re(MPC)\n"
"                      -1 implies use mpfr_round\n"
"    imag_round:     rounding mode for Im(MPC)\n"
"                      -1 implies use real_round\n"
"    e_max:          maximum allowed exponent\n"
"    e_min:          minimum allowed exponent\n"
"    subnormalize:   if True, subnormalized results can be returned\n"
"    trap_underflow: if True, raise exception for underflow\n"
"                    if False, set underflow flag\n"
"    trap_overflow:  if True, raise exception for overflow\n"
"                    if False, set overflow flag and return Inf or -Inf\n"
"    trap_inexact:   if True, raise exception for inexact result\n"
"                    if False, set inexact flag\n"
"    trap_invalid:   if True, raise exception for invalid operation\n"
"                    if False, set invalid flag and return NaN\n"
"    trap_erange:    if True, raise exception for range error\n"
"                    if False, set erange flag\n"
"    trap_divzero:   if True, raise exception for division by zero\n"
"                    if False, set divzero flag and return Inf or -Inf\n"
"    trap_expbound:  if True, raise exception when mpfr/mpc exponent\n"
"                        no longer valid in current context\n"
"                    if False, mpfr/mpc with exponent out-of-bounds\n"
"                        will be coerced to either 0 or Infinity\n"
"    allow_complex:  if True, allow mpfr functions to return mpc\n"
"                    if False, mpfr functions cannot return an mpc\n");

#else

PyDoc_STRVAR(doc_context,
"context() -> context\n\n"
"Return a new context for controlling MPFR arithmetic. To load the\n"
"new context, use set_context(). Options can only be specified as\n"
"keyword arguments. \n\n"
"    precision:      precision, in bits, of an MPFR result\n"
"    round:          rounding mode for MPFR\n"
"    e_max:          maximum allowed exponent\n"
"    e_min:          minimum allowed exponent\n"
"    subnormalize:   if True, subnormalized results can be returned\n"
"    trap_underflow: if True, raise exception for underflow\n"
"                    if False, set underflow flag\n"
"    trap_overflow:  if True, raise exception for overflow\n"
"                    if False, set overflow flag and return Inf or -Inf\n"
"    trap_inexact:   if True, raise exception for inexact result\n"
"                    if False, set inexact flag\n"
"    trap_invalid:   if True, raise exception for invalid operation\n"
"                    if False, set invalid flag and return NaN\n"
"    trap_erange:    if True, raise exception for range error\n"
"                    if False, set erange flag\n"
"    trap_divzero:   if True, raise exception for division by zero\n"
"                    if False, set divzero flag and return Inf or -Inf\n"
"    trap_expbound:  if True, raise exception when mpfr/mpc exponent\n"
"                        no longer valid in current context\n"
"                    if False, mpfr/mpc with exponent out-of-bounds will\n"
"                        coerced to either 0 or Infinity\n");

#endif

static PyObject *
GMPyContext_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPyContextObject *result;

#ifdef WITHMPC
    static char *kwlist[] = {
        "precision", "real_prec", "imag_prec", "round",
        "real_round", "imag_round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero", "trap_expbound",
        "allow_complex", NULL };
#else
    static char *kwlist[] = {
        "precision", "round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero",
        "trap_expbound", NULL };
#endif

    if (PyTuple_GET_SIZE(args)) {
        VALUE_ERROR("context() only supports keyword arguments");
        return NULL;
    }

    if (!(result = (GMPyContextObject*)GMPyContext_new()))
        return NULL;

#ifdef WITHMPC
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|llliiilliiiiiiiii", kwlist,
            &result->ctx.mpfr_prec,
            &result->ctx.real_prec,
            &result->ctx.imag_prec,
            &result->ctx.mpfr_round,
            &result->ctx.real_round,
            &result->ctx.imag_round,
            &result->ctx.emax,
            &result->ctx.emin,
            &result->ctx.subnormalize,
            &result->ctx.trap_underflow,
            &result->ctx.trap_overflow,
            &result->ctx.trap_inexact,
            &result->ctx.trap_invalid,
            &result->ctx.trap_erange,
            &result->ctx.trap_divzero,
            &result->ctx.trap_expbound,
            &result->ctx.allow_complex))) {
#else
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|lilliiiiiiii", kwlist,
            &result->ctx.mpfr_prec,
            &result->ctx.mpfr_round,
            &result->ctx.emax,
            &result->ctx.emin,
            &result->ctx.subnormalize,
            &result->ctx.trap_underflow,
            &result->ctx.trap_overflow,
            &result->ctx.trap_inexact,
            &result->ctx.trap_invalid,
            &result->ctx.trap_erange,
            &result->ctx.trap_divzero,
            &result->ctx.trap_expbound))) {
#endif
        VALUE_ERROR("invalid keyword arguments in context()");
        return NULL;
    }

    /* Sanity check for values. */
    if (result->ctx.mpfr_prec < MPFR_PREC_MIN ||
        result->ctx.mpfr_prec > MPFR_PREC_MAX) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }

#ifdef WITHMPC
    if (!(result->ctx.real_prec == GMPY_DEFAULT ||
        (result->ctx.real_prec >= MPFR_PREC_MIN &&
        result->ctx.real_prec <= MPFR_PREC_MAX))) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for real_prec");
        return NULL;
    }
    if (!(result->ctx.imag_prec == GMPY_DEFAULT ||
        (result->ctx.imag_prec >= MPFR_PREC_MIN &&
        result->ctx.imag_prec <= MPFR_PREC_MAX))) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for imag_prec");
        return NULL;
    }
#endif

    if (!(result->ctx.mpfr_round == MPFR_RNDN ||
        result->ctx.mpfr_round == MPFR_RNDZ ||
        result->ctx.mpfr_round == MPFR_RNDU ||
        result->ctx.mpfr_round == MPFR_RNDD ||
        result->ctx.mpfr_round == MPFR_RNDA)) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for round");
        return NULL;
    }

#ifdef WITHMPC
    if (result->ctx.mpfr_round == MPFR_RNDA) {
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
           to MPFR_RNDN. */
        result->ctx.real_round = MPFR_RNDN;
        result->ctx.imag_round = MPFR_RNDN;
    }
    if (!(result->ctx.real_round == MPFR_RNDN ||
        result->ctx.real_round == MPFR_RNDZ ||
        result->ctx.real_round == MPFR_RNDU ||
        result->ctx.real_round == MPFR_RNDD ||
        result->ctx.real_round == GMPY_DEFAULT)) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for real_round");
        return NULL;
    }
    if (!(result->ctx.imag_round == MPFR_RNDN ||
        result->ctx.imag_round == MPFR_RNDZ ||
        result->ctx.imag_round == MPFR_RNDU ||
        result->ctx.imag_round == MPFR_RNDD ||
        result->ctx.imag_round == GMPY_DEFAULT)) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for imag_round");
        return NULL;
    }
#endif

    if (!(result->ctx.emin < 0 && result->ctx.emax > 0)) {
        VALUE_ERROR("invalid values for emin and/or emax");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    if (mpfr_set_emin(result->ctx.emin)) {
        VALUE_ERROR("invalid value for emin");
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    if (mpfr_set_emax(result->ctx.emax)) {
        VALUE_ERROR("invalid value for emax");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    result->ctx.underflow = 0;
    result->ctx.overflow = 0;
    result->ctx.inexact = 0;
    result->ctx.invalid = 0;
    result->ctx.erange = 0;
    result->ctx.divzero = 0;
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_set_context,
"set_context(context)\n\n"
"Activate a context object controlling MPFR and MPC arithmetic.\n");

static PyObject *
GMPyContext_set_context(PyObject *self, PyObject *other)
{
    if (GMPyContext_Check(other)) {
        Py_INCREF((PyObject*)other);
        Py_DECREF((PyObject*)context);
        context = (GMPyContextObject*)other;
        mpfr_set_emin(context->ctx.emin);
        mpfr_set_emax(context->ctx.emax);
        Py_RETURN_NONE;
    }
    else {
        VALUE_ERROR("set_context() requires a context argument");
        return NULL;
    }
}

static PyObject *
GMPyContextManager_enter(PyObject *self, PyObject *args)
{
    GMPyContextObject *result;

    if (!(result = (GMPyContextObject*)GMPyContext_new()))
        return NULL;

    result->ctx = ((GMPyContextManagerObject*)self)->new_ctx;
    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)result;
    Py_INCREF((PyObject*)context);
    mpfr_set_emin(context->ctx.emin);
    mpfr_set_emax(context->ctx.emax);
    return (PyObject*)result;
}

static PyObject *
GMPyContextManager_exit(PyObject *self, PyObject *args)
{
    GMPyContextObject *result;

    if (!(result = (GMPyContextObject*)GMPyContext_new()))
        return NULL;

    result->ctx = ((GMPyContextManagerObject*)self)->old_ctx;
    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)result;
    mpfr_set_emin(context->ctx.emin);
    mpfr_set_emax(context->ctx.emax);
    Py_RETURN_NONE;
}

static PyObject *
GMPyContext_enter(PyObject *self, PyObject *args)
{
    GMPyContextObject *result;

    if (!(result = (GMPyContextObject*)GMPyContext_new()))
        return NULL;

    result->ctx = ((GMPyContextObject*)self)->ctx;
    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)result;
    Py_INCREF((PyObject*)context);
    mpfr_set_emin(context->ctx.emin);
    mpfr_set_emax(context->ctx.emax);
    return (PyObject*)result;
}

static PyObject *
GMPyContext_exit(PyObject *self, PyObject *args)
{
    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)self;
    Py_INCREF((PyObject*)context);
    mpfr_set_emin(context->ctx.emin);
    mpfr_set_emax(context->ctx.emax);
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_context_clear_flags,
"clear_flags()\n\n"
"Clear all MPFR exception flags.");
static PyObject *
GMPyContext_clear_flags(PyObject *self, PyObject *args)
{
    ((GMPyContextObject*)self)->ctx.underflow = 0;
    ((GMPyContextObject*)self)->ctx.overflow = 0;
    ((GMPyContextObject*)self)->ctx.inexact = 0;
    ((GMPyContextObject*)self)->ctx.invalid = 0;
    ((GMPyContextObject*)self)->ctx.erange = 0;
    ((GMPyContextObject*)self)->ctx.divzero = 0;
    Py_RETURN_NONE;
}

/* Define the get/set functions. */

#define GETSET_BOOLEAN(NAME) \
static PyObject * \
GMPyContext_get_##NAME(GMPyContextObject *self, void *closure) \
{ \
    return PyBool_FromLong(self->ctx.NAME); \
}; \
static int \
GMPyContext_set_##NAME(GMPyContextObject *self, PyObject *value, void *closure) \
{ \
    if (!(PyBool_Check(value))) { \
        TYPE_ERROR(#NAME " must be True or False"); \
        return -1; \
    } \
    self->ctx.NAME = (value == Py_True) ? 1 : 0; \
    return 0; \
}

GETSET_BOOLEAN(subnormalize);
GETSET_BOOLEAN(underflow);
GETSET_BOOLEAN(overflow);
GETSET_BOOLEAN(inexact);
GETSET_BOOLEAN(invalid);
GETSET_BOOLEAN(erange);
GETSET_BOOLEAN(divzero);
GETSET_BOOLEAN(trap_underflow);
GETSET_BOOLEAN(trap_overflow);
GETSET_BOOLEAN(trap_inexact);
GETSET_BOOLEAN(trap_invalid);
GETSET_BOOLEAN(trap_erange);
GETSET_BOOLEAN(trap_divzero);
GETSET_BOOLEAN(trap_expbound);
#ifdef WITHMPC
GETSET_BOOLEAN(allow_complex)
#endif

static PyObject *
GMPyContext_get_precision(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->ctx.mpfr_prec));
}

static int
GMPyContext_set_precision(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("precision must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX || PyErr_Occurred()) {
        VALUE_ERROR("invalid value for precision");
        return -1;
    }
    self->ctx.mpfr_prec = (mpfr_prec_t)temp;
    return 0;
}

#ifdef WITHMPC
static PyObject *
GMPyContext_get_real_prec(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_REAL_PREC(self)));
}

static int
GMPyContext_set_real_prec(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("real_prec must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    if (temp == -1) {
        if (PyErr_Occurred()) {
            VALUE_ERROR("invalid value for real_prec");
            return -1;
        }
    }
    else if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for real_prec");
        return -1;
    }
    self->ctx.real_prec = (mpfr_prec_t)temp;
    return 0;
}

static PyObject *
GMPyContext_get_imag_prec(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_IMAG_PREC(self)));
}

static int
GMPyContext_set_imag_prec(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("imag_prec must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    if (temp == -1) {
        if (PyErr_Occurred()) {
            VALUE_ERROR("invalid value for imag_prec");
            return -1;
        }
    }
    else if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for imag_prec");
        return -1;
    }
    self->ctx.imag_prec = (mpfr_prec_t)temp;
    return 0;
}
#endif

static PyObject *
GMPyContext_get_round(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)(self->ctx.mpfr_round));
}

static int
GMPyContext_set_round(GMPyContextObject *self, PyObject *value, void *closure)
{
    long temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("round mode must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsLong(value);
    if (temp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    if (temp == MPFR_RNDN)
        self->ctx.mpfr_round = MPFR_RNDN;
    else if (temp == MPFR_RNDZ)
        self->ctx.mpfr_round = MPFR_RNDZ;
    else if (temp == MPFR_RNDU)
        self->ctx.mpfr_round = MPFR_RNDU;
    else if (temp == MPFR_RNDD)
        self->ctx.mpfr_round = MPFR_RNDD;
    else if (temp == MPFR_RNDA) {
        self->ctx.mpfr_round = MPFR_RNDA;
#ifdef WITHMPC
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
           to MPFR_RNDN. */
        self->ctx.real_round = MPFR_RNDN;
        self->ctx.imag_round = MPFR_RNDN;
#endif
    }
    else {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    return 0;
}

#ifdef WITHMPC
static PyObject *
GMPyContext_get_real_round(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_REAL_ROUND(self));
}

static int
GMPyContext_set_real_round(GMPyContextObject *self, PyObject *value, void *closure)
{
    long temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("round mode must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsLong(value);
    if (temp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    if (temp == GMPY_DEFAULT || temp == MPFR_RNDN || temp == MPFR_RNDZ ||
        temp == MPFR_RNDU || temp == MPFR_RNDD) {
        self->ctx.real_round = (int)temp;
    }
    else {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    return 0;
}

static PyObject *
GMPyContext_get_imag_round(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_IMAG_ROUND(self));
}

static int
GMPyContext_set_imag_round(GMPyContextObject *self, PyObject *value, void *closure)
{
    long temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("round mode must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsLong(value);
    if (temp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    if (temp == GMPY_DEFAULT || temp == MPFR_RNDN || temp == MPFR_RNDZ ||
        temp == MPFR_RNDU || temp == MPFR_RNDD) {
        self->ctx.imag_round = (int)temp;
    }
    else {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    return 0;
}
#endif

static PyObject *
GMPyContext_get_emin(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong(self->ctx.emin);
}

static int
GMPyContext_set_emin(GMPyContextObject *self, PyObject *value, void *closure)
{
    long exp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("emin must be Python integer");
        return -1;
    }
    exp = PyIntOrLong_AsLong(value);
    if (exp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("requested minimum exponent is invalid");
        return -1;
    }
    if (mpfr_set_emin((mpfr_prec_t)exp)) {
        VALUE_ERROR("requested minimum exponent is invalid");
        return -1;
    }
    self->ctx.emin = exp;
    mpfr_set_emin(exp);
    return 0;
}

static PyObject *
GMPyContext_get_emax(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong(self->ctx.emax);
}

static int
GMPyContext_set_emax(GMPyContextObject *self, PyObject *value, void *closure)
{
    long exp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("emax must be Python integer");
        return -1;
    }
    exp = PyIntOrLong_AsLong(value);
    if (exp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("requested maximum exponent is invalid");
        return -1;
    }
    if (mpfr_set_emax((mpfr_prec_t)exp)) {
        VALUE_ERROR("requested maximum exponent is invalid");
        return -1;
    }
    self->ctx.emax = exp;
    mpfr_set_emax(exp);
    return 0;
}

#define ADD_GETSET(NAME) \
    {#NAME, \
        (getter)GMPyContext_get_##NAME, \
        (setter)GMPyContext_set_##NAME, NULL, NULL}

static PyGetSetDef GMPyContext_getseters[] = {
    ADD_GETSET(precision),
#ifdef WITHMPC
    ADD_GETSET(real_prec),
    ADD_GETSET(imag_prec),
#endif
    ADD_GETSET(round),
#ifdef WITHMPC
    ADD_GETSET(real_round),
    ADD_GETSET(imag_round),
#endif
    ADD_GETSET(emax),
    ADD_GETSET(emin),
    ADD_GETSET(subnormalize),
    ADD_GETSET(underflow),
    ADD_GETSET(overflow),
    ADD_GETSET(inexact),
    ADD_GETSET(invalid),
    ADD_GETSET(erange),
    ADD_GETSET(divzero),
    ADD_GETSET(trap_underflow),
    ADD_GETSET(trap_overflow),
    ADD_GETSET(trap_inexact),
    ADD_GETSET(trap_invalid),
    ADD_GETSET(trap_erange),
    ADD_GETSET(trap_divzero),
    ADD_GETSET(trap_expbound),
#ifdef WITHMPC
    ADD_GETSET(allow_complex),
#endif
    {NULL}
};

static PyMethodDef GMPyContext_methods[] =
{
    { "clear_flags", GMPyContext_clear_flags, METH_NOARGS,
            doc_context_clear_flags },
    { "__enter__", GMPyContext_enter, METH_NOARGS, NULL },
    { "__exit__", GMPyContext_exit, METH_VARARGS, NULL },
    { NULL, NULL, 1 }
};

static PyTypeObject GMPyContext_Type =
{
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "gmpy2 context",                        /* tp_name          */
    sizeof(GMPyContextObject),              /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    (destructor) GMPyContext_dealloc,       /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPyContext_repr,            /* tp_repr          */
        0,                                  /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
        0,                                  /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
    "GMPY2 Context Object",                 /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
        0,                                  /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    GMPyContext_methods,                    /* tp_methods       */
        0,                                  /* tp_members       */
    GMPyContext_getseters,                  /* tp_getset        */
};

static PyMethodDef GMPyContextManager_methods[] =
{
    { "__enter__", GMPyContextManager_enter, METH_NOARGS, NULL },
    { "__exit__", GMPyContextManager_exit, METH_VARARGS, NULL },
    { NULL, NULL, 1 }
};

static PyTypeObject GMPyContextManager_Type =
{
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                   /* ob_size          */
#endif
    "gmpy2 context",                         /* tp_name          */
    sizeof(GMPyContextManagerObject),        /* tp_basicsize     */
        0,                                   /* tp_itemsize      */
    (destructor) GMPyContextManager_dealloc, /* tp_dealloc       */
        0,                                   /* tp_print         */
        0,                                   /* tp_getattr       */
        0,                                   /* tp_setattr       */
        0,                                   /* tp_reserved      */
    (reprfunc) GMPyContextManager_repr,      /* tp_repr          */
        0,                                   /* tp_as_number     */
        0,                                   /* tp_as_sequence   */
        0,                                   /* tp_as_mapping    */
        0,                                   /* tp_hash          */
        0,                                   /* tp_call          */
        0,                                   /* tp_str           */
        0,                                   /* tp_getattro      */
        0,                                   /* tp_setattro      */
        0,                                   /* tp_as_buffer     */
    Py_TPFLAGS_DEFAULT,                      /* tp_flags         */
    "GMPY2 Context manager",                 /* tp_doc           */
        0,                                   /* tp_traverse      */
        0,                                   /* tp_clear         */
        0,                                   /* tp_richcompare   */
        0,                                   /* tp_weaklistoffset*/
        0,                                   /* tp_iter          */
        0,                                   /* tp_iternext      */
    GMPyContextManager_methods,              /* tp_methods       */
        0,                                   /* tp_members       */
        0,                                   /* tp_getset        */
};
