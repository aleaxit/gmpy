/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_context.c                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
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

/* This file implements contexts and context managers.
 *
 * Public API
 * ==========
 *   TBD
 *
 * Private API
 * ===========
 *   GMPy_CTXT_New
 *   GMPy_CTXT_Dealloc
 *   GMPy_CTXT_Set
 *   GMPy_CTXT_Get
 *   GMPy_CTXT_Copy
 *   GMPy_CTXT_ieee
 *   GMPy_CTXT_Local
 *   GMPy_CTXT_Context
 *   GMPy_CTXT_Repr_Slot
 *   GMPy_CTXT_Enter
 *   GMPy_CTXT_Exit
 *   GMPy_CTXT_Clear_Flags
 *   GMPy_CTXT_Manager_New
 *   GMPy_CTXT_Manager_Dealloc
 *   GMPy_CTXT_Manager_Repr_Slot
 *   GMPy_CTXT_Manager_Enter
 *   GMPy_CTXT_Manager_Exit
 *     plus getters & setters....
 *
 * Internal functions
 * ==================
 *   GMPy_current_context
 */

/* Create and delete Context objects. */

static PyObject *
GMPy_CTXT_New(void)
{
    CTXT_Object *result;

    if ((result = PyObject_New(CTXT_Object, &CTXT_Type))) {
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
        result->ctx.traps = TRAP_NONE;
        result->ctx.real_prec = -1;
        result->ctx.imag_prec = -1;
        result->ctx.real_round = -1;
        result->ctx.imag_round = -1;
        result->ctx.allow_complex = 0;
        result->ctx.rational_division = 0;
        result->ctx.allow_release_gil = 0;

#ifndef WITHOUT_THREADS
        result->tstate = NULL;
#endif

    }
    return (PyObject*)result;
};

static void
GMPy_CTXT_Dealloc(CTXT_Object *self)
{
    PyObject_Del(self);
};

/* Support for global and thread local contexts. */

/* Doc-string, alternate definitions below. */

PyDoc_STRVAR(GMPy_doc_set_context,
"set_context(context)\n\n"
"Activate a context object controlling gmpy2 arithmetic.\n");

#ifdef WITHOUT_THREADS

static PyObject *
GMPy_CTXT_Set(PyObject *self, PyObject *other)
{
    if (!CTXT_Check(other)) {
        VALUE_ERROR("set_context() requires a context argument");
        return NULL;
    }

    Py_DECREF((PyObject*)module_context);
    Py_INCREF((PyObject*)other);
    module_context = (CTXT_Object*)other;

    Py_RETURN_NONE;
}

#else

/* Begin support for thread local contexts. */

/* Get the context from the thread state dictionary. */
static CTXT_Object *
current_context_from_dict(void)
{
    PyObject *dict;
    PyObject *tl_context;
    PyThreadState *tstate;

    dict = PyThreadState_GetDict();
    if (dict == NULL) {
        RUNTIME_ERROR("cannot get thread state");
        return NULL;
    }

#ifdef PY3
    tl_context = PyDict_GetItemWithError(dict, tls_context_key);
#else
    tl_context = PyDict_GetItem(dict, tls_context_key);
#endif
    if (!tl_context) {
#ifdef PY3
        if (PyErr_Occurred()) {
            return NULL;
        }
#endif

        /* Set up a new thread local context. */
        tl_context = GMPy_CTXT_New();
        if (!tl_context) {
            return NULL;
        }

        if (PyDict_SetItem(dict, tls_context_key, tl_context) < 0) {
            Py_DECREF(tl_context);
            return NULL;
        }
        Py_DECREF(tl_context);
    }

    /* Cache the context of the current thread, assuming that it
     * will be accessed several times before a thread switch. */
    tstate = PyThreadState_GET();
    if (tstate) {
        cached_context = (CTXT_Object*)tl_context;
        cached_context->tstate = tstate;
    }

    /* Borrowed reference with refcount==1 */
    return (CTXT_Object*)tl_context;
}

/* Return borrowed reference to thread local context. */
static CTXT_Object *
GMPy_current_context(void)
{
    PyThreadState *tstate = PyThreadState_GET();

    if (cached_context && cached_context->tstate == tstate) {
        return (CTXT_Object*)cached_context;
    }

    return current_context_from_dict();
}

/* Set the thread local context to a new context, decrement old reference */
static PyObject *
GMPy_CTXT_Set(PyObject *self, PyObject *other)
{
    PyObject *dict;
    PyThreadState *tstate;

    if (!CTXT_Check(other)) {
        VALUE_ERROR("set_context() requires a context argument");
        return NULL;
    }

    dict = PyThreadState_GetDict();
    if (dict == NULL) {
        RUNTIME_ERROR("cannot get thread state");
        return NULL;
    }

    if (PyDict_SetItem(dict, tls_context_key, other) < 0) {
        return NULL;
    }

    /* Cache the context of the current thread, assuming that it
     * will be accessed several times before a thread switch. */
    cached_context = NULL;
    tstate = PyThreadState_GET();
    if (tstate) {
        cached_context = (CTXT_Object*)other;
        cached_context->tstate = tstate;
    }

    Py_RETURN_NONE;
}
#endif

PyDoc_STRVAR(GMPy_doc_context_ieee,
"ieee(size[,subnormalize=True]) -> context\n\n"
"Return a new context corresponding to a standard IEEE floating point\n"
"format. The supported sizes are 16, 32, 64, 128, and multiples of\n"
"32 greater than 128.");

static PyObject *
GMPy_CTXT_ieee(PyObject *self, PyObject *args, PyObject *kwargs)
{
    long bitwidth;
    double bitlog2;
    int sub_mode=1;
    PyObject *temp;
    CTXT_Object *result;
    static char *kwlist[] = {"subnormalize", NULL};

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("ieee() requires 'int' argument");
        return NULL;
    }

    bitwidth = PyIntOrLong_AsLong(PyTuple_GET_ITEM(args, 0));
    if (bitwidth == -1 && PyErr_Occurred()) {
        TYPE_ERROR("ieee() requires 'int' argument");
        return NULL;
    }

    if (bitwidth <= 0) {
        VALUE_ERROR("ieee() requires positive value for size");
        return NULL;
    }

    /* Process just the two keyword arguments. */

    if (!(temp = PyTuple_New(0))) {
        return NULL;
    }

    if (!(PyArg_ParseTupleAndKeywords(temp, kwargs,
            "|ii", kwlist, &sub_mode))) {
        VALUE_ERROR("invalid keyword arguments for ieee()");
        Py_DECREF(temp);
        return NULL;
    }
    Py_DECREF(temp);

    if (sub_mode)
        sub_mode = 1;

    if (!(result = (CTXT_Object*)GMPy_CTXT_New())) {
        return NULL;
    }

    if (bitwidth == 16) {
        result->ctx.mpfr_prec = 11;
        result->ctx.emax = 16;
    }
    else if (bitwidth == 32) {
        result->ctx.mpfr_prec = 24;
        result->ctx.emax = 128;
    }
    else if (bitwidth == 64) {
        result->ctx.mpfr_prec = 53;
        result->ctx.emax = 1024;
    }
    else if (bitwidth == 128) {
        result->ctx.mpfr_prec = 113;
        result->ctx.emax = 16384;
    }
    else {
        if ((bitwidth < 128) && (bitwidth & 31)) {
            VALUE_ERROR("bitwidth must be 16, 32, 64, 128; or must be greater than 128 and divisible by 32.");
            Py_DECREF((PyObject*)result);
            return NULL;
        }
        bitlog2 = floor((4 * log(bitwidth) / log(2.0)) + 0.5);
        result->ctx.mpfr_prec = bitwidth - (long)(bitlog2) + 13;
        result->ctx.emax = 1 << (bitwidth - result->ctx.mpfr_prec - 1);
    }

    result->ctx.subnormalize = sub_mode;
    result->ctx.emin = 4 - result->ctx.emax - result->ctx.mpfr_prec;
    return (PyObject*)result;
}

/* Create and delete ContextManager objects. */

static PyObject *
GMPy_CTXT_Manager_New(void)
{
    PyObject *result;

    result = (PyObject*)PyObject_New(CTXT_Manager_Object, &CTXT_Manager_Type);
    ((CTXT_Manager_Object*)(result))->new_context = NULL;
    ((CTXT_Manager_Object*)(result))->old_context = NULL;
    return result;
}

static void
GMPy_CTXT_Manager_Dealloc(CTXT_Manager_Object *self)
{
    Py_XDECREF(self->new_context);
    Py_XDECREF(self->old_context);
    PyObject_Del(self);
}

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
}

static PyObject *
GMPy_CTXT_Repr_Slot(CTXT_Object *self)
{
    PyObject *format;
    PyObject *tuple;
    PyObject *result = NULL;
    int i = 0;

    tuple = PyTuple_New(24);
    if (!tuple)
        return NULL;

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
            "        allow_complex=%s,\n"
            "        rational_division=%s,\n"
            "        allow_release_gil=%s)"
            );
    if (!format) {
        Py_DECREF(tuple);
        return NULL;
    }

    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.mpfr_prec));
    if (self->ctx.real_prec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, i++, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.real_prec));
    if (self->ctx.imag_prec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, i++, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.imag_prec));
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->ctx.mpfr_round));
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->ctx.real_round));
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->ctx.imag_round));
    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.emax));
    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->ctx.emin));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.subnormalize));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.traps & TRAP_UNDERFLOW));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.underflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.traps & TRAP_OVERFLOW));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.overflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.traps & TRAP_INEXACT));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.inexact));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.traps & TRAP_INVALID));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.invalid));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.traps & TRAP_ERANGE));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.erange));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.traps & TRAP_DIVZERO));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.divzero));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.allow_complex));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.rational_division));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->ctx.allow_release_gil));

    if (!PyErr_Occurred())
        result = Py2or3String_Format(format, tuple);
    else
        SYSTEM_ERROR("internal error in GMPy_CTXT_Repr");

    Py_DECREF(format);
    Py_DECREF(tuple);
    return result;
}

static PyObject *
GMPy_CTXT_Manager_Repr_Slot(CTXT_Manager_Object *self)
{
    return Py_BuildValue("s", "<gmpy2.ContextManagerObject>");
}

PyDoc_STRVAR(GMPy_doc_get_context,
"get_context() -> gmpy2 context\n\n"
"Return a reference to the current context.");

static PyObject *
GMPy_CTXT_Get(PyObject *self, PyObject *args)
{
    CTXT_Object *context;

    CURRENT_CONTEXT(context);
    Py_XINCREF((PyObject*)context);
    return (PyObject*)context;
}

PyDoc_STRVAR(GMPy_doc_context_copy,
"context.copy() -> gmpy2 context\n\n"
"Return a copy of a context.");

static PyObject *
GMPy_CTXT_Copy(PyObject *self, PyObject *other)
{
    CTXT_Object *result;

    result = (CTXT_Object*)GMPy_CTXT_New();
    result->ctx = ((CTXT_Object*)self)->ctx;
    return (PyObject*)result;
}

/* Parse the keyword arguments available to a context. Returns 1 if no
 * error occurred; returns 0 is an error occurred.
 */

static int
_parse_context_args(CTXT_Object *ctxt, PyObject *kwargs)
{
    PyObject *args;
    int x_trap_underflow = 0, x_trap_overflow = 0, x_trap_inexact = 0;
    int x_trap_invalid = 0, x_trap_erange = 0, x_trap_divzero = 0;

    static char *kwlist[] = {
        "precision", "real_prec", "imag_prec", "round",
        "real_round", "imag_round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero", "allow_complex",
        "rational_division", "allow_release_gil", NULL };

    /* Create an empty dummy tuple to use for args. */

    if (!(args = PyTuple_New(0))) {
        return 0;
    }

    /* Convert the trap bit positions into ints for the benefit of
     * PyArg_ParseTupleAndKeywords().
     */
    x_trap_underflow = ctxt->ctx.traps & TRAP_UNDERFLOW;
    x_trap_overflow = ctxt->ctx.traps & TRAP_OVERFLOW;
    x_trap_inexact = ctxt->ctx.traps & TRAP_INEXACT;
    x_trap_invalid = ctxt->ctx.traps & TRAP_INVALID;
    x_trap_erange = ctxt->ctx.traps & TRAP_ERANGE;
    x_trap_divzero = ctxt->ctx.traps & TRAP_DIVZERO;

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|llliiilliiiiiiiiii", kwlist,
            &ctxt->ctx.mpfr_prec,
            &ctxt->ctx.real_prec,
            &ctxt->ctx.imag_prec,
            &ctxt->ctx.mpfr_round,
            &ctxt->ctx.real_round,
            &ctxt->ctx.imag_round,
            &ctxt->ctx.emax,
            &ctxt->ctx.emin,
            &ctxt->ctx.subnormalize,
            &x_trap_underflow,
            &x_trap_overflow,
            &x_trap_inexact,
            &x_trap_invalid,
            &x_trap_erange,
            &x_trap_divzero,
            &ctxt->ctx.allow_complex,
            &ctxt->ctx.rational_division,
            &ctxt->ctx.allow_release_gil))) {
        VALUE_ERROR("invalid keyword arguments for context");
        Py_DECREF(args);
        return 0;
    }
    Py_DECREF(args);

    ctxt->ctx.traps = TRAP_NONE;
    if (x_trap_underflow)
        ctxt->ctx.traps |= TRAP_UNDERFLOW;
    if (x_trap_overflow)
        ctxt->ctx.traps |= TRAP_OVERFLOW;
    if (x_trap_inexact)
        ctxt->ctx.traps |= TRAP_INEXACT;
    if (x_trap_invalid)
        ctxt->ctx.traps |= TRAP_INVALID;
    if (x_trap_erange)
        ctxt->ctx.traps |= TRAP_ERANGE;
    if (x_trap_divzero)
        ctxt->ctx.traps |= TRAP_DIVZERO;

    if (ctxt->ctx.subnormalize)
        ctxt->ctx.subnormalize = 1;

    /* Sanity check for values. */
    if (ctxt->ctx.mpfr_prec < MPFR_PREC_MIN ||
        ctxt->ctx.mpfr_prec > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return 0;
    }

    if (!(ctxt->ctx.real_prec == GMPY_DEFAULT ||
        (ctxt->ctx.real_prec >= MPFR_PREC_MIN &&
        ctxt->ctx.real_prec <= MPFR_PREC_MAX))) {
        VALUE_ERROR("invalid value for real_prec");
        return 0;
    }

    if (!(ctxt->ctx.imag_prec == GMPY_DEFAULT ||
        (ctxt->ctx.imag_prec >= MPFR_PREC_MIN &&
        ctxt->ctx.imag_prec <= MPFR_PREC_MAX))) {
        VALUE_ERROR("invalid value for imag_prec");
        return 0;
    }

    if (!(ctxt->ctx.mpfr_round == MPFR_RNDN ||
        ctxt->ctx.mpfr_round == MPFR_RNDZ ||
        ctxt->ctx.mpfr_round == MPFR_RNDU ||
        ctxt->ctx.mpfr_round == MPFR_RNDD ||
        ctxt->ctx.mpfr_round == MPFR_RNDA)) {
        VALUE_ERROR("invalid value for round");
        return 0;
    }

    if (ctxt->ctx.mpfr_round == MPFR_RNDA) {
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
         * to MPFR_RNDN.
         */
        ctxt->ctx.real_round = MPFR_RNDN;
        ctxt->ctx.imag_round = MPFR_RNDN;
    }

    if (!(ctxt->ctx.real_round == MPFR_RNDN ||
        ctxt->ctx.real_round == MPFR_RNDZ ||
        ctxt->ctx.real_round == MPFR_RNDU ||
        ctxt->ctx.real_round == MPFR_RNDD ||
        ctxt->ctx.real_round == GMPY_DEFAULT)) {
        VALUE_ERROR("invalid value for real_round");
        return 0;
    }

    if (!(ctxt->ctx.imag_round == MPFR_RNDN ||
        ctxt->ctx.imag_round == MPFR_RNDZ ||
        ctxt->ctx.imag_round == MPFR_RNDU ||
        ctxt->ctx.imag_round == MPFR_RNDD ||
        ctxt->ctx.imag_round == GMPY_DEFAULT)) {
        VALUE_ERROR("invalid value for imag_round");
        return 0;
    }

    if (ctxt->ctx.emin < mpfr_get_emin_min() ||
        ctxt->ctx.emin > mpfr_get_emin_max()) {
        VALUE_ERROR("invalid value for emin");
        return 0;
    }

    if (ctxt->ctx.emax < mpfr_get_emax_min() ||
        ctxt->ctx.emax > mpfr_get_emax_max()) {
        VALUE_ERROR("invalid value for emax");
        return 0;
    }

    return 1;
}

PyDoc_STRVAR(GMPy_doc_local_context,
"local_context([context[,keywords]]) -> context manager\n\n"
"Create a context manager object that will restore the current context\n"
"when the 'with ...' block terminates. The temporary context for the\n"
"'with ...' block is based on the current context if no context is\n"
"specified. Keyword arguments are supported and will modify the\n"
"temporary new context.");

static PyObject *
GMPy_CTXT_Local(PyObject *self, PyObject *args, PyObject *kwargs)
{
    CTXT_Manager_Object *result;
    int arg_context = 0;
    CTXT_Object *context, *temp;

    CURRENT_CONTEXT(context);

    if (PyTuple_GET_SIZE(args) == 1 && CTXT_Check(PyTuple_GET_ITEM(args, 0))) {
        arg_context = 1;
    }
    else if (PyTuple_GET_SIZE(args)) {
        VALUE_ERROR("local_context() only supports [context[,keyword]] arguments");
        return NULL;
    }

    if (!(result = (CTXT_Manager_Object*)GMPy_CTXT_Manager_New()))
        return NULL;

    if (arg_context) {
        temp = (CTXT_Object*)PyTuple_GET_ITEM(args, 0);
        result->new_context = temp;
        Py_INCREF((PyObject*)(result->new_context));
    }
    else {
        result->new_context = context;
        Py_INCREF((PyObject*)(result->new_context));
    }

    result->old_context = (CTXT_Object*)GMPy_CTXT_Copy((PyObject*)context, NULL);
    if (!(result->old_context)) {
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    if (!_parse_context_args(result->new_context, kwargs)) {
        /* There was an error parsing the keyword arguments. */
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    else {
        /* Parsing was successful. */
        return (PyObject*)result;
    }
}

PyDoc_STRVAR(GMPy_doc_context,
"context() -> context manager\n\n"
"Return a new context for controlling MPFR and MPC arithmetic. To load\n"
"the new context, use set_context(). Options can only be specified as\n"
"keyword arguments. \n"
"\nOptions\n"
"    precision:         precision, in bits, of an MPFR result\n"
"    real_prec:         precision, in bits, of Re(MPC)\n"
"                         -1 implies use mpfr_prec\n"
"    imag_prec:         precision, in bits, of Im(MPC)\n"
"                         -1 implies use real_prec\n"
"    round:             rounding mode for MPFR\n"
"    real_round:        rounding mode for Re(MPC)\n"
"                         -1 implies use mpfr_round\n"
"    imag_round:        rounding mode for Im(MPC)\n"
"                         -1 implies use real_round\n"
"    e_max:             maximum allowed exponent\n"
"    e_min:             minimum allowed exponent\n"
"    subnormalize:      if True, subnormalized results can be returned\n"
"    trap_underflow:    if True, raise exception for underflow\n"
"                       if False, set underflow flag\n"
"    trap_overflow:     if True, raise exception for overflow\n"
"                       if False, set overflow flag and return Inf or -Inf\n"
"    trap_inexact:      if True, raise exception for inexact result\n"
"                       if False, set inexact flag\n"
"    trap_invalid:      if True, raise exception for invalid operation\n"
"                       if False, set invalid flag and return NaN\n"
"    trap_erange:       if True, raise exception for range error\n"
"                       if False, set erange flag\n"
"    trap_divzero:      if True, raise exception for division by zero\n"
"                       if False, set divzero flag and return Inf or -Inf\n"
"    allow_complex:     if True, allow mpfr functions to return mpc\n"
"                       if False, mpfr functions cannot return an mpc\n"
"    rational_division: if True, mpz/mpz returns an mpq\n"
"                       if False, mpz/mpz follows default behavior\n"
"    allow_release_gil: if True, mpq operations may release the GIL\n"
"                       if False, mpq operations may not release the GIL\n");
#if 0
"\nMethods\n"
"    abs(x)          return absolute value of x\n"
"    acos(x)         return inverse cosine of x\n"
"    acosh(x)        return inverse hyperbolic cosine of x\n"
"    add(x,y)        return x + y\n"
"    agm(x,y)        return arthimetic-geometric mean of x and y\n"
"    ai(x)           return the Airy function of x\n"
"    asin(x)         return inverse sine of x\n"
"    asinh(x)        return inverse hyperbolic sine of x\n"
"    atan(x)         return inverse tangent of x\n"
"    atan2(y,x)      return inverse tangent of (y / x)\n"
"    atanh(x)        return inverse hyperbolic tangent of x\n"
"    cbrt(x)         return cube root of x\n"
"    ceil(x)         return ceiling of x\n"
"    check_range(x)  return value with exponents within current range\n"
"    clear_flags()   clear all exception flags\n"
"    const_catalan() return Catalan constant (0.91596559...)\n"
"    const_euler()   return Euler contstant (0.57721566...)\n"
"    const_log()     return natural log of 2 (0.69314718...)\n"
"    const_pi()      return Pi (3.14159265...)\n"
"    copy()          return a copy of the context\n"
"    cos(x)          return cosine of x\n"
"    cosh(x)         return hyperbolic cosine of x\n"
"    cot(x)          return cotangent of x\n"
"    coth(x)         return hyperbolic cotangent of x\n"
"    csc(x)          return cosecant of x\n"
"    csch(x)         return hyperbolic cosecant of x\n"
"    degrees(x)      convert value in radians to degrees\n"
"    digamma(x)      return the digamma of x\n"
"    div(x,y)        return x / y\n"
"    divmod(x,y)     return integer quotient and remainder\n"
"    div_2exp(x,n)   return x / 2**n)\n"
"    eint(x)         return exponential integral of x\n"
"    erf(x)          return error function of x\n"
"    erfc(x)         return complementary error function of x\n"
"    exp(x)          return e**x\n"
"    exp10(x)        return 10**x\n"
"    exp2(x)         return 2**x\n"
"    expm1(x)        return e**x - 1\n"
"    factorial(n)    return floating-point approximation to n!\n"
"    floor(x)        return floor of x\n"
"    fma(x,y,z)      return correctly rounded (x * y) + z\n"
"    fmod(x,y)       return x - int(x / y) * y, rounding to 0\n"
"    fms(x,y,z)      return correctly rounded (x * y) - z\n"
"    fsum(i)         return accurate sum of iterable i\n"
"    gamma(x)        return gamma of x\n"
"    hypot(y,x)      return square root of (x**2 + y**2)\n"
"    is_finite(x)    return True if x is finite\n"
"    is_infinite(x)\n"
"    is_nan(x)\n"
"    is_zero(x)\n"
"    j0(x)           return Bessel of first kind of order 0 of x\n"
"    j1(x)           return Bessel of first kind of order 1 of x\n"
"    jn(x,n)         return Bessel of first kind of order n of x\n"
"    lgamma(x)       return tuple (log(abs(gamma(x)), sign(gamma(x)))\n"
"    li2(x)          return real part of dilogarithm of x\n"
"    lngamma(x)      return logarithm of gamma of x\n"
"    log(x)          return natural logarithm of x\n"
"    log10(x)        return base-10 logarithm of x\n"
"    log2(x)         return base-2 logarithm of x\n"
"    max2(x,y)       return maximum of x and y, rounded to context\n"
"    mpc(...)        create a new instance of an mpc\n"
"    mpfr(...)       create a new instance of an mpfr\n"
"    minus(x)        return -x\n"
"    min2(x,y)       return minimum of x and y, rounded to context\n"
"    mul(x,y)        return x * y\n"
"    mul_2exp(x,n)   return x * 2**n\n"
"    next_above(x)   return next mpfr towards +Infinity\n"
"    next_below(x)   return next mpfr towards -Infinity\n"
"    plus(x)         return +x\n"
"    pow(x,y)        return x ** y\n"
"    radians(x)      convert value in degrees to radians\n"
"    rec_sqrt(x)     return 1 / sqrt(x)\n"
"    rel_diff(x,y)   return abs(x - y) / x\n"
"    remainder(x,y)  return x - int(x / y) * y, rounding to even\n"
"    remquo(x,y)     return tuple of remainder(x,y) and low bits of\n"
"                    the quotient\n"
"    rint(x)         return x rounded to integer with current rounding\n"
"    rint_ceil(x)    ...\n"
"    rint_floor(x)   ...\n"
"    rint_round(x)   ...\n"
"    rint_trunc(x)   ...\n"
"    root(x,n)       return the n-th of x\n"
#ifdef MPC_110
"    root_of_unity() return the k-th power of the n-th root of mpc(1)\n"
#endif
"    round2(x,n)     return x rounded to n bits.\n"
"    round_away(x)   return x rounded to integer, ties away from 0\n"
"    sec(x)          return secant of x\n"
"    sech(x)         return hyperbolic secant of x\n"
"    sin(x)          return sine of x\n"
"    sin_cos(x)      return tuple (sin(x), cos(x))\n"
"    sinh(x)         return hyperbolic sine of x\n"
"    sinh_cosh(x)    return tuple (sinh(x), cosh(x))\n"
"    sqrt(x)         return square root of x\n"
"    square(x)       return x * x\n"
"    sub(x)          return x - y\n"
"    tan(x)          return tangent of x\n"
"    tanh(x)         return hyperbolic tangent of x\n"
"    trunc(x)        return x rounded towards 0\n"
"    y0(x)           return Bessel of second kind of order 0 of x\n"
"    y1(x)           return Bessel of second kind of order 1 of x\n"
"    yn(x,n)         return Bessel of second kind of order n of x\n"
"    zeta(x)         return Riemann zeta of x"
#endif

static PyObject *
GMPy_CTXT_Context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    CTXT_Object *result;

    if (PyTuple_GET_SIZE(args)) {
        VALUE_ERROR("context() only supports keyword arguments");
        return NULL;
    }

    if (!(result = (CTXT_Object*)GMPy_CTXT_New()))
        return NULL;

    if (!_parse_context_args(result, kwargs)) {
        /* There was an error parsing the keyword arguments. */
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    else {
        /* Parsing was successful. */
        return (PyObject*)result;
    }
}

static PyObject *
GMPy_CTXT_Manager_Enter(PyObject *self, PyObject *args)
{
    PyObject *temp;

    temp = GMPy_CTXT_Set(NULL, (PyObject*)((CTXT_Manager_Object*)self)->new_context);
    if (!temp)
        return NULL;
    Py_DECREF(temp);

    Py_INCREF((PyObject*)(((CTXT_Manager_Object*)self)->new_context));
    return (PyObject*)(((CTXT_Manager_Object*)self)->new_context);
}

static PyObject *
GMPy_CTXT_Manager_Exit(PyObject *self, PyObject *args)
{
    PyObject *temp;

    temp = GMPy_CTXT_Set(NULL, (PyObject*)((CTXT_Manager_Object*)self)->old_context);
    if (!temp)
        return NULL;
    Py_DECREF(temp);

    Py_RETURN_NONE;
}

static PyObject *
GMPy_CTXT_Enter(PyObject *self, PyObject *args)
{
    PyObject *temp;
    PyObject *result;

    result = GMPy_CTXT_Copy(self, NULL);
    if (!result)
        return NULL;

    temp = GMPy_CTXT_Set(NULL, result);
    if (!temp)
        return NULL;
    Py_DECREF(temp);

    return result;
}

static PyObject *
GMPy_CTXT_Exit(PyObject *self, PyObject *args)
{
    PyObject *temp;

    temp = GMPy_CTXT_Set(NULL, self);
    if (!temp)
        return NULL;
    Py_DECREF(temp);

    Py_RETURN_NONE;
}

PyDoc_STRVAR(GMPy_doc_context_clear_flags,
"clear_flags()\n\n"
"Clear all MPFR exception flags.");

static PyObject *
GMPy_CTXT_Clear_Flags(PyObject *self, PyObject *args)
{
    ((CTXT_Object*)self)->ctx.underflow = 0;
    ((CTXT_Object*)self)->ctx.overflow = 0;
    ((CTXT_Object*)self)->ctx.inexact = 0;
    ((CTXT_Object*)self)->ctx.invalid = 0;
    ((CTXT_Object*)self)->ctx.erange = 0;
    ((CTXT_Object*)self)->ctx.divzero = 0;
    Py_RETURN_NONE;
}

/* Define the get/set functions. */

#define GETSET_BOOLEAN(NAME) \
static PyObject * \
GMPy_CTXT_Get_##NAME(CTXT_Object *self, void *closure) \
{ \
    return PyBool_FromLong(self->ctx.NAME); \
}; \
static int \
GMPy_CTXT_Set_##NAME(CTXT_Object *self, PyObject *value, void *closure) \
{ \
    if (!(PyBool_Check(value))) { \
        TYPE_ERROR(#NAME " must be True or False"); \
        return -1; \
    } \
    self->ctx.NAME = (value == Py_True) ? 1 : 0; \
    return 0; \
}

/* Define the get/set functions. This version works with the individual
 * bits in the traps field.
 */

#define GETSET_BOOLEAN_BIT(NAME, TRAP) \
static PyObject * \
GMPy_CTXT_Get_##NAME(CTXT_Object *self, void *closure) \
{ \
    return PyBool_FromLong(self->ctx.traps & TRAP); \
}; \
static int \
GMPy_CTXT_Set_##NAME(CTXT_Object *self, PyObject *value, void *closure) \
{ \
    if (!(PyBool_Check(value))) { \
        TYPE_ERROR(#NAME " must be True or False"); \
        return -1; \
    } \
    if (value == Py_True) \
        self->ctx.traps |= TRAP; \
    else \
        self->ctx.traps &= ~(TRAP); \
    return 0; \
}

GETSET_BOOLEAN(subnormalize);
GETSET_BOOLEAN(underflow);
GETSET_BOOLEAN(overflow);
GETSET_BOOLEAN(inexact);
GETSET_BOOLEAN(invalid);
GETSET_BOOLEAN(erange);
GETSET_BOOLEAN(divzero);
GETSET_BOOLEAN_BIT(trap_underflow, TRAP_UNDERFLOW);
GETSET_BOOLEAN_BIT(trap_overflow, TRAP_OVERFLOW);
GETSET_BOOLEAN_BIT(trap_inexact, TRAP_INEXACT);
GETSET_BOOLEAN_BIT(trap_invalid, TRAP_INVALID);
GETSET_BOOLEAN_BIT(trap_erange, TRAP_ERANGE);
GETSET_BOOLEAN_BIT(trap_divzero, TRAP_DIVZERO);
GETSET_BOOLEAN(allow_complex)
GETSET_BOOLEAN(rational_division)
GETSET_BOOLEAN(allow_release_gil)

static PyObject *
GMPy_CTXT_Get_precision(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->ctx.mpfr_prec));
}

static int
GMPy_CTXT_Set_precision(CTXT_Object *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("precision must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    /* A return value of -1 indicates an error has occurred. Since -1 is not
     * a legal value, we don't specifically check for an error condition.
     */
    if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for precision");
        return -1;
    }
    self->ctx.mpfr_prec = (mpfr_prec_t)temp;
    return 0;
}

static PyObject *
GMPy_CTXT_Get_real_prec(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_REAL_PREC(self)));
}

static int
GMPy_CTXT_Set_real_prec(CTXT_Object *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("real_prec must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    /* A return value of -1 indicates an error has occurred. Since -1 is not
     * a legal value, we don't specifically check for an error condition.
     */
    if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for real_prec");
        return -1;
    }
    self->ctx.real_prec = (mpfr_prec_t)temp;
    return 0;
}

static PyObject *
GMPy_CTXT_Get_imag_prec(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_IMAG_PREC(self)));
}

static int
GMPy_CTXT_Set_imag_prec(CTXT_Object *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("imag_prec must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    /* A return value of -1 indicates an error has occurred. Since -1 is not
     * a legal value, we don't specifically check for an error condition.
     */
    if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for imag_prec");
        return -1;
    }
    self->ctx.imag_prec = (mpfr_prec_t)temp;
    return 0;
}

static PyObject *
GMPy_CTXT_Get_round(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromLong((long)(self->ctx.mpfr_round));
}

static int
GMPy_CTXT_Set_round(CTXT_Object *self, PyObject *value, void *closure)
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
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
           to MPFR_RNDN. */
        self->ctx.real_round = MPFR_RNDN;
        self->ctx.imag_round = MPFR_RNDN;
    }
    else {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    return 0;
}

static PyObject *
GMPy_CTXT_Get_real_round(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_REAL_ROUND(self));
}

static int
GMPy_CTXT_Set_real_round(CTXT_Object *self, PyObject *value, void *closure)
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
GMPy_CTXT_Get_imag_round(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_IMAG_ROUND(self));
}

static int
GMPy_CTXT_Set_imag_round(CTXT_Object *self, PyObject *value, void *closure)
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

static PyObject *
GMPy_CTXT_Get_emin(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromLong(self->ctx.emin);
}

static int
GMPy_CTXT_Set_emin(CTXT_Object *self, PyObject *value, void *closure)
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
    if (!(exp >= MPFR_EMIN_MIN && exp <= MPFR_EMIN_MAX)) {
        VALUE_ERROR("requested minimum exponent is invalid");
        return -1;
    }
    self->ctx.emin = exp;
    return 0;
}

static PyObject *
GMPy_CTXT_Get_emax(CTXT_Object *self, void *closure)
{
    return PyIntOrLong_FromLong(self->ctx.emax);
}

static int
GMPy_CTXT_Set_emax(CTXT_Object *self, PyObject *value, void *closure)
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
    if (!(exp >= MPFR_EMAX_MIN && exp <= MPFR_EMAX_MAX)) {
        VALUE_ERROR("requested maximum exponent is invalid");
        return -1;
    }
    self->ctx.emax = exp;
    return 0;
}

#define ADD_GETSET(NAME) \
    {#NAME, \
        (getter)GMPy_CTXT_Get_##NAME, \
        (setter)GMPy_CTXT_Set_##NAME, NULL, NULL}

static PyGetSetDef GMPyContext_getseters[] = {
    ADD_GETSET(precision),
    ADD_GETSET(real_prec),
    ADD_GETSET(imag_prec),
    ADD_GETSET(round),
    ADD_GETSET(real_round),
    ADD_GETSET(imag_round),
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
    ADD_GETSET(allow_complex),
    ADD_GETSET(rational_division),
    ADD_GETSET(allow_release_gil),
    {NULL}
};

static PyMethodDef GMPyContext_methods[] =
{
    { "abs", GMPy_Context_Abs, METH_O, GMPy_doc_context_abs },
    { "acos", GMPy_Context_Acos, METH_O, GMPy_doc_context_acos },
    { "acosh", GMPy_Context_Acosh, METH_O, GMPy_doc_context_acosh },
    { "add", GMPy_Context_Add, METH_VARARGS, GMPy_doc_context_add },
    { "agm", GMPy_Context_AGM, METH_VARARGS, GMPy_doc_context_agm },
    { "ai", GMPy_Context_Ai, METH_O, GMPy_doc_context_ai },
    { "asin", GMPy_Context_Asin, METH_O, GMPy_doc_context_asin },
    { "asinh", GMPy_Context_Asinh, METH_O, GMPy_doc_context_asinh },
    { "atan", GMPy_Context_Atan, METH_O, GMPy_doc_context_atan },
    { "atanh", GMPy_Context_Atanh, METH_O, GMPy_doc_context_atanh },
    { "atan2", GMPy_Context_Atan2, METH_VARARGS, GMPy_doc_context_atan2 },
    { "clear_flags", GMPy_CTXT_Clear_Flags, METH_NOARGS, GMPy_doc_context_clear_flags },
    { "cbrt", GMPy_Context_Cbrt, METH_O, GMPy_doc_context_cbrt },
    { "ceil", GMPy_Context_Ceil, METH_O, GMPy_doc_context_ceil },
    { "check_range", GMPy_Context_CheckRange, METH_O, GMPy_doc_context_check_range },
    { "const_catalan", GMPy_Context_Const_Catalan, METH_NOARGS, GMPy_doc_context_const_catalan },
    { "const_euler", GMPy_Context_Const_Euler, METH_NOARGS, GMPy_doc_context_const_euler },
    { "const_log2", GMPy_Context_Const_Log2, METH_NOARGS, GMPy_doc_context_const_log2 },
    { "const_pi", GMPy_Context_Const_Pi, METH_NOARGS, GMPy_doc_context_const_pi },
    { "cos", GMPy_Context_Cos, METH_O, GMPy_doc_context_cos },
    { "cosh", GMPy_Context_Cosh, METH_O, GMPy_doc_context_cosh },
    { "cot", GMPy_Context_Cot, METH_O, GMPy_doc_context_cot },
    { "coth", GMPy_Context_Coth, METH_O, GMPy_doc_context_coth },
    { "copy", GMPy_CTXT_Copy, METH_NOARGS, GMPy_doc_context_copy },
    { "csc", GMPy_Context_Csc, METH_O, GMPy_doc_context_csc },
    { "csch", GMPy_Context_Csch, METH_O, GMPy_doc_context_csch },
    { "degrees", GMPy_Context_Degrees, METH_O, GMPy_doc_context_degrees },
    { "digamma", GMPy_Context_Digamma, METH_O, GMPy_doc_context_digamma },
    { "div", GMPy_Context_TrueDiv, METH_VARARGS, GMPy_doc_context_truediv },
    { "divmod", GMPy_Context_DivMod, METH_VARARGS, GMPy_doc_context_divmod },
    { "div_2exp", GMPy_Context_Div_2exp, METH_VARARGS, GMPy_doc_context_div_2exp },
    { "eint", GMPy_Context_Eint, METH_O, GMPy_doc_context_eint },
    { "erf", GMPy_Context_Erf, METH_O, GMPy_doc_context_erf },
    { "erfc", GMPy_Context_Erfc, METH_O, GMPy_doc_context_erfc },
    { "exp", GMPy_Context_Exp, METH_O, GMPy_doc_context_exp },
    { "expm1", GMPy_Context_Expm1, METH_O, GMPy_doc_context_expm1 },
    { "exp10", GMPy_Context_Exp10, METH_O, GMPy_doc_context_exp10 },
    { "exp2", GMPy_Context_Exp2, METH_O, GMPy_doc_context_exp2 },
    { "factorial", GMPy_Context_Factorial, METH_O, GMPy_doc_context_factorial },
    { "floor", GMPy_Context_Floor, METH_O, GMPy_doc_context_floor },
    { "floor_div", GMPy_Context_FloorDiv, METH_VARARGS, GMPy_doc_context_floordiv },
    { "fma", GMPy_Context_FMA, METH_VARARGS, GMPy_doc_context_fma },
    { "fms", GMPy_Context_FMS, METH_VARARGS, GMPy_doc_context_fms },
#if MPFR_VERSION_MAJOR > 3
    { "fmma", GMPy_Context_FMMA, METH_VARARGS, GMPy_doc_context_fmma },
    { "fmms", GMPy_Context_FMMS, METH_VARARGS, GMPy_doc_context_fmms },
#endif
    { "fmod", GMPy_Context_Fmod, METH_VARARGS, GMPy_doc_context_fmod },
    { "frac", GMPy_Context_Frac, METH_O, GMPy_doc_context_frac },
    { "frexp", GMPy_Context_Frexp, METH_O, GMPy_doc_context_frexp },
    { "fsum", GMPy_Context_Fsum, METH_O, GMPy_doc_context_fsum },
    { "gamma", GMPy_Context_Gamma, METH_O, GMPy_doc_context_gamma },
    { "hypot", GMPy_Context_Hypot, METH_VARARGS, GMPy_doc_context_hypot },
    { "is_finite", GMPy_Context_Is_Finite, METH_O, GMPy_doc_context_is_finite },
    { "is_infinite", GMPy_Context_Is_Infinite, METH_O, GMPy_doc_context_is_infinite },
    { "is_integer", GMPy_Context_Is_Integer, METH_O, GMPy_doc_context_is_integer },
    { "is_nan", GMPy_Context_Is_NAN, METH_O, GMPy_doc_context_is_nan },
    { "is_regular", GMPy_Context_Is_Regular, METH_O, GMPy_doc_context_is_regular },
    { "is_signed", GMPy_Context_Is_Signed, METH_O, GMPy_doc_context_is_signed },
    { "is_zero", GMPy_Context_Is_Zero, METH_O, GMPy_doc_context_is_zero },
    { "jn", GMPy_Context_Jn, METH_VARARGS, GMPy_doc_context_jn },
    { "j0", GMPy_Context_J0, METH_O, GMPy_doc_context_j0 },
    { "j1", GMPy_Context_J1, METH_O, GMPy_doc_context_j1 },
    { "li2", GMPy_Context_Li2, METH_O, GMPy_doc_context_li2 },
    { "lgamma", GMPy_Context_Lgamma, METH_O, GMPy_doc_context_lgamma },
    { "lngamma", GMPy_Context_Lngamma, METH_O, GMPy_doc_context_lngamma },
    { "log", GMPy_Context_Log, METH_O, GMPy_doc_context_log },
    { "log10", GMPy_Context_Log10, METH_O, GMPy_doc_context_log10 },
    { "log1p", GMPy_Context_Log1p, METH_O, GMPy_doc_context_log1p },
    { "log2", GMPy_Context_Log2, METH_O, GMPy_doc_context_log2 },
    { "maxnum", GMPy_Context_Maxnum, METH_VARARGS, GMPy_doc_context_maxnum },
    { "minnum", GMPy_Context_Minnum, METH_VARARGS, GMPy_doc_context_minnum },
    { "minus", GMPy_Context_Minus, METH_VARARGS, GMPy_doc_context_minus },
    { "mod", GMPy_Context_Mod, METH_VARARGS, GMPy_doc_context_mod },
    { "modf", GMPy_Context_Modf, METH_O, GMPy_doc_context_modf },
    { "mul", GMPy_Context_Mul, METH_VARARGS, GMPy_doc_context_mul },
    { "mul_2exp", GMPy_Context_Mul_2exp, METH_VARARGS, GMPy_doc_context_mul_2exp },
    { "next_above", GMPy_Context_NextAbove, METH_O, GMPy_doc_context_next_above },
    { "next_below", GMPy_Context_NextBelow, METH_O, GMPy_doc_context_next_below },
    { "next_toward", GMPy_Context_NextToward, METH_VARARGS, GMPy_doc_context_next_toward },
    { "norm", GMPy_Context_Norm, METH_O, GMPy_doc_context_norm },
    { "phase", GMPy_Context_Phase, METH_O, GMPy_doc_context_phase },
    { "plus", GMPy_Context_Plus, METH_VARARGS, GMPy_doc_context_plus },
    { "polar", GMPy_Context_Polar, METH_O, GMPy_doc_context_polar },
    { "proj", GMPy_Context_Proj, METH_O, GMPy_doc_context_proj },
    { "pow", GMPy_Context_Pow, METH_VARARGS, GMPy_doc_context_pow },
    { "radians", GMPy_Context_Radians, METH_O, GMPy_doc_context_radians },
    { "rect", GMPy_Context_Rect, METH_VARARGS, GMPy_doc_context_rect },
    { "rec_sqrt", GMPy_Context_RecSqrt, METH_O, GMPy_doc_context_rec_sqrt },
    { "reldiff", GMPy_Context_RelDiff, METH_VARARGS, GMPy_doc_context_reldiff },
    { "remainder", GMPy_Context_Remainder, METH_VARARGS, GMPy_doc_context_remainder },
    { "remquo", GMPy_Context_RemQuo, METH_VARARGS, GMPy_doc_context_remquo },
    { "rint", GMPy_Context_Rint, METH_O, GMPy_doc_context_rint },
    { "rint_ceil", GMPy_Context_RintCeil, METH_O, GMPy_doc_context_rint_ceil },
    { "rint_floor", GMPy_Context_RintFloor, METH_O, GMPy_doc_context_rint_floor },
    { "rint_round", GMPy_Context_RintRound, METH_O, GMPy_doc_context_rint_round },
    { "rint_trunc", GMPy_Context_RintTrunc, METH_O, GMPy_doc_context_rint_trunc },
    { "root", GMPy_Context_Root, METH_VARARGS, GMPy_doc_context_root },
    { "rootn", GMPy_Context_Rootn, METH_VARARGS, GMPy_doc_context_rootn },
#ifdef MPC_110
    { "root_of_unity", GMPy_Context_Root_Of_Unity, METH_VARARGS, GMPy_doc_context_root_of_unity },
#endif
    { "round2", GMPy_Context_Round2, METH_VARARGS, GMPy_doc_context_round2 },
    { "round_away", GMPy_Context_RoundAway, METH_O, GMPy_doc_context_round_away },
    { "sec", GMPy_Context_Sec, METH_O, GMPy_doc_context_sec },
    { "sech", GMPy_Context_Sech, METH_O, GMPy_doc_context_sech },
    { "sin", GMPy_Context_Sin, METH_O, GMPy_doc_context_sin },
    { "sin_cos", GMPy_Context_Sin_Cos, METH_O, GMPy_doc_context_sin_cos },
    { "sinh", GMPy_Context_Sinh, METH_O, GMPy_doc_context_sinh },
    { "sinh_cosh", GMPy_Context_Sinh_Cosh, METH_O, GMPy_doc_context_sinh_cosh },
    { "sqrt", GMPy_Context_Sqrt, METH_O, GMPy_doc_context_sqrt },
    { "square", GMPy_Context_Square, METH_O, GMPy_doc_context_square },
    { "sub", GMPy_Context_Sub, METH_VARARGS, GMPy_doc_context_sub },
    { "tan", GMPy_Context_Tan, METH_O, GMPy_doc_context_tan },
    { "tanh", GMPy_Context_Tanh, METH_O, GMPy_doc_context_tanh },
    { "trunc", GMPy_Context_Trunc, METH_O, GMPy_doc_context_trunc },
#ifdef VECTOR
    { "vector", GMPy_Context_Vector, METH_O, GMPy_doc_context_vector },
    { "vector2", GMPy_Context_Vector2, METH_VARARGS, GMPy_doc_context_vector2 },
#endif
    { "yn", GMPy_Context_Yn, METH_VARARGS, GMPy_doc_context_yn },
    { "y0", GMPy_Context_Y0, METH_O, GMPy_doc_context_y0 },
    { "y1", GMPy_Context_Y1, METH_O, GMPy_doc_context_y1 },
    { "zeta", GMPy_Context_Zeta, METH_O, GMPy_doc_context_zeta },
    { "__enter__", GMPy_CTXT_Enter, METH_NOARGS, NULL },
    { "__exit__", GMPy_CTXT_Exit, METH_VARARGS, NULL },
    { NULL, NULL, 1 }
};

static PyTypeObject CTXT_Type =
{
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "gmpy2 context",                        /* tp_name          */
    sizeof(CTXT_Object),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    (destructor) GMPy_CTXT_Dealloc,         /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_CTXT_Repr_Slot,         /* tp_repr          */
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
    { "__enter__", GMPy_CTXT_Manager_Enter, METH_NOARGS, NULL },
    { "__exit__", GMPy_CTXT_Manager_Exit, METH_VARARGS, NULL },
    { NULL, NULL, 1 }
};

static PyTypeObject CTXT_Manager_Type =
{
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                   /* ob_size          */
#endif
    "gmpy2 context",                         /* tp_name          */
    sizeof(CTXT_Manager_Object),             /* tp_basicsize     */
        0,                                   /* tp_itemsize      */
    (destructor) GMPy_CTXT_Manager_Dealloc,  /* tp_dealloc       */
        0,                                   /* tp_print         */
        0,                                   /* tp_getattr       */
        0,                                   /* tp_setattr       */
        0,                                   /* tp_reserved      */
    (reprfunc) GMPy_CTXT_Manager_Repr_Slot,  /* tp_repr          */
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
