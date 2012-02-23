/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_context.c                                                          *
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

static GMPyContextObject *
GMPyContext_new(void)
{
    GMPyContextObject *result;

    if ((result = PyObject_New(GMPyContextObject, &GMPyContext_Type))) {
        result->now.mpfr_prec = DBL_MANT_DIG;
#ifdef WITHMPC
        result->now.mpc_rprec = -1;
        result->now.mpc_iprec = -1;
#endif
        result->now.mpfr_round = MPFR_RNDN;
#ifdef WITHMPC
        result->now.mpc_rround = -1;
        result->now.mpc_iround = -1;
#endif
        result->now.emax = mpfr_get_emax();
        result->now.emin = mpfr_get_emin();
        result->now.subnormalize = 0;
        result->now.underflow = 0;
        result->now.overflow = 0;
        result->now.inexact = 0;
        result->now.invalid = 0;
        result->now.erange = 0;
        result->now.divzero = 0;
        result->now.trap_underflow = 0;
        result->now.trap_overflow = 0;
        result->now.trap_inexact = 0;
        result->now.trap_invalid = 0;
        result->now.trap_erange = 0;
        result->now.trap_divzero = 0;
        result->now.trap_expbound = 0;
#ifdef WITHMPC
        result->now.allow_complex = 0;
#endif
        result->orig = NULL;
    }
    return result;
};

static void
GMPyContext_dealloc(GMPyContextObject *self)
{
    Py_XDECREF((PyObject*)self->orig);
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
            "context(precision=%s, mpc_rprec=%s, mpc_iprec=%s,\n"
            "        round=%s, mpc_rround=%s, mpc_iround=%s,\n"
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

    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->now.mpfr_prec));
#ifdef WITHMPC
    if (self->now.mpc_rprec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, i++, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->now.mpc_rprec));
    if (self->now.mpc_iprec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, i++, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->now.mpc_iprec));
#endif
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->now.mpfr_round));
#ifdef WITHMPC
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->now.mpc_rround));
    PyTuple_SET_ITEM(tuple, i++, _round_to_name(self->now.mpc_iround));
#endif
    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->now.emax));
    PyTuple_SET_ITEM(tuple, i++, PyIntOrLong_FromLong(self->now.emin));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.subnormalize));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_underflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.underflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_overflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.overflow));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_inexact));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.inexact));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_invalid));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.invalid));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_erange));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.erange));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_divzero));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.divzero));
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.trap_expbound));
#ifdef WITHMPC
    PyTuple_SET_ITEM(tuple, i++, PyBool_FromLong(self->now.allow_complex));
#endif

    if (!PyErr_Occurred())
        result = Py2or3String_Format(format, tuple);
    else
        SYSTEM_ERROR("internal error in GMPyContext_repr");

    Py_DECREF(format);
    Py_DECREF(tuple);
    return result;
};

PyDoc_STRVAR(doc_get_context,
"get_context() -> context manager\n\n"
"Return a copy of the current context manager.");

static PyObject *
Pygmpy_get_context(PyObject *self, PyObject *args)
{
    GMPyContextObject *result;

    if ((result = GMPyContext_new()))
        result->now = context->now;

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_local_context,
"local_context() -> context manager\n\n"
"Return a reference to the current context manager controlling\n"
"MPFR and MPC arithmetic. The returned object no longer refers to\n"
"the current context after a new context is loaded using\n"
"set_context() or 'with <context>'. See 'context()' for the\n"
"description of context options.");

static PyObject *
Pygmpy_local_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    gmpy_context old;
#ifdef WITHMPC
    static char *kwlist[] = {
        "precision", "mpc_rprec", "mpc_iprec", "round",
        "mpc_rround", "mpc_iround", "emax", "emin", "subnormalize",
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
        VALUE_ERROR("local_context() only supports keyword arguments");
        return NULL;
    }

    /* Temporarily save the options of the global context manager in case
     * there is an error while processing the arguments.
     */
     
    old = context->now;

#ifdef WITHMPC
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|llliiilliiiiiiiii", kwlist,
            &context->now.mpfr_prec,
            &context->now.mpc_rprec,
            &context->now.mpc_iprec,
            &context->now.mpfr_round,
            &context->now.mpc_rround,
            &context->now.mpc_iround,
            &context->now.emax,
            &context->now.emin,
            &context->now.subnormalize,
            &context->now.trap_underflow,
            &context->now.trap_overflow,
            &context->now.trap_inexact,
            &context->now.trap_invalid,
            &context->now.trap_erange,
            &context->now.trap_divzero,
            &context->now.trap_expbound,
            &context->now.allow_complex))) {
#else
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|lilliiiiiiii", kwlist,
            &context->now.mpfr_prec,
            &context->now.mpfr_round,
            &context->now.emax,
            &context->now.emin,
            &context->now.subnormalize,
            &context->now.trap_underflow,
            &context->now.trap_overflow,
            &context->now.trap_inexact,
            &context->now.trap_invalid,
            &context->now.trap_erange,
            &context->now.trap_divzero,
            &context->now.trap_expbound))) {
#endif
        VALUE_ERROR("invalid keyword arguments in local_context()");
        return NULL;
    }

    /* Sanity check for values. */
    if (context->now.mpfr_prec < MPFR_PREC_MIN ||
        context->now.mpfr_prec > MPFR_PREC_MAX) {
        context->now = old;
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }

#ifdef WITHMPC
    if (!(context->now.mpc_rprec == GMPY_DEFAULT ||
        (context->now.mpc_rprec >= MPFR_PREC_MIN &&
        context->now.mpc_rprec <= MPFR_PREC_MAX))) {
        context->now = old;
        VALUE_ERROR("invalid value for mpc_rprec");
        return NULL;
    }
    if (!(context->now.mpc_iprec == GMPY_DEFAULT ||
        (context->now.mpc_iprec >= MPFR_PREC_MIN &&
        context->now.mpc_iprec <= MPFR_PREC_MAX))) {
        context->now = old;
        VALUE_ERROR("invalid value for mpc_iprec");
        return NULL;
    }
#endif

    if (!(context->now.mpfr_round == MPFR_RNDN ||
        context->now.mpfr_round == MPFR_RNDZ ||
        context->now.mpfr_round == MPFR_RNDU ||
        context->now.mpfr_round == MPFR_RNDD ||
        context->now.mpfr_round == MPFR_RNDA)) {
        context->now = old;
        VALUE_ERROR("invalid value for mpfr_round");
        return NULL;
    }

#ifdef WITHMPC
    if (context->now.mpfr_round == MPFR_RNDA) {
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
         * to MPFR_RNDN.
         */
        context->now.mpc_rround = MPFR_RNDN;
        context->now.mpc_iround = MPFR_RNDN;
    }
    if (!(context->now.mpc_rround == MPFR_RNDN ||
        context->now.mpc_rround == MPFR_RNDZ ||
        context->now.mpc_rround == MPFR_RNDU ||
        context->now.mpc_rround == MPFR_RNDD ||
        context->now.mpc_rround == GMPY_DEFAULT)) {
        context->now = old;
        VALUE_ERROR("invalid value for mpc_rround");
        return NULL;
    }
    if (!(context->now.mpc_iround == MPFR_RNDN ||
        context->now.mpc_iround == MPFR_RNDZ ||
        context->now.mpc_iround == MPFR_RNDU ||
        context->now.mpc_iround == MPFR_RNDD ||
        context->now.mpc_iround == GMPY_DEFAULT)) {
        context->now = old;
        VALUE_ERROR("invalid value for mpc_iround");
        return NULL;
    }
#endif

    if (!(context->now.emin < 0 && context->now.emax > 0)) {
        VALUE_ERROR("invalid values for emin and/or emax");
        context->now = old;
        return NULL;
    }

    if (mpfr_set_emin(context->now.emin)) {
        VALUE_ERROR("invalid value for emin");
        context->now = old;
        return NULL;
    }
    if (mpfr_set_emax(context->now.emax)) {
        VALUE_ERROR("invalid value for emax");
        context->now = old;
        return NULL;
    }

    /* The values in the global context manager have been changed. Also return
     * another reference to that context manager.
     */
     
    Py_INCREF((PyObject*)context);
    return (PyObject*)context;
}

#ifdef WITHMPC

PyDoc_STRVAR(doc_context,
"context() -> context manager\n\n"
"Return a new context manager controlling MPFR and MPC arithmetic. Options\n"
"can only be specified as keyword arguments. \n\n"
"    precision:      precision, in bits, of an MPFR result\n"
"    mpc_rprec:      precision, in bits, of Re(MPC)\n"
"                      -1 implies use mpfr_prec\n"
"    mpc_iprec:      precision, in bits, of Im(MPC)\n"
"                      -1 implies use mpc_rprec\n"
"    round:          rounding mode for MPFR\n"
"    mpc_rround:     rounding mode for Re(MPC)\n"
"                      -1 implies use mpfr_round\n"
"    mpc_iround:     rounding mode for Im(MPC)\n"
"                      -1 implies use mpc_rround\n"
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
"                        coerced to either 0 or Infinity\n"
"    allow_complex:  if True, allow mpfr functions to return mpc\n"
"                    if False, mpfr functions cannot return an mpc\n");

#else

PyDoc_STRVAR(doc_context,
"context() -> context\n\n"
"Return a new context manager controlling MPFR arithmetic. Options can\n"
"only be specified as keyword arguments. Options are also available as\n"
"instance attributes.\n\n"
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
Pygmpy_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPyContextObject *result;

#ifdef WITHMPC
    static char *kwlist[] = {
        "precision", "mpc_rprec", "mpc_iprec", "round",
        "mpc_rround", "mpc_iround", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero",
        "allow_complex", NULL };
#else
    static char *kwlist[] = {
        "precision", "round", "emax", "emin", "subnormalize",
        "trap_underflow", "trap_overflow", "trap_inexact",
        "trap_invalid", "trap_erange", "trap_divzero", NULL };
#endif

    if (PyTuple_GET_SIZE(args)) {
        VALUE_ERROR("context() only supports keyword arguments");
        return NULL;
    }

    if (!(result = GMPyContext_new()))
        return NULL;

#ifdef WITHMPC
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|llliiilliiiiiiiii", kwlist,
            &result->now.mpfr_prec,
            &result->now.mpc_rprec,
            &result->now.mpc_iprec,
            &result->now.mpfr_round,
            &result->now.mpc_rround,
            &result->now.mpc_iround,
            &result->now.emax,
            &result->now.emin,
            &result->now.subnormalize,
            &result->now.trap_underflow,
            &result->now.trap_overflow,
            &result->now.trap_inexact,
            &result->now.trap_invalid,
            &result->now.trap_erange,
            &result->now.trap_divzero,
            &result->now.trap_expbound,
            &result->now.allow_complex))) {
#else
    if (!(PyArg_ParseTupleAndKeywords(args, kwargs,
            "|lilliiiiiiii", kwlist,
            &result->now.mpfr_prec,
            &result->now.mpfr_round,
            &result->now.emax,
            &result->now.emin,
            &result->now.subnormalize,
            &result->now.trap_underflow,
            &result->now.trap_overflow,
            &result->now.trap_inexact,
            &result->now.trap_invalid,
            &result->now.trap_erange,
            &result->now.trap_divzero,
            &result->now.trap_expbound))) {
#endif
        VALUE_ERROR("invalid keyword arguments in context()");
        return NULL;
    }

    /* Sanity check for values. */
    if (result->now.mpfr_prec < MPFR_PREC_MIN ||
        result->now.mpfr_prec > MPFR_PREC_MAX) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for precision");
        return NULL;
    }

#ifdef WITHMPC
    if (!(result->now.mpc_rprec == GMPY_DEFAULT ||
        (result->now.mpc_rprec >= MPFR_PREC_MIN &&
        result->now.mpc_rprec <= MPFR_PREC_MAX))) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for mpc_rprec");
        return NULL;
    }
    if (!(result->now.mpc_iprec == GMPY_DEFAULT ||
        (result->now.mpc_iprec >= MPFR_PREC_MIN &&
        result->now.mpc_iprec <= MPFR_PREC_MAX))) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for mpc_iprec");
        return NULL;
    }
#endif

    if (!(result->now.mpfr_round == MPFR_RNDN ||
        result->now.mpfr_round == MPFR_RNDZ ||
        result->now.mpfr_round == MPFR_RNDU ||
        result->now.mpfr_round == MPFR_RNDD ||
        result->now.mpfr_round == MPFR_RNDA)) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for mpfr_round");
        return NULL;
    }

#ifdef WITHMPC
    if (result->now.mpfr_round == MPFR_RNDA) {
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
           to MPFR_RNDN. */
        result->now.mpc_rround = MPFR_RNDN;
        result->now.mpc_iround = MPFR_RNDN;
    }
    if (!(result->now.mpc_rround == MPFR_RNDN ||
        result->now.mpc_rround == MPFR_RNDZ ||
        result->now.mpc_rround == MPFR_RNDU ||
        result->now.mpc_rround == MPFR_RNDD ||
        result->now.mpc_rround == GMPY_DEFAULT)) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for mpc_rround");
        return NULL;
    }
    if (!(result->now.mpc_iround == MPFR_RNDN ||
        result->now.mpc_iround == MPFR_RNDZ ||
        result->now.mpc_iround == MPFR_RNDU ||
        result->now.mpc_iround == MPFR_RNDD ||
        result->now.mpc_iround == GMPY_DEFAULT)) {
        Py_DECREF((PyObject*)result);
        VALUE_ERROR("invalid value for mpc_iround");
        return NULL;
    }
#endif

    if (!(result->now.emin < 0 && result->now.emax > 0)) {
        VALUE_ERROR("invalid values for emin and/or emax");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    if (mpfr_set_emin(result->now.emin)) {
        VALUE_ERROR("invalid value for emin");
        Py_DECREF((PyObject*)result);
        return NULL;
    }
    if (mpfr_set_emax(result->now.emax)) {
        VALUE_ERROR("invalid value for emax");
        Py_DECREF((PyObject*)result);
        return NULL;
    }

    result->now.underflow = 0;
    result->now.overflow = 0;
    result->now.inexact = 0;
    result->now.invalid = 0;
    result->now.erange = 0;
    result->now.divzero = 0;
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_set_context,
"set_context(context)\n\n"
"Activate a context manager controlling MPFR and MPC arithmetic.\n");

static PyObject *
Pygmpy_set_context(PyObject *self, PyObject *other)
{
    if (GMPyContext_Check(other)) {
        Py_INCREF((PyObject*)other);
        Py_DECREF((PyObject*)context);
        context = (GMPyContextObject*)other;
        mpfr_set_emin(context->now.emin);
        mpfr_set_emax(context->now.emax);
        Py_RETURN_NONE;
    }
    else {
        VALUE_ERROR("set_context() requires a context argument");
        return NULL;
    }
}

static PyObject *
GMPyContext_enter(PyObject *self, PyObject *args)
{
    GMPyContextObject *result, *save;

    if (((GMPyContextObject*)self)->orig != NULL) {
        SYSTEM_ERROR("Internal error in GMPyContext_enter");
        return NULL;
    }

    if ((save = GMPyContext_new())) {
        save->now = context->now;
        save->orig = NULL;
    }
    else {
        SYSTEM_ERROR("Internal error in GMPyContext_enter");
        return NULL;
    }
        

    ((GMPyContextObject*)self)->orig = (PyObject*)save;
    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)self;
    mpfr_set_emin(context->now.emin);
    mpfr_set_emax(context->now.emax);
    Py_INCREF((PyObject*)context);
    result = context;
    Py_INCREF((PyObject*)result);
    return (PyObject*)result;
}

static PyObject *
GMPyContext_exit(PyObject *self, PyObject *args)
{
    if (((GMPyContextObject*)self)->orig == NULL) {
        SYSTEM_ERROR("Internal error in GMPyContext_exit");
        return NULL;
    }

    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)(((GMPyContextObject*)self)->orig);
    mpfr_set_emin(context->now.emin);
    mpfr_set_emax(context->now.emax);
    ((GMPyContextObject*)self)->orig = NULL;
    Py_RETURN_NONE;
}

PyDoc_STRVAR(doc_context_clear_flags,
"clear_flags()\n\n"
"Clear all MPFR exception flags.");
static PyObject *
GMPyContext_clear_flags(PyObject *self, PyObject *args)
{
    ((GMPyContextObject*)self)->now.underflow = 0;
    ((GMPyContextObject*)self)->now.overflow = 0;
    ((GMPyContextObject*)self)->now.inexact = 0;
    ((GMPyContextObject*)self)->now.invalid = 0;
    ((GMPyContextObject*)self)->now.erange = 0;
    ((GMPyContextObject*)self)->now.divzero = 0;
    Py_RETURN_NONE;
}

/* Define the get/set functions. */

#define GETSET_BOOLEAN(NAME) \
static PyObject * \
GMPyContext_get_##NAME(GMPyContextObject *self, void *closure) \
{ \
    return PyBool_FromLong(self->now.NAME); \
}; \
static int \
GMPyContext_set_##NAME(GMPyContextObject *self, PyObject *value, void *closure) \
{ \
    if (!(PyBool_Check(value))) { \
        TYPE_ERROR(#NAME " must be True or False"); \
        return -1; \
    } \
    self->now.NAME = (value == Py_True) ? 1 : 0; \
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
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.mpfr_prec));
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
    self->now.mpfr_prec = (mpfr_prec_t)temp;
    return 0;
}

#ifdef WITHMPC
static PyObject *
GMPyContext_get_mpc_rprec(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_MPC_RPREC(self)));
}

static int
GMPyContext_set_mpc_rprec(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("mpc_rprec must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    if (temp == -1) {
        if (PyErr_Occurred()) {
            VALUE_ERROR("invalid value for mpc_rprec");
            return -1;
        }
    }
    else if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for mpc_rprec");
        return -1;
    }
    self->now.mpc_rprec = (mpfr_prec_t)temp;
    return 0;
}

static PyObject *
GMPyContext_get_mpc_iprec(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_MPC_IPREC(self)));
}

static int
GMPyContext_set_mpc_iprec(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t temp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("mpc_iprec must be Python integer");
        return -1;
    }
    temp = PyIntOrLong_AsSsize_t(value);
    if (temp == -1) {
        if (PyErr_Occurred()) {
            VALUE_ERROR("invalid value for mpc_iprec");
            return -1;
        }
    }
    else if (temp < MPFR_PREC_MIN || temp > MPFR_PREC_MAX) {
        VALUE_ERROR("invalid value for mpc_iprec");
        return -1;
    }
    self->now.mpc_iprec = (mpfr_prec_t)temp;
    return 0;
}
#endif

static PyObject *
GMPyContext_get_round(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)(self->now.mpfr_round));
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
        self->now.mpfr_round = temp;
    else if (temp == MPFR_RNDZ)
        self->now.mpfr_round = temp;
    else if (temp == MPFR_RNDU)
        self->now.mpfr_round = temp;
    else if (temp == MPFR_RNDD)
        self->now.mpfr_round = temp;
    else if (temp == MPFR_RNDA) {
        self->now.mpfr_round = temp;
#ifdef WITHMPC
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
           to MPFR_RNDN. */
        self->now.mpc_rround = MPFR_RNDN;
        self->now.mpc_iround = MPFR_RNDN;
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
GMPyContext_get_mpc_rround(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_MPC_RROUND(self));
}

static int
GMPyContext_set_mpc_rround(GMPyContextObject *self, PyObject *value, void *closure)
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
    if (temp == GMPY_DEFAULT)
        self->now.mpc_rround = temp;
    else if (temp == MPFR_RNDN)
        self->now.mpc_rround = temp;
    else if (temp == MPFR_RNDZ)
        self->now.mpc_rround = temp;
    else if (temp == MPFR_RNDU)
        self->now.mpc_rround = temp;
    else if (temp == MPFR_RNDD)
        self->now.mpc_rround = temp;
    else {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    return 0;
}

static PyObject *
GMPyContext_get_mpc_iround(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_MPC_IROUND(self));
}

static int
GMPyContext_set_mpc_iround(GMPyContextObject *self, PyObject *value, void *closure)
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
    if (temp == GMPY_DEFAULT)
        self->now.mpc_iround = temp;
    else if (temp == MPFR_RNDN)
        self->now.mpc_iround = temp;
    else if (temp == MPFR_RNDZ)
        self->now.mpc_iround = temp;
    else if (temp == MPFR_RNDU)
        self->now.mpc_iround = temp;
    else if (temp == MPFR_RNDD)
        self->now.mpc_iround = temp;
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
    return PyIntOrLong_FromLong(self->now.emin);
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
    self->now.emin = exp;
    mpfr_set_emin(exp);
    return 0;
}

static PyObject *
GMPyContext_get_emax(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromLong(self->now.emax);
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
    self->now.emax = exp;
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
    ADD_GETSET(mpc_rprec),
    ADD_GETSET(mpc_iprec),
#endif
    ADD_GETSET(round),
#ifdef WITHMPC
    ADD_GETSET(mpc_rround),
    ADD_GETSET(mpc_iround),
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
    { "__enter__", GMPyContext_enter, METH_NOARGS, NULL },
    { "__exit__", GMPyContext_exit, METH_VARARGS, NULL },
    { "clear_flags", GMPyContext_clear_flags, METH_NOARGS,
            doc_context_clear_flags },
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
    "GMPY2 Context manager",                /* tp_doc           */
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
