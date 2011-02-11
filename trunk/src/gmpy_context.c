/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_context.c                                                          *
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

/* Create a context manager type. */

#define GET_MPC_RPREC(c) ((c->now.mpc_rprec==GMPY_DEFAULT)?c->now.mpfr_prec:c->now.mpc_rprec)
#define GET_MPC_IPREC(c) ((c->now.mpc_iprec==GMPY_DEFAULT)?GET_MPC_RPREC(c):c->now.mpc_iprec)

#define GET_MPC_RROUND(c) ((c->now.mpc_rround==GMPY_DEFAULT)?c->now.mpfr_round:c->now.mpc_rround)
#define GET_MPC_IROUND(c) ((c->now.mpc_iround==GMPY_DEFAULT)?GET_MPC_RROUND(c):c->now.mpc_iround)
#define GET_MPC_ROUND(c) (RNDC(GET_MPC_RROUND(c), GET_MPC_IROUND(c)))

static GMPyContextObject *
GMPyContext_new(void)
{
    GMPyContextObject *self;

    if ((self = PyObject_New(GMPyContextObject, &GMPyContext_Type))) {
        self->now.nonstop = 0;
        self->now.subnormalize = 0;
        self->now.mpfr_prec = DBL_MANT_DIG;
        self->now.mpc_rprec = -1;
        self->now.mpc_iprec = -1;
        self->now.mpfr_round = MPFR_RNDN;
        self->now.mpc_rround = -1;
        self->now.mpc_iround = -1;
        self->now.mpc_round = MPC_RNDNN;
        self->now.emax = mpfr_get_emax();
        self->now.emin = mpfr_get_emin();
        self->now.underflow = 0;
        self->now.overflow = 0;
        self->now.inexact = 0;
        self->now.invalid = 0;
        self->now.erange = 0;
        self->now.divzero = 0;
        self->now.raise_underflow = 0;
        self->now.raise_overflow = 0;
        self->now.raise_inexact = 0;
        self->now.raise_invalid = 0;
        self->now.raise_erange = 0;
        self->now.raise_divzero = 0;
        self->orig = NULL;
    }
    return self;
};

static void
GMPyContext_dealloc(GMPyContextObject *self)
{
    Py_XDECREF((PyObject*)self->orig);
    PyObject_Del(self);
};

/* Helper value to convert to convert a rounding mode to a string. */

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

    tuple = PyTuple_New(22);
    if (!tuple) return NULL;

    format = Py2or3String_FromString(
            "context(nonstop=%s,\n"
            "        subnormalize=%s,\n"
            "        precision=%s, mpc_rprec=%s, mpc_iprec=%s,\n"
            "        round=%s, mpc_rround=%s, mpc_iround=%s,\n"
            "        emax=%s, emin=%s,\n"
            "        raise_underflow=%s, underflow=%s,\n"
            "        raise_overflow=%s, overflow=%s,\n"
            "        raise_inexact=%s, inexact=%s,\n"
            "        raise_invalid=%s, invalid=%s,\n"
            "        raise_erange=%s, erange=%s,\n"
            "        raise_divzero=%s, divzero=%s)"
            );
    if (!format) return NULL;

    PyTuple_SET_ITEM(tuple, 0, PyBool_FromLong(self->now.nonstop));
    PyTuple_SET_ITEM(tuple, 1, PyBool_FromLong(self->now.subnormalize));
    PyTuple_SET_ITEM(tuple, 2, PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.mpfr_prec)));

    if (self->now.mpc_rprec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, 3, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, 3, PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.mpc_rprec)));

    if (self->now.mpc_iprec == GMPY_DEFAULT)
        PyTuple_SET_ITEM(tuple, 4, Py2or3String_FromString("Default"));
    else
        PyTuple_SET_ITEM(tuple, 4, PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.mpc_iprec)));

    PyTuple_SET_ITEM(tuple, 5, _round_to_name(self->now.mpfr_round));
    PyTuple_SET_ITEM(tuple, 6, _round_to_name(self->now.mpc_rround));
    PyTuple_SET_ITEM(tuple, 7, _round_to_name(self->now.mpc_iround));
    PyTuple_SET_ITEM(tuple, 8, PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.emax)));
    PyTuple_SET_ITEM(tuple, 9, PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.emin)));

    PyTuple_SET_ITEM(tuple, 10, PyBool_FromLong(self->now.raise_underflow));
    PyTuple_SET_ITEM(tuple, 11, PyBool_FromLong(self->now.underflow));
    PyTuple_SET_ITEM(tuple, 12, PyBool_FromLong(self->now.raise_overflow));
    PyTuple_SET_ITEM(tuple, 13, PyBool_FromLong(self->now.overflow));
    PyTuple_SET_ITEM(tuple, 14, PyBool_FromLong(self->now.raise_inexact));
    PyTuple_SET_ITEM(tuple, 15, PyBool_FromLong(self->now.inexact));
    PyTuple_SET_ITEM(tuple, 16, PyBool_FromLong(self->now.raise_invalid));
    PyTuple_SET_ITEM(tuple, 17, PyBool_FromLong(self->now.invalid));
    PyTuple_SET_ITEM(tuple, 18, PyBool_FromLong(self->now.raise_erange));
    PyTuple_SET_ITEM(tuple, 19, PyBool_FromLong(self->now.erange));
    PyTuple_SET_ITEM(tuple, 20, PyBool_FromLong(self->now.raise_divzero));
    PyTuple_SET_ITEM(tuple, 21, PyBool_FromLong(self->now.divzero));

    if (!PyErr_Occurred())
        result = Py2or3String_Format(format, tuple);
    else
        SYSTEM_ERROR("internal error in GMPyContext_repr");

    Py_DECREF(format);
    Py_DECREF(tuple);
    return result;
};

/* Return a reference to the current context. */

PyDoc_STRVAR(doc_current,
"current() -> context\n\n"
"Return a reference to the current context manager controlling\n"
"MPFR and MPC arithmetic. The returned object no longer refers to\n"
"the current context after a new context is loaded using \n"
"set_context() or with ...\n");

static PyObject *
Pygmpy_current(PyObject *self, PyObject *args)
{
    GMPyContextObject *result;

    Py_INCREF((PyObject*)context);
    result = context;
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_get_context,
"get_context() -> context\n\n"
"Return a copy of the current context manager controlling\n"
"MPFR and MPC arithmetic.\n\n"
"    nonstop:      if True, return nan or inf\n"
"                  if False, raise exception\n"
"    subnormalize: if True, subnormalized results can be returned\n"
"    precision:    precision, in bits, of an MPFR result\n"
"    mpc_rprec:    precision, in bits, of Re(MPC)\n"
"                    -1 implies use mpfr_prec\n"
"    mpc_iprec:    precision, in bits, of Im(MPC)\n"
"                    -1 implies use mpc_rprec\n"
"    round:        rounding mode for MPFR\n"
"    mpc_rround:   rounding mode for Re(MPC)\n"
"                    -1 implies use mpfr_round\n"
"    mpc_iround:   rounding mode for Im(MPC)\n"
"                    -1 implies use mpc_rround\n"
"    e_max:        maximum allowed exponent\n"
"    e_min:        minimum allowed exponent\n");

static PyObject *
Pygmpy_get_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    GMPyContextObject *result;

    if ((result = GMPyContext_new())) {
        result->now = context->now;
        result->orig = NULL;
    }
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

    ((GMPyContextObject*)self)->orig = (PyObject*)save;
    Py_DECREF((PyObject*)context);
    context = (GMPyContextObject*)self;
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
    ((GMPyContextObject*)self)->orig = NULL;
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

GETSET_BOOLEAN(nonstop);
GETSET_BOOLEAN(subnormalize);
GETSET_BOOLEAN(underflow);
GETSET_BOOLEAN(overflow);
GETSET_BOOLEAN(inexact);
GETSET_BOOLEAN(invalid);
GETSET_BOOLEAN(erange);
GETSET_BOOLEAN(divzero);
GETSET_BOOLEAN(raise_underflow);
GETSET_BOOLEAN(raise_overflow);
GETSET_BOOLEAN(raise_inexact);
GETSET_BOOLEAN(raise_invalid);
GETSET_BOOLEAN(raise_erange);
GETSET_BOOLEAN(raise_divzero);

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
        /* Since RNDA is not supported for MPC, set the MPC rounding modes
           to MPFR_RNDN. */
        self->now.mpc_rround = MPFR_RNDN;
        self->now.mpc_iround = MPFR_RNDN;
    }
    else {
        VALUE_ERROR("invalid value for round mode");
        return -1;
    }
    return 0;
}

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

static PyObject *
GMPyContext_get_emin(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.emin));
}

static int
GMPyContext_set_emin(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t exp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("emin must be Python integer");
        return -1;
    }
    exp = PyIntOrLong_AsSsize_t(value);
    if (exp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("requested minimum exponent is invalid");
        return -1;
    }
    if (mpfr_set_emin((mpfr_prec_t)exp)) {
        VALUE_ERROR("requested minimum exponent is invalid");
        return -1;
    }
    self->now.emin = exp;
    return 0;
}

static PyObject *
GMPyContext_get_emax(GMPyContextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.emax));
}

static int
GMPyContext_set_emax(GMPyContextObject *self, PyObject *value, void *closure)
{
    Py_ssize_t exp;

    if (!(PyIntOrLong_Check(value))) {
        TYPE_ERROR("emax must be Python integer");
        return -1;
    }
    exp = PyIntOrLong_AsSsize_t(value);
    if (exp == -1 && PyErr_Occurred()) {
        VALUE_ERROR("requested maximum exponent is invalid");
        return -1;
    }
    if (mpfr_set_emax((mpfr_prec_t)exp)) {
        VALUE_ERROR("requested maximum exponent is invalid");
        return -1;
    }
    self->now.emax = exp;
    return 0;
}

#define ADD_GETSET(NAME) \
    {#NAME, \
        (getter)GMPyContext_get_##NAME, \
        (setter)GMPyContext_set_##NAME, NULL, NULL}

static PyGetSetDef GMPyContext_getseters[] = {
    ADD_GETSET(nonstop),
    ADD_GETSET(subnormalize),
    ADD_GETSET(precision),
    ADD_GETSET(mpc_rprec),
    ADD_GETSET(mpc_iprec),
    ADD_GETSET(round),
    ADD_GETSET(mpc_rround),
    ADD_GETSET(mpc_iround),
    ADD_GETSET(emax),
    ADD_GETSET(emin),
    ADD_GETSET(underflow),
    ADD_GETSET(overflow),
    ADD_GETSET(inexact),
    ADD_GETSET(invalid),
    ADD_GETSET(erange),
    ADD_GETSET(divzero),
    ADD_GETSET(raise_underflow),
    ADD_GETSET(raise_overflow),
    ADD_GETSET(raise_inexact),
    ADD_GETSET(raise_invalid),
    ADD_GETSET(raise_erange),
    ADD_GETSET(raise_divzero),
    {NULL}
};

static PyMethodDef GMPyContext_methods[] =
{
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
