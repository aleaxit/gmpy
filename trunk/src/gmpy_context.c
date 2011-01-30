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

PyDoc_STRVAR(doc_current,
"current()\n\n"
"Access attributes in the current context manager controlling\n"
"MPFR and MPC arithmetic.\n\n"
"    nonstop:      if True, return nan or inf\n"
"                  if False, raise exception\n"
"    subnormalize: if True, subnormalized results can be returned\n"
"    mpfr_prec:    precision, in bits, of an MPFR result\n"
"    mpc_rprec:    precision, in bits, of Re(MPC)\n"
"                    -1 implies use mpfr_prec\n"
"    mpc_iprec:    precision, in bits, of Im(MPC)\n"
"                    -1 implies use mpc_rprec\n"
"    mpfr_round:   rounding mode for MPFR\n"
"    mpc_rround:   rounding mdoe for Re(MPC)\n"
"                    -1 implies use mpfr_round\n"
"    mpc_iround:   rounding mode for Im(MPC)\n"
"                    -1 implies use mpc_rround\n"
"    e_max:        maximum allowed exponent\n"
"    e_min:        minimum allowed exponent\n");

#define GET_MPC_RPREC(c) ((c->now.mpc_rprec==GMPY_DEFAULT)?c->now.mpfr_prec:c->now.mpc_rprec)
#define GET_MPC_IPREC(c) ((c->now.mpc_iprec==GMPY_DEFAULT)?GET_MPC_RPREC(c):c->now.mpc_iprec)

#define GET_MPC_RROUND(c) ((c->now.mpc_rround==GMPY_DEFAULT)?c->now.mpfr_round:c->now.mpc_rround)
#define GET_MPC_IROUND(c) ((c->now.mpc_iround==GMPY_DEFAULT)?GET_MPC_RROUND(c):c->now.mpc_iround)
#define GET_MPC_ROUND(c) (RNDC(GET_MPC_RROUND(c), GET_MPC_IROUND(c)))

static PycontextObject *
Pycontext_new(void)
{
    PycontextObject *self;

    TRACE("Entering Pycontext_new\n");

    if (!(self = PyObject_New(PycontextObject, &Pycontext_Type)))
        return NULL;
    if (!context) {
        self->orig.nonstop = 0;
        self->orig.subnormalize = 0;
        self->orig.mpfr_prec = DBL_MANT_DIG;
        self->orig.mpc_rprec = -1;
        self->orig.mpc_iprec = -1;
        self->orig.mpfr_round = MPFR_RNDN;
        self->orig.mpc_rround = -1;
        self->orig.mpc_iround = -1;
        self->orig.mpc_round = MPC_RNDNN;
        self->orig.emax = mpfr_get_emax();
        self->orig.emin = mpfr_get_emin();
        self->now = self->orig;
    }
    else {
        self->orig = context->now;
        self->now = context->now;
    }
    return self;
};

static void
Pycontext_dealloc(PycontextObject *self)
{
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
Pycontext_repr(PycontextObject *self)
{
    PyObject *format;
    PyObject *tuple;
    PyObject *result = NULL;

    tuple = PyTuple_New(10);
    if (!tuple) return NULL;

    format = Py2or3String_FromString(
            "context(nonstop=%s,\n"
            "        subnormalize=%s,\n"
            "        precision=%s,\n"
            "        mpc_rprec=%s,\n"
            "        mpc_iprec=%s,\n"
            "        round=%s,\n"
            "        mpc_rround=%s,\n"
            "        mpc_iround=%s,\n"
            "        emax=%s,\n"
            "        emin=%s)"
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

    if (!PyErr_Occurred())
        result = Py2or3String_Format(format, tuple);
    else
        SYSTEM_ERROR("internal error in Pycontext_repr");

    Py_DECREF(format);
    Py_DECREF(tuple);
    return result;
};

/* Return a reference to the current context. */

static PycontextObject *
Pygmpy_current(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PycontextObject *result;

    result = context;
    Py_INCREF((PyObject*)context);
    return result;
}

/* Define the get/set functions. */

static PyObject *
Pycontext_get_nonstop(PycontextObject *self, void *closure)
{
    return PyBool_FromLong(self->now.nonstop);
}

static int
Pycontext_set_nonstop(PycontextObject *self, PyObject *value, void *closure)
{
    if (!(PyBool_Check(value))) {
        TYPE_ERROR("nonstop must be True or False");
        return -1;
    }
    if (value == Py_True)
        self->now.nonstop = 1;
    else
        self->now.nonstop = 0;
    return 0;
}

static PyObject *
Pycontext_get_subnormalize(PycontextObject *self, void *closure)
{
    return PyBool_FromLong(self->now.nonstop);
}

static int
Pycontext_set_subnormalize(PycontextObject *self, PyObject *value, void *closure)
{
    if (!(PyBool_Check(value))) {
        TYPE_ERROR("subnormalize must be True or False");
        return -1;
    }
    if (value == Py_True)
        self->now.subnormalize = 1;
    else
        self->now.subnormalize = 0;
    return 0;
}

static PyObject *
Pycontext_get_mpfr_prec(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.mpfr_prec));
}

static int
Pycontext_set_mpfr_prec(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_mpc_rprec(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_MPC_RPREC(self)));
}

static int
Pycontext_set_mpc_rprec(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_mpc_iprec(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(GET_MPC_IPREC(self)));
}

static int
Pycontext_set_mpc_iprec(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_mpfr_round(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)(self->now.mpfr_round));
}

static int
Pycontext_set_mpfr_round(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_mpc_rround(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_MPC_RROUND(self));
}

static int
Pycontext_set_mpc_rround(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_mpc_iround(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromLong((long)GET_MPC_IROUND(self));
}

static int
Pycontext_set_mpc_iround(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_emin(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.emin));
}

static int
Pycontext_set_emin(PycontextObject *self, PyObject *value, void *closure)
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
Pycontext_get_emax(PycontextObject *self, void *closure)
{
    return PyIntOrLong_FromSsize_t((Py_ssize_t)(self->now.emax));
}

static int
Pycontext_set_emax(PycontextObject *self, PyObject *value, void *closure)
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

static PyGetSetDef Pycontext_getseters[] = {
    {"nonstop",
        (getter)Pycontext_get_nonstop,
        (setter)Pycontext_set_nonstop, NULL, NULL},
    {"subnormalize",
        (getter)Pycontext_get_subnormalize,
        (setter)Pycontext_set_subnormalize, NULL, NULL},
    {"precision",
        (getter)Pycontext_get_mpfr_prec,
        (setter)Pycontext_set_mpfr_prec, NULL, NULL},
    {"mpc_rprec",
        (getter)Pycontext_get_mpc_rprec,
        (setter)Pycontext_set_mpc_rprec, NULL, NULL},
    {"mpc_iprec",
        (getter)Pycontext_get_mpc_iprec,
        (setter)Pycontext_set_mpc_iprec, NULL, NULL},
    {"round",
        (getter)Pycontext_get_mpfr_round,
        (setter)Pycontext_set_mpfr_round, NULL, NULL},
    {"mpc_rround",
        (getter)Pycontext_get_mpc_rround,
        (setter)Pycontext_set_mpc_rround, NULL, NULL},
    {"mpc_iround",
        (getter)Pycontext_get_mpc_iround,
        (setter)Pycontext_set_mpc_iround, NULL, NULL},
    {"emax",
        (getter)Pycontext_get_emax,
        (setter)Pycontext_set_emax, NULL, NULL},
    {"emin",
        (getter)Pycontext_get_emin,
        (setter)Pycontext_set_emin, NULL, NULL},
    {NULL}
};

static PyTypeObject Pycontext_Type =
{
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "context",                              /* tp_name          */
    sizeof(PycontextObject),                /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    (destructor) Pycontext_dealloc,         /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) Pycontext_repr,              /* tp_repr          */
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
        0,                                  /* tp_methods       */
        0,                                  /* tp_members       */
    Pycontext_getseters,                    /* tp_getset        */
};
