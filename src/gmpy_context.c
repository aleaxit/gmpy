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

PyDoc_STRVAR(doc_context,
"context()\n\n"
"Context manager controlling MPFR and MPC arithmetic.\n\n"
"    nonstop:      if True, return nan or inf\n"
"                  if False, raise exception\n"
"    subnormalize: if True, subnormalized results can be returned\n"
"    mpfr_prec:    precision, in bits, of an MPFR result\n"
"    mpc_rprec:    precision, in bits, of Re(MPC)\n"
"                  -1 implies use mpfr_prec\n"
"    mpc_iprec:    precision, in bits, of Im(MPC)\n"
"                  -1 implies use mpc_rprec\n"
"    mpfr_round:   rounding mode for MPFR\n"
"    mpc_rround:   rounding mdoe for Re(MPC)\n"
"                  -1 implies use mpfr_round\n"
"    mpc_iround:   rounding mode for Im(MPC)\n"
"                  -1 implies use mpc_rround\n"
"    e_max:        maximum allowed exponent\n"
"    e_min:        minimum allowed exponent\n");

static PycontextObject *
Pycontext_new(void)
{
    PycontextObject *self;

    TRACE("Entering Pycontext_new\n");

    if (!(self = PyObject_New(PycontextObject, &Pycontext_Type)))
        return NULL;
    self->orig = context;
    self->now = context;
    return self;
};

static void
Pycontext_dealloc(PycontextObject *self)
{
    PyObject_Del(self);
};

static PyObject *
Pycontext_repr(PycontextObject *self)
{
    PyObject *format;
    PyObject *tuple;
    PyObject *result;

    tuple = PyTuple_New(1);
    if (!tuple) return NULL;

    format = Py2or3String_FromString("context(nonstop=%s)");
    if (!format) return NULL;

    PyTuple_SET_ITEM(tuple, 0, (PyObject*)PyBool_FromLong(self->now.raise));

    if (PyErr_Occurred()) {
        result = NULL;
        goto cleanup;
    };

    result = Py2or3String_Format(format, tuple);

  cleanup:
    Py_XDECREF(PyTuple_GET_ITEM(tuple,0));
    Py_DECREF(tuple);
    return result;
};

static PycontextObject *
Pygmpy_context(PyObject *self, PyObject *args, PyObject *kwargs)
{
    return Pycontext_new();
}

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
};
