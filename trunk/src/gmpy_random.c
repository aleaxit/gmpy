/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_random.c                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen             *
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

static GMPYRandomState *
GMPYRandomState_New(void)
{
    GMPYRandomState *result;

    if ((result = PyObject_New(GMPYRandomState, &GMPYRandomState_Type))) {
        gmp_randinit_default(result->state);
    }
    return result;
};

static void
GMPYRandomState_Dealloc(GMPYRandomState *self)
{
    gmp_randclear(self->state);
    PyObject_Del(self);
};

static PyObject *
GMPYRandomState_Repr(GMPyContextObject *self)
{
    return Py_BuildValue("s", "<random number state>");
};

PyDoc_STRVAR(doc_random_state,
"random_state([seed]) -> object\n\n"
"Return new object containing state information for the random number\n"
"generator. An optional integer can be specified as the seed value.");

static PyObject *
GMPY_random_state(PyObject *self, PyObject *args)
{
    GMPYRandomState *result;
    PympzObject *temp;
    
    if (!(result = GMPYRandomState_New()))
        return NULL;

    if (PyTuple_GET_SIZE(args) == 0) {
        gmp_randseed_ui(result->state, 0);
    }
    else if (PyTuple_GET_SIZE(args) == 1) {
        if (!(temp = Pympz_From_Integer(other))) {
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("seed must be an integer");
            return NULL;
        }
        gmp_randseed(result->state, temp->z);
        Py_DECREF((PyObject*)temp);
    }
    else {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("random_state() require 0 or 1 integer arguments");
        return NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_urandomm,
"mpz_urandomm(random state, n) -> mpz\n\n"
"Return uniformly distributed integer between 0 and n-1.");

static PyObject *
GMPY_mpz_urandomm(PyObject *self, PyObject *args)
{
    PympzObject *result, *temp;
    
    if (!(temp = Pympz_From_Integer(other))) {
        TYPE_ERROR("seed must be an integer");
        return NULL;
    }
    if (!(result = GMPYRandomState_New())) {
        Py_DECREF((PyObject*)temp);
        return NULL;
    }
    gmp_randseed(result->state, temp->z);
    Py_DECREF((PyObject*)temp);
    return (PyObject*)result;
}

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
