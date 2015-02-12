/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_random.c                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015 Case Van Horsen                                          *
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

static GMPYRandomStateObject *
GMPYRandomState_New(void)
{
    GMPYRandomStateObject *result;

    if ((result = PyObject_New(GMPYRandomStateObject,
                               &GMPYRandomState_Type))) {
        gmp_randinit_default(result->state);
    }
    return result;
};

static void
GMPYRandomState_Dealloc(GMPYRandomStateObject *self)
{
    gmp_randclear(self->state);
    PyObject_Del(self);
};

static PyObject *
GMPYRandomState_Repr(GMPYRandomStateObject *self)
{
    return Py_BuildValue("s", "<gmpy2.RandomState>");
};

PyDoc_STRVAR(doc_random_state,
"random_state([seed]) -> object\n\n"
"Return new object containing state information for the random number\n"
"generator. An optional integer can be specified as the seed value.");

static PyObject *
GMPY_random_state(PyObject *self, PyObject *args)
{
    GMPYRandomStateObject *result;
    PympzObject *temp;

    if (!(result = GMPYRandomState_New()))
        return NULL;

    if (PyTuple_GET_SIZE(args) == 0) {
        gmp_randseed_ui(result->state, 0);
    }
    else if (PyTuple_GET_SIZE(args) == 1) {
        if (!(temp = Pympz_From_Integer(PyTuple_GET_ITEM(args,0)))) {
            Py_DECREF((PyObject*)result);
            TYPE_ERROR("seed must be an integer");
            return NULL;
        }
        gmp_randseed(result->state, temp->z);
        Py_DECREF((PyObject*)temp);
    }
    else {
        Py_DECREF((PyObject*)result);
        TYPE_ERROR("random_state() requires 0 or 1 integer arguments");
        return NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_urandomb,
"mpz_urandomb(random_state, bit_count) -> mpz\n\n"
"Return uniformly distributed random integer between 0 and\n"
"2**bit_count-1.");

static PyObject *
GMPY_mpz_urandomb(PyObject *self, PyObject *args)
{
    PympzObject *result;
    mp_bitcnt_t len;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mpz_urandomb() requires 2 arguments");
        return NULL;
    }

    if (!GMPYRandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpz_urandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    len = MP_BITCNT_FROM_INTEGER(PyTuple_GET_ITEM(args, 1));
    if (len == (mp_bitcnt_t)-1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz_urandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    if ((result = (PympzObject*)Pympz_new())) {
        mpz_urandomb(Pympz_AS_MPZ(result),
                     PyObj_AS_STATE(PyTuple_GET_ITEM(args, 0)),
                     len);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_rrandomb,
"mpz_rrandomb(random_state, bit_count) -> mpz\n\n"
"Return a random integer between 0 and 2**bit_count-1 with long\n"
"sequences of zeros and one in its binary representation.");

static PyObject *
GMPY_mpz_rrandomb(PyObject *self, PyObject *args)
{
    PympzObject *result;
    mp_bitcnt_t len;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mpz_rrandomb() requires 2 arguments");
        return NULL;
    }

    if (!GMPYRandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpz_rrandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    len = MP_BITCNT_FROM_INTEGER(PyTuple_GET_ITEM(args, 1));
    if (len == (mp_bitcnt_t)-1 && PyErr_Occurred()) {
        TYPE_ERROR("mpz_rrandomb() requires 'random_state' and 'bit_count' arguments");
        return NULL;
    }

    if ((result = (PympzObject*)Pympz_new())) {
        mpz_rrandomb(Pympz_AS_MPZ(result),
                     PyObj_AS_STATE(PyTuple_GET_ITEM(args, 0)),
                     len);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpz_random,
"mpz_random(random_state, int) -> mpz\n\n"
"Return uniformly distributed random integer between 0 and n-1.");

static PyObject *
GMPY_mpz_random(PyObject *self, PyObject *args)
{
    PympzObject *result, *temp;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("mpz_random() requires 2 arguments");
        return NULL;
    }

    if (!GMPYRandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpz_random() requires 'random_state' and 'int' arguments");
        return NULL;
    }

    if (!(temp = Pympz_From_Integer(PyTuple_GET_ITEM(args, 1)))) {
        TYPE_ERROR("mpz_random() requires 'random_state' and 'int' arguments");
        return NULL;
    }

    if ((result = (PympzObject*)Pympz_new())) {
        mpz_urandomm(Pympz_AS_MPZ(result),
                     PyObj_AS_STATE(PyTuple_GET_ITEM(args, 0)),
                     Pympz_AS_MPZ(temp));
    }

    Py_DECREF((PyObject*)temp);
    return (PyObject*)result;
}

#ifdef WITHMPFR
PyDoc_STRVAR(doc_mpfr_random,
"mpfr_random(random_state) -> mpfr\n\n"
"Return uniformly distributed number between [0,1].");

static PyObject *
GMPY_mpfr_random(PyObject *self, PyObject *args)
{
    PympfrObject *result;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfr_random() requires 1 argument");
        return NULL;
    }

    if (!GMPYRandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpfr_random() requires 'random_state' argument");
        return NULL;
    }

    if ((result = (PympfrObject*)Pympfr_new(0))) {
        mpfr_urandom(Pympfr_AS_MPFR(result),
                     PyObj_AS_STATE(PyTuple_GET_ITEM(args, 0)),
                     context->ctx.mpfr_round);
    }

    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpfr_grandom,
"mpfr_grandom(random_state) -> (mpfr, mpfr)\n\n"
"Return two random numbers with gaussian distribution.");

static PyObject *
GMPY_mpfr_grandom(PyObject *self, PyObject *args)
{
    PympfrObject *result1, *result2;
    PyObject *result;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfr_grandom() requires 1 argument");
        return NULL;
    }

    if (!GMPYRandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpfr_grandom() requires 'random_state' argument");
        return NULL;
    }

    result1 = (PympfrObject*)Pympfr_new(0);
    result2 = (PympfrObject*)Pympfr_new(0);
    if (!result1 || !result2) {
        Py_XDECREF((PyObject*)result1);
        Py_XDECREF((PyObject*)result2);
        return NULL;
    }

    mpfr_grandom(Pympfr_AS_MPFR(result1), Pympfr_AS_MPFR(result2),
                 PyObj_AS_STATE(PyTuple_GET_ITEM(args, 0)),
                 context->ctx.mpfr_round);

    result = Py_BuildValue("(NN)", (PyObject*)result1, (PyObject*)result2);
    if (!result) {
        Py_DECREF((PyObject*)result1);
        Py_DECREF((PyObject*)result2);
    }
    return result;
}
#endif

#ifdef WITHMPC
PyDoc_STRVAR(doc_mpc_random,
"mpfc_random(random_state) -> mpc\n\n"
"Return uniformly distributed number in the unit square [0,1]x[0,1].");

static PyObject *
GMPY_mpc_random(PyObject *self, PyObject *args)
{
    PympcObject *result;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("mpfc_random() requires 1 argument");
        return NULL;
    }

    if (!GMPYRandomState_Check(PyTuple_GET_ITEM(args, 0))) {
        TYPE_ERROR("mpc_random() requires 'random_state' argument");
        return NULL;
    }

    if ((result = (PympcObject*)Pympc_new(0,0))) {
        mpc_urandom(Pympc_AS_MPC(result),
                     PyObj_AS_STATE(PyTuple_GET_ITEM(args, 0)));
    }

    return (PyObject*)result;
}
#endif

static PyTypeObject GMPYRandomState_Type =
{
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "gmpy2 random status",                  /* tp_name          */
    sizeof(GMPYRandomStateObject),          /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    (destructor) GMPYRandomState_Dealloc,   /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPYRandomState_Repr,        /* tp_repr          */
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
    "GMPY2 Random number generator state",  /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
        0,                                  /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
        0,                                  /* tp_methods       */
        0,                                  /* tp_members       */
        0,                                  /* tp_getset        */
};
