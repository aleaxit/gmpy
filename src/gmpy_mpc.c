/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.c                                                              *
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

/* Verify that a valid rounding mode is specified for complex arithmetic.
 * Returns 0 (false) if the rounding mode is not valid else returns 1 (true).
 */

static int
Pymisc_verify_mpc_round(int rmode)
{
    if ( rmode == MPC_RNDNN || rmode == MPC_RNDNZ ||
         rmode == MPC_RNDNU || rmode == MPC_RNDND ||
         rmode == MPC_RNDZN || rmode == MPC_RNDZZ ||
         rmode == MPC_RNDZU || rmode == MPC_RNDZD ||
         rmode == MPC_RNDUN || rmode == MPC_RNDUZ ||
         rmode == MPC_RNDUU || rmode == MPC_RNDUD ||
         rmode == MPC_RNDDN || rmode == MPC_RNDDZ ||
         rmode == MPC_RNDDU || rmode == MPC_RNDDD )
        return 1;
    else
        return 0;
}

/* Verify that valid precisions are requested for complex arithmetic.
 * Returns 0 if the precisions are not valid else returns 1.
 */

static int
Pymisc_verify_mpc_precision(Py_ssize_t rprec, Py_ssize_t iprec)
{
    if ( rprec < MPFR_PREC_MIN || rprec > MPFR_PREC_MAX ||
         iprec < MPFR_PREC_MIN || iprec > MPFR_PREC_MAX )
        return 0;
    else
        return 1;
}




static PympcObject *
Pympc2Pympc(PyObject *self, mpfr_prec_t rprec, mpfr_prec_t iprec)
{
    PympcObject *newob;

    assert(Pympc_Check(self));
    if (rprec == 0 || iprec == 0)
        mpc_get_prec2(&rprec, &iprec, Pympc_AS_MPC(self));
    if ((newob = Pympc_new(rprec, iprec)))
        mpc_set(newob->c, Pympc_AS_MPC(self), context->now.mpc_round);
    return newob;
}

PyDoc_STRVAR(doc_g_mpc,
"mpc(n, [prec=(i,i) | rnd=(i,i]) -> mpc object\n\n"
"Return an mpc object by converting a numeric value 'n' into a\n"
"complex number. If 'prec' is omitted, then get_mpc_precision() is\n"
"used. If 'rnd' is omitted, then get_mpc_round() is used.\n\n"
"mpc(s, [base = 10 | prec=(i,i) | rnd=(i,i)]) -> mpc object\n\n"
"Return an mpc object by converting a string 's' into a complex\n"
"number. If 'base' is omitted, then a base 10 representation is\n"
"assumed otherwise a base between 2 and 36 can be specified.\n"
"If 'prec' is omitted, then get_mpc_precision() is used. If 'rnd'\n"
"is omitted, then get_mpc_round() is used.");
static PyObject *
Pygmpy_mpc(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PympcObject *result;
    PyObject *arg0 = NULL, *prec_obj = NULL;
    Py_ssize_t dummy, rprec, iprec;
    int base, rmode;
    static char *kwlist[] = {"n", "base", "prec", "rnd", NULL};

    dummy = 0;
    base = 10;
    rprec = context->now.mpc_rprec;
    iprec = context->now.mpc_iprec;
    rmode = context->now.mpc_round;

    if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "|OiOi", kwlist,
                                      &arg0, &base, &prec_obj, &rmode)))
        return NULL;

    if (base < 2 || base > 36) {
        VALUE_ERROR("base for gmpy2.mpc() must be in the interval 2..36.");
        return NULL;
    }

    if (rmode != context->now.mpc_round && Pymisc_verify_mpc_round(rmode) == -1) {
        VALUE_ERROR("invalid rounding mode for complex arithmetic.");
        return NULL;
    }

    if (prec_obj && PyIntOrLong_Check(prec_obj)) {
        rprec = PyLong_AsSsize_t(prec_obj);
        iprec = rprec;

        /* Also catches return value of -1. */
        if (!Pymisc_verify_mpc_precision(rprec, iprec)) {
            VALUE_ERROR("invalid precision for complex arithmetic");
            return NULL;
        }
    }
    else if (prec_obj && PyTuple_Check(prec_obj)) {
        if (!(PyArg_ParseTuple(prec_obj,
                    "nn;invalid precision for complex arithmetic.",
                    &rprec, &iprec)))
            return NULL;

        if (!Pymisc_verify_mpc_precision(rprec, iprec)) {
            VALUE_ERROR("invalid precision for complex arithmetic");
            return NULL;
        }
    }
    else {
        VALUE_ERROR("invalid precision for complex arithmetic");
        return NULL;
    }

    if (arg0 == NULL) {
        if (!(result = Pympc_new(rprec, iprec)))
            return NULL;
        /* May want to return (NaN NaN) ?? */
        mpc_set_ui(result->c, 0, rmode);
        return (PyObject*)result;
    }

    Py_INCREF(arg0);
    return Py_BuildValue("OiOi", arg0, base,
                         Py_BuildValue("nn", rprec, iprec),
                         rmode);
}

#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
    0,            /* nb_add                  */
    0,            /* nb_subtract             */
    0,            /* nb_multiply             */
    0,            /* nb_remaider             */
    0,         /* nb_divmod               */
    0,           /* nb_power                */
    0,              /* nb_negative             */
    0,              /* nb_positive             */
    0,              /* nb_absolute             */
    0,            /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    0,           /* nb_int                  */
        0,                               /* nb_reserved             */
    0,          /* nb_float                */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    0,       /* nb_floor_divide         */
    0,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
    0,            /* nb_add                  */
    0,            /* nb_subtract             */
    0,            /* nb_multiply             */
    0,           /* nb_divide               */
    0,            /* nb_remaider             */
    0,         /* nb_divmod               */
    0,           /* nb_power                */
    0,              /* nb_negative             */
    0,               /* nb_positive             */
    0,              /* nb_absolute             */
    0,             /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    0,            /* nb_int                  */
    0,           /* nb_long                 */
    0,          /* nb_float                */
        0,                               /* nb_oct                  */
        0,                               /* nb_hex                  */
        0,                               /* nb_inplace_add          */
        0,                               /* nb_inplace_subtract     */
        0,                               /* nb_inplace_multiply     */
        0,                               /* nb_inplace_divide       */
        0,                               /* nb_inplace_remainder    */
        0,                               /* nb_inplace_power        */
        0,                               /* nb_inplace_lshift       */
        0,                               /* nb_inplace_rshift       */
        0,                               /* nb_inplace_and          */
        0,                               /* nb_inplace_xor          */
        0,                               /* nb_inplace_or           */
    0,       /* nb_floor_divide         */
    0,        /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyMethodDef Pympc_methods[] =
{
    { NULL, NULL, 1 }
};

static PyTypeObject Pympc_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpc",                                  /* tp_name          */
    sizeof(PympcObject),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) Pympc_dealloc,             /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
        0,                                  /* tp_repr          */
    &mpc_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) 0,                           /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "MPC-based complex number",             /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
        0,                                  /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympc_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
        0,                                  /* tp_getset        */
};


