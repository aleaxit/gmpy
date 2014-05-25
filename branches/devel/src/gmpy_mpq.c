/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpq.c                                                              *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Case Van Horsen      *
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

PyDoc_STRVAR(GMPy_doc_mpq_factory,
"mpq() -> mpq(0,1)\n\n"
"     If no argument is given, return mpq(0,1).\n\n"
"mpq(n) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n. Decimal and\n"
"     Fraction values are converted exactly.\n\n"
"mpq(n,m) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n/m.\n\n"
"mpq(s[, base=10]) -> mpq\n\n"
"     Return an 'mpq' object from a string s made up of digits in\n"
"     the given base. s may be made up of two numbers in the same\n"
"     base separated by a '/' character.\n");

static PyObject *
GMPy_MPQ_Factory(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPQ_Object *result, *temp;
    PyObject *n, *m;
    int base = 10;
    Py_ssize_t argc, keywdc = 0;
    static char *kwlist[] = {"s", "base", NULL };
    CTXT_Object *context = NULL;

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    argc = PyTuple_Size(args);
    if (keywds) {
        keywdc = PyDict_Size(keywds);
    }

    if (argc + keywdc > 2) {
        TYPE_ERROR("mpq() takes at most 2 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPQ_New(context))) {
            mpq_set_ui(result->q, 0, 1);
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpq() requires at least one non-keyword argument");
        return NULL;
    }

    n = PyTuple_GetItem(args, 0);
    
    /* Handle the case where the first argument is a string. */
    if (PyStrOrUnicode_Check(n)) {
        /* keyword base is legal */
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, &n, &base))) {
                return NULL;
            }
        }

        if ((base != 0) && ((base < 2) || (base > 62))) {
            VALUE_ERROR("base for mpq() must be 0 or in the interval [2, 62]");
            return NULL;
        }

        return (PyObject*)GMPy_MPQ_From_PyStr(n, base, context);
    }

    /* Handle 1 argument. It must be non-complex number. */
    if (argc == 1) {
        if (IS_REAL(n)) {
            return (PyObject*)GMPy_MPQ_From_Number(n, context);
        }
    }

    /* Handle 2 arguments. Both arguments must be integer or rational. */
    if (argc == 2) {
        m = PyTuple_GetItem(args, 1);

        if (IS_RATIONAL(n) && IS_RATIONAL(m)) {
           result = GMPy_MPQ_From_Rational(n, context);
           temp = GMPy_MPQ_From_Rational(m, context);
           if (!result || !temp) {
               Py_XDECREF((PyObject*)result);
               Py_XDECREF((PyObject*)temp);
               return NULL;
            }

            if (mpq_sgn(temp->q) == 0) {
                ZERO_ERROR("zero denominator in mpq()");
                Py_DECREF((PyObject*)result);
                Py_DECREF((PyObject*)temp);
                return NULL;
            }

            mpq_div(result->q, result->q, temp->q);
            Py_DECREF((PyObject*)temp);
            return (PyObject*)result;
        }
    }

    TYPE_ERROR("mpq() requires numeric or string argument");
    return NULL;
}

#ifdef PY3
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) GMPy_MPQ_Add_Slot,         /* nb_add                  */
    (binaryfunc) GMPy_MPQ_Sub_Slot,         /* nb_subtract             */
    (binaryfunc) GMPy_MPQ_Mul_Slot,         /* nb_multiply             */
    (binaryfunc) GMPy_MPQ_Mod_Slot,         /* nb_remainder            */
    (binaryfunc) GMPy_MPQ_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,      /* nb_power                */
    (unaryfunc) GMPy_MPQ_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPQ_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPQ_Abs_Slot,          /* nb_absolute             */
    (inquiry) Pympq_nonzero,                /* nb_bool                 */
        0,                                  /* nb_invert               */
        0,                                  /* nb_lshift               */
        0,                                  /* nb_rshift               */
        0,                                  /* nb_and                  */
        0,                                  /* nb_xor                  */
        0,                                  /* nb_or                   */
    (unaryfunc) GMPy_MPQ_Int_Slot,          /* nb_int                  */
        0,                                  /* nb_reserved             */
    (unaryfunc) GMPy_MPQ_Float_Slot,        /* nb_float                */
        0,                                  /* nb_inplace_add          */
        0,                                  /* nb_inplace_subtract     */
        0,                                  /* nb_inplace_multiply     */
        0,                                  /* nb_inplace_remainder    */
        0,                                  /* nb_inplace_power        */
        0,                                  /* nb_inplace_lshift       */
        0,                                  /* nb_inplace_rshift       */
        0,                                  /* nb_inplace_and          */
        0,                                  /* nb_inplace_xor          */
        0,                                  /* nb_inplace_or           */
    (binaryfunc) GMPy_MPQ_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_MPQ_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                  /* nb_inplace_floor_divide */
        0,                                  /* nb_inplace_true_divide  */
        0,                                  /* nb_index                */
};
#else
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) GMPy_MPQ_Add_Slot,         /* nb_add                  */
    (binaryfunc) GMPy_MPQ_Sub_Slot,         /* nb_subtract             */
    (binaryfunc) GMPy_MPQ_Mul_Slot,         /* nb_multiply             */
    (binaryfunc) GMPy_MPQ_TrueDiv_Slot,     /* nb_divide               */
    (binaryfunc) GMPy_MPQ_Mod_Slot,         /* nb_remainder            */
    (binaryfunc) GMPy_MPQ_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,      /* nb_power                */
    (unaryfunc) GMPy_MPQ_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPQ_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPQ_Abs_Slot,          /* nb_absolute             */
    (inquiry) Pympq_nonzero,                /* nb_bool                 */
        0,                                  /* nb_invert               */
        0,                                  /* nb_lshift               */
        0,                                  /* nb_rshift               */
        0,                                  /* nb_and                  */
        0,                                  /* nb_xor                  */
        0,                                  /* nb_or                   */
        0,                                  /* nb_coerce               */
    (unaryfunc) GMPy_MPQ_Int_Slot,          /* nb_int                  */
    (unaryfunc) GMPy_MPQ_Long_Slot,         /* nb_long                 */
    (unaryfunc) GMPy_MPQ_Float_Slot,        /* nb_float                */
        0,                                  /* nb_oct                  */
        0,                                  /* nb_hex                  */
        0,                                  /* nb_inplace_add;         */
        0,                                  /* nb_inplace_subtract     */
        0,                                  /* nb_inplace_multiply     */
        0,                                  /* nb_inplace_divide       */
        0,                                  /* nb_inplace_remainder    */
        0,                                  /* nb_inplace_power        */
        0,                                  /* nb_inplace_lshift       */
        0,                                  /* nb_inplace_rshift       */
        0,                                  /* nb_inplace_and          */
        0,                                  /* nb_inplace_xor          */
        0,                                  /* nb_inplace_or           */
    (binaryfunc) GMPy_MPQ_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_MPQ_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                  /* nb_inplace_floor_divide */
        0,                                  /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympq_getseters[] =
{
    { "numerator", (getter)Pympq_getnumer, NULL, "numerator", NULL },
    { "denominator", (getter)Pympq_getdenom, NULL, "denominator", NULL },
    {NULL}
};

static PyMethodDef Pympq_methods [] =
{
    { "__ceil__", Pympq_ceil, METH_NOARGS, doc_mpq_ceil },
    { "__floor__", Pympq_floor, METH_NOARGS, doc_mpq_floor },
    { "__round__", Pympq_round, METH_VARARGS, doc_mpq_round },
    { "__sizeof__", Pympq_sizeof, METH_NOARGS, doc_mpq_sizeof },
    { "__trunc__", Pympq_trunc, METH_NOARGS, doc_mpq_trunc },
    { "digits", GMPy_MPQ_Digits_Method, METH_VARARGS, GMPy_doc_mpq_digits_method },
    { NULL, NULL, 1 }
};

static PyTypeObject MPQ_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpq",                                  /* tp_name          */
    sizeof(MPQ_Object),                     /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPQ_Dealloc,          /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPQ_Repr_Slot,          /* tp_repr          */
    &mpq_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) GMPy_MPQ_Hash_Slot,          /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPQ_Str_Slot,           /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE |
        Py_TPFLAGS_CHECKTYPES,              /* tp_flags         */
#endif
    "Multiple precision rational",          /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&GMPy_RichCompare_Slot,    /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympq_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympq_getseters,                        /* tp_getset        */
};

