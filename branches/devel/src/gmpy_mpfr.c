/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpfr.c                                                             *
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

PyDoc_STRVAR(GMPy_doc_mpfr_factory,
"mpfr() -> mpfr(0.0)\n\n"
"      If no argument is given, return mpfr(0.0).\n\n"
"mpfr(n [, precision=0]) -> mpfr\n\n"
"      Return an 'mpfr' object after converting a numeric value. See\n"
"      below for the interpretation of precision.\n\n"
"mpfr(s [, precision=0 [, base=0]]) -> mpfr\n\n"
"      Return a new 'mpfr' object by converting a string s made of\n"
"      digits in the given base, possibly with fraction-part (with a\n"
"      period as a separator) and/or exponent-part (with an exponent\n"
"      marker 'e' for base<=10, else '@'). The base of the string\n"
"      representation must be 0 or in the interval [2,62]. If the base\n"
"      is 0, the leading digits of the string are used to identify the\n"
"      base: 0b implies base=2, 0x implies base=16, otherwise base=10\n"
"      is assumed.\n\n"
"Note: If a precision greater than or equal to 2 is specified, then it\n"
"      is used.\n\n"
"      A precision of 0 (the default) implies the precision of the\n"
"      current context is used.\n\n"
"      A precision of 1 minimizes the loss of precision by following\n"
"      these rules:\n"
"        1) If n is a radix-2 floating point number, then the full\n"
"           precision of n is retained.\n"
"        2) For all other n, the precision of the result is the context\n"
"           precision + guard_bits.\n" );

static PyObject *
GMPy_MPFR_Factory(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPFR_Object *result = NULL;
    PyObject *arg0 = NULL;
    int base = 10;
    Py_ssize_t argc, keywdc = 0;
    CTXT_Object *context = NULL;

    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t prec = 0;

    static char *kwlist_s[] = {"s", "precision", "base", NULL};
    static char *kwlist_n[] = {"n", "precision", NULL};

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

    if (argc + keywdc > 3) {
        TYPE_ERROR("mpfr() takes at most 3 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPFR_New(0, context))) {
            mpfr_set_ui(result->f, 0, MPFR_RNDN);
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpfr() requires at least one non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    /* A string can have both precision and base additional arguments. */
    if (PyStrOrUnicode_Check(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|li", kwlist_s,
                                              &arg0, &prec, &base)))
                return NULL;
        }

        if (prec < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
            return NULL;
        }

        if (base != 0 && (base < 2 || base > 62)) {
            VALUE_ERROR("base for mpfr() must be 0 or in the interval [2, 62]");
            return NULL;
        }

        return (PyObject*)GMPy_MPFR_From_PyStr(arg0, base, prec, context);
    }

    /* A number can only have precision additional argument. */
    if (IS_REAL(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, keywds, "O|l", kwlist_n,
                                              &arg0, &prec)))
                return NULL;
        }

        if (prec < 0) {
            VALUE_ERROR("precision for mpfr() must be >= 0");
            return NULL;
        }

        return (PyObject*)GMPy_MPFR_From_Real(arg0, prec, context);
    }

    TYPE_ERROR("mpfr() requires numeric or string argument");
    return NULL;
}

#ifdef PY3
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) GMPy_MPFR_Add_Slot,         /* nb_add                  */
    (binaryfunc) GMPy_MPFR_Sub_Slot,         /* nb_subtract             */
    (binaryfunc) GMPy_MPFR_Mul_Slot,         /* nb_multiply             */
    (binaryfunc) GMPy_MPFR_Mod_Slot,         /* nb_remainder            */
    (binaryfunc) GMPy_MPFR_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,       /* nb_power                */
    (unaryfunc) GMPy_MPFR_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPFR_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPFR_Abs_Slot,          /* nb_absolute             */
    (inquiry) GMPy_MPFR_NonZero_Slot,        /* nb_bool                 */
        0,                                   /* nb_invert               */
        0,                                   /* nb_lshift               */
        0,                                   /* nb_rshift               */
        0,                                   /* nb_and                  */
        0,                                   /* nb_xor                  */
        0,                                   /* nb_or                   */
    (unaryfunc) GMPy_MPFR_Int_Slot,          /* nb_int                  */
        0,                                   /* nb_reserved             */
    (unaryfunc) GMPy_MPFR_Float_Slot,        /* nb_float                */
        0,                                   /* nb_inplace_add          */
        0,                                   /* nb_inplace_subtract     */
        0,                                   /* nb_inplace_multiply     */
        0,                                   /* nb_inplace_remainder    */
        0,                                   /* nb_inplace_power        */
        0,                                   /* nb_inplace_lshift       */
        0,                                   /* nb_inplace_rshift       */
        0,                                   /* nb_inplace_and          */
        0,                                   /* nb_inplace_xor          */
        0,                                   /* nb_inplace_or           */
    (binaryfunc) GMPy_MPFR_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_MPFR_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                   /* nb_inplace_floor_divide */
        0,                                   /* nb_inplace_true_divide  */
        0,                                   /* nb_index                */
};
#else
static PyNumberMethods mpfr_number_methods =
{
    (binaryfunc) GMPy_MPFR_Add_Slot,         /* nb_add                  */
    (binaryfunc) GMPy_MPFR_Sub_Slot,         /* nb_subtract             */
    (binaryfunc) GMPy_MPFR_Mul_Slot,         /* nb_multiply             */
    (binaryfunc) GMPy_MPFR_TrueDiv_Slot,     /* nb_divide               */
    (binaryfunc) GMPy_MPFR_Mod_Slot,         /* nb_remainder            */
    (binaryfunc) GMPy_MPFR_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,       /* nb_power                */
    (unaryfunc) GMPy_MPFR_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPFR_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPFR_Abs_Slot,          /* nb_absolute             */
    (inquiry) GMPy_MPFR_NonZero_Slot,        /* nb_bool                 */
        0,                                   /* nb_invert               */
        0,                                   /* nb_lshift               */
        0,                                   /* nb_rshift               */
        0,                                   /* nb_and                  */
        0,                                   /* nb_xor                  */
        0,                                   /* nb_or                   */
        0,                                   /* nb_coerce               */
    (unaryfunc) GMPy_MPFR_Int_Slot,          /* nb_int                  */
    (unaryfunc) GMPy_MPFR_Long_Slot,         /* nb_long                 */
    (unaryfunc) GMPy_MPFR_Float_Slot,        /* nb_float                */
        0,                                   /* nb_oct                  */
        0,                                   /* nb_hex                  */
        0,                                   /* nb_inplace_add          */
        0,                                   /* nb_inplace_subtract     */
        0,                                   /* nb_inplace_multiply     */
        0,                                   /* nb_inplace_divide       */
        0,                                   /* nb_inplace_remainder    */
        0,                                   /* nb_inplace_power        */
        0,                                   /* nb_inplace_lshift       */
        0,                                   /* nb_inplace_rshift       */
        0,                                   /* nb_inplace_and          */
        0,                                   /* nb_inplace_xor          */
        0,                                   /* nb_inplace_or           */
    (binaryfunc) GMPy_MPFR_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_MPFR_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                   /* nb_inplace_floor_divide */
        0,                                   /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympfr_getseters[] =
{
    {"precision", (getter)GMPy_MPFR_GetPrec_Attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)GMPy_MPFR_GetRc_Attrib, NULL, "return code", NULL},
    {"imag", (getter)GMPy_MPFR_GetImag_Attrib, NULL, "imaginary component", NULL},
    {"real", (getter)GMPy_MPFR_GetReal_Attrib, NULL, "real component", NULL},
    {NULL}
};

static PyMethodDef Pympfr_methods [] =
{
    { "__ceil__", GMPy_MPFR_Method_Ceil, METH_NOARGS, GMPy_doc_mpfr_ceil_method },
    { "__floor__", GMPy_MPFR_Method_Floor, METH_NOARGS, GMPy_doc_mpfr_floor_method },
    { "__format__", GMPy_MPFR_Format, METH_VARARGS, GMPy_doc_mpfr_format },
    { "__round__", GMPy_MPFR_Method_Round10, METH_VARARGS, GMPy_doc_method_round10 },
    { "__sizeof__", GMPy_MPFR_SizeOf_Method, METH_NOARGS, GMPy_doc_mpfr_sizeof_method },
    { "__trunc__", GMPy_MPFR_Method_Trunc, METH_NOARGS, GMPy_doc_mpfr_trunc_method },
    { "as_integer_ratio", GMPy_MPFR_Integer_Ratio_Method, METH_NOARGS, GMPy_doc_method_integer_ratio },
    { "as_mantissa_exp", GMPy_MPFR_Mantissa_Exp_Method, METH_NOARGS, GMPy_doc_method_mantissa_exp },
    { "as_simple_fraction", (PyCFunction)GMPy_MPFR_Simple_Fraction_Method, METH_VARARGS | METH_KEYWORDS, GMPy_doc_method_simple_fraction },
    { "conjugate", GMPy_MPFR_Conjugate_Method, METH_NOARGS, GMPy_doc_mpfr_conjugate_method },
    { "digits", GMPy_MPFR_Digits_Method, METH_VARARGS, GMPy_doc_mpfr_digits_method },
    { "is_finite", GMPy_MPFR_Is_Finite_Method, METH_NOARGS, GMPy_doc_method_is_finite },
    { "is_infinite", GMPy_MPFR_Is_Infinite_Method, METH_NOARGS, GMPy_doc_method_is_infinite },
    { "is_integer", GMPy_MPFR_Is_Integer_Method, METH_NOARGS, GMPy_doc_method_is_integer },
    { "is_nan", GMPy_MPFR_Is_NAN_Method, METH_NOARGS, GMPy_doc_method_is_nan },
    { "is_signed", GMPy_MPFR_Is_Regular_Method, METH_NOARGS, GMPy_doc_method_is_regular },
    { "is_signed", GMPy_MPFR_Is_Signed_Method, METH_NOARGS, GMPy_doc_method_is_signed },
    { "is_zero", GMPy_MPFR_Is_Zero_Method, METH_NOARGS, GMPy_doc_method_is_zero },
    { NULL, NULL, 1 }
};

static PyTypeObject MPFR_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpfr",                                 /* tp_name          */
    sizeof(MPFR_Object),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPFR_Dealloc,         /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPFR_Repr_Slot,         /* tp_repr          */
    &mpfr_number_methods,                   /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) GMPy_MPFR_Hash_Slot,         /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPFR_Str_Slot,          /* tp_str           */
    (getattrofunc) 0,                       /* tp_getattro      */
    (setattrofunc) 0,                       /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_RICHCOMPARE|Py_TPFLAGS_CHECKTYPES,  /* tp_flags */
#endif
    "Multiple precision real",              /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&GMPy_RichCompare_Slot,    /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympfr_methods,                         /* tp_methods       */
        0,                                  /* tp_members       */
    Pympfr_getseters,                       /* tp_getset        */
};

