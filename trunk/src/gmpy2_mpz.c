/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpz.c                                                             *
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

PyDoc_STRVAR(GMPy_doc_mpz_factory,
"mpz() -> mpz(0)\n\n"
"     If no argument is given, return mpz(0).\n\n"
"mpz(n) -> mpz\n\n"
"     Return an 'mpz' object with a numeric value 'n' (truncating n\n"
"     to its integer part if it's a Fraction, 'mpq', Decimal, float\n"
"     or 'mpfr').\n\n"
"mpz(s[, base=0]):\n\n"
"     Return an 'mpz' object from a string 's' made of digits in the\n"
"     given base.  If base=0, binary, octal, or hex Python strings\n"
"     are recognized by leading 0b, 0o, or 0x characters, otherwise\n"
"     the string is assumed to be decimal. Values for base can range\n"
"     between 2 and 62.");

static PyObject *
GMPy_MPZ_Factory(PyObject *self, PyObject *args, PyObject *keywds)
{
    MPZ_Object *result = NULL;
    PyObject *n;
    int base = 0;
    Py_ssize_t argc;
    static char *kwlist[] = {"s", "base", NULL };
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context)
    
    /* Optimize the most common use cases first; either 0 or 1 argument */
    
    argc = PyTuple_GET_SIZE(args);
    
    if (argc == 0) {
        if ((result = GMPy_MPZ_New(context))) {
            mpz_set_ui(result->z, 0);
        }
        return (PyObject*)result;
    }

    if (argc == 1 && !keywds) {
        n = PyTuple_GET_ITEM(args, 0);
        if (IS_REAL(n)) {
            result = GMPy_MPZ_From_Number(n, context);
        }
        else if (PyStrOrUnicode_Check(n)) {
            result = GMPy_MPZ_From_PyStr(n, base, context);
        }
        else {
            TYPE_ERROR("mpz() requires numeric or string argument");
        }
        return (PyObject*)result;
    }

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|i", kwlist, &n, &base)) {
        return NULL;
    }

    if ((base != 0) && ((base < 2)|| (base > 62))) {
        VALUE_ERROR("base for mpz() must be 0 or in the interval [2, 62]");
        return NULL;
    }

    if (PyStrOrUnicode_Check(n)) {
        result = GMPy_MPZ_From_PyStr(n, base, context);
    }
    else if (IS_REAL(n)) {
        TYPE_ERROR("mpz() with number argument only takes 1 argument");
    }
    else {
        TYPE_ERROR("mpz() requires numeric or string (and optional base) arguments");
    }
    return (PyObject*)result;
}

#ifdef PY3
static PyNumberMethods GMPy_MPZ_number_methods =
{
    (binaryfunc) GMPy_MPZ_Add_Slot,        /* nb_add                  */
    (binaryfunc) GMPy_MPZ_Sub_Slot,        /* nb_subtract             */
    (binaryfunc) GMPy_MPZ_Mul_Slot,        /* nb_multiply             */
    (binaryfunc) GMPy_MPZ_Mod_Slot,        /* nb_remainder            */
    (binaryfunc) GMPy_MPZ_DivMod_Slot,     /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,     /* nb_power                */
    (unaryfunc) GMPy_MPZ_Minus_Slot,       /* nb_negative             */
    (unaryfunc) GMPy_MPZ_Plus_Slot,        /* nb_positive             */
    (unaryfunc) GMPy_MPZ_Abs_Slot,         /* nb_absolute             */
    (inquiry) GMPy_MPZ_NonZero_Slot,       /* nb_bool                 */
    (unaryfunc) GMPy_MPZ_Invert_Slot,      /* nb_invert               */
    (binaryfunc) GMPy_MPZ_Lshift_Slot,     /* nb_lshift               */
    (binaryfunc) GMPy_MPZ_Rshift_Slot,     /* nb_rshift               */
    (binaryfunc) GMPy_MPZ_And_Slot,        /* nb_and                  */
    (binaryfunc) GMPy_MPZ_Xor_Slot,        /* nb_xor                  */
    (binaryfunc) GMPy_MPZ_Ior_Slot,        /* nb_or                   */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_int                  */
        0,                                 /* nb_reserved             */
    (unaryfunc) GMPy_MPZ_Float_Slot,       /* nb_float                */
    (binaryfunc) GMPy_MPZ_IAdd_Slot,       /* nb_inplace_add          */
    (binaryfunc) GMPy_MPZ_ISub_Slot,       /* nb_inplace_subtract     */
    (binaryfunc) GMPy_MPZ_IMul_Slot,       /* nb_inplace_multiply     */
    (binaryfunc) GMPy_MPZ_IRem_Slot,       /* nb_inplace_remainder    */
    (ternaryfunc) GMPy_MPZ_IPow_Slot,      /* nb_inplace_power        */
    (binaryfunc) GMPy_MPZ_ILshift_Slot,    /* nb_inplace_lshift       */
    (binaryfunc) GMPy_MPZ_IRshift_Slot,    /* nb_inplace_rshift       */
        0,                                 /* nb_inplace_and          */
        0,                                 /* nb_inplace_xor          */
        0,                                 /* nb_inplace_or           */
    (binaryfunc) GMPy_MPZ_FloorDiv_Slot,   /* nb_floor_divide         */
    (binaryfunc) GMPy_MPZ_TrueDiv_Slot,    /* nb_true_divide          */
    (binaryfunc) GMPy_MPZ_IFloorDiv_Slot,  /* nb_inplace_floor_divide */
        0,                                 /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_index                */
};

#else
static PyNumberMethods GMPy_MPZ_number_methods =
{
    (binaryfunc) GMPy_MPZ_Add_Slot,        /* nb_add                  */
    (binaryfunc) GMPy_MPZ_Sub_Slot,        /* nb_subtract             */
    (binaryfunc) GMPy_MPZ_Mul_Slot,        /* nb_multiply             */
    (binaryfunc) GMPy_MPZ_Div2_Slot,       /* nb_divide               */
    (binaryfunc) GMPy_MPZ_Mod_Slot,        /* nb_remainder            */
    (binaryfunc) GMPy_MPZ_DivMod_Slot,     /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,     /* nb_power                */
    (unaryfunc) GMPy_MPZ_Minus_Slot,       /* nb_negative             */
    (unaryfunc) GMPy_MPZ_Plus_Slot,        /* nb_positive             */
    (unaryfunc) GMPy_MPZ_Abs_Slot,         /* nb_absolute             */
    (inquiry) GMPy_MPZ_NonZero_Slot,       /* nb_bool                 */
    (unaryfunc) GMPy_MPZ_Invert_Slot,      /* nb_invert               */
    (binaryfunc) GMPy_MPZ_Lshift_Slot,     /* nb_lshift               */
    (binaryfunc) GMPy_MPZ_Rshift_Slot,     /* nb_rshift               */
    (binaryfunc) GMPy_MPZ_And_Slot,        /* nb_and                  */
    (binaryfunc) GMPy_MPZ_Xor_Slot,        /* nb_xor                  */
    (binaryfunc) GMPy_MPZ_Ior_Slot,        /* nb_or                   */
        0,                                 /* nb_coerce               */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_int                  */
    (unaryfunc) GMPy_MPZ_Long_Slot,        /* nb_long                 */
    (unaryfunc) GMPy_MPZ_Float_Slot,       /* nb_float                */
    (unaryfunc) GMPy_MPZ_Oct_Slot,         /* nb_oct                  */
    (unaryfunc) GMPy_MPZ_Hex_Slot,         /* nb_hex                  */
    (binaryfunc) GMPy_MPZ_IAdd_Slot,       /* nb_inplace_add          */
    (binaryfunc) GMPy_MPZ_ISub_Slot,       /* nb_inplace_subtract     */
    (binaryfunc) GMPy_MPZ_IMul_Slot,       /* nb_inplace_multiply     */
        0,                                 /* nb_inplace_divide       */
    (binaryfunc) GMPy_MPZ_IRem_Slot,       /* nb_inplace_remainder    */
    (ternaryfunc) GMPy_MPZ_IPow_Slot,      /* nb_inplace_power        */
    (binaryfunc) GMPy_MPZ_ILshift_Slot,    /* nb_inplace_lshift       */
    (binaryfunc) GMPy_MPZ_IRshift_Slot,    /* nb_inplace_rshift       */
        0,                                 /* nb_inplace_and          */
        0,                                 /* nb_inplace_xor          */
        0,                                 /* nb_inplace_or           */
    (binaryfunc) GMPy_MPZ_FloorDiv_Slot,   /* nb_floor_divide         */
    (binaryfunc) GMPy_MPZ_TrueDiv_Slot,    /* nb_true_divide          */
    (binaryfunc) GMPy_MPZ_IFloorDiv_Slot,  /* nb_inplace_floor_divide */
        0,                                 /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_MPZ_Int_Slot,         /* nb_index                */
};
#endif

static PyMappingMethods GMPy_MPZ_mapping_methods = {
    (lenfunc)GMPy_MPZ_Method_Length,
    (binaryfunc)GMPy_MPZ_Method_SubScript,
    NULL
};

static PyGetSetDef GMPy_MPZ_getseters[] =
{
    { "numerator", (getter)GMPy_MPZ_Attrib_GetNumer, NULL, "numerator", NULL },
    { "denominator", (getter)GMPy_MPZ_Attrib_GetDenom, NULL, "denominator", NULL },
    {NULL}
};

static PyMethodDef GMPy_MPZ_methods [] =
{
    { "__format__", GMPy_MPZ_Format, METH_VARARGS, GMPy_doc_mpz_format },
    { "__ceil__", GMPy_MPZ_Method_Ceil, METH_NOARGS, GMPy_doc_mpz_method_ceil },
    { "__floor__", GMPy_MPZ_Method_Floor, METH_NOARGS, GMPy_doc_mpz_method_floor },
    { "__round__", GMPy_MPZ_Method_Round, METH_VARARGS, GMPy_doc_mpz_method_round },
    { "__sizeof__", GMPy_MPZ_Method_SizeOf, METH_NOARGS, GMPy_doc_mpz_method_sizeof },
    { "__trunc__", GMPy_MPZ_Method_Trunc, METH_NOARGS, GMPy_doc_mpz_method_trunc },
    { "bit_clear", GMPy_MPZ_bit_clear_method, METH_O, doc_bit_clear_method },
    { "bit_flip", GMPy_MPZ_bit_flip_method, METH_O, doc_bit_flip_method },
    { "bit_length", GMPy_MPZ_bit_length_method, METH_NOARGS, doc_bit_length_method },
    { "bit_scan0", GMPy_MPZ_bit_scan0_method, METH_VARARGS, doc_bit_scan0_method },
    { "bit_scan1", GMPy_MPZ_bit_scan1_method, METH_VARARGS, doc_bit_scan1_method },
    { "bit_set", GMPy_MPZ_bit_set_method, METH_O, doc_bit_set_method },
    { "bit_test", GMPy_MPZ_bit_test_method, METH_O, doc_bit_test_method },
    { "digits", GMPy_MPZ_Digits_Method, METH_VARARGS, GMPy_doc_mpz_digits_method },
    { "is_congruent", GMPy_MPZ_Method_IsCongruent, METH_VARARGS, GMPy_doc_mpz_method_is_congruent },
    { "is_divisible", GMPy_MPZ_Method_IsDivisible, METH_O, GMPy_doc_mpz_method_is_divisible },
    { "num_digits", GMPy_MPZ_Method_NumDigits, METH_VARARGS, GMPy_doc_mpz_method_num_digits },
    { NULL, NULL, 1 }
};

static PyTypeObject MPZ_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(0, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "mpz",                                  /* tp_name          */
    sizeof(MPZ_Object),                     /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPZ_Dealloc,          /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPZ_Repr_Slot,          /* tp_repr          */
    &GMPy_MPZ_number_methods,               /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &GMPy_MPZ_mapping_methods,              /* tp_as_mapping    */
    (hashfunc) GMPy_MPZ_Hash_Slot,          /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPZ_Str_Slot,           /* tp_str           */
        0,                                  /* tp_getattro      */
        0,                                  /* tp_setattro      */
        0,                                  /* tp_as_buffer     */
#ifdef PY3
    Py_TPFLAGS_DEFAULT,                     /* tp_flags         */
#else
    Py_TPFLAGS_HAVE_INDEX|Py_TPFLAGS_HAVE_RICHCOMPARE| \
    Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_HAVE_CLASS| \
    Py_TPFLAGS_HAVE_INPLACEOPS,
#endif
    "Multiple precision integer",           /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&GMPy_RichCompare_Slot,    /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    GMPy_MPZ_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    GMPy_MPZ_getseters,                        /* tp_getset        */
};

