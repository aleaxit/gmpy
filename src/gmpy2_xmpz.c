/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz.c                                                            *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
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

PyDoc_STRVAR(GMPy_doc_xmpz,
"xmpz() -> xmpz(0)\n\n"
"     If no argument is given, return xmpz(0).\n\n"
"xmpz(n) -> xmpz\n\n"
"     Return an 'xmpz' object with a numeric value 'n' (truncating n\n"
"     to its integer part if it's a Fraction, 'mpq', float or 'mpfr').\n\n"
"xmpz(s[, base=0]):\n\n"
"     Return an 'xmpz' object from a string 's' made of digits in the\n"
"     given base.  If base=0, binary, octal, or hex Python strings\n"
"     are recognized by leading 0b, 0o, or 0x characters, otherwise\n"
"     the string is assumed to be decimal. Values for base can range\n"
"     between 2 and 62.\n\n"
"     Note: 'xmpz' is a mutable integer. It can be faster when used\n"
"     for augmented assignment (+=, *=, etc.). 'xmpz' objects cannot\n"
"     be used as dictionary keys. The use of 'mpz' objects is recommended\n"
"     in most cases.");


#ifdef PY3
static PyNumberMethods GMPy_XMPZ_number_methods =
{
    (binaryfunc) GMPy_Number_Add_Slot,       /* nb_add                  */
    (binaryfunc) GMPy_Number_Sub_Slot,       /* nb_subtract             */
    (binaryfunc) GMPy_Number_Mul_Slot,       /* nb_multiply             */
    (binaryfunc) GMPy_Number_Mod_Slot,       /* nb_remainder            */
    (binaryfunc) GMPy_Number_DivMod_Slot,    /* nb_divmod               */
    (ternaryfunc) GMPy_Number_Pow_Slot,      /* nb_power                */
    (unaryfunc) GMPy_XMPZ_Neg_Slot,          /* nb_negative             */
    (unaryfunc) GMPy_XMPZ_Pos_Slot,          /* nb_positive             */
    (unaryfunc) GMPy_XMPZ_Abs_Slot,          /* nb_absolute             */
    (inquiry) GMPy_XMPZ_NonZero_Slot,        /* nb_bool                 */
    (unaryfunc) GMPy_XMPZ_Com_Slot,          /* nb_invert               */
    (binaryfunc) GMPy_MPZ_Lshift_Slot,       /* nb_lshift               */
    (binaryfunc) GMPy_MPZ_Rshift_Slot,       /* nb_rshift               */
    (binaryfunc) GMPy_MPZ_And_Slot,          /* nb_and                  */
    (binaryfunc) GMPy_MPZ_Xor_Slot,          /* nb_xor                  */
    (binaryfunc) GMPy_MPZ_Ior_Slot,          /* nb_or                   */
    (unaryfunc) GMPy_MPZ_Int_Slot,           /* nb_int                  */
        0,                                   /* nb_reserved             */
    (unaryfunc) GMPy_MPZ_Float_Slot,         /* nb_float                */
    (binaryfunc) GMPy_XMPZ_IAdd_Slot,        /* nb_inplace_add          */
    (binaryfunc) GMPy_XMPZ_ISub_Slot,        /* nb_inplace_subtract     */
    (binaryfunc) GMPy_XMPZ_IMul_Slot,        /* nb_inplace_multiply     */
    (binaryfunc) GMPy_XMPZ_IRem_Slot,        /* nb_inplace_remainder    */
    (ternaryfunc) GMPy_XMPZ_IPow_Slot,       /* nb_inplace_power        */
    (binaryfunc) GMPy_XMPZ_ILshift_Slot,     /* nb_inplace_lshift       */
    (binaryfunc) GMPy_XMPZ_IRshift_Slot,     /* nb_inplace_rshift       */
    (binaryfunc) GMPy_XMPZ_IAnd_Slot,        /* nb_inplace_and          */
    (binaryfunc) GMPy_XMPZ_IXor_Slot,        /* nb_inplace_xor          */
    (binaryfunc) GMPy_XMPZ_IIor_Slot,        /* nb_inplace_or           */
    (binaryfunc) GMPy_Number_FloorDiv_Slot,  /* nb_floor_divide         */
    (binaryfunc) GMPy_Number_TrueDiv_Slot,   /* nb_true_divide          */
    (binaryfunc) GMPy_XMPZ_IFloorDiv_Slot,   /* nb_inplace_floor_divide */
        0,                                   /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_MPZ_Int_Slot,           /* nb_index                */
};

#else
static PyNumberMethods GMPy_XMPZ_number_methods =
{
    (binaryfunc) GMPy_Number_Add_Slot,       /* nb_add                  */
    (binaryfunc) GMPy_Number_Sub_Slot,       /* nb_subtract             */
    (binaryfunc) GMPy_Number_Mul_Slot,       /* nb_multiply             */
    (binaryfunc) GMPy_Number_Div2_Slot,      /* nb_divide               */
    (binaryfunc) GMPy_Number_Mod_Slot,       /* nb_remainder            */
    (binaryfunc) GMPy_Number_DivMod_Slot,    /* nb_divmod               */
    (ternaryfunc) GMPy_Number_Pow_Slot,      /* nb_power                */
    (unaryfunc) GMPy_XMPZ_Neg_Slot,          /* nb_negative             */
    (unaryfunc) GMPy_XMPZ_Pos_Slot,          /* nb_positive             */
    (unaryfunc) GMPy_XMPZ_Abs_Slot,          /* nb_absolute             */
    (inquiry) GMPy_XMPZ_NonZero_Slot,        /* nb_bool                 */
    (unaryfunc) GMPy_XMPZ_Com_Slot,          /* nb_invert               */
    (binaryfunc) GMPy_MPZ_Lshift_Slot,       /* nb_lshift               */
    (binaryfunc) GMPy_MPZ_Rshift_Slot,       /* nb_rshift               */
    (binaryfunc) GMPy_MPZ_And_Slot,          /* nb_and                  */
    (binaryfunc) GMPy_MPZ_Xor_Slot,          /* nb_xor                  */
    (binaryfunc) GMPy_MPZ_Ior_Slot,          /* nb_or                   */
        0,                                   /* nb_coerce               */
    (unaryfunc) GMPy_MPZ_Int_Slot,           /* nb_int                  */
    (unaryfunc) GMPy_MPZ_Long_Slot,          /* nb_long                 */
    (unaryfunc) GMPy_MPZ_Float_Slot,         /* nb_float                */
    (unaryfunc) GMPy_XMPZ_Oct_Slot,          /* nb_oct                  */
    (unaryfunc) GMPy_XMPZ_Hex_Slot,          /* nb_hex                  */
    (binaryfunc) GMPy_XMPZ_IAdd_Slot,        /* nb_inplace_add          */
    (binaryfunc) GMPy_XMPZ_ISub_Slot,        /* nb_inplace_subtract     */
    (binaryfunc) GMPy_XMPZ_IMul_Slot,        /* nb_inplace_multiply     */
        0,                                   /* nb_inplace_divide       */
    (binaryfunc) GMPy_XMPZ_IRem_Slot,        /* nb_inplace_remainder    */
    (ternaryfunc) GMPy_XMPZ_IPow_Slot,       /* nb_inplace_power        */
    (binaryfunc) GMPy_XMPZ_ILshift_Slot,     /* nb_inplace_lshift       */
    (binaryfunc) GMPy_XMPZ_IRshift_Slot,     /* nb_inplace_rshift       */
    (binaryfunc) GMPy_XMPZ_IAnd_Slot,        /* nb_inplace_and          */
    (binaryfunc) GMPy_XMPZ_IXor_Slot,        /* nb_inplace_xor          */
    (binaryfunc) GMPy_XMPZ_IIor_Slot,        /* nb_inplace_or           */
    (binaryfunc) GMPy_Number_FloorDiv_Slot,  /* nb_floor_divide         */
    (binaryfunc) GMPy_Number_TrueDiv_Slot,   /* nb_true_divide          */
    (binaryfunc) GMPy_XMPZ_IFloorDiv_Slot,   /* nb_inplace_floor_divide */
        0,                                   /* nb_inplace_true_divide  */
    (unaryfunc) GMPy_MPZ_Int_Slot,           /* nb_index                */
};
#endif

static PyMappingMethods GMPy_XMPZ_mapping_methods = {
    (lenfunc)GMPy_XMPZ_Method_Length,
    (binaryfunc)GMPy_XMPZ_Method_SubScript,
    (objobjargproc)GMPy_XMPZ_Method_AssignSubScript
};

static PyGetSetDef GMPy_XMPZ_getseters[] =
{
    { "numerator", (getter)GMPy_XMPZ_Attrib_GetNumer, NULL,
        "the numerator of a rational number in lowest terms", NULL },
    { "denominator", (getter)GMPy_XMPZ_Attrib_GetDenom, NULL,
        "the denominator of a rational number in lowest terms", NULL },
    { "real", (getter)GMPy_XMPZ_Attrib_GetReal, NULL,
        "the real part of a complex number", NULL },
    { "denominator", (getter)GMPy_XMPZ_Attrib_GetImag, NULL,
        "the imaginary part of a complex number", NULL },
    {NULL}
};

static PyMethodDef GMPy_XMPZ_methods [] =
{
    { "__format__", GMPy_MPZ_Format, METH_VARARGS, GMPy_doc_mpz_format },
    { "__sizeof__", GMPy_XMPZ_Method_SizeOf, METH_NOARGS, GMPy_doc_xmpz_method_sizeof },
    { "bit_clear", GMPy_MPZ_bit_clear_method, METH_O, doc_bit_clear_method },
    { "bit_flip", GMPy_MPZ_bit_flip_method, METH_O, doc_bit_flip_method },
    { "bit_length", GMPy_MPZ_bit_length_method, METH_NOARGS, doc_bit_length_method },
    { "bit_scan0", GMPy_MPZ_bit_scan0_method, METH_VARARGS, doc_bit_scan0_method },
    { "bit_scan1", GMPy_MPZ_bit_scan1_method, METH_VARARGS, doc_bit_scan1_method },
    { "bit_set", GMPy_MPZ_bit_set_method, METH_O, doc_bit_set_method },
    { "bit_test", GMPy_MPZ_bit_test_method, METH_O, doc_bit_test_method },
    { "conjugate", GMPy_MP_Method_Conjugate, METH_NOARGS, GMPy_doc_mp_method_conjugate },
    { "copy", GMPy_XMPZ_Method_Copy, METH_NOARGS, GMPy_doc_xmpz_method_copy },
    { "digits", GMPy_XMPZ_Digits_Method, METH_VARARGS, GMPy_doc_mpz_digits_method },
    { "iter_bits", (PyCFunction)GMPy_XMPZ_Method_IterBits, METH_VARARGS | METH_KEYWORDS, GMPy_doc_xmpz_method_iter_bits },
    { "iter_clear", (PyCFunction)GMPy_XMPZ_Method_IterClear, METH_VARARGS | METH_KEYWORDS, GMPy_doc_xmpz_method_iter_clear },
    { "iter_set", (PyCFunction)GMPy_XMPZ_Method_IterSet, METH_VARARGS | METH_KEYWORDS, GMPy_doc_xmpz_method_iter_set },
    { "make_mpz", GMPy_XMPZ_Method_MakeMPZ, METH_NOARGS, GMPy_doc_xmpz_method_make_mpz },
    { "num_digits", GMPy_MPZ_Method_NumDigits, METH_VARARGS, GMPy_doc_mpz_method_num_digits },
    { "num_limbs", GMPy_XMPZ_Method_NumLimbs, METH_NOARGS, GMPy_doc_xmpz_method_num_limbs },
    { "limbs_read", GMPy_XMPZ_Method_LimbsRead, METH_NOARGS, GMPy_doc_xmpz_method_limbs_read },
    { "limbs_write", GMPy_XMPZ_Method_LimbsWrite, METH_O, GMPy_doc_xmpz_method_limbs_write },
    { "limbs_modify", GMPy_XMPZ_Method_LimbsModify, METH_O, GMPy_doc_xmpz_method_limbs_modify },
    { "limbs_finish", GMPy_XMPZ_Method_LimbsFinish, METH_O, GMPy_doc_xmpz_method_limbs_finish },
    { NULL, NULL, 1 }
};

static PyTypeObject XMPZ_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
        0,                                  /* ob_size          */
#endif
    "xmpz",                                 /* tp_name          */
    sizeof(XMPZ_Object),                    /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_XMPZ_Dealloc,         /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_XMPZ_Repr_Slot,         /* tp_repr          */
    &GMPy_XMPZ_number_methods,              /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
    &GMPy_XMPZ_mapping_methods,             /* tp_as_mapping    */
        0,                                  /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_XMPZ_Str_Slot,          /* tp_str           */
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
    GMPy_doc_xmpz,                          /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&GMPy_RichCompare_Slot,    /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    GMPy_XMPZ_methods,                      /* tp_methods       */
        0,                                  /* tp_members       */
    GMPy_XMPZ_getseters,                    /* tp_getset        */
        0,                                  /* tp_base          */
        0,                                  /* tp_dict          */
        0,                                  /* tp_descr_get     */
        0,                                  /* tp_descr_set     */
        0,                                  /* tp_dictoffset    */
        0,                                  /* tp_init          */
        0,                                  /* tp_alloc         */
    GMPy_XMPZ_NewInit,                       /* tp_new           */
        0,                                  /* tp_free          */
};



