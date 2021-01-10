/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpq.c                                                             *
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

PyDoc_STRVAR(GMPy_doc_mpq,
"mpq() -> mpq(0,1)\n\n"
"     If no argument is given, return mpq(0,1).\n\n"
"mpq(n) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n. Fraction values\n"
"     are converted exactly.\n\n"
"mpq(n,m) -> mpq\n\n"
"     Return an 'mpq' object with a numeric value n/m.\n\n"
"mpq(s[, base=10]) -> mpq\n\n"
"     Return an 'mpq' object from a string s made up of digits in\n"
"     the given base. s may be made up of two numbers in the same\n"
"     base separated by a '/' character.\n");

/* Since `gmpy2.mpq` is now a type and no longer a factory function, see
 * gmpy2_cache.c/GMPy_MPQ_NewInit for details on creation.
 */

#ifdef PY3
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) GMPy_Number_Add_Slot,      /* nb_add                  */
    (binaryfunc) GMPy_Number_Sub_Slot,      /* nb_subtract             */
    (binaryfunc) GMPy_Number_Mul_Slot,      /* nb_multiply             */
    (binaryfunc) GMPy_Number_Mod_Slot,      /* nb_remainder            */
    (binaryfunc) GMPy_Number_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_Number_Pow_Slot,      /* nb_power                */
    (unaryfunc) GMPy_MPQ_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPQ_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPQ_Abs_Slot,          /* nb_absolute             */
    (inquiry) GMPy_MPQ_NonZero_Slot,        /* nb_bool                 */
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
    (binaryfunc) GMPy_Number_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_Number_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                  /* nb_inplace_floor_divide */
        0,                                  /* nb_inplace_true_divide  */
        0,                                  /* nb_index                */
};
#else
static PyNumberMethods mpq_number_methods =
{
    (binaryfunc) GMPy_Number_Add_Slot,      /* nb_add                  */
    (binaryfunc) GMPy_Number_Sub_Slot,      /* nb_subtract             */
    (binaryfunc) GMPy_Number_Mul_Slot,      /* nb_multiply             */
    (binaryfunc) GMPy_Number_TrueDiv_Slot,     /* nb_divide               */
    (binaryfunc) GMPy_Number_Mod_Slot,      /* nb_remainder            */
    (binaryfunc) GMPy_Number_DivMod_Slot,      /* nb_divmod               */
    (ternaryfunc) GMPy_Number_Pow_Slot,      /* nb_power                */
    (unaryfunc) GMPy_MPQ_Minus_Slot,        /* nb_negative             */
    (unaryfunc) GMPy_MPQ_Plus_Slot,         /* nb_positive             */
    (unaryfunc) GMPy_MPQ_Abs_Slot,          /* nb_absolute             */
    (inquiry) GMPy_MPQ_NonZero_Slot,        /* nb_bool                 */
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
    (binaryfunc) GMPy_Number_FloorDiv_Slot,    /* nb_floor_divide         */
    (binaryfunc) GMPy_Number_TrueDiv_Slot,     /* nb_true_divide          */
        0,                                  /* nb_inplace_floor_divide */
        0,                                  /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef GMPy_MPQ_getseters[] =
{
    { "numerator", (getter)GMPy_MPQ_Attrib_GetNumer, NULL,
        "the numerator of a rational number in lowest terms", NULL },
    { "denominator", (getter)GMPy_MPQ_Attrib_GetDenom, NULL,
        "the denominator of a rational number in lowest terms", NULL },
    { "real", (getter)GMPy_MPQ_Attrib_GetReal, NULL,
        "the real part of a complex number", NULL },
    { "imag", (getter)GMPy_MPQ_Attrib_GetImag, NULL,
        "the imaginary part of a complex number", NULL },
    {NULL}
};

static PyMethodDef GMPy_MPQ_methods [] =
{
    { "__ceil__", GMPy_MPQ_Method_Ceil, METH_NOARGS, GMPy_doc_mpq_method_ceil },
    { "__floor__", GMPy_MPQ_Method_Floor, METH_NOARGS, GMPy_doc_mpq_method_floor },
    { "__round__", GMPy_MPQ_Method_Round, METH_VARARGS, GMPy_doc_mpq_method_round },
    { "__sizeof__", GMPy_MPQ_Method_Sizeof, METH_NOARGS, GMPy_doc_mpq_method_sizeof },
    { "__trunc__", GMPy_MPQ_Method_Trunc, METH_NOARGS, GMPy_doc_mpq_method_trunc },
    { "conjugate", GMPy_MP_Method_Conjugate, METH_NOARGS, GMPy_doc_mp_method_conjugate },
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
    GMPy_doc_mpq,                           /* tp_doc           */
        0,                                  /* tp_traverse      */
        0,                                  /* tp_clear         */
    (richcmpfunc)&GMPy_RichCompare_Slot,    /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    GMPy_MPQ_methods,                       /* tp_methods       */
        0,                                  /* tp_members       */
    GMPy_MPQ_getseters,                     /* tp_getset        */
        0,                                  /* tp_base          */
        0,                                  /* tp_dict          */
        0,                                  /* tp_descr_get     */
        0,                                  /* tp_descr_set     */
        0,                                  /* tp_dictoffset    */
        0,                                  /* tp_init          */
        0,                                  /* tp_alloc         */
    GMPy_MPQ_NewInit,                       /* tp_new           */
        0,                                  /* tp_free          */
};

