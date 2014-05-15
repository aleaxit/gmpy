/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpc.c                                                              *
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

PyDoc_STRVAR(GMPy_doc_mpc_factory,
"mpc() -> mpc(0.0+0.0j)\n\n"
"      If no argument is given, return mpc(0.0+0.0j).\n\n"
"mpc(c [, precision=0]) -> mpc\n\n"
"      Return a new 'mpc' object from an existing complex number (either\n"
"      a Python complex object or another 'mpc' object).\n\n"
"mpc(real [,imag=0 [, precision=0]]) -> mpc\n\n"
"      Return a new 'mpc' object by converting two non-complex numbers\n"
"      into the real and imaginary components of an 'mpc' object.\n\n"
"mpc(s [, precision=0 [, base=10]]) -> mpc\n\n"
"      Return a new 'mpc' object by converting a string s into a complex\n"
"      number. If base is omitted, then a base-10 representation is\n"
"      assumed otherwise the base must be in the interval [2,36].\n\n"
"Note: The precision can be specified either a single number that\n"
"      is used for both the real and imaginary components, or as a\n"
"      tuple that can specify different precisions for the real\n"
"      and imaginary components.\n\n"
"      If a precision greater than or equal to 2 is specified, then it\n"
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
GMPy_MPC_Factory(PyObject *self, PyObject *args, PyObject *kwargs)
{
    MPC_Object *result = NULL;
    MPFR_Object *tempreal = NULL, *tempimag = NULL;
    PyObject *arg0 = NULL, *arg1 = NULL, *prec = NULL;
    int base = 10;
    Py_ssize_t argc = 0, keywdc = 0;
    CTXT_Object *context = NULL;

    /* Assumes mpfr_prec_t is the same as a long. */
    mpfr_prec_t rprec = 0, iprec = 0;

    static char *kwlist_c[] = {"c", "precision", NULL};
    static char *kwlist_r[] = {"real", "imag", "precision", NULL};
    static char *kwlist_s[] = {"s", "precision", "base", NULL};

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    argc = PyTuple_Size(args);
    if (kwargs) {
        keywdc = PyDict_Size(kwargs);
    }

    if (argc + keywdc > 3) {
        TYPE_ERROR("mpc() takes at most 3 arguments");
        return NULL;
    }

    if (argc + keywdc == 0) {
        if ((result = GMPy_MPC_New(0, 0, context))) {
            mpc_set_ui(result->c, 0, GET_MPC_ROUND(context));
        }
        return (PyObject*)result;
    }

    if (argc == 0) {
        TYPE_ERROR("mpc() requires at least one non-keyword argument");
        return NULL;
    }

    arg0 = PyTuple_GET_ITEM(args, 0);

    /* If building an mpc from a string, there can be upto two additional
     * arguments. Note that precision can be either a single integer or
     * a tuple.
     */
    if (PyStrOrUnicode_Check(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|Oi", kwlist_s,
                                              &arg0, &prec, &base)))
                return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            else {
                TYPE_ERROR("precision for mpc() must be integer or tuple");
                return NULL;
            }

            if (rprec < 0 || iprec < 0) {
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc()");
                }
                else {
                    VALUE_ERROR("precision for mpc() must be >= 0");
                }
                return NULL;
            }
        }

        if (base < 2 || base > 36) {
            VALUE_ERROR("base for mpc() must be in the interval [2,36]");
            return NULL;
        }

        return (PyObject*)GMPy_MPC_From_PyStr(arg0, base, rprec, iprec, context);
    }

    if (IS_REAL(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|OO", kwlist_r,
                                            &arg0, &arg1, &prec)))
                return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            else {
                TYPE_ERROR("precision for mpc() must be integer or tuple");
                return NULL;
            }

            if (rprec < 0 || iprec < 0) {
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc()");
                }
                else {
                    VALUE_ERROR("precision for mpc() must be >= 0");
                }
                return NULL;
            }
        }

        if (arg1 && !IS_REAL(arg1)) {
            TYPE_ERROR("invalid type for imaginary component in mpc()");
            return NULL;
        }

        tempreal = GMPy_MPFR_From_Real(arg0, rprec, context);
        if (arg1) {
            tempimag = GMPy_MPFR_From_Real(arg1, iprec, context);
        }
        else {
            if ((tempimag = GMPy_MPFR_New(iprec, context))) {
                mpfr_set_ui(tempimag->f, 0, MPFR_RNDN);
            }
        }

        result = GMPy_MPC_New(rprec, iprec, context);
        if (!tempreal || !tempimag || !result) {
            Py_XDECREF(tempreal);
            Py_XDECREF(tempimag);
            Py_XDECREF(result);
            TYPE_ERROR("mpc() requires string or numeric argument.");
            return NULL;
        }

        mpc_set_fr_fr(result->c, tempreal->f, tempimag->f, GET_MPC_ROUND(context));
        Py_DECREF(tempreal);
        Py_DECREF(tempimag);
        return (PyObject*)result;
    }

    if (IS_COMPLEX_ONLY(arg0)) {
        if (keywdc || argc > 1) {
            if (!(PyArg_ParseTupleAndKeywords(args, kwargs, "O|O", kwlist_c,
                                          &arg0, &prec)))
            return NULL;
        }

        if (prec) {
            if (PyIntOrLong_Check(prec)) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(prec);
                iprec = rprec;
            }
            else if (PyTuple_Check(prec) && PyTuple_Size(prec) == 2) {
                rprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 0));
                iprec = (mpfr_prec_t)PyIntOrLong_AsLong(PyTuple_GET_ITEM(prec, 1));
            }
            else {
                TYPE_ERROR("precision for mpc() must be integer or tuple");
                return NULL;
            }

            if (rprec < 0 || iprec < 0) {
                if (PyErr_Occurred()) {
                    VALUE_ERROR("invalid value for precision in mpc()");
                }
                else {
                    VALUE_ERROR("precision for mpc() must be >= 0");
                }
                return NULL;
            }
        }

        if (PyComplex_Check(arg0)) {
            result = GMPy_MPC_From_PyComplex(arg0, rprec, iprec, context);
        }
        else {
            result = GMPy_MPC_From_MPC((MPC_Object*)arg0, rprec, iprec, context);
        }
        return (PyObject*)result;
    }

    TYPE_ERROR("mpc() requires numeric or string argument");
    return NULL;
}

/* Implement the conjugate() method. */

PyDoc_STRVAR(GMPy_doc_mpc_conjugate_method,
"x.conjugate() -> mpc\n\n"
"Returns the conjugate of x.");

static PyObject *
GMPy_MPC_Conjugate_Method(PyObject *self, PyObject *args)
{
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        return NULL;
    }

    result->rc = mpc_conj(result->c, MPC(self), GET_MPC_ROUND(context));

    GMPY_MPC_CLEANUP(result, context, "conjugate()");
    return (PyObject*)result;
}

/* Implement the .precision attribute of an mpfr. */

static PyObject *
GMPy_MPC_GetPrec_Attrib(MPC_Object *self, void *closure)
{
    mpfr_prec_t rprec = 0, iprec = 0;

    mpc_get_prec2(&rprec, &iprec, self->c);
    return Py_BuildValue("(nn)", rprec, iprec);
}

/* Implement the .rc attribute of an mpfr. */

static PyObject *
GMPy_MPC_GetRc_Attrib(MPC_Object *self, void *closure)
{
    return Py_BuildValue("(ii)", MPC_INEX_RE(self->rc), MPC_INEX_IM(self->rc));
}

/* Implement the .imag attribute of an mpfr. */

static PyObject *
GMPy_MPC_GetImag_Attrib(MPC_Object *self, void *closure)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(0, context))) {
        result->rc = mpc_imag(result->f, self->c, GET_MPFR_ROUND(context));
        GMPY_MPFR_CLEANUP(result, context, "imag()");
    }
    return (PyObject*)result;
}

/* Implement the .real attribute of an mpfr. */

static PyObject *
GMPy_MPC_GetReal_Attrib(MPC_Object *self, void *closure)
{
    MPFR_Object *result = NULL;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT(context);

    if ((result = GMPy_MPFR_New(0, context))) {
        result->rc = mpc_real(result->f, self->c, context->ctx.mpfr_round);
        GMPY_MPFR_CLEANUP(result, context, "real()");
    }
    return (PyObject*)result;
}

/* Implement the nb_bool slot. */

static int
GMPy_MPC_NonZero_Slot(MPC_Object *self)
{
    return !MPC_IS_ZERO_P(self->c);
}

PyDoc_STRVAR(doc_mpc_phase,
"phase(x) -> mpfr\n\n"
"Return the phase angle, also known as argument, of a complex x.");

static PyObject *
Pympc_phase(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("phase() requires 'mpc' argument");

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_arg(result->f, MPC(self),
                         context->ctx.mpfr_round);
    Py_DECREF((PyObject*)self);

    MPFR_SUBNORMALIZE(result);
    MPFR_CHECK_OVERFLOW(result, "overflow in 'mpc' phase()");
    MPFR_CHECK_INVALID(result, "invalid operation 'mpc' phase()");
    MPFR_CHECK_UNDERFLOW(result, "underflow in 'mpc' phase()");
    MPFR_CHECK_INEXACT(result, "inexact operation in 'mpc' phase()");
  done:
    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpc_norm,
"norm(x) -> mpfr\n\n"
"Return the norm of a complex x. The norm(x) is defined as\n"
"x.real**2 + x.imag**2. abs(x) is the square root of norm(x).\n");

static PyObject *
Pympc_norm(PyObject *self, PyObject *other)
{
    MPFR_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("norm() requires 'mpc' argument");

    if (!(result = GMPy_MPFR_New(0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_norm(result->f, MPC(self),
                          context->ctx.mpfr_round);
    Py_DECREF((PyObject*)self);

    MPFR_SUBNORMALIZE(result);
    MPFR_CHECK_OVERFLOW(result, "overflow in 'mpc' norm()");
    MPFR_CHECK_INVALID(result, "invalid operation 'mpc' norm()");
    MPFR_CHECK_UNDERFLOW(result, "underflow in 'mpc' norm()");
    MPFR_CHECK_INEXACT(result, "inexact operation in 'mpc' norm()");
  done:
    if (PyErr_Occurred()) {
        Py_DECREF((PyObject*)result);
        result = NULL;
    }
    return (PyObject*)result;
}

PyDoc_STRVAR(doc_mpc_polar,
"polar(x) -> (abs(x), phase(x))\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

static PyObject *
Pympc_polar(PyObject *self, PyObject *other)
{
    PyObject *abs, *phase, *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("norm() requires 'mpc' argument");

    if (!(abs = GMPy_Complex_Abs(self, context))) {
        Py_DECREF(self);
        return NULL;
    }
    if (!(phase = Pympc_phase(self, other))) {
        Py_DECREF(abs);
        Py_DECREF(self);
        return NULL;
    }

    result = Py_BuildValue("(NN)", abs, phase);
    if (!result) {
        Py_DECREF(abs);
        Py_DECREF(phase);
    }
    Py_DECREF(self);
    return result;
}

PyDoc_STRVAR(doc_mpc_rect,
"rect(x) -> mpc\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

/* Note: does not properly check for inexact or underflow */

static PyObject *
Pympc_rect(PyObject *self, PyObject *args)
{
    PyObject *other;
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_TWO_MPFR_ARGS(other, "rect() requires 'mpfr','mpfr' arguments");

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        Py_DECREF(other);
        return NULL;
    }

    mpfr_cos(mpc_realref(result->c), MPFR(other),
             GET_REAL_ROUND(context));
    mpfr_mul(mpc_realref(result->c), mpc_realref(result->c),
             MPFR(self), GET_REAL_ROUND(context));
    mpfr_sin(mpc_imagref(result->c), MPFR(other),
             GET_IMAG_ROUND(context));
    mpfr_mul(mpc_imagref(result->c), mpc_imagref(result->c),
             MPFR(self), GET_IMAG_ROUND(context));
    Py_DECREF(self);
    Py_DECREF(other);

    MPC_CLEANUP(result, "rect()");
}

PyDoc_STRVAR(doc_mpc_proj,
"proj(x) -> mpc\n\n"
"Returns the projection of a complex x on to the Riemann sphere.");

static PyObject *
Pympc_proj(PyObject *self, PyObject *other)
{
    MPC_Object *result;
    CTXT_Object *context = NULL;

    CHECK_CONTEXT_SET_EXPONENT(context);

    PARSE_ONE_MPC_OTHER("proj() requires 'mpc' argument");

    if (!(result = GMPy_MPC_New(0, 0, context))) {
        Py_DECREF(self);
        return NULL;
    }

    result->rc = mpc_proj(result->c, MPC(self),
                          GET_MPC_ROUND(context));
    Py_DECREF(self);

    MPC_CLEANUP(result, "proj()");
}

PyDoc_STRVAR(GMPy_doc_mpc_sizeof_method,
"x.__sizeof__()\n\n"
"Returns the amount of memory consumed by x.");

static PyObject *
GMPy_MPC_SizeOf_Method(PyObject *self, PyObject *other)
{
    return PyIntOrLong_FromSize_t(sizeof(MPC_Object) + \
        (((mpc_realref(MPC(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t)) + \
        (((mpc_imagref(MPC(self))->_mpfr_prec + mp_bits_per_limb - 1) / \
        mp_bits_per_limb) * sizeof(mp_limb_t)));
}

static PyMethodDef Pympc_methods[] =
{
    { "__complex__", GMPy_PyComplex_From_MPC, METH_O, GMPy_doc_mpc_complex },
    { "__format__", GMPy_MPC_Format, METH_VARARGS, GMPy_doc_mpc_format },
    { "__sizeof__", GMPy_MPC_SizeOf_Method, METH_NOARGS, GMPy_doc_mpc_sizeof_method },
    { "conjugate", GMPy_MPC_Conjugate_Method, METH_NOARGS, GMPy_doc_mpc_conjugate_method },
    { "digits", GMPy_MPC_Digits_Method, METH_VARARGS, GMPy_doc_mpc_digits_method },
    { "is_finite", GMPy_MPC_Is_Finite_Method, METH_NOARGS, GMPy_doc_method_is_finite },
    { "is_infinite", GMPy_MPC_Is_Infinite_Method, METH_NOARGS, GMPy_doc_method_is_infinite },
    { "is_nan", GMPy_MPC_Is_NAN_Method, METH_NOARGS, GMPy_doc_method_is_nan },
    { "is_zero", GMPy_MPC_Is_Zero_Method, METH_NOARGS, GMPy_doc_method_is_zero },
    { NULL, NULL, 1 }
};


#ifdef PY3
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) GMPy_MPC_Add_Slot,      /* nb_add                  */
    (binaryfunc) GMPy_MPC_Sub_Slot,      /* nb_subtract             */
    (binaryfunc) GMPy_MPC_Mul_Slot,      /* nb_multiply             */
    (binaryfunc) GMPy_MPC_Mod_Slot,      /* nb_remainder            */
    (binaryfunc) GMPy_MPC_DivMod_Slot,   /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,   /* nb_power                */
    (unaryfunc) GMPy_MPC_Minus_Slot,     /* nb_negative             */
    (unaryfunc) GMPy_MPC_Plus_Slot,      /* nb_positive             */
    (unaryfunc) GMPy_MPC_Abs_Slot,       /* nb_absolute             */
    (inquiry) GMPy_MPC_NonZero_Slot,     /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
    (unaryfunc) GMPy_MPC_Int_Slot,       /* nb_int                  */
        0,                               /* nb_reserved             */
    (unaryfunc) GMPy_MPC_Float_Slot,     /* nb_float                */
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
    (binaryfunc) GMPy_MPC_FloorDiv_Slot, /* nb_floor_divide         */
    (binaryfunc) GMPy_MPC_TrueDiv_Slot,  /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
        0,                               /* nb_index                */
};
#else
static PyNumberMethods mpc_number_methods =
{
    (binaryfunc) GMPy_MPC_Add_Slot,      /* nb_add                  */
    (binaryfunc) GMPy_MPC_Sub_Slot,      /* nb_subtract             */
    (binaryfunc) GMPy_MPC_Mul_Slot,      /* nb_multiply             */
    (binaryfunc) GMPy_MPC_TrueDiv_Slot,  /* nb_divide               */
    (binaryfunc) GMPy_MPC_Mod_Slot,      /* nb_remainder            */
    (binaryfunc) GMPy_MPC_DivMod_Slot,   /* nb_divmod               */
    (ternaryfunc) GMPy_MPANY_Pow_Slot,   /* nb_power                */
    (unaryfunc) GMPy_MPC_Minus_Slot,     /* nb_negative             */
    (unaryfunc) GMPy_MPC_Plus_Slot,      /* nb_positive             */
    (unaryfunc) GMPy_MPC_Abs_Slot,       /* nb_absolute             */
    (inquiry) GMPy_MPC_NonZero_Slot,     /* nb_bool                 */
        0,                               /* nb_invert               */
        0,                               /* nb_lshift               */
        0,                               /* nb_rshift               */
        0,                               /* nb_and                  */
        0,                               /* nb_xor                  */
        0,                               /* nb_or                   */
        0,                               /* nb_coerce               */
    (unaryfunc) GMPy_MPC_Int_Slot,       /* nb_int                  */
    (unaryfunc) GMPy_MPC_Long_Slot,      /* nb_long                 */
    (unaryfunc) GMPy_MPC_Float_Slot,     /* nb_float                */
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
    (binaryfunc) GMPy_MPC_FloorDiv_Slot, /* nb_floor_divide         */
    (binaryfunc) GMPy_MPC_TrueDiv_Slot,  /* nb_true_divide          */
        0,                               /* nb_inplace_floor_divide */
        0,                               /* nb_inplace_true_divide  */
};
#endif

static PyGetSetDef Pympc_getseters[] =
{
    {"precision", (getter)GMPy_MPC_GetPrec_Attrib, NULL, "precision in bits", NULL},
    {"rc", (getter)GMPy_MPC_GetRc_Attrib, NULL, "return code", NULL},
    {"imag", (getter)GMPy_MPC_GetImag_Attrib, NULL, "imaginary component", NULL},
    {"real", (getter)GMPy_MPC_GetReal_Attrib, NULL, "real component", NULL},
    {NULL}
};

static PyTypeObject MPC_Type =
{
    /* PyObject_HEAD_INIT(&PyType_Type) */
#ifdef PY3
    PyVarObject_HEAD_INIT(NULL, 0)
#else
    PyObject_HEAD_INIT(0)
    0,                                      /* ob_size          */
#endif
    "mpc",                                  /* tp_name          */
    sizeof(MPC_Object),                     /* tp_basicsize     */
        0,                                  /* tp_itemsize      */
    /* methods */
    (destructor) GMPy_MPC_Dealloc,          /* tp_dealloc       */
        0,                                  /* tp_print         */
        0,                                  /* tp_getattr       */
        0,                                  /* tp_setattr       */
        0,                                  /* tp_reserved      */
    (reprfunc) GMPy_MPC_Repr_Slot,          /* tp_repr          */
    &mpc_number_methods,                    /* tp_as_number     */
        0,                                  /* tp_as_sequence   */
        0,                                  /* tp_as_mapping    */
    (hashfunc) GMPy_MPC_Hash_Slot,          /* tp_hash          */
        0,                                  /* tp_call          */
    (reprfunc) GMPy_MPC_Str_Slot,           /* tp_str           */
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
    (richcmpfunc)&GMPy_RichCompare_Slot,    /* tp_richcompare   */
        0,                                  /* tp_weaklistoffset*/
        0,                                  /* tp_iter          */
        0,                                  /* tp_iternext      */
    Pympc_methods,                          /* tp_methods       */
        0,                                  /* tp_members       */
    Pympc_getseters,                        /* tp_getset        */
};


