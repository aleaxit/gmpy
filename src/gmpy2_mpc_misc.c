/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_mpc_misc.c                                                        *
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

PyDoc_STRVAR(GMPy_doc_context_phase,
"context.phase(x) -> mpfr\n\n"
"Return the phase angle, also known as argument, of a complex x.");

PyDoc_STRVAR(GMPy_doc_function_phase,
"phase(x) -> mpfr\n\n"
"Return the phase angle, also known as argument, of a complex x.");

static PyObject *
GMPy_Complex_Phase(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;
    MPC_Object *tempx;

    CHECK_CONTEXT(context)

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        return NULL;
    }

    result->rc = mpc_arg(result->f, tempx->c, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    GMPY_MPFR_CLEANUP(result, context, "phase()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_Phase(PyObject *x, CTXT_Object *context)
{
    if (IS_COMPLEX_ONLY(x))
        return GMPy_Complex_Phase(x, context);

    TYPE_ERROR("phase() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Phase(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("phase() requires 1 argument");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Phase(PyTuple_GET_ITEM(args, 0), context);
}

PyDoc_STRVAR(GMPy_doc_context_norm,
"context.norm(x) -> mpfr\n\n"
"Return the norm of a complex x. The norm(x) is defined as\n"
"x.real**2 + x.imag**2. abs(x) is the square root of norm(x).\n");

PyDoc_STRVAR(GMPy_doc_function_norm,
"norm(x) -> mpfr\n\n"
"Return the norm of a complex x. The norm(x) is defined as\n"
"x.real**2 + x.imag**2. abs(x) is the square root of norm(x).\n");

static PyObject *
GMPy_Complex_Norm(PyObject *x, CTXT_Object *context)
{
    MPFR_Object *result;
    MPC_Object *tempx;

    CHECK_CONTEXT(context);

    result = GMPy_MPFR_New(0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        return NULL;
    }

    mpfr_clear_flags();
    result->rc = mpc_norm(result->f, tempx->c, GET_MPFR_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    GMPY_MPFR_CLEANUP(result, context, "norm()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_Norm(PyObject *x, CTXT_Object *context)
{
    if (IS_COMPLEX_ONLY(x))
        return GMPy_Complex_Norm(x, context);

    TYPE_ERROR("norm() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Norm(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("norm() requires 1 argument");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Norm(PyTuple_GET_ITEM(args, 0), context);
}

PyDoc_STRVAR(GMPy_doc_context_polar,
"context.polar(x) -> (abs(x), phase(x))\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

PyDoc_STRVAR(GMPy_doc_function_polar,
"polar(x) -> (abs(x), phase(x))\n\n"
"Return the polar coordinate form of a complex x that is in\n"
"rectangular form.");

static PyObject *
GMPy_Complex_Polar(PyObject *x, CTXT_Object *context)
{
    PyObject *tempx, *abs, *phase, *result;

    CHECK_CONTEXT(context);

    if (!(tempx = (PyObject*)GMPy_MPC_From_Complex(x, 1, 1, context))) {
        return NULL;
    }

    abs = GMPy_Complex_Abs(tempx, context);
    phase = GMPy_Complex_Phase(tempx, context);
    Py_DECREF(tempx);
    result = PyTuple_New(2);
    if (!abs || !phase || !result) {
        Py_XDECREF(abs);
        Py_XDECREF(phase);
        Py_XDECREF(result);
        return NULL;
    }
    PyTuple_SET_ITEM(result, 0, abs);
    PyTuple_SET_ITEM(result, 1, phase);
    return result;
}

static PyObject *
GMPy_Number_Polar(PyObject *x, CTXT_Object *context)
{
    if (IS_COMPLEX_ONLY(x))
        return GMPy_Complex_Polar(x, context);

    TYPE_ERROR("polar() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Polar(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("polar() requires 1 argument");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Polar(PyTuple_GET_ITEM(args, 0), context);
}

PyDoc_STRVAR(GMPy_doc_context_rect,
"context.rect(r, phi) -> mpc\n\n"
"Return the rectangular coordinate form of a complex number that is\n"
"given in polar form.");

PyDoc_STRVAR(GMPy_doc_function_rect,
"rect(r, phi) -> mpc\n\n"
"Return the rectangular coordinate form of a complex number that is\n"
"given in polar form.");

static PyObject *
GMPy_Complex_Rect(PyObject *x, PyObject *y, CTXT_Object *context)
{
    MPFR_Object *tempx, *tempy;
    MPC_Object *result;

    CHECK_CONTEXT(context);

    tempx = GMPy_MPFR_From_Real(x, 1, context);
    tempy = GMPy_MPFR_From_Real(y, 1, context);
    result = GMPy_MPC_New(0, 0, context);
    if (!tempx || !tempy || !result) {
        Py_XDECREF((PyObject*)tempx);
        Py_XDECREF((PyObject*)tempy);
        Py_XDECREF((PyObject*)result);
        return NULL;
    }

    mpfr_cos(mpc_realref(result->c), tempy->f, GET_REAL_ROUND(context));
    mpfr_mul(mpc_realref(result->c), mpc_realref(result->c), tempx->f, GET_REAL_ROUND(context));
    mpfr_sin(mpc_imagref(result->c), tempy->f, GET_IMAG_ROUND(context));
    mpfr_mul(mpc_imagref(result->c), mpc_imagref(result->c), tempx->f, GET_IMAG_ROUND(context));
    Py_DECREF((PyObject*)tempx);
    Py_DECREF((PyObject*)tempy);

    GMPY_MPC_CLEANUP(result, context, "rect()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_Rect(PyObject *x, PyObject *y, CTXT_Object *context)
{
    if (IS_REAL(x) && IS_REAL(y))
        return GMPy_Complex_Rect(x, y, context);

    TYPE_ERROR("rect() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Rect(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 2) {
        TYPE_ERROR("rect() requires 2 arguments");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Rect(PyTuple_GET_ITEM(args, 0), PyTuple_GET_ITEM(args, 1), context);
}

PyDoc_STRVAR(GMPy_doc_context_proj,
"context.proj(x) -> mpc\n\n"
"Returns the projection of a complex x on to the Riemann sphere.");

PyDoc_STRVAR(GMPy_doc_function_proj,
"proj(x) -> mpc\n\n"
"Returns the projection of a complex x on to the Riemann sphere.");

static PyObject *
GMPy_Complex_Proj(PyObject *x, CTXT_Object *context)
{
    MPC_Object *result, *tempx;

    CHECK_CONTEXT(context);

    result = GMPy_MPC_New(0, 0, context);
    tempx = GMPy_MPC_From_Complex(x, 1, 1, context);
    if (!result || !tempx) {
        return NULL;
    }

    result->rc = mpc_proj(result->c, tempx->c, GET_MPC_ROUND(context));
    Py_DECREF((PyObject*)tempx);

    GMPY_MPC_CLEANUP(result, context, "proj()");
    return (PyObject*)result;
}

static PyObject *
GMPy_Number_Proj(PyObject *x, CTXT_Object *context)
{
    if (IS_COMPLEX_ONLY(x))
        return GMPy_Complex_Proj(x, context);

    TYPE_ERROR("proj() argument type not supported");
    return NULL;
}

static PyObject *
GMPy_Context_Proj(PyObject *self, PyObject *args)
{
    CTXT_Object *context = NULL;

    if (PyTuple_GET_SIZE(args) != 1) {
        TYPE_ERROR("proj() requires 1 argument");
        return NULL;
    }

    if (self && CTXT_Check(self)) {
        context = (CTXT_Object*)self;
    }
    else {
        CHECK_CONTEXT(context);
    }

    return GMPy_Number_Proj(PyTuple_GET_ITEM(args, 0), context);
}
