/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_context.h                                                         *
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

#ifndef GMPY_CONTEXT_H
#define GMPY_CONTEXT_H

#ifdef __cplusplus
extern "C" {
#endif

#define TRAP_NONE      0
#define TRAP_UNDERFLOW 1
#define TRAP_OVERFLOW  2
#define TRAP_INEXACT   4
#define TRAP_INVALID   8
#define TRAP_ERANGE    16
#define TRAP_DIVZERO   32

/* The actual typedefs have been moved to gmpy2_types.h. */

static PyTypeObject CTXT_Type;
static PyTypeObject CTXT_Manager_Type;

#ifdef WITHOUT_THREADS
#define CURRENT_CONTEXT(obj) obj = module_context;
#define GMPY_MAYBE_BEGIN_ALLOW_THREADS /* NOTHING */
#define GMPY_MAYBE_ALLOW_THREADS do {} while(0)
#else
#define CURRENT_CONTEXT(obj) obj = GMPy_current_context();
#define GMPY_MAYBE_BEGIN_ALLOW_THREADS(context) { \
        PyThreadState *_save; \
        _save = GET_THREAD_MODE(context) ? PyEval_SaveThread() : NULL;
#define GMPY_MAYBE_END_ALLOW_THREADS(context) \
        if (_save) PyEval_RestoreThread(_save); \
    }
#endif

#define CHECK_CONTEXT(context) \
    if (!context) CURRENT_CONTEXT(context);

#define CTXT_Check(v) (((PyObject*)v)->ob_type == &CTXT_Type)
#define CTXT_Manager_Check(v) (((PyObject*)v)->ob_type == &CTXT_Manager_Type)

#define GET_MPFR_PREC(c) (c->ctx.mpfr_prec)
#define GET_REAL_PREC(c) ((c->ctx.real_prec==GMPY_DEFAULT)?GET_MPFR_PREC(c):c->ctx.real_prec)
#define GET_IMAG_PREC(c) ((c->ctx.imag_prec==GMPY_DEFAULT)?GET_REAL_PREC(c):c->ctx.imag_prec)
#define GET_MPFR_ROUND(c) (c->ctx.mpfr_round)
#define GET_REAL_ROUND(c) ((c->ctx.real_round==GMPY_DEFAULT)?GET_MPFR_ROUND(c):c->ctx.real_round)
#define GET_IMAG_ROUND(c) ((c->ctx.imag_round==GMPY_DEFAULT)?GET_REAL_ROUND(c):c->ctx.imag_round)
#define GET_MPC_ROUND(c) (MPC_RND(GET_REAL_ROUND(c), GET_IMAG_ROUND(c)))

#define GET_DIV_MODE(c) (c->ctx.rational_division)

#define GET_THREAD_MODE(c) (c->ctx.allow_release_gil)


static PyObject *    GMPy_CTXT_Manager_New(void);
static void          GMPy_CTXT_Manager_Dealloc(CTXT_Manager_Object *self);
static PyObject *    GMPy_CTXT_Manager_Repr_Slot(CTXT_Manager_Object *self);
static PyObject *    GMPy_CTXT_Manager_Enter(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Manager_Exit(PyObject *self, PyObject *args);

static PyObject *    GMPy_CTXT_New(void);
static void          GMPy_CTXT_Dealloc(CTXT_Object *self);
static PyObject *    GMPy_CTXT_Repr_Slot(CTXT_Object *self);
static PyObject *    GMPy_CTXT_Get(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Local(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *    GMPy_CTXT_Context(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *    GMPy_CTXT_Set(PyObject *self, PyObject *other);
static PyObject *    GMPy_CTXT_Clear_Flags(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Copy(PyObject *self, PyObject *other);
static PyObject *    GMPy_CTXT_ieee(PyObject *self, PyObject *args, PyObject *kwargs);
static PyObject *    GMPy_CTXT_Enter(PyObject *self, PyObject *args);
static PyObject *    GMPy_CTXT_Exit(PyObject *self, PyObject *args);

#ifndef WITHOUT_THREADS
static CTXT_Object * GMPy_current_context(void);
#endif

#ifdef __cplusplus
}
#endif
#endif

