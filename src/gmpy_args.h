/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_args.h                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
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

/* Various macros for parsing arguments. */

#ifndef GMPY_ARG_H
#define GMPY_ARG_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Create two 'mpz' and a 2-tuple.
 */

#define CREATE_TWO_MPZ_TUPLE(q, r, t) \
    q = GMPy_MPZ_New(context); \
    r = GMPy_MPZ_New(context); \
    t = PyTuple_New(2); \
    if (!q || !r || !t) { \
        Py_XDECREF(t); \
        Py_XDECREF((PyObject*)q); \
        Py_XDECREF((PyObject*)r); \
        return NULL; \
    }

/*
 * Parses one, and only one, argument into "self" and converts it to an
 * mpz. Is faster, but not as generic, as using PyArg_ParseTuple. It
 * supports either gmpy.fname(z) or z.fname(). "self" must be decref'ed.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. Replaces SELF_MPZ_NO_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ(msg) \
    if (self && CHECK_MPZANY(self)) { \
        if (PyTuple_GET_SIZE(args) != 0) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        Py_INCREF(self); \
    } else { \
        if (PyTuple_GET_SIZE(args) != 1) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        self = PyTuple_GET_ITEM(args, 0); \
        self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        if (!self) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
    }

/*
 * Parses one argument into "self" and an optional second argument into
 * "var". The second argument is converted into a C long. If there is not a
 * second argument, "var" is unchanged. Is faster, but not as generic, as
 * using PyArg_ParseTuple with "|l". It supports either gmpy.fname(z,l) or
 * z.fname(l). "self" must be decref'ed. "var" must be a pointer to a long.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. Replaces some uses of SELF_MPZ_ONE_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_OPT_CLONG(var, msg) \
    if (self && CHECK_MPZANY(self)) { \
        if (PyTuple_GET_SIZE(args) == 1) { \
            var = clong_From_Integer(PyTuple_GET_ITEM(args, 0)); \
            if (var == -1 && PyErr_Occurred()) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
        } \
        else if (PyTuple_GET_SIZE(args) > 1) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        Py_INCREF(self); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) == 2) { \
            var = clong_From_Integer(PyTuple_GET_ITEM(args, 1)); \
            if (var == -1 && PyErr_Occurred()) { \
                TYPE_ERROR(msg); \
                return NULL; \
            } \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        } \
        else if (PyTuple_GET_SIZE(args) == 1) { \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        } \
        else { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
        if (!self) { \
            TYPE_ERROR(msg); \
            return NULL; \
        } \
    }

/* Parses one argument into "self" and an optional second argument into
 * "var". The second argument is converted into an mpir_si. If there is not a
 * second argument, "var" is unchanged. It supports either gmpy.fname(z,l) or
 * z.fname(l). "self" must be decref'ed. "var" must be a pointer to an mpir_si.
 * "msg" should be an error message that includes the function name and
 * describes the required arguments. Replaces some uses of SELF_MPZ_ONE_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_OPT_SI(var, msg) \
    if (self && CHECK_MPZANY(self)) { \
        if (PyTuple_GET_SIZE(args) == 1) { \
            *var = SI_From_Integer(PyTuple_GET_ITEM(args, 0)); \
            if (*var == -1 && PyErr_Occurred()) { \
                PyErr_SetString(PyExc_TypeError, msg); \
                return NULL; \
            } \
        } \
        else if (PyTuple_GET_SIZE(args) > 1) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        Py_INCREF(self); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) == 2) { \
            *var = SI_From_Integer(PyTuple_GET_ITEM(args, 1)); \
            if (*var == -1 && PyErr_Occurred()) { \
                PyErr_SetString(PyExc_TypeError, msg); \
                return NULL; \
            } \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        } \
        else if (PyTuple_GET_SIZE(args) == 1) { \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        } \
        else { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        if (!self) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
    }

/*
 * Parses one argument into "self" and an optional second argument into "var".
 * The second argument is converted into a Py_ssize_t. If there is not a
 * second argument, "var" is unchanged. Is faster, but not as generic, as
 * using PyArg_ParseTuple with "|l". It supports either gmpy.fname(z,l) or
 * z.fname(l). "self" must be decref'ed. "var" must be a pointer to a
 * Py_ssize_t. "msg" should be an error message that includes the function
 * name and describes the required arguments. Replaces some uses of
 * SELF_MPZ_ONE_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_OPT_SSIZE_T(var, msg) \
    if (self && CHECK_MPZANY(self)) { \
        if (PyTuple_GET_SIZE(args) == 1) { \
            *var = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 0)); \
            if(*var == -1 && PyErr_Occurred()) { \
                PyErr_SetString(PyExc_TypeError, msg); \
                return NULL; \
            } \
        } \
        else if (PyTuple_GET_SIZE(args) > 1) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        Py_INCREF(self); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) == 2) { \
            *var = ssize_t_From_Integer(PyTuple_GET_ITEM(args, 1)); \
            if (*var == -1 && PyErr_Occurred()) { \
                PyErr_SetString(PyExc_TypeError, msg); \
                return NULL; \
            } \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        } \
        else if (PyTuple_GET_SIZE(args) == 1) { \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        } \
        else { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        if (!self) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
    }

/*
 * Parses one argument into "self" and a required second argument into
 * "var". The second argument is converted into a C long. Is faster, but not
 * as generic, as using PyArg_ParseTuple with "l". It supports either
 * gmpy.fname(z,l) or z.fname(l). "self" must be decref'ed. "var" must be a
 * pointer to a long. "msg" should be an error message that includes the
 * function name and describes the required arguments. Replaces some uses of
 * SELF_MPZ_ONE_ARG.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_REQ_CLONG(var, msg) \
    if (self && CHECK_MPZANY(self)) { \
        if (PyTuple_GET_SIZE(args) != 1) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        else { \
            *var = clong_From_Integer(PyTuple_GET_ITEM(args, 0)); \
            if (*var == -1 && PyErr_Occurred()) { \
                PyErr_SetString(PyExc_TypeError, msg); \
                return NULL; \
            } \
        } \
        Py_INCREF(self); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) != 2) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        else {\
            *var = clong_From_Integer(PyTuple_GET_ITEM(args, 1)); \
            if (*var == -1 && PyErr_Occurred()) { \
                PyErr_SetString(PyExc_TypeError, msg); \
                return NULL; \
            } \
            self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0)); \
        } \
        if (!self) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
    }

/* Parses one argument into "RES" and a required second argument into
 * "VAR". RES is assumed to be an MPZ_Object. The second argument is converted
 * into an mpir_si. Only use in a  function (gmpy.fname(z,l)) is supported. If
 * required, the context is in "CTX". "MSG" should be an error message that
 * includes the function name and describes the required arguments.
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_ONE_MPZ_REQ_SI_FUNCTION(RES, VAR, CTX, MSG) \
    if (PyTuple_GET_SIZE(args) != 2) { \
        TYPE_ERROR(MSG); \
        return NULL; \
    } \
    else { \
        VAR = SI_From_Integer(PyTuple_GET_ITEM(args, 1)); \
        if (VAR == -1 && PyErr_Occurred()) { \
            return NULL; \
        } \
        if (!(RES = GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), CTX))) { \
            return NULL; \
        } \
    } \

/*
 * Parses two, and only two, arguments into "self" and "var" and converts
 * them both to mpz. Is faster, but not as generic, as using PyArg_ParseTuple.
 * It supports either gmpy.fname(z,z) or z.fname(z). "self" & "var" must be
 * decref'ed after use. "msg" should be an error message that includes the
 * function name and describes the required arguments. Replaces
 * SELF_MPZ_ONE_ARG_CONVERTED(var).
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_TWO_MPZ(var, msg) \
    if (self && CHECK_MPZANY(self)) { \
        if (PyTuple_GET_SIZE(args) != 1) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        var = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        if (!var) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        Py_INCREF(self); \
    } \
    else { \
        if (PyTuple_GET_SIZE(args) != 2) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            return NULL; \
        } \
        self = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 0), context); \
        var = (PyObject*)GMPy_MPZ_From_Integer(PyTuple_GET_ITEM(args, 1), context); \
        if (!self || !var) { \
            PyErr_SetString(PyExc_TypeError, msg); \
            Py_XDECREF((PyObject*)self); \
            Py_XDECREF((PyObject*)var); \
            return NULL; \
        } \
    }

#ifdef __cplusplus
}
#endif
#endif /* !defined(Py_GMPYMODULE_H */
