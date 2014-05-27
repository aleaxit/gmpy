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
 * Parses two, and only two, arguments into "self" and "var" and converts
 * them both to mpz. Is faster, but not as generic, as using PyArg_ParseTuple.
 * It supports either gmpy.fname(z,z) or z.fname(z). "self" & "var" must be
 * decref'ed after use. "msg" should be an error message that includes the
 * function name and describes the required arguments. Replaces
 * SELF_MPZ_ONE_ARG_CONVERTED(var).
 *
 * Also considers an 'xmpz' to be equivalent to an 'mpz'.
 */

#define PARSE_xxxTWO_MPZ(var, msg) \
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
