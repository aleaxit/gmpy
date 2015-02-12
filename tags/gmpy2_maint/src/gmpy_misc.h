/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_misc.c                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014,                     *
 *           2015 Case Van Horsen                                          *
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

#ifndef GMPY_MISC_H
#define GMPY_MISC_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * Pygmpy_get_license(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_version(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_cvsid(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_mp_version(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_mpfr_version(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_mpc_version(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_mp_limbsize(PyObject *self, PyObject *args);
static PyObject * Pygmpy_get_cache(PyObject *self, PyObject *args);
static PyObject * Pygmpy_set_cache(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
