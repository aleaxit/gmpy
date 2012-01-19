/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_binary.c                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *      Copyright 2000 - 2009 Alex Martelli                                *
 *      Copyright 2008, 2009, 2010, 2011, 2012 Case Van Horsen             *
 *                                                                         *
 * This library is free software; you can redistribute it and/or modify it *
 * under the terms of the GNU Lesser General Public License as published   *
 * by the Free Software Foundation; either version 2.1 of the License, or  *
 * (at your option) any later version.                                     *
 *                                                                         *
 * This library is distributed in the hope that it will be useful,         *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU       *
 * Lesser General Public License for more details.                         *
 *                                                                         *
 * You should have received a copy of the GNU Lesser General Public        *
 * License along with this library; if not, write to the Free Software     *
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA           *
 * 02110-1301  USA                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef GMPY_BINARY_H
#define GMPY_BINARY_H

#ifdef __cplusplus
extern "C" {
#endif

/* Conversion routines between GMPY2 objects and a compact, portable
 * binary representation. The binary format of GMPY2 is not compatible
 * with GMPY 1.x. Methods to read the old format are provided.
 */

static PyObject * Pympz_From_Old_Binary(PyObject *self, PyObject *other);
static PyObject * Pympq_From_Old_Binary(PyObject *self, PyObject *other);
#ifdef WITHMPFR
static PyObject * Pympfr_From_Old_Binary(PyObject *self, PyObject *other);
#endif

#ifdef __cplusplus
}
#endif
#endif
