/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_divmod2exp.h                                                   *
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

#ifndef GMPY_MPZ_DIVMOD2EXP_H
#define GMPY_MPZ_DIVMOD2EXP_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * Pygmpy_c_divmod_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_c_div_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_c_mod_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_f_divmod_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_f_div_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_f_mod_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_t_divmod_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_t_div_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_t_mod_2exp(PyObject *self, PyObject *args);
static PyObject * Pygmpy_pack(PyObject *self, PyObject *args);
static PyObject * Pygmpy_unpack(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
