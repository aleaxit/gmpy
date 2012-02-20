/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_mpz_prp.h                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2011 David Cleaver                                            *
 *                                                                         *
 * Copyright 2012 Case Van Horsen                                          *
 *                                                                         *
 * The original file is available at:                                      *
 *   <http://sourceforge.net/projects/mpzprp/files/>                       *
 *                                                                         *
 * Modified by Case Van Horsen for inclusion into GMPY2.                   *
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

#ifndef GMPY_PRP_H
#define GMPY_PRP_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject * GMPY_mpz_is_fermat_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_euler_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_strong_prp(PyObject *self, PyObject *args);
static PyObject * GMPY_mpz_is_fibonacci_prp(PyObject *self, PyObject *args);

#ifdef __cplusplus
}
#endif
#endif
