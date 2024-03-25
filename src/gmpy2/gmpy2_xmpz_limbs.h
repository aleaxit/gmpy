/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_xmpz_limbs.h                                                      *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP, MPFR, and MPC multiple precision           *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2020 Tyler Lanphear                                           *
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

#ifndef GMPY_XMPZ_LIMBS_H
#define GMPY_XMPZ_LIMBS_H

#ifdef __cplusplus
extern "C" {
#endif

static PyObject* GMPy_XMPZ_Method_NumLimbs(PyObject* obj, PyObject* other);
static PyObject* GMPy_XMPZ_Method_LimbsRead(PyObject* obj, PyObject* other);
static PyObject* GMPy_XMPZ_Method_LimbsWrite(PyObject* obj, PyObject* other);
static PyObject* GMPy_XMPZ_Method_LimbsModify(PyObject* obj, PyObject* other);
static PyObject* GMPy_XMPZ_Method_LimbsFinish(PyObject* obj, PyObject* other);

#ifdef __cplusplus
}
#endif
#endif
