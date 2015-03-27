/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_const.c                                                           *
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

PyDoc_STRVAR(GMPy_doc_function_const_pi,
"const_pi([precision=0]) -> number\n\n"
"Return the constant pi using the specified precision. If no\n"
"precision is specified, the default precision is used.");

PyDoc_STRVAR(GMPy_doc_context_const_pi,
"context.const_pi() -> number\n\n"
"Return the constant pi using the context's precision.");

GMPY_MPFR_CONST(Const_Pi, const_pi)
GMPY_MPFR_NOOP(Const_Pi, const_pi)

PyDoc_STRVAR(GMPy_doc_function_const_euler,
"const_euler([precision=0]) -> number\n\n"
"Return the euler constant using the specified precision. If no\n"
"precision is specified, the default precision is used.");

PyDoc_STRVAR(GMPy_doc_context_const_euler,
"context.const_euler() -> number\n\n"
"Return the euler constant using the context's precision.");

GMPY_MPFR_CONST(Const_Euler, const_euler)
GMPY_MPFR_NOOP(Const_Euler, const_euler)

PyDoc_STRVAR(GMPy_doc_function_const_log2,
"const_log2([precision=0]) -> number\n\n"
"Return the log2 constant using the specified precision. If no\n"
"precision is specified, the default precision is used.");

PyDoc_STRVAR(GMPy_doc_context_const_log2,
"context.const_log2() -> number\n\n"
"Return the log2 constant using the context's precision.");

GMPY_MPFR_CONST(Const_Log2, const_log2)
GMPY_MPFR_NOOP(Const_Log2, const_log2)

PyDoc_STRVAR(GMPy_doc_function_const_catalan,
"const_catalan([precision=0]) -> number\n\n"
"Return the catalan constant using the specified precision. If no\n"
"precision is specified, the default precision is used.");

PyDoc_STRVAR(GMPy_doc_context_const_catalan,
"context.const_catalan() -> number\n\n"
"Return the catalan constant using the context's precision.");

GMPY_MPFR_CONST(Const_Catalan, const_catalan)
GMPY_MPFR_NOOP(Const_Catalan, const_catalan)

