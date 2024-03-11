/* Copyright 2011-2015 David Cleaver
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * This file has been extensively modified by Case Van Horsen to simplify use
 * with gmpy2.
 *
 * Summary of changes:
 *  - gmpy2 already includes modified copies of the PRP functions that were
 *    included in this file. They have been removed.
 *  - The 64 bit integer type is a "long".
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

#ifndef __MPZ_APRCL__
#define __MPZ_APRCL__

#ifndef HAVE_U64_T
#define HAVE_U64_T
typedef long s64_t;
typedef unsigned long u64_t;
#endif

#include "jacobi_sum.h"

/* **********************************************************************************
 * APR-CL (also known as APRT-CLE) is a prime proving algorithm developed by:
 * L. Adleman, C. Pomerance, R. Rumely, H. Cohen, and H. W. Lenstra
 * APRT-CLE = Adleman-Pomerance-Rumely Test Cohen-Lenstra Extended version
 * You can find all the details of this implementation in the Cohen & Lenstra paper:
 *    H. Cohen and A. K. Lenstra, "Implementation of a new primality test",
 *    Math. Comp., 48 (1987) 103--121
 *
 * ----------------------------------------------------------------------------------
 *
 * This C/GMP version is a conversion of Dario Alpern's Java based APRT-CLE code
 * His code was based on Yuji Kida's UBASIC code
 *
 * Based on APRT-CLE Written by Dario Alejandro Alpern (Buenos Aires - Argentina)
 * From: Last updated September 10th, 2011. See http://www.alpertron.com.ar/ECM.HTM
 *
 * On 2012/11/12 Dario Alpern has approved the conversion, from Java to C/GMP, of
 * his implementation of the APR-CL algorithm, and that it be licensed under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * With improvements based on Jason Moxham's APRCL v1.15 code, from 2003/01/01
 *
 * On 2013/04/14 Toby Moxham has approved the APR-CL code and data tables,
 * originally written by his brother Jason Moxham on 2003/01/01, to be released
 * under the LGPL.
 *
 * ----------------------------------------------------------------------------------
 *
 * v1.0 to v1.1 improvements:
 *
 * [The following fix was recommended by Dana Jacobsen and verified by Jon Grantham]
 *      - Bug fix: Removed unnecessary vl==0 check in mpz_extrastronglucas_prp
 * [The following improvements/fixes were recommended by Laurent Desnogues in 2013/08]
 *      - Speed improvement 1: Removed extraneous NormalizeJS calls in ARPCL
 *      - Speed improvement 2: Removed/consolidated calls to mpz_mod in APRCL
 *        (these improvements make the APRCL code about 1.5-2.2x faster)
 *      - Bug fix: Final test in APRCL routine is now correct
 *
 *
 * *********************************************************************************/

PyDoc_STRVAR(doc_mpz_is_aprcl_prime,
"is_aprcl_prime(n) -> boolean\n\n"
"Return True if n is a proven prime number verified by the APRCL test.\n"
"Return False if n is composite. An exception is be raised if n is too\n"
"large for this implementation of APRCL test.");

static PyObject *
GMPY_mpz_is_aprcl_prime(PyObject *self, PyObject *other)

#endif
