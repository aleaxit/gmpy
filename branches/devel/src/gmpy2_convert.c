/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy_convert.c                                                          *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,               *
 *           2008, 2009 Alex Martelli                                      *
 *                                                                         *
 * Copyright 2008, 2009, 2010, 2011, 2012, 2013 Case Van Horsen            *
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

/* This file contains all the conversion functions for gmpy2.
 *
 * Overview
 * --------
 * gmpy2 tries to optimize the performance and accuracy of conversions from
 * other numeric types. gmpy2 uses a LBYL (Look Before You Leap) approach and
 * identifies the numeric type before conversion before conversion to a gmpy2
 * type. The basic operations (+, -, *, /) are optimized to directly work with
 * some basic types such as C longs or doubles.
 *
 * Support for the Decimal type is a challenge. For the basic operations, it
 * is most accurate to convert a Decimal instance into an mpq and then use
 * MPFR's functions to accurately operate on an mpfr and mpq. This approach is
 * challenging because (1) a large exponent can create a very large mpq and
 * (2) the changes made to C-coded version of Decimal in Python 3.3.
 *
 */

static int isInteger(PyObject* obj)
{
    if (MPZ_Check(obj))         return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (XMPZ_Check(obj))      return 1;

    return 0;
}

static int isFraction(PyObject* obj)
{
    if (!strcmp(Py_TYPE(obj)->tp_name, "Fraction")) return 1;

    return 0;
}

static int isDecimal(PyObject* obj)
{
#if PY_VERSION_HEX < 0x03030000
    if (!strcmp(Py_TYPE(obj)->tp_name, "Decimal")) return 1;
#else
    if (!strcmp(Py_TYPE(obj)->tp_name, "decimal.Decimal")) return 1;
#endif

    return 0;
}

static int isRational(PyObject* obj)
{
    if (MPZ_Check(obj))         return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (MPQ_Check(obj))         return 1;
    if (XMPZ_Check(obj))        return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

static int isReal(PyObject* obj)
{
    if (MPZ_Check(obj))         return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (MPQ_Check(obj))         return 1;
    if (MPFR_Check(obj))        return 1;
    if (XMPZ_Check(obj))        return 1;
    if (PyFloat_Check(obj))     return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}

static int isComplex(PyObject* obj)
{
    if (MPZ_Check(obj))         return 1;
    if (PyIntOrLong_Check(obj)) return 1;
    if (MPQ_Check(obj))         return 1;
    if (MPFR_Check(obj))        return 1;
    if (XMPZ_Check(obj))        return 1;
    if (MPC_Check(obj))         return 1;
    if (PyFloat_Check(obj))     return 1;
    if (PyComplex_Check(obj))   return 1;
    if (isDecimal(obj))         return 1;
    if (isFraction(obj))        return 1;

    return 0;
}



