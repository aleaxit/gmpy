/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * gmpy2_misc.c                                                            *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Python interface to the GMP or MPIR, MPFR, and MPC multiple precision   *
 * libraries.                                                              *
 *                                                                         *
 * Copyright 2000 - 2009 Alex Martelli                                     *
 *                                                                         *
 * Copyright 2008 - 2021 Case Van Horsen                                   *
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

/* Miscellaneous module-level functions and helper functions. */

PyDoc_STRVAR(GMPy_doc_license,
"license() -> string\n\n"
"Return string giving license information.");

static PyObject *
GMPy_get_license(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_license);
}

PyDoc_STRVAR(GMPy_doc_version,
"version() -> string\n\n"
"Return string giving current GMPY2 version.");

static PyObject *
GMPy_get_version(PyObject *self, PyObject *args)
{
    return Py_BuildValue("s", gmpy_version);
}

PyDoc_STRVAR(GMPy_doc_mp_version,
"mp_version() -> string\n\n"
"Return string giving the name and version of the multiple precision\n"
"library used.");

static PyObject *
GMPy_get_mp_version(PyObject *self, PyObject *args)
{
#ifndef __MPIR_VERSION
    return Py2or3String_FromFormat("%s %s", "GMP", gmp_version);
#else
    return Py2or3String_FromFormat("%s %s", "MPIR", mpir_version);
#endif
}

PyDoc_STRVAR(GMPy_doc_mpfr_version,
"mpfr_version() -> string\n\n"
"Return string giving current MPFR version.");

static PyObject *
GMPy_get_mpfr_version(PyObject *self, PyObject *args)
{
    return Py2or3String_FromFormat("%s %s", "MPFR", MPFR_VERSION_STRING);
}

PyDoc_STRVAR(GMPy_doc_mpc_version,
"mpc_version() -> string\n\n"
"Return string giving current MPC version.");

static PyObject *
GMPy_get_mpc_version(PyObject *self, PyObject *args)
{
    return Py2or3String_FromFormat("%s %s", "MPC", MPC_VERSION_STRING);
}

PyDoc_STRVAR(GMPy_doc_mp_limbsize,
"mp_limbsize() -> integer\n\n\
Return the number of bits per limb.");

static PyObject *
GMPy_get_mp_limbsize(PyObject *self, PyObject *args)
{
    return Py_BuildValue("i", mp_bits_per_limb);
}

/*
 * access cache options
 */

PyDoc_STRVAR(GMPy_doc_get_cache,
"get_cache() -> (cache_size, object_size)\n\n\
Return the current cache size (number of objects) and maximum size\n\
per object (number of limbs) for all GMPY2 objects.");

static PyObject *
GMPy_get_cache(PyObject *self, PyObject *args)
{
    return Py_BuildValue("(ii)", global.cache_size, global.cache_obsize);
}

PyDoc_STRVAR(GMPy_doc_set_cache,
"set_cache(cache_size, object_size)\n\n\
Set the current cache size (number of objects) and the maximum size\n\
per object (number of limbs). Raises ValueError if cache size exceeds\n\
1000 or object size exceeds 16384.");

static PyObject *
GMPy_set_cache(PyObject *self, PyObject *args)
{
    int newcache = -1, newsize = -1;

    if (!PyArg_ParseTuple(args, "ii", &newcache, &newsize))
        return NULL;
    if (newcache<0 || newcache>MAX_CACHE) {
        VALUE_ERROR("cache size must between 0 and 1000");
        return NULL;
    }
    if (newsize<0 || newsize>MAX_CACHE_LIMBS) {
        VALUE_ERROR("object size must between 0 and 16384");
        return NULL;
    }

    global.cache_size = newcache;
    global.cache_obsize = newsize;
    set_gmpympzcache();
    set_gmpympqcache();
    set_gmpyxmpzcache();
    set_gmpympfrcache();
    set_gmpympccache();
    Py_RETURN_NONE;
}

PyDoc_STRVAR(GMPy_doc_function_printf,
"_printf(fmt, x) -> string\n\n"
"Return a Python string by formatting 'x' using the format string\n"
"'fmt'.\n\n"
"WARNING: Invalid format strings will cause a crash. Please see the\n"
"         GMP and MPFR manuals for details on the format code. 'mpc'\n"
"         objects are not supported.");

static PyObject *
GMPy_printf(PyObject *self, PyObject *args)
{
    PyObject *result = NULL, *x = NULL;
    char *buffer = NULL, *fmtcode = NULL;
    void *generic;
    int buflen;

    if (!PyArg_ParseTuple(args, "sO", &fmtcode, &x))
        return NULL;

    if (CHECK_MPZANY(x) || MPQ_Check(x)) {
        if (CHECK_MPZANY(x))
            generic = MPZ(x);
        else
            generic = MPQ(x);
        buflen = mpfr_asprintf(&buffer, fmtcode, generic);
        if (buflen < 0) {
            VALUE_ERROR("_printf() could not format the 'mpz' or 'mpq' object");
        }
        else {
            result = Py_BuildValue("s", buffer);
            mpfr_free_str(buffer);
        }
        return result;
    }
    else if(MPFR_Check(x)) {
        generic = MPFR(x);
        buflen = mpfr_asprintf(&buffer, fmtcode, generic);
        if (buflen < 0) {
            VALUE_ERROR("_printf() could not format the 'mpfr' object");
        }
        else {
            result = Py_BuildValue("s", buffer);
            mpfr_free_str(buffer);
        }
        return result;
    }
    else if(MPC_Check(x)) {
        TYPE_ERROR("_printf() does not support 'mpc'");
        return NULL;
    }
    else {
        TYPE_ERROR("_printf() argument type not supported");
        return NULL;
    }
}

