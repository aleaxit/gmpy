/*
dnl  Copyright 2009 Jason Moxham

dnl  This file is part of the MPIR Library.

dnl  The MPIR Library is free software; you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as published
dnl  by the Free Software Foundation; either version 2.1 of the License, or (at
dnl  your option) any later version.

dnl  The MPIR Library is distributed in the hope that it will be useful, but
dnl  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
dnl  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
dnl  License for more details.

dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with the MPIR Library; see the file COPYING.LIB.  If not, write
dnl  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
dnl  Boston, MA 02110-1301, USA.
*/

long cpuid(char *p,int i)
{int a[4];

__cpuid(a,i);
memcpy(p,&(a[1]),4);
memcpy(p+4,&(a[3]),4);
memcpy(p+8,&(a[2]),4);
return a[0];}

#define CONFIG_GUESS		1
// safe to always try 32bit
#define CONFIG_GUESS_32BIT	1
#define CONFIG_GUESS_64BIT	0
#define FAT32			0
#define FAT64			0
#define INFAT			0
#include "cpuid.c"

main (int argc,char *argv[])
{char *p,*modelstr;int a=sizeof(p)*8;

modelstr=__gmpn_cpu(0);
if(argc==2){printf("set GCPU=%s\nset GBITS=%d\n",modelstr,a);return 0;}
printf("%s-pc-Win%d\n",modelstr,a);return 0;}
