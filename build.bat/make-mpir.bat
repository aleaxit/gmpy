@echo off
::  Copyright 2011 Jason Moxham
::
::  This file is part of the MPIR Library.
::
::  The MPIR Library is free software; you can redistribute it and/or modify
::  it under the terms of the GNU Lesser General Public License as published
::  by the Free Software Foundation; either version 2.1 of the License, or (at
::  your option) any later version.
::
::  The MPIR Library is distributed in the hope that it will be useful, but
::  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
::  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
::  License for more details.
::
::  You should have received a copy of the GNU Lesser General Public License
::  along with the MPIR Library; see the file COPYING.LIB.  If not, write
::  to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
::  Boston, MA 02110-1301, USA.

::  Modified by Case Van Horsen to optimize build for gmpy2.


if "%1" == "clean" goto :clean

if not exist config.params.bat (
	echo run configure first
	exit /b 1
)
call config.params.bat
if "%1" == "" goto :make
if "%1" == "check" goto :check
if "%1" == "install" goto :install

echo Usage: make [clean^|check^|install]
exit /b 1

:make
cd ..\mpir\win
md mpn fft mpz mpq mpf printf scanf cxx > nul 2>&1

copy ..\build.vc10\gen_mpir_h.bat .
copy ..\build.vc10\out_copy_rename.bat .
copy ..\build.vc10\gen_config_h.bat .
copy ..\build.vc10\cfg.h .

if %ABI% == 64 (set LOCALABI=x64)
if %ABI% == 32 (set LOCALABI=win32)
set YASMFLAG=
if %LIBTYPE% == dll (set YASMFLAG=-D DLL)
if %LIBTYPE% == dll (set FLAGS1=%FLAGS1% /D "__GMP_WITHIN_GMP" /D "__GMP_WITHIN_GMPXX")
call gen_mpir_h %LOCALABI%
copy ..\mpn\generic\gmp-mparam.h .. > nul 2>&1
copy ..\longlong_pre.h+..\mpn\generic\longlong_inc.h+..\longlong_post.h ..\longlong.h > nul 2>&1
for %%X in ( %MPNPATH% ) do (
	copy ..\mpn\%%X\gmp-mparam.h .. > nul 2>&1
	copy ..\longlong_pre.h+..\mpn\%%X\longlong_inc.h+..\longlong_post.h ..\longlong.h > nul 2>&1
	call gen_config_h ..\mpn\%%X\ > nul 2>&1
)

set HAVE_STDINT=yes
echo #include ^<stdint.h^> > comptest.c
echo int main(void){return 0;} >> comptest.c
cl /nologo comptest.c > nul 2>&1
if errorlevel 1 (
	echo #undef HAVE_STDINT_H >> ..\config.h
	echo #undef HAVE_INTMAX_T >> ..\config.h
	echo #undef HAVE_UINTMAX_T >> ..\config.h
	echo #undef HAVE_PTRDIFF_T >> ..\config.h
	echo #undef HAVE_UINT_LEAST32_T >> ..\config.h
	echo #undef SIZEOF_UINTMAX_T >> ..\config.h
	set HAVE_STDINT=no
)
del comptest.*

if %ABI% == 64 (set OPT=/D "_WIN64" %FLAGS% %FLAGS1% /c)
if %ABI% == 32 (set OPT=%FLAGS% %FLAGS1% /c)

if %ABI% == 64 (set LOCALDIR=x86_64w)
if %ABI% == 32 (set LOCALDIR=x86w)
:: just compile all generic and just overwrite with asm
cd mpn
for %%X in ( ..\..\mpn\generic\*.c) do (
	:: exclude udiv_w_sdiv.c from all builds
	if not "%%X" == "..\..\mpn\generic\udiv_w_sdiv.c" (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
	)
)
for %%X in ( %MPNPATH% ) do (
	for %%i in ( ..\..\mpn\%%X\*.asm ) do (
		%YASMEXE% %YASMFLAG% -I ..\..\mpn\%LOCALDIR% -f %LOCALABI% %%i
		if errorlevel 1 goto :err
		echo assembling %%i
	)
)
:: dont know what the asm version have so delete them
del preinv_divrem_1.obj preinv_mod_1.obj divrem_1.obj mod_1.obj divrem_euclidean_qr_1.obj > nul 2>&1
cl %OPT% -I..\.. ..\..\mpn\generic\divrem_1.c
if errorlevel 1 goto :err
cl %OPT% -I..\.. ..\..\mpn\generic\mod_1.c
if errorlevel 1 goto :err
cl %OPT% -I..\.. ..\..\mpn\generic\divrem_euclidean_qr_1.c
if errorlevel 1 goto :err
cl /D "USE_PREINV_DIVREM_1" %OPT% -I..\.. ..\..\mpn\generic\preinv_divrem_1.c
if errorlevel 1 goto :err
cl %OPT% -I..\.. ..\..\mpn\generic\preinv_mod_1.c
if errorlevel 1 goto :err
cd ..

cd mpz
for %%X in ( ..\..\mpz\*.c) do (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
if %HAVE_STDINT% == no (
	del *_sx.obj *_ux.obj > nul 2>&1
)
cd ..

cd fft
for %%X in ( ..\..\fft\*.c) do (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
cd ..

cd mpf
for %%X in ( ..\..\mpf\*.c) do (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
cd ..

cd mpq
for %%X in ( ..\..\mpq\*.c) do (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
cd ..

cd printf
for %%X in ( ..\..\printf\*.c) do (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
cd ..

cd scanf
for %%X in ( ..\..\scanf\*.c) do (
	cl %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
cd ..

for %%X in ( ..\assert.c ..\errno.c ..\extract-dbl.c ..\invalid.c ..\memory.c ..\mp_bpl.c ..\mp_clz_tab.c ..\mp_dv_tab.c ..\mp_get_fns.c ..\mp_minv_tab.c ..\mp_set_fns.c ..\randbui.c ..\randclr.c ..\randdef.c ..\randiset.c ..\randlc2s.c ..\randlc2x.c ..\randmt.c ..\randmts.c ..\randmui.c ..\rands.c ..\randsd.c ..\randsdui.c ..\tal-reent.c ..\version.c ) do (
	cl %OPT% -I.. %%X
	if errorlevel 1 goto :err
)

cd cxx
for %%X in ( ..\..\cxx\*.cc) do (
	cl /EHsc %OPT% -I..\.. %%X
	if errorlevel 1 goto :err
)
cd ..

if %LIBTYPE% == dll (
	link /DLL /NODEFAULTLIB:LIBCMT.lib /nologo scanf\*.obj printf\*.obj mpz\*.obj fft\*.obj mpq\*.obj mpf\*.obj mpn\*.obj *.obj /out:mpir.%LIBTYPE%
	if errorlevel 1 goto :err
	link /DLL /NODEFAULTLIB:LIBCMT.lib /nologo scanf\*.obj printf\*.obj mpz\*.obj fft\*.obj mpq\*.obj mpf\*.obj mpn\*.obj *.obj /out:mpirxx.%LIBTYPE%
	if errorlevel 1 goto :err
)
if %LIBTYPE% == lib (
	lib /nologo scanf\*.obj printf\*.obj mpz\*.obj fft\*.obj mpq\*.obj mpf\*.obj mpn\*.obj *.obj /out:mpir.%LIBTYPE%
	if errorlevel 1 goto :err
	lib /nologo cxx/*.obj /out:mpirxx.%LIBTYPE%
	if errorlevel 1 goto :err
)
cd ..\..\batch
exit /b 0


:check
cd ..\mpir\win

md tests
cd tests
md mpn fft mpz mpq mpf rand misc cxx devel > nul 2>&1
if %LIBTYPE% == dll (
	for %%X in ( mpn fft mpz mpq mpf rand misc cxx devel .) do (
		copy ..\mpir.dll %%X
	)
	copy ..\mpirxx.dll cxx
)
cd ..

set OPT=%FLAGS% %FLAGS1%

set MPIRLIB=/link ..\..\mpir.lib
if %LIBTYPE% == lib (set MPIRLIB=..\..\mpir.lib)
if %LIBTYPE% == lib (set MPIRXXLIB=..\..\mpirxx.lib)

set MPIRLIB1=/link ..\mpir.lib
if %LIBTYPE% == lib (set MPIRLIB1=..\mpir.lib)

cd tests

cl %OPT% /c ..\..\tests\memory.c /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\tests\misc.c   /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\tests\trace.c  /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\tests\refmpn.c /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\tests\refmpz.c /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\tests\refmpq.c /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\tests\refmpf.c /I..\..
if errorlevel 1 goto :err
:: these only needed for try.exe
cl %OPT% /c ..\..\tests\spinner.c /I..\..
if errorlevel 1 goto :err
cl %OPT% /c ..\..\build.vc10\getopt.c /I..\..
if errorlevel 1 goto :err

for %%X in ( ..\..\tests\t-*.c) do (
	cl %OPT% /I..\.. /I..\..\tests %%X misc.obj memory.obj trace.obj refmpn.obj %MPIRLIB1%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing %%X
	%%X
	if errorlevel 1 goto :err
)

cd mpn
for %%X in ( ..\..\..\tests\mpn\t-*.c) do (
	if not %%X == ..\..\..\tests\mpn\t-fat.c (
		cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj %MPIRLIB%
		if errorlevel 1 goto :err
	)
)
for %%X in ( *.exe) do (
	echo testing mpn_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd fft
for %%X in ( ..\..\..\tests\fft\t-*.c) do (
	if not %%X == ..\..\..\tests\fft\t-fat.c (
		cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj %MPIRLIB%
		if errorlevel 1 goto :err
	)
)
for %%X in ( *.exe) do (
	echo testing fft_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd mpz
for %%X in ( ..\..\..\tests\mpz\*.c) do (
	cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj ..\refmpz.obj %MPIRLIB%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing mpz_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd mpq
for %%X in ( ..\..\..\tests\mpq\t-*.c) do (
	cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj ..\refmpz.obj ..\refmpq.obj %MPIRLIB%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing mpq_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd mpf
for %%X in ( ..\..\..\tests\mpf\*.c) do (
	cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj ..\refmpf.obj %MPIRLIB%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing mpf_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd rand
for %%X in ( ..\..\..\tests\rand\t-*.c) do (
	cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj %MPIRLIB%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing rand_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd misc
for %%X in ( ..\..\..\tests\misc\t-*.c) do (
	cl %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj %MPIRLIB%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing misc_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd cxx
for %%X in ( ..\..\..\tests\cxx\t-*.cc) do (
	cl /EHsc %OPT% /I..\..\.. /I..\..\..\tests %%X ..\misc.obj ..\memory.obj ..\trace.obj ..\refmpn.obj %MPIRLIB% %MPIRXXLIB%
	if errorlevel 1 goto :err
)
for %%X in ( *.exe) do (
	echo testing cxx_%%X
	%%X
	if errorlevel 1 goto :err
)
cd ..

cd devel
cl %OPT% ..\..\..\tests\devel\try.c /I..\..\..\ /I..\..\..\tests ..\refmpn.obj ..\refmpz.obj ..\trace.obj ..\spinner.obj ..\misc.obj ..\memory.obj ..\getopt.obj %MPIRLIB%
if errorlevel 1 goto :err
cd ..

cd ..
cd ..\..\batch
exit /b 0


:err
cd C:\src\batch
echo ERROR
exit /b 1

:clean
cd ..\mpir\win

del *.obj ..\mpir.h ..\config.h ..\gmp-mparam.h mpir.lib *.idb *.pdb > nul 2>&1
rmdir /S/Q mpn fft mpz mpq mpf scanf printf tests cxx tune speed > nul 2>&1
del gen_mpir_h.bat out_copy_rename.bat gen_config_h.bat cfg.h > nul 2>&1
del getopt.h getrusage.h gettimeofday.h unistd.h win_timing.h > nul 2>&1
del mpir.dll mpir.dll.manifest mpir.exp > nul 2>&1
del mpirxx.lib mpirxx.dll mpirxx.exp mpirxx.dll.manifest > nul 2>&1
cd ..\..\batch
del config.guess.bat config.guess.exe config.guess.obj config.params.bat > nul 2>&1
exit /b 0

:install
mkdir ..\lib > nul 2>&1
mkdir ..\include > nul 2>&1
xcopy /y ..\mpir\mpir.h ..\include > nul 2>&1
xcopy /y ..\mpir\gmp.h ..\include > nul 2>&1
xcopy /y ..\mpir\gmp-mparam.h ..\include > nul 2>&1
xcopy /y ..\mpir\gmp-impl.h ..\include > nul 2>&1
xcopy /y ..\mpir\win\mpir.lib ..\lib > nul 2>&1
exit /b 0
