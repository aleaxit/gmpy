@echo off
::  Copyright 2009,2011 Jason Moxham
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
::
::  * * * * * * * NOTE * * * * * * *
::
:: The following directory structure is assumed:
:: C:\src\mpir
::       \mpfr
::       \mpc
::       \gmpy2
::       \batch
::
:: The contents of gmpy2\build.bat should be copied to C:\src\batch.
::
:: The batch files should be run from the "C:\src\batch" directory.
:: For example, C:\src\batch>configure2.bat

:: set default lib type as static
set LIBTYPE=lib
:: set default abi to ?
set ABI=?
:: set default cpu to ?
set CPU=?
:: set default C compiler to ?
set CC=?
:: set defaults for non-debug build
set DFLAGS=/Ox /Ot /Oi /D "NDEBUG"
set DFLAGS1=/MT

:: parse params
:lp
shift
if "%0" == "" goto :exitlp
if "%0" == "--enable-shared" ( set LIBTYPE=dll)
if "%0" == "--enable-static" ( set LIBTYPE=lib)
if "%0" == "--disable-shared" ( set LIBTYPE=lib)
if "%0" == "--disable-static" ( set LIBTYPE=dll)
if "%0" == "--debug" (
   shift
   set DFLAGS=/Z7 /D "DEBUG"
   set DFLAGS1=/MTd
)
if "%0" == "ABI" (
	shift
	set ABI=%1
)
if "%0" == "--cpu" (
	shift
	set CPU=%1
)
if "%0" == "CPU" (
	shift
	set CPU=%1
)
if "%0" == "CC" (
    shift
    set CC=%1
)
goto :lp
:exitlp

:: ARCH is native ABI
set ARCH=64
if %PROCESSOR_ARCHITECTURE% == x86 ( set ARCH=32)
set VCTARGET=
if %ARCH% == 64 (
	if %ABI% == 64 (set VCTARGET=amd64)
	if %ABI% == 32 (set VCTARGET=x86)
)
if %ARCH% == 32 (
	if %ABI% == 64 (set VCTARGET=x86_amd64)
	if %ABI% == 32 (set VCTARGET=x86)
)

:: If a compiler is already present on the path, CC is ignored. If no compiler
:: present, then CC is checked to see if it requests a specific compiler. If
:: CC is empty, then look for VS 2010 or VS 2008, in that order.

:: Check if a compiler is already present....
:: Copy config.guess.c and cpuid.c to batch directory and modify config.guess.c to
:: change hard-coded include of cpuid.c.
copy /y ..\mpir\cpuid.c . > nul 2>&1
cl /nologo config.guess.c > nul 2>&1
if errorlevel 1 goto :findcc

config.guess.exe print > config.guess.bat
call config.guess.bat
if %ABI% == ? goto :gotcc
if %GBITS% == %ABI% goto :gotcc

:findcc
:: For simplicity, we assume if no ABI is specified, we want to use the platform architecture.
:: This may cause an issue if an x86 compiler is the only compiler present on an x64 host.
:: In that case, just specify the desired ABI.

if %ABI% == ? (set ABI=%ARCH%)

if %CC% == SDK70 (
    set YASMPATH=C:\Program Files (x86^)\Microsoft Visual Studio 9.0\VC\bin
    set MSSdk=1
    set DISTUTILS_USE_SDK=1
    if %ARCH% == 64 (
        if %ABI% == 32 (
            if exist "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars32.bat" (
                call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars32.bat"
                goto :checkcc
            )
        )
        if %ABI% == 64 (
            if exist "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars64.bat" (
                call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars64.bat"
                goto :checkcc
            )
        )
    )
    if %ARCH% == 32 (
        if %ABI% == 32 (
            if exist "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars32.bat" (
                call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvars32.bat"
                goto :checkcc
            )
        )
        if %ABI% == 64 (
            if exist "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvarx86_amd64.bat" (
                call "C:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\bin\vcvarx86_amd64.bat"
                goto :checkcc
            )
        )
    )
)

if %CC% == SDK71 (
    set YASMPATH=C:\Program Files (x86^)\Microsoft Visual Studio 10.0\VC\bin
    set MSSdk=1
    set DISTUTILS_USE_SDK=1
    if %ABI% == 32 (
        if exist "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" (
            :: Fix cl.exe fails with mspdb100.dll error.
            if not exist "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\mspdb100.dll" (
                xcopy /Y "C:\Program Files (x86)\Microsoft Visual Studio 10.0\Common7\IDE\ms*100.dll" "C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin"
            )
	        call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /Release /x86
            color
	        goto :checkcc
        )
    )
    if %ABI% == 64 (
        if exist "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" (
	        call "C:\Program Files\Microsoft SDKs\Windows\v7.1\Bin\SetEnv.cmd" /Release /x64
            color
	        goto :checkcc
        )
    )
)

if %CC% == VS2008 (
    if exist "%VS90COMNTOOLS%\..\..\VC\vcvarsall.bat" (
	    call "%VS90COMNTOOLS%\..\..\VC\vcvarsall.bat" %VCTARGET%
	    goto :checkcc
    )
)

if %CC% == VS2010 (
    if exist "%VS100COMNTOOLS%\..\..\VC\vcvarsall.bat" (
	    call "%VS100COMNTOOLS%\..\..\VC\vcvarsall.bat" %VCTARGET%
	    goto :checkcc
    )
)

:: No compiler was specified so just try to find one.
if exist "%VS100COMNTOOLS%\..\..\VC\vcvarsall.bat" (
	call "%VS100COMNTOOLS%\..\..\VC\vcvarsall.bat" %VCTARGET%
	goto :checkcc
)

if exist "%VS90COMNTOOLS%\..\..\VC\vcvarsall.bat" (
	call "%VS90COMNTOOLS%\..\..\VC\vcvarsall.bat" %VCTARGET%
	goto :checkcc
)

:nocc
echo Could not find a compiler.
exit /b 1

:checkcc
cl /nologo config.guess.c > nul 2>&1
if errorlevel 1 goto :nocc
config.guess.exe print > config.guess.bat
call config.guess.bat
if %ABI% == ? goto :gotcc
if %GBITS% == %ABI% goto :gotcc
goto :nocc

:gotcc
if %ABI% == ? ( set ABI=%GBITS%)
if %CPU% == ? ( set CPU=%GCPU%)

:: now find yasm , not needed for C build?
if %CPU% == none goto :gotyasm
set YASMEXE=yasm.exe
yasm.exe --version > nul 2>&1
if not errorlevel 1 goto :gotyasm
set YASMEXE=vsyasm.exe
vsyasm.exe --version > nul 2>&1
if not errorlevel 1 goto :gotyasm
echo testing %YASMPATH%
if exist "%YASMPATH%\yasm.exe" (
	set YASMEXE="%YASMPATH%\yasm.exe"
	goto :gotyasm
)
if exist "%YASMPATH%\vsyasm.exe" (
	set YASMEXE="%YASMPATH%\vsyasm.exe"
	goto :gotyasm
)
if exist "%VS100COMNTOOLS%\..\..\VC\bin\vsyasm.exe" (
	set YASMEXE="%VS100COMNTOOLS%\..\..\VC\bin\vsyasm.exe"
	goto :gotyasm
)
if exist "c:\Program Files (x86)\MSBuild\Microsoft.Cpp\v4.0\BuildCustomizations\vsyasm.exe" (
	set YASMEXE="c:\Program Files (x86)\MSBuild\Microsoft.Cpp\v4.0\BuildCustomizations\vsyasm.exe"
	goto :gotyasm
)
if exist "%VS90COMNTOOLS%\..\..\VC\bin\yasm.exe" (
	set YASMEXE="%VS90COMNTOOLS%\..\..\VC\bin\yasm.exe"
	goto :gotyasm
)
echo cant find yasm
exit /b 1
:gotyasm

:: set config.params.bat to the settings needed by make etc
echo (set LIBTYPE=%LIBTYPE%) > config.params.bat
:: Added /W3 and /WX- to produce more warinings
:: echo (set FLAGS=/W2 /WX- /Ox /Ot /D "NDEBUG" /D "HAVE_CONFIG_H" /nologo /D "_MBCS" /GS-) >> config.params.bat
echo (set FLAGS=%DFLAGS% /D "HAVE_CONFIG_H" /nologo /D "_MBCS" /GS-) >> config.params.bat
:: Changed /MT to /MD for GMPY2 compatibility
:: if %LIBTYPE% == lib (set FLAGS1=/Oi /D "_LIB" /D "PIC" /MT)
if %LIBTYPE% == lib (set FLAGS1=/D "_LIB" /D "PIC" %DFLAGS1%)
if %LIBTYPE% == dll (set FLAGS1=/D "__GMP_LIBGMP_DLL" /D "_WINDLL" /GF /EHsc /MD)
echo (set FLAGS1=%FLAGS1%) >> config.params.bat
set MPNPATH=.
if %ABI% == 32 (
	if %CPU% == x86		(set MPNPATH=x86w)
	if %CPU% == i386	(set MPNPATH=x86w)
	if %CPU% == i486	(set MPNPATH=x86w)
	if %CPU% == i586	(set MPNPATH=x86w)
	if %CPU% == pentium	(set MPNPATH=x86w)
	if %CPU% == pentiummmx	(set MPNPATH=x86w)
	if %CPU% == pentiumpro	(set MPNPATH=x86w x86w\p6)
	if %CPU% == i686	(set MPNPATH=x86w x86w\p6)
	if %CPU% == pentium2	(set MPNPATH=x86w x86w\p6 x86w\p6\mmx)
	if %CPU% == pentium3	(set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == k6		(set MPNPATH=x86w)
	if %CPU% == k62		(set MPNPATH=x86w)
	if %CPU% == k63		(set MPNPATH=x86w)
	if %CPU% == k7		(set MPNPATH=x86w)
	if %CPU% == athlon	(set MPNPATH=x86w)
	if %CPU% == pentium4	(set MPNPATH=x86w x86w\pentium4 x86w\pentium4\mmx x86w\pentium4\sse2)
	if %CPU% == prescott	(set MPNPATH=x86w x86w\pentium4 x86w\pentium4\mmx x86w\pentium4\sse2)
	if %CPU% == core	(set MPNPATH=x86w x86w\pentium4 x86w\pentium4\mmx x86w\pentium4\sse2)
	if %CPU% == viac3	(set MPNPATH=x86w)
	if %CPU% == viac32	(set MPNPATH=x86w)
	if %CPU% == x86_64	(set MPNPATH=x86w)
	if %CPU% == netburst	(set MPNPATH=x86w x86w\pentium4 x86w\pentium4\mmx x86w\pentium4\sse2)
	if %CPU% == netburstlahf	(set MPNPATH=x86w x86w\pentium4 x86w\pentium4\mmx x86w\pentium4\sse2)
	if %CPU% == k8           (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == k10          (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == k102         (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == k103         (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == bulldozer    (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == bobcat       (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == core2        (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == penryn       (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == nehalem      (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == westmere     (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == sandybridge  (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == atom         (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
	if %CPU% == nano         (set MPNPATH=x86w x86w\p6 x86w\p6\mmx x86w\p6\p3mmx)
)
if %ABI% == 64 (
	if %CPU% == x86_64       (set MPNPATH=x86_64w)
	if %CPU% == netburst     (set MPNPATH=x86_64w x86_64w\netburst)
	if %CPU% == netburstlahf (set MPNPATH=x86_64w x86_64w\netburst)
	if %CPU% == k8           (set MPNPATH=x86_64w x86_64w\k8 x86_64w\k8\k8only)
	if %CPU% == k10          (set MPNPATH=x86_64w x86_64w\k8 x86_64w\k8\k10)
	if %CPU% == k102         (set MPNPATH=x86_64w x86_64w\k8 x86_64w\k8\k10 x86_64w\k8\k10\k102)
	if %CPU% == k103         (set MPNPATH=x86_64w x86_64w\k8 x86_64w\k8\k10 x86_64w\k8\k10\k102)
	if %CPU% == bulldozer    (set MPNPATH=x86_64w x86_64w\k8 x86_64w\k8\k10 x86_64w\k8\k10\k102)
	if %CPU% == bobcat       (set MPNPATH=x86_64w x86_64w\bobcat)
	if %CPU% == core2        (set MPNPATH=x86_64w x86_64w\core2)
	if %CPU% == penryn       (set MPNPATH=x86_64w x86_64w\core2 x86_64w\core2\penryn)
	if %CPU% == nehalem      (set MPNPATH=x86_64w x86_64w\nehalem)
	if %CPU% == westmere     (set MPNPATH=x86_64w x86_64w\nehalem x86_64w\nehalem\westmere)
	if %CPU% == sandybridge  (set MPNPATH=x86_64w x86_64w\sandybridge)
	if %CPU% == atom         (set MPNPATH=x86_64w x86_64w\atom)
	if %CPU% == nano         (set MPNPATH=x86_64w x86_64w\k8 x86_64w\k8\k8only)
)
echo (set MPNPATH=%MPNPATH%) >> config.params.bat
echo (set ABI=%ABI%) >> config.params.bat
echo .
echo CPU identified as %CPU%.
echo .
echo Setting parameters to:
type config.params.bat
exit /b 0
