Windows Build from Source
=========================

Building gmpy2 and the related libraries from scratch requires several steps. The following
instructions create a dedicated build environment using standalone

Install MSYS2
-------------

Browse to https://www.msys2.org and download and run the current installer. Accept the default
path of C:\msys64

Update MSYS2
------------

pacman -Syuu

pacman -S patch m4 lzip wget tar make diffutils

Install the MinGW-w64 Compiler Suite
------------------------------------

MSYS2 follows a "rolling-release" model. New versions of applications, including the MinGW-w64
compiler suite, are released frequently. The continual release cycles make it difficult to 
provide a reproducible build environment. We will be using MinGW-w64 as distributed by 
http://winlibs.com/. This documentation was based on the following version:

https://github.com/brechtsanders/winlibs_mingw/releases/download/10.3.0-12.0.0-9.0.0-r2/winlibs-x86_64-posix-seh-gcc-10.3.0-mingw-w64-9.0.0-r2.zip

Unzip the downloaded file and copy the "mingw64" directory to the desired location.
The following documentation assumes the files are copied to C:\mingw64.

The Build Environment
---------------------

The build process is split into two phases. During phase 1 we use the MSYS2 shell environment
and the winlibs.com MinGW-w64 compiler suite to build GMP, MPFR, and MPC. Once those libraries
are built, we will use the Windows command prompt and the same compiler to build the actual
gmpy2 extension (DLL). 

The MSYS2 environment provides different command line operating environments. We will just
use the generic "MSYS2 MSYS" environment.

 * MSYS2 MSYS
   This is the general-purpose MSYS shell but it does not provide
   access to the MinGW compiler suite. We will use this shell and
   manually add the location of the winlibs.com version.

MSYS2 does include versions of GMP, MPFR, and MPC but we will compile our own version directly
from the source. We do this so we can create reproducible builds.

Compiling GMP, MPFR, and MPC
----------------------------

Start the appropriate shell: MSYS2 MSYS

Add the MinGW-w64 compiler to the PATH.

PATH="/c/mingw64/bin:"$PATH

In your home directory, create directories for the various files that are created.

mkdir $HOME/64
mkdir $HOME/64/src
mkdir $HOME/64/static
mkdir $HOME/64/shared

Download GMP, MPFR, and MPC
---------------------------

Execute the file