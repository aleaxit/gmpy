#!/bin/bash

set -e -x

GMP_VERSION=6.3.0
MPFR_VERSION=4.2.1
MPC_VERSION=1.3.1

unset CFLAGS

PREFIX="$(pwd)/.local/"

# -- build GMP --
curl -s -O https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz
tar -xf gmp-${GMP_VERSION}.tar.xz
cd gmp-${GMP_VERSION}
# Patch the mp_bitcnt_t to "unsigned long long int" on WINDOWS AMD64:
patch -N -Z -p0 < ../scripts/mp_bitcnt_t.diff
patch -N -Z -p0 < ../scripts/fat_build_fix.diff
if [ "$OSTYPE" = "msys" ] || [ "$OSTYPE" = "cygwin" ]
then
  patch -N -Z -p0 < ../scripts/dll-importexport.diff
fi
# config.guess uses microarchitecture and configfsf.guess doesn't
# We replace config.guess with configfsf.guess to avoid microarchitecture
# specific code in common code.
rm config.guess && mv configfsf.guess config.guess && chmod +x config.guess
./configure --enable-fat \
            --enable-shared \
            --disable-static \
            --with-pic \
            --prefix=$PREFIX
make -j6
make install
cd ../

# -- build MPFR --
curl -s -O https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz
tar -xf mpfr-${MPFR_VERSION}.tar.gz
cd mpfr-${MPFR_VERSION}
./configure --enable-shared \
            --disable-static \
            --with-pic \
            --with-gmp=$PREFIX \
            --prefix=$PREFIX
make -j6
make install
cd ../
# -- build MPC --
curl -s -O https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz
tar -xf mpc-${MPC_VERSION}.tar.gz
cd mpc-${MPC_VERSION}

## add pkg-config data
#cat > mpc.pc.in <<'EOF'
#prefix=@prefix@
#exec_prefix=@exec_prefix@
#includedir=@includedir@
#libdir=@libdir@
#
#Name: @PACKAGE_NAME@
#Description: GNU MPC is a complex floating-point library with exact rounding.
#URL: https://www.multiprecision.org/
#Version: @PACKAGE_VERSION@
#Cflags: -I${includedir}
#Libs: -L${libdir} -lmpc
#EOF
#patch -N -Z -p0 < ../scripts/mpc-pkg-config.diff
#autoreconf -vfi

./configure --enable-shared \
            --disable-static \
            --with-pic \
            --with-gmp=$PREFIX \
            --with-mpfr=$PREFIX \
            --prefix=$PREFIX
make -j6
make install

mkdir -p ${PREFIX}/lib/pkgconfig
cat > ${PREFIX}/lib/pkgconfig/mpc.pc <<EOF
prefix=${PREFIX}
exec_prefix=${PREFIX}/bin/
includedir=${PREFIX}/include/
libdir=${PREFIX}/lib/

Name: mpc
Description: GNU MPC is a complex floating-point library with exact rounding.
URL: https://www.multiprecision.org/
Version: 1.3.1
Cflags: -I${PREFIX}/include/
Libs: -L${PREFIX}/lib/ -lmpc
EOF

cd ../
