#!/bin/bash

set -e -x

GMP_VERSION=6.3.0
MPFR_VERSION=4.2.1
MPC_VERSION=1.3.1

PREFIX="$(pwd)/.local/"

# -- build GMP --
curl -s -O https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz
tar -xf gmp-${GMP_VERSION}.tar.xz
cd gmp-${GMP_VERSION}
# config.guess uses microarchitecture and configfsf.guess doesn't
# We replace config.guess with configfsf.guess to avoid microarchitecture
# specific code in common code.
rm config.guess && mv configfsf.guess config.guess && chmod +x config.guess
./configure --enable-fat \
            --enable-shared \
            --disable-static \
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
            --with-gmp=$PREFIX \
            --prefix=$PREFIX
make -j6
make install
cd ../
# -- build MPC --
curl -s -O https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz
tar -xf mpc-${MPC_VERSION}.tar.gz
cd mpc-${MPC_VERSION}
./configure --enable-shared \
            --disable-static \
            --with-gmp=$PREFIX \
            --with-mpfr=$PREFIX \
            --prefix=$PREFIX
make -j6
make install
cd ../

# -- copy headers --
cp $PREFIX/include/{gmp,mpfr,mpc}.h gmpy2/
