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
./configure --enable-shared \
            --disable-static \
            --with-pic \
            --with-gmp=$PREFIX \
            --with-mpfr=$PREFIX \
            --prefix=$PREFIX
make -j6
make install
cd ../

# -- copy headers --
cp $PREFIX/include/{gmp,mpfr,mpc}.h gmpy2/

# -- generate *.lib files from *.dll on M$ Windows --
if [ "$OSTYPE" = "msys" ] || [ "$OSTYPE" = "cygwin" ]
then
  # Set path to dumpbin & lib
  PATH="$PATH:$(find "/c/Program Files/Microsoft Visual Studio/2022/" -name "Hostx86")/x64/"

  # See http://stackoverflow.com/questions/9946322/
  cd .local/bin
  for dll_file in libgmp-10.dll libmpfr-6.dll libmpc-3.dll
  do
    lib_name=$(basename -s .dll ${dll_file})
    exports_file=${lib_name}-exports.txt
    def_file=${lib_name}.def
    lib_file=${lib_name}.lib
    name=$(echo ${lib_name}|sed 's/^lib//;s/-[0-9]\+//')

    dumpbin //exports ${dll_file} > ${exports_file}

    echo LIBRARY ${lib_name} > ${def_file}
    echo EXPORTS >> ${def_file}
    cat ${exports_file} | awk 'NR>19 && $4 != "" {print $4 " @"$1}' >> ${def_file}
    sed -i 's/$/\r/' ${def_file}

    lib //def:${def_file} //out:${lib_file} //machine:x64

    rm ${exports_file} ${def_file} ${lib_name}.exp
    mv ${lib_file} ${name}.lib
  done
fi
