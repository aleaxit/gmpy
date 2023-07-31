PATH="/c/mingw64/bin:"$PATH

GMP_VERSION=6.2.1
MPFR_VERSION=4.2.0
MPC_VERSION=1.3.1

cd /c/64

mkdir src
mkdir shared

cd /c/64/src

if [ ! -f gmp-${GMP_VERSION}.tar.lz ]; then
    wget https://gmplib.org/download/gmp/gmp-${GMP_VERSION}.tar.lz
fi
tar xf gmp-${GMP_VERSION}.tar.lz

if [ ! -f mpfr-${MPFR_VERSION}.tar.bz2 ]; then
    wget https://www.mpfr.org/mpfr-current/mpfr-${MPFR_VERSION}.tar.bz2
fi
tar xf mpfr-${MPFR_VERSION}.tar.bz2
cd mpfr-${MPFR_VERSION}/
wget https://www.mpfr.org/mpfr-current/allpatches
patch -N -Z -p1 < allpatches
cd ..

if [ ! -f mpc-${MPC_VERSION}.tar.gz ]; then
    wget ftp://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz
fi
tar xf mpc-${MPC_VERSION}.tar.gz

cd gmp-${GMP_VERSION}/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/c/64/shared --enable-shared --disable-static --enable-fat --with-pic
mv gmp.h gmp.original
sed 's/typedef\s*unsigned\s*long\s*int\s*mp_bitcnt_t/typedef unsigned long long int  mp_bitcnt_t\n/g' gmp.original > gmp.h
make -j8
make check
make install
cd ..

cd mpfr-${MPFR_VERSION}/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/c/64/shared --enable-shared --disable-static --disable-decimal-float --disable-float128 --with-pic --with-gmp=/c/64/shared
make -j8
make check
make install
cd ..

cd mpc-${MPC_VERSION}/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/c/64/shared --enable-shared --disable-static --with-pic --with-gmp=/c/64/shared --with-mpfr=/c/64/shared
make -j8
make check
make install
cd ..
