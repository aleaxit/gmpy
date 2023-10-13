GMP_VERSION=6.3.0
MPFR_VERSION=4.2.1
MPC_VERSION=1.3.1

cd $HOME/temp

mkdir src
mkdir shared

cd $HOME/temp/src

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
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=$HOME/temp/shared --enable-shared --disable-static --enable-fat --with-pic

# Uncomment the following lines to change the mp_bitcnt_t type to "unsigned long long int"
mv gmp.h gmp.original
sed 's/typedef\s*unsigned\s*long\s*int\s*mp_bitcnt_t/typedef unsigned long long int  mp_bitcnt_t\n/g' gmp.original > gmp.h
# Note: also need to fix mis-matched function signatures in millerrabin.c

make -j32
make -j16 check
make install
cd ..

cd mpfr-${MPFR_VERSION}/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=$HOME/temp/shared --enable-shared --disable-static --disable-decimal-float --disable-float128 --with-pic --with-gmp=$HOME/temp/shared
make -j32
make check
make install
cd ..

cd mpc-${MPC_VERSION}/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=$HOME/temp/shared --enable-shared --disable-static --with-pic --with-gmp=$HOME/temp/shared --with-mpfr=$HOME/temp/shared
make -j16
make check
make install
cd ..
