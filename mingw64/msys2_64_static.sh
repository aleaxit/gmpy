PATH="/d/mingw64/bin:"$PATH

cd /d/64

mkdir src
mkdir static

cd /d/64/src

if [ ! -f gmp-6.2.1.tar.lz ]; then
    wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
fi
tar xf gmp-6.2.1.tar.lz

if [ ! -f mpfr-4.1.0.tar.bz2 ]; then
    wget https://www.mpfr.org/mpfr-current/mpfr-4.1.0.tar.bz2
fi
tar xf mpfr-4.1.0.tar.bz2
cd mpfr-4.1.0/
wget https://www.mpfr.org/mpfr-current/allpatches
patch -N -Z -p1 < allpatches
cd ..

if [ ! -f mpc-1.2.1.tar.gz ]; then
    wget ftp://ftp.gnu.org/gnu/mpc/mpc-1.2.1.tar.gz
fi
tar xf mpc-1.2.1.tar.gz

cd gmp-6.2.1/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/d/64/static --enable-static --disable-shared --enable-fat --with-pic
mv gmp.h gmp.original
sed 's/typedef\s*unsigned\s*long\s*int\s*mp_bitcnt_t/typedef unsigned long long int  mp_bitcnt_t\n/g' gmp.original > gmp.h
make -j4
make check
make install
cd ..

cd mpfr-4.1.0/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/d/64/static --enable-static --disable-shared --disable-decimal-float --disable-float128 --with-pic --with-gmp=/d/64/static
make -j4
make check
make install
cd ..

cd mpc-1.2.1/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/d/64/static --enable-static --disable-shared --with-pic --with-gmp=/d/64/static --with-mpfr=/d/64/static
make -j4
make check
make install
cd ..
