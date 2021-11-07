PATH="/d/mingw64/bin:"$PATH

cd /d/64

mkdir src
mkdir shared

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
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/d/64/shared --enable-shared --disable-static --enable-fat --with-pic
mv gmp.h gmp.original
make -j4
make check
make install
cd ..

cd mpfr-4.1.0/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/d/64/shared --enable-shared --disable-static --disable-decimal-float --disable-float128 --with-pic --with-gmp=/d/64/shared
make -j4
make check
make install
cd ..

cd mpc-1.2.1/
make distclean
./configure --build=x86_64-pc-mingw64 --host=x86_64-pc-mingw64 --prefix=/d/64/shared --enable-shared --disable-static --with-pic --with-gmp=/d/64/shared --with-mpfr=/d/64/shared
make -j4
make check
make install
cd ..
