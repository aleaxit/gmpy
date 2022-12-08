GMP_VERSION=6.2.1
MPFR_VERSION=4.1.1
MPC_VERSION=1.2.1
export CPPFLAGS=" --target=arm64-apple-macos11"
export LDFLAGS=" -arch arm64"
EXTRA="--build=x86_64-apple-darwin --host=aarch64-apple-darwin --target=aarch64-apple-darwin"
if [ ! -f finish_before_ci_build ]; then
  if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux-musl" || "$OSTYPE" == "darwin"* ]]; then
    curl -O https://gmplib.org/download/gmp/gmp-${GMP_VERSION}.tar.xz
    tar -xf gmp-${GMP_VERSION}.tar.xz
    # need to set host to the oldest triple to avoid building binaries
    # that use build machine micro-architecure. configfsf.guess is the one that
    # comes with autotools which is micro-architecture agnostic.
    # config.guess is a custom gmp script which knows about micro-architectures.
    cd gmp-${GMP_VERSION} && ./configure $EXTRA --enable-fat && make -j4 && make install && cd ../
    curl -O -k https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz
    tar -xf mpfr-${MPFR_VERSION}.tar.gz
    cd mpfr-${MPFR_VERSION} && ./configure $EXTRA && make -j4 && make install && cd ../
    curl -O https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz
    tar -xf mpc-${MPC_VERSION}.tar.gz
    cd mpc-${MPC_VERSION} && ./configure $EXTRA && make -j4 && make install && cd ../
    pip install Cython
  fi
  touch finish_before_ci_build
else
  echo "has finished before ci build"
fi
