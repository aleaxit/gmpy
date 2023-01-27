set -e

GMP_VERSION=6.2.1
MPFR_VERSION=4.1.1
MPC_VERSION=1.2.1
if [ ! -f finish_before_ci_build ]; then
  if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux-musl" || "$OSTYPE" == "darwin"* ]]; then
    curl -O https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz
    tar -xvf gmp-${GMP_VERSION}.tar.xz
    # need to set host to the oldest triple to avoid building binaries
    # that use build machine micro-architecure. configfsf.guess is the one that
    # comes with autotools which is micro-architecture agnostic.
    # config.guess is a custom gmp script which knows about micro-architectures.
    cd gmp-${GMP_VERSION}
    patch -N -Z -p0 < ../scripts/patch-arm64.diff && cd ..
    cd gmp-${GMP_VERSION} && ./configure --enable-fat --host=$(./configfsf.guess) && make -j4 && make install && cd ../
    curl -O -k https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz
    tar -xvf mpfr-${MPFR_VERSION}.tar.gz
    cd mpfr-${MPFR_VERSION} && ./configure && make -j4 && make install && cd ../
    curl -O https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz
    tar -xvf mpc-${MPC_VERSION}.tar.gz
    cd mpc-${MPC_VERSION} && ./configure && make -j4 && make install && cd ../
    pip install Cython
  fi
  touch finish_before_ci_build
else
  echo "has finished before ci build"
fi
