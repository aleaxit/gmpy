set -e -x

GMP_VERSION=6.3.0
MPFR_VERSION=4.2.1
MPC_VERSION=1.3.1
if [ ! -f finish_before_ci_build ]; then
  if [[ "$OSTYPE" == "linux-gnu" || "$OSTYPE" == "linux-musl" || "$OSTYPE" == "darwin"* ]]; then
    curl -s -O https://ftp.gnu.org/gnu/gmp/gmp-${GMP_VERSION}.tar.xz
    tar -xf gmp-${GMP_VERSION}.tar.xz
    cd gmp-${GMP_VERSION}
    # config.guess uses microarchitecture and configfsf.guess doesn't
    # We replace config.guess with configfsf.guess to avoid microarchitecture
    # specific code in common code.
    rm config.guess && mv configfsf.guess config.guess && chmod +x config.guess
    ./configure --enable-fat && make -j4 && make install && cd ../
    curl -s -O https://ftp.gnu.org/gnu/mpfr/mpfr-${MPFR_VERSION}.tar.gz
    tar -xf mpfr-${MPFR_VERSION}.tar.gz
    cd mpfr-${MPFR_VERSION} && ./configure && make -j4 && make install && cd ../
    curl -s -O https://ftp.gnu.org/gnu/mpc/mpc-${MPC_VERSION}.tar.gz
    tar -xf mpc-${MPC_VERSION}.tar.gz
    cd mpc-${MPC_VERSION} && ./configure && make -j4 && make install && cd ../
  fi
  touch finish_before_ci_build
else
  (cd gmp-*[0-9] && make install)
  (cd mpfr-*[0-9] && make install)
  (cd mpc-*[0-9] && make install)
fi
