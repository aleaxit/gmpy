if [ ! -f finish_before_ci_build ]; then
  if [[ "$OSTYPE" == "linux-gnu" ]]; then
    yum install -y wget lzip
    wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.lz
    tar -xvf gmp-6.2.0.tar.lz
    cd gmp-6.2.0 && ./configure && make -j4 && make install && cd ../
    wget https://www.mpfr.org/mpfr-current/mpfr-4.0.2.tar.gz
    tar -xvf mpfr-4.0.2.tar.gz
    cd mpfr-4.0.2 && ./configure && make -j4 && make install && cd ../
    wget https://ftp.gnu.org/gnu/mpc/mpc-1.1.0.tar.gz
    tar -xvf mpc-1.1.0.tar.gz
    cd mpc-1.1.0 && ./configure && make -j4 && make install && cd ../
    pip install Cython
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew install gmp mpfr libmpc
  fi
  touch finish_before_ci_build
else
  echo "has finished before ci build"
fi
