if [ ! -f finish_before_ci_build ]; then
  if [[ "$OSTYPE" == "linux-gnu" ]]; then
    apt update
    apt install -y libmpc-dev
    pip install Cython
  elif [[ "$OSTYPE" == "darwin"* ]]; then
    brew install gmp mpfr libmpc
  fi
  touch finish_before_ci_build
else
  echo "has finished before ci build"
fi
