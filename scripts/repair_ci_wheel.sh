dest_dir=$1
wheel=$2
wheel unpack --dest ${dest_dir} ${wheel}
cp /usr/local/include/{gmp,mpfr,mpc}.h ${dest_dir}/gmpy2-*/gmpy2/
(cd ${dest_dir} && wheel pack gmpy2-*[0-9])
cp ${dest_dir}/gmpy2-*.whl ${wheel}
rm -rf ${dest_dir}/*
if [[ "$OSTYPE" == "darwin"* ]]
then
  delocate-wheel --require-archs $3 --lib-sdir ../gmpy2.libs -w ${dest_dir} -v ${wheel}
else
  auditwheel repair -w ${dest_dir} ${wheel}
fi
# cleanup libs from before_ci_build*.sh
(cd gmp-*[0-9] && make uninstall)
(cd mpfr-*[0-9] && make uninstall)
(cd mpc-*[0-9] && make uninstall)
