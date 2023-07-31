@echo off
 doskey py37=py -3.7 $*
 doskey py38=py -3.8 $*
 doskey py39=py -3.9 $*
 doskey py310=py -3.10 $*
 doskey py311=py -3.11 $*
 doskey py312=py -3.12 $*

 doskey py37_build=py -3.7 setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py38_build=py -3.8 setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py39_build=py -3.9 setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py310_build=py -3.10 setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py311_build=py -3.11 setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py312_build=py -3.12 setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel

 doskey py37_install=cmd /c "for %f in (dist\gmpy2*37*) do py -3.7 -m pip install %f --force-reinstall"
 doskey py38_install=cmd /c "for %f in (dist\gmpy2*38*) do py -3.8 -m pip install %f --force-reinstall"
 doskey py39_install=cmd /c "for %f in (dist\gmpy2*39*) do py -3.9 -m pip install %f --force-reinstall"
 doskey py310_install=cmd /c "for %f in (dist\gmpy2*310*) do py -3.10 -m pip install %f --force-reinstall"
 doskey py311_install=cmd /c "for %f in (dist\gmpy2*311*) do py -3.11 -m pip install %f --force-reinstall"
 doskey py312_install=cmd /c "for %f in (dist\gmpy2*312*) do py -3.12 -m pip install %f --force-reinstall"

set Path=C:\mingw64\bin;%Path%
