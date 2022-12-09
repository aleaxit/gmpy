@echo off
 doskey py27=D:\64\Python27\python.exe $*
 doskey py35=D:\64\Python35\python.exe $*
 doskey py36=D:\64\Python36\python.exe $*
 doskey py37=D:\64\Python37\python.exe $*
 doskey py38=D:\64\Python38\python.exe $*
 doskey py39=D:\64\Python39\python.exe $*
 doskey py310=D:\64\Python310\python.exe $*
 doskey py311=D:\64\Python311\python.exe $*

 doskey py27_build=D:\64\Python27\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py35_build=D:\64\Python35\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py36_build=D:\64\Python36\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py37_build=D:\64\Python37\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py38_build=D:\64\Python38\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py39_build=D:\64\Python39\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py310_build=D:\64\Python310\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py311_build=D:\64\Python311\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel

 doskey py27_install=D:\64\Python27\python.exe -m pip install dist\gmpy2-2.1.4-cp27-cp27m-win_amd64.whl --force-reinstall
 doskey py35_install=D:\64\Python35\python.exe -m pip install dist\gmpy2-2.1.4-cp35-cp35m-win_amd64.whl --force-reinstall
 doskey py36_install=D:\64\Python36\python.exe -m pip install dist\gmpy2-2.1.4-cp36-cp36m-win_amd64.whl --force-reinstall
 doskey py37_install=D:\64\Python37\python.exe -m pip install dist\gmpy2-2.1.4-cp37-cp37m-win_amd64.whl --force-reinstall
 doskey py38_install=D:\64\Python38\python.exe -m pip install dist\gmpy2-2.1.4-cp38-cp38-win_amd64.whl --force-reinstall
 doskey py39_install=D:\64\Python39\python.exe -m pip install dist\gmpy2-2.1.4-cp39-cp39-win_amd64.whl --force-reinstall
 doskey py310_install=D:\64\Python310\python.exe -m pip install dist\gmpy2-2.1.4-cp310-cp310-win_amd64.whl --force-reinstall
 doskey py311_install=D:\64\Python311\python.exe -m pip install dist\gmpy2-2.1.4-cp311-cp311-win_amd64.whl --force-reinstall

 set Path=D:\mingw64\bin;%Path%
