@echo off
 doskey py37=D:\64\Python37\python.exe $*
 doskey py38=D:\64\Python38\python.exe $*
 doskey py39=D:\64\Python39\python.exe $*
 doskey py310=D:\64\Python310\python.exe $*
 doskey py311=D:\64\Python311\python.exe $*

 doskey py37_build=D:\64\Python37\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py38_build=D:\64\Python38\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py39_build=D:\64\Python39\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py310_build=D:\64\Python310\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel
 doskey py311_build=D:\64\Python311\python.exe setup.py build_ext -cmingw32 -Imingw64\shared\include -Lmingw64\shared\lib -f bdist_wheel

 doskey py37_install=cmd /c "for %f in (dist\gmpy2*37*) do D:\64\Python37\python.exe -m pip install %f --force-reinstall"
 doskey py38_install=cmd /c "for %f in (dist\gmpy2*38*) do D:\64\Python38\python.exe -m pip install %f --force-reinstall"
 doskey py39_install=cmd /c "for %f in (dist\gmpy2*39*) do D:\64\Python39\python.exe -m pip install %f --force-reinstall"
 doskey py310_install=cmd /c "for %f in (dist\gmpy2*310*) do D:\64\Python310\python.exe -m pip install %f --force-reinstall"
 doskey py311_install=cmd /c "for %f in (dist\gmpy2*311*) do D:\64\Python311\python.exe -m pip install %f --force-reinstall"

set Path=D:\mingw64\bin;%Path%
