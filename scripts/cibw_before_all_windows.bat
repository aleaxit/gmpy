@echo on
cd .local\bin
set "PATH=C:\Program Files\Microsoft Visual Studio\2022\Enterprise\VC\Tools\MSVC\14.38.33130\bin\HostX86\x64\"
call D:\a\gmpy\gmpy\scripts\dll2lib.bat 64 libgmp-10.dll
call D:\a\gmpy\gmpy\scripts\dll2lib.bat 64 libmpfr-6.dll
call D:\a\gmpy\gmpy\scripts\dll2lib.bat 64 libmpc-3.dll
ren libgmp-10.lib gmp.lib
ren libmpfr-6.lib mpfr.lib
ren libmpc-3.lib mpc.lib
