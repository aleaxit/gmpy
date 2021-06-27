cd C:\64\Python27\libs
gendef c:\Windows\System32\python27.dll
dlltool --dllname c:\Windows\System32\python27.dll --def python27.def --output-lib libpython27.a

cd C:\64\Python35\libs
gendef c:\64\Python35\python35.dll
dlltool --dllname c:\Windows\System32\python35.dll --def python35.def --output-lib libpython35.a

cd C:\64\Python36\libs
gendef c:\64\Python36\python36.dll
dlltool --dllname c:\Windows\System32\python36.dll --def python36.def --output-lib libpython36.a

cd C:\64\Python37\libs
gendef c:\64\Python37\python37.dll
dlltool --dllname c:\Windows\System32\python37.dll --def python37.def --output-lib libpython37.a

cd C:\64\Python38\libs
gendef c:\64\Python38\python38.dll
dlltool --dllname c:\Windows\System32\python38.dll --def python38.def --output-lib libpython38.a

cd C:\64\Python39\libs
gendef c:\64\Python39\python39.dll
dlltool --dllname c:\Windows\System32\python39.dll --def python39.def --output-lib libpython39.a

cd C:\64\Python310\libs
gendef c:\64\Python310\python310.dll
dlltool --dllname c:\Windows\System32\python310.dll --def python310.def --output-lib libpython310.a
