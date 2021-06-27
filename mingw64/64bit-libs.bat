cd C:\64\Python27\libs
gendef c:\Windows\System32\python27.dll
dlltool --dllname c:\Windows\System32\python27.dll --def python27.def --output-lib libpython27.a

cd C:\64\Python35\libs
gendef ..\python35.dll
dlltool --dllname ..\python35.dll --def python35.def --output-lib libpython35.a

cd C:\64\Python36\libs
gendef ..\python36.dll
dlltool --dllname ..\python36.dll --def python36.def --output-lib libpython36.a

cd C:\64\Python37\libs
gendef ..\python37.dll
dlltool --dllname ..\python37.dll --def python37.def --output-lib libpython37.a

cd C:\64\Python38\libs
gendef ..\python38.dll
dlltool --dllname ..\python38.dll --def python38.def --output-lib libpython38.a

cd C:\64\Python39\libs
gendef ..\python39.dll
dlltool --dllname ..\python39.dll --def python39.def --output-lib libpython39.a

cd C:\64\Python310\libs
gendef ..\python310.dll
dlltool --dllname ..\python310.dll --def python310.def --output-lib libpython310.a
