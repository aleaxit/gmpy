pushd D:\64\Python27\libs
gendef C:\Windows\System32\python27.dll
popd dlltool --dllname C:\Windows\System32\python27.dll --def python27.def --output-lib libpython27.a

pushd D:\64\Python35\libs
gendef D:\64\Python35\python35.dll
dlltool --dllname D:\64\Python35\python35.dll --def python35.def --output-lib libpython35.a
popd

pushd D:\64\Python36\libs
gendef D:\64\Python36\python36.dll
dlltool --dllname D:\64\Python36\python36.dll --def python36.def --output-lib libpython36.a
popd

pushd D:\64\Python37\libs
gendef D:\64\Python37\python37.dll
dlltool --dllname D:\64\Python37\python37.dll --def python37.def --output-lib libpython37.a
popd

pushd D:\64\Python38\libs
gendef D:\64\Python38\python38.dll
dlltool --dllname D:\64\Python38\python38.dll --def python38.def --output-lib libpython38.a
popd

pushd D:\64\Python39\libs
gendef D:\64\Python39\python39.dll
dlltool --dllname D:\64\Python39\python39.dll --def python39.def --output-lib libpython39.a
popd

pushd D:\64\Python310\libs
gendef D:\64\Python310\python310.dll
dlltool --dllname D:\64\Python310\python310.dll --def python310.def --output-lib libpython310.a
popd
