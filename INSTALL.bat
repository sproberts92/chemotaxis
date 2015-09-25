@echo off

REM echo Cloning submodules...
REM git submodule update --init
REM echo Cloning submodules... Done.

echo Setting compiler environment variables...
call vcvarsall amd64
echo Done.

REM echo Compiling [...] library... 
REM pushd ext\[...]
REM nmake /f [...]
REM popd
REM echo Compiling [...] library... Done.

echo removing .gitignore...

pushd src
del /f /q .gitignore
popd

pushd ext
del /f /q .gitignore
popd

echo Done.

echo Install complete.