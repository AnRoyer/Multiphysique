@echo off
:: Ligne de commande qui va bien pour compiler mes brols
::
:: Utilisation:
::    - installer les libs (dans c:\local p. expl.)
::    - definir la variable d'env MYLOCAL  (=c:\local p. expl. )
::    - double cliquer sur myenv.bat
::    - cd build
::    - cmake -G "Visual Studio 11 Win64" -DMP_USE_MUMPS=ON ..
::    - cmake --build . --config Release
::    - ctest -C Release
::

set PATH=%PATH%;%MYLOCAL%\swigwin-2.0.12
set INCLUDE=%MYLOCAL%\include;%MYLOCAL%\MUMPS\include
set LIB=%MYLOCAL%\MUMPS\lib

:: marche pas...
:: C:\Windows\SysWOW64\cmd.exe /E:ON /V:ON /K ""C:\Program Files (x86)\Intel\Composer XE 2013 SP1\mkl\bin\mklvars.bat" intel64"

:: ok (MKL only)
::C:\Windows\SysWOW64\cmd.exe /E:ON /V:ON /K "C:\Program Files (x86)\Intel\Composer XE\mkl\bin\intel64\mklvars_intel64.bat"

:: fait freezer cmake...
:: C:\Windows\SysWOW64\cmd.exe /E:ON /V:ON /K ""C:\Program Files (x86)\Intel\Composer XE 2013 SP1\bin\ipsxe-comp-vars.bat" intel64 vs2012"
:: C:\Windows\SysWOW64\cmd.exe /E:ON /V:ON /K ""C:\Program Files (x86)\Intel\Composer XE 2013 SP1\bin\compilervars.bat" intel64 vs2012"

:: ok (MKL + VS en ligne de commande)
call "C:\Program Files (x86)\Intel\Composer XE\mkl\bin\intel64\mklvars_intel64.bat"
%comspec% /K ""C:\Program Files (x86)\Microsoft Visual Studio 11.0\VC\vcvarsall.bat" amd64"
