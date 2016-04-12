@echo off
::
:: compilation mingw - codeblocks
::
:: cmake -G "CodeBlocks - MinGW Makefiles" -DBLA_VENDOR=OpenBlas -DMP_USE_MUMPS=ON ..
:: cmake --build . 
::   ou
:: codeblocks Multiphysique.cbp
::
:: MP.exe ..\geo\carreSimple.msh ..\geo\carreSimple.phy

echo setting CodeBlocks environment

set INCLUDE=C:\Users\Boman\devgcc\OpenBLAS-v0.2.14-Win32\include
set INCLUDE=%INCLUDE%;C:\Users\Boman\devgcc\MUMPS\include

set PATH=C:\Program Files (x86)\CodeBlocks\MinGW\bin
set PATH=%PATH%;C:\Users\Boman\devgcc\OpenBLAS-v0.2.14-Win32\bin
set PATH=%PATH%;C:\Program Files (x86)\CMake\bin
set PATH=%PATH%;C:\Users\Boman\devgcc\MUMPS\bin
set PATH=%PATH%;C:\Program Files (x86)\CodeBlocks

%comspec%
