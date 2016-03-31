# Multiphysique

## build on Ubuntu (with basic linear solver)
```
mkdir build
cd build
cmake ..
make -j 6
```
## build on Ubuntu (with MUMPS)

*install MUMPS*
```
sudo apt-get install libmumps-seq-dev
```
*build project with MUMPS*
```
cmake -DMP_USE_MUMPS=ON ..
```

## run a test
```
time ./MP ../geo/carreSimple.msh ../geo/carreSimple.phy
```
Example: carreSimple (square meshed with 25x25x2 elements)
* CPU time (without MUMPS): 1 min 27s 
* CPU time (with MUMPS): 1s


