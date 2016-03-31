# Multiphysique

## build on Ubuntu (with basic linear solver)
```
mkdir build
cd build
cmake ..
make
```
## build on Ubuntu (with MUMPS)

### install MUMPS
```
sudo apt-get install libmumps-seq-dev
```

### build project with MUMPS
```
cmake -DMP_USE_MUMPS=ON ..
```


## run a test
```
./MP ../geo/carreSimple.msh ../geo/carreSimple.phy
```
