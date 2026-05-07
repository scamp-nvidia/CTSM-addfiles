nvfortran sam-original.f90 -O3 -c 
nvfortran driver.f90 -O3 -c 
nvfortran sam-original.o driver.o -O3 -o test
./test
