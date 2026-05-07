nvfortran sam-openacc.f90 -O3 -c -acc=gpu 
nvfortran driver.f90 -O3 -c -acc=gpu
nvfortran sam-openacc.o driver.o -O3 -o test -acc=gpu
./test
