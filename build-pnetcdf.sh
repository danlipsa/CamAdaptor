#! /bin/sh

CC=gcc;export CC
FC=gfortran;export FC
F90=gfortran;export F90
CXX=g++;export CXX;
MPICC=mpicc;export MPICC
MPIF77=mpif77;export MPIF77
~/src/parallel-netcdf-1.5.0/configure
make -j8

