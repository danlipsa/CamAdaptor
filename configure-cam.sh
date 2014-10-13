#! /bin/sh

INC_NETCDF=/usr/include;export INC_NETCDF
LIB_NETCDF=/usr/lib;export LIB_NETCDF

CAM_ROOT=/home/danlipsa/src/cesm1_2_1
camcfg=$CAM_ROOT/models/atm/cam/bld

# serial
#$camcfg/configure -fc gfortran -dyn fv -hgrid 10x15 -nospmd -nosmp -debug -test
# SPMD (mpi)
$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn fv -hgrid 10x15 -ntasks 6 -debug -nosmp -catalyst -test
# SMP (multi-threading)
#$camcfg/configure -fc=gfortran -dyn fv -hgrid 10x15 -nospmd -nthreads 6 -test

# make -j8 | tee make.out
