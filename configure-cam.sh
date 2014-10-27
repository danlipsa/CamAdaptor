#! /bin/sh

INC_NETCDF=/usr/include;export INC_NETCDF
LIB_NETCDF=/usr/lib;export LIB_NETCDF

CAM_ROOT=/home/danlipsa/src/cesm1_2_1
camcfg=$CAM_ROOT/models/atm/cam/bld

# serial
# $camcfg/configure -fc gfortran -dyn fv -hgrid 10x15 -nospmd -nosmp -debug -catalyst -test

# SPMD (mpi) Finite Volume
# -dyn (dynamical core) (default): fv
# - phys (physics package) (default): cam5
# - chem (prognostic chemistry package) (default): trop_mam3
#$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn fv -hgrid 10x15 -ntasks 2 -debug -nosmp -catalyst -test

# SPMD (mpi) Spectral Element
$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn se -hgrid ne30np4 -ntasks 2 -debug -nosmp -catalyst -test

# SMP (multi-threading)
#$camcfg/configure -fc=gfortran -dyn fv -hgrid 10x15 -nospmd -nthreads 6 -test
