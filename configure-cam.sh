#! /bin/sh

# Configures the simulation
INC_NETCDF=/usr/include;export INC_NETCDF
LIB_NETCDF=/usr/lib;export LIB_NETCDF

CAM_ROOT=/home/danlipsa/src/cesm1_2_1
camcfg=$CAM_ROOT/models/atm/cam/bld

# serial
# $camcfg/configure -fc gfortran -dyn fv -hgrid 10x15 -nospmd -nosmp -debug -catalyst -test

# SPMD (mpi) Finite Volume
# fv, cam5, trop_mam3
#$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn fv -hgrid 10x15 -ntasks 2 -debug -nosmp -catalyst -test

# SPMD (mpi)
# fv, cam5, trop_mam3
$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn fv -hgrid 4x5 -ntasks 4 -debug -nosmp -catalyst -test


# SPMD (mpi) Finite Volume
# fv, cam5, trop_mam3
#$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn fv -hgrid 2.5x3.33 -ntasks 8 -debug -nosmp -catalyst -test


# SPMD (mpi) Spectral Element
#$camcfg/configure -fc mpif90 -fc_type gnu -cc mpicc -dyn se -hgrid ne30np4 -ntasks 8 -debug -nosmp -catalyst -test

# SMP (multi-threading)
#$camcfg/configure -fc=gfortran -dyn fv -hgrid 10x15 -nospmd -nthreads 6 -test
