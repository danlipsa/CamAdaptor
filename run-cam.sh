#! /bin/sh

# run cam5 on one processor
# LD_LIBRARY_PATH=~/build/cam-5.3/cam_catalyst_adapter ~/build/cam-5.3/cam | tee cam.log

# runs the cam5 simulation on xx MPI processors. 
# WARNING
# Make sure the value for -np matches the value for -ntasks in configure-cam.sh
pwd=`pwd`
DIR_NAME=`basename $pwd`
CAM_BUILD=~/build/$DIR_NAME
LD_LIBRARY_PATH=${CAM_BUILD}/cam_catalyst_adapter mpiexec -np 4 ${CAM_BUILD}/cam > cam.log 2>&1
