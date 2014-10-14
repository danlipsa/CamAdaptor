#! /bin/sh

# run cam5 on one processor
LD_LIBRARY_PATH=~/build/cam-5.3/cam_catalyst_adapter ~/build/cam-5.3/cam | tee cam.log

# runs the cam5 simulation on xx MPI processors. Make sure the value
# for -np matches the value for -ntasks in configure-cam.sh
#LD_LIBRARY_PATH=~/build/cam-5.3/cam_catalyst_adapter mpiexec -np 2 ~/build/cam-5.3/cam | tee cam.log
