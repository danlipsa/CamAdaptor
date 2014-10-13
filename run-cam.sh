#! /bin/sh

# runs the cam5 simulation using MPI on 6 processors
LD_LIBRARY_PATH=~/build/cam-5.3/cam_catalyst_adapter mpiexec -np 6 ~/build/cam-5.3/cam | tee cam.log
