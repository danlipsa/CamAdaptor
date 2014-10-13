#! /bin/sh

CAM_ROOT=/home/danlipsa/src/cesm1_2_1
camcfg=$CAM_ROOT/models/atm/cam/bld
CAM_BUILD=/home/danlipsa/build/cam-5.3
CSMDATA=/home/danlipsa/src/cesm-data;export CSMDATA

# should be run inside the run directory
$camcfg/build-namelist -test -config $CAM_BUILD/config_cache.xml -namelist "&seq_timemgr_inparm stop_n=3 stop_option='ndays' / &cam_inparm  nhtfrq=-12, -24, -24, -24, -24, -24 hfilename_spec='%c.cam.h%t.%y-%m-%d-%s.nc' print_step_cost=.true. / "

