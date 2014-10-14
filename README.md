This project provides a ParaView Catalyst adaptor for CAM5 simulation.

To build CAM5 and the adapter:

1. Download CESM (Community Earth System Model) using instructions in:

    <http://www.cesm.ucar.edu/models/cesm1.2/cesm/doc/usersguide/book1.html>
        
   and CAM5 (Community Atmosphere Model) data using instructions in:
   
   <http://www.cesm.ucar.edu/models/cesm1.2/cam/docs/ug5_3/ug.html>

2. Create a CAM5 build directory, cd into it and configure it:

        mkdir ~/build/cam-5.3
        cd ~/build/cam-5.3
        ~/src/cesm1_2_1/models/atm/cam/cam_catalyst_adapter/configure-cam.sh

3. Try building CAM5. This step creates mod files used by the adapter.

        ~/src/cesm1_2_1/models/atm/cam/cam_catalyst_adapter/build-cam.sh
   You'll get errors as the adapter is not built yet.

3. Create adapter's build directory and configure it. Make sure you
specify ParaView's build (or installation) directory. Build the
adapter:

        mkdir ~/build/cam-5.3/cam_catalyst_adapter
        cd ~/build/cam-5.3/cam_catalyst_adapter
        ccmake ~/src/cesm1_2_1/models/atm/cam/cam_catalyst_adapter
        make

4. Build CAM5:

        cd ..
        ~/src/cesm1_2_1/models/atm/cam/cam_catalyst_adapter/build-cam.sh

5. Create a run directory for the simulation, create the namelist
   (parameters for the simulation) and run the simulation.

        mkdir ~/run/cam-5.3
        cd ~/run/cam-5.3
        ~/src/cesm1_2_1/models/atm/cam/cam_catalyst_adapter/build-namelist.sh
        ~/src/cesm1_2_1/models/atm/cam/cam_catalyst_adapter/run-cam.sh
   Namelist for the simulation are listed in:

   <http://www.cesm.ucar.edu/cgi-bin/eaton/namelist/nldef2html-cam5_3>
