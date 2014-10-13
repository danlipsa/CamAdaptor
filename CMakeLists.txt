cmake_minimum_required(VERSION 2.6)
PROJECT(cam_catalyst_adapter CXX Fortran)

find_package(ParaView 4.1 REQUIRED COMPONENTS vtkPVPythonCatalyst)

INCLUDE(${PARAVIEW_USE_FILE})
if(NOT PARAVIEW_USE_MPI)
  message(SEND_ERROR "ParaView must be built with MPI enabled")
endif()
if(NOT MPI_Fortran_LIBRARIES)
  find_package(MPI)
endif()

set(CAM_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/../" CACHE PATH "Path to CAM build directory")
include_directories(${CAM_BUILD_DIR})

ADD_LIBRARY(cam_catalyst_adapter SHARED cam_catalyst_adapter.F90 create_data.cxx)
TARGET_LINK_LIBRARIES(cam_catalyst_adapter vtkPVPythonCatalyst ${MPI_Fortran_LIBRARIES})
