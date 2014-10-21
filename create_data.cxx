// Adaptor for getting Fortran simulation code into ParaView CoProcessor.

// CoProcessor specific headers
#include <iostream>
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkPointData.h"
#include "vtkImageData.h"

// Fortran specific header
#include "vtkCPPythonAdaptorAPI.h"


// These will be called from the Fortran "glue" code"
// Completely dependent on data layout, structured vs. unstructured, etc.
// since VTK/ParaView uses different internal layouts for each.

// Creates the data container for the CoProcessor.
extern "C" void create_grid_(
  int* dim, 
  double* lonCoord, double* latCoord, double* levCoord)
{
  std::cerr << "\nDimensions: ";
  std::cerr << dim[0] << ", " << dim[1] << ", " << dim[2] << std::endl;
  std::cerr << "lon: ";
  for (int i = 0; i < dim[0]; ++i)
    std::cerr << lonCoord[i] << " ";
  std::cerr << "\nlat: ";
  for (int i = 0; i < dim[1]; ++i)
    std::cerr << latCoord[i] << " ";
  std::cerr << "\nlev: ";
  for (int i = 0; i < dim[2]; ++i)
    std::cerr << levCoord[i] << " ";
  std::cerr << endl;
}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void add_field_(const char* name)
{
  std::cerr << name << endl;
}
