// Adaptor for getting Fortran simulation code into ParaView CoProcessor.

// CoProcessor specific headers
#include <iostream>
#include <sstream>
#include <algorithm>
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellType.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "vtkCPPythonScriptPipeline.h"
#include "vtkSmartPointer.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include "vtkUnstructuredGrid.h"

// Fortran specific header
#include "vtkCPPythonAdaptorAPI.h"

namespace
{
  const double PI = atan(1.0) * 4.0;
  int g_rank = -1;
  int g_nstep = -1;
  int g_maxNcols;     // maximum number of (vertical) columns in a chunk

  // 2d grid
  vtkIdType g_currentPoint2d;
  int g_rankArrayIndex2d = -1;
  int g_psArrayIndex = -1;

  // 3d grid
  int g_levelSize;
  vtkIdType g_currentPoint3d;
  int g_rankArrayIndex3d = -1;  
  int g_tArrayIndex = -1;
  int g_uArrayIndex = -1;
  int g_vArrayIndex = -1;
  double* g_level;

  template<typename T>
  int addAttribute(vtkUnstructuredGrid* grid, const char* name, vtkIdType size)
  {
    vtkSmartPointer<T> a = vtkSmartPointer<T>::New();
    a->SetNumberOfComponents(1);
    a->SetNumberOfValues(size);
    a->SetName(name);
    return grid->GetCellData()->AddArray(a);    
  }
};

extern "C" void cxx_finalize_()
{
  if (g_level)
    {
    delete[] g_level;
    }
}

extern "C" void cxx_create_grid_(
  int* dim, double* lonCoord, double* latCoord, double* levCoord,
  int* nPoints2d, int* maxNcols,
  int* myRank)
{
  // debugging
  std::ostringstream ostr;
  ostr << "\nDimensions: ";
  ostr << dim[0] << ", " << dim[1] << ", " << dim[2] << std::endl;
  ostr << "lon: ";
  for (int i = 0; i < dim[0]; ++i)
    ostr << lonCoord[i] << " ";
  ostr << "\nlat: ";
  for (int i = 0; i < dim[1]; ++i)
    ostr << latCoord[i] << " ";
  ostr << "\nlev: ";
  for (int i = 0; i < dim[2]; ++i)
    ostr << levCoord[i] << " ";
  ostr << endl;
  ostr << "\nnPoints2d: " << *nPoints2d 
       << "\nmyRank: " << *myRank << endl;
  std::cerr << ostr.str();


  g_rank = *myRank;
  g_nstep = -1;
  g_maxNcols = *maxNcols;
  if (!vtkCPPythonAdaptorAPI::GetCoProcessorData())
    {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
    }

  // 2d grid
  vtkSmartPointer<vtkUnstructuredGrid> grid2d = 
    vtkSmartPointer<vtkUnstructuredGrid>::New();  
  vtkSmartPointer<vtkPoints> points2d = vtkSmartPointer<vtkPoints>::New();
  points2d->SetNumberOfPoints(*nPoints2d);  
  grid2d->SetPoints(points2d);
  vtkCPInputDataDescription* idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->
    GetInputDescriptionByName("input");
  if (!idd)
    {
    std::ostringstream ostr;
    ostr << "No vtkCPInputDataDescription for 'input'";
    vtkGenericWarningMacro(<< ostr.str());
    return;
    }
  idd->SetGrid(grid2d);

  // 3D grid
  g_levelSize = dim[2];
  g_level = new double[g_levelSize];
  std::copy(levCoord, levCoord + g_levelSize, g_level);
  vtkSmartPointer<vtkUnstructuredGrid> grid3d = 
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points3d = vtkSmartPointer<vtkPoints>::New();
  points3d->SetNumberOfPoints(*nPoints2d * g_levelSize);  
  grid3d->SetPoints(points3d);
  idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->
    GetInputDescriptionByName("input3D");
  if (!idd)
    {
    std::ostringstream ostr;
    ostr << "No vtkCPInputDataDescription for 'input3D'";
    vtkGenericWarningMacro(<< ostr.str());
    return;
    }
  idd->SetGrid(grid3d);
}


// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void cxx_add_chunk_(
  int* nstep, int* ncols,
  double* lonIndexes, double* latIndexes, 
  double* psScalar, double *tScalar, double* uScalar, double* vScalar)
{
  // get the grids from Catalyst
  vtkCPInputDataDescription* idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->
    GetInputDescriptionByName("input");
  if (!idd)
    {
    return;
    }
  vtkUnstructuredGrid* grid2d = 
    vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  idd = vtkCPPythonAdaptorAPI::GetCoProcessorData()->
    GetInputDescriptionByName("input3D");
  if (!idd)
    {
    return;
    }
  vtkUnstructuredGrid* grid3d = 
    vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  int nPoints2d = grid2d->GetPoints()->GetNumberOfPoints();

  //debugging
  std::ostringstream ostr;
  ostr << "add_chunk:"
       << "\nncols: " << *ncols
       << "\nnPoints: " << nPoints2d
       << "\nmyRank: " << g_rank
       << "\nlon: ";
  for (int i = 0; i < *ncols; ++i)
    {
    ostr << lonIndexes[i] / PI * 180.0 << " ";
    }
  ostr << "\nlat: ";  
  for (int i = 0; i < *ncols; ++i)
    {
    ostr << latIndexes[i] / PI * 180.0 << " ";
    }
  ostr << "\nPS: ";
  for (int i = 0; i < *ncols; ++i)
    {
    ostr << psScalar[i] << " ";
    }
  ostr << "\nt: ";
  for (int i = 0; i < 2*(*ncols); ++i)
    {
    ostr << tScalar[i] << " ";
    }
  ostr << endl;
  std::cerr << ostr.str();

  if (*nstep > g_nstep)
    {
    g_nstep = *nstep;

    // 2d grid
    g_currentPoint2d = 0;
    grid2d->GetCells()->Initialize();
    grid2d->Allocate(nPoints2d);
    g_rankArrayIndex2d = addAttribute<vtkIntArray>(grid2d, "Rank", nPoints2d);
    g_psArrayIndex = addAttribute<vtkDoubleArray>(grid2d, "PS", nPoints2d);

    // 3d grid
    g_currentPoint3d = 0;
    grid3d->GetCells()->Initialize();
    grid3d->Allocate(nPoints2d * g_levelSize);
    g_rankArrayIndex3d = addAttribute<vtkIntArray>(
      grid3d, "Rank", nPoints2d * g_levelSize);
    g_tArrayIndex = addAttribute<vtkDoubleArray>(
      grid3d, "T", nPoints2d * g_levelSize);
    g_uArrayIndex = addAttribute<vtkDoubleArray>(
      grid3d, "U", nPoints2d * g_levelSize);
    g_vArrayIndex = addAttribute<vtkDoubleArray>(
      grid3d, "V", nPoints2d * g_levelSize);
    }

  for (int i = 0; i < *ncols; ++i)
    {
    // 2d grid
    grid2d->GetPoints()->SetPoint(g_currentPoint2d, 
                                  lonIndexes[i] / PI * 180.0, 
                                  latIndexes[i] / PI * 180.0, 0);
    grid2d->InsertNextCell(VTK_VERTEX, 1, &g_currentPoint2d);
    vtkDoubleArray::SafeDownCast(
      grid2d->GetCellData()->GetArray(g_psArrayIndex))
      ->SetValue(g_currentPoint2d, psScalar[i]);
    vtkIntArray::SafeDownCast(
      grid2d->GetCellData()->GetArray(g_rankArrayIndex2d))
      ->SetValue(g_currentPoint2d, g_rank);
    ++g_currentPoint2d;
    
    // 3d grid
    for (int j = 0; j < g_levelSize; ++j)
      {
      grid3d->GetPoints()->SetPoint(g_currentPoint3d, 
                                    lonIndexes[i] / PI * 180.0, 
                                    latIndexes[i] / PI * 180.0, g_level[j]);
      grid3d->InsertNextCell(VTK_VERTEX, 1, &g_currentPoint3d);
      vtkIntArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_rankArrayIndex3d))
        ->SetValue(g_currentPoint3d, g_rank);
      vtkDoubleArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_tArrayIndex))
        ->SetValue(g_currentPoint3d, tScalar[i + j*g_maxNcols]);
      vtkDoubleArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_uArrayIndex))
        ->SetValue(g_currentPoint3d, uScalar[i + j*g_maxNcols]);
      vtkDoubleArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_vArrayIndex))
        ->SetValue(g_currentPoint3d, vScalar[i + j*g_maxNcols]);
      ++g_currentPoint3d;
      }
    }
}
