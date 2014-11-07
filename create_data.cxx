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

namespace
{
  const double PI = atan(1.0) * 4.0;

  // coprocessor
  vtkCPProcessor* g_coprocessor;
  vtkCPDataDescription* g_coprocessorData;
  bool g_isTimeDataSet;

  int g_rank = -1;
  int g_chunkCapacity;  // maximum number of (vertical) columns in a chunk
  int g_nCells2d;       // total number of 2D cells in a MPI processor
  int g_dim[3];         // lon x lat x lev
  double g_lonStep;
  double g_latStep;

  // 2d grid
  int g_rankArrayIndex2d = -1;
  int g_psArrayIndex = -1;
  vtkIdType* g_cellId2d;      // 2d array with cellId at lon x lat index
  vtkIdType* g_pointId2d;     // 2d array with pointIds at lon x lat index

  // 3d grid
  int g_rankArrayIndex3d = -1;
  int g_tArrayIndex = -1;
  int g_uArrayIndex = -1;
  int g_vArrayIndex = -1;
  double* g_level;            // level values
  vtkIdType* g_cellId3d;
  vtkIdType* g_pointId3d;


enum GridType
{
  RECTILINEAR,
  SPHERE
};

double toDegrees(double rad)
{
  return rad / PI * 180.0;
}

double toRadians(double deg)
{
  return deg * PI / 180.0;
}

template<typename T>
int addAttribute(vtkUnstructuredGrid* grid, const char* name, vtkIdType size)
{
  vtkSmartPointer<T> a = vtkSmartPointer<T>::New();
  a->SetNumberOfComponents(1);
  a->SetNumberOfValues(size);
  a->SetName(name);
  return grid->GetCellData()->AddArray(a);
}

void printCreateGrid(int* dim, double* lonCoord, double* latCoord,
                     double* levCoord, int* nCells2d, int* maxNcols,
                     int* myRank)
{

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
  ostr << "\nnCells2d: " << *nCells2d
       << "\nmyRank: " << *myRank << endl;
  std::cerr << ostr.str();
}

void printAddChunk(
  int* nstep, int nCells2d, int nPoints2d, int nCellArray, int* chunkCols,
  double* lonRad, double* latRad, double* psScalar, double *tScalar)
{
  std::ostringstream ostr;
  ostr << "add_chunk: " << *nstep
       << "\nchunkCols: " << *chunkCols
       << "\nnCells: " << nCells2d
       << "\nnPoints: " << nPoints2d
       << "\nnCellArrays: " << nCellArray
       << "\nmyRank: " << g_rank
       << "\nlon: ";
  for (int i = 0; i < *chunkCols; ++i)
    {
    ostr << toDegrees(lonRad[i]) << " ";
    }
  ostr << "\nlat: ";
  for (int i = 0; i < *chunkCols; ++i)
    {
    ostr << toDegrees(latRad[i]) << " ";
    }
  ostr << "\nPS: ";
  for (int i = 0; i < *chunkCols; ++i)
    {
    ostr << psScalar[i] << " ";
    }
  ostr << "\nt: ";
  for (int i = 0; i < 2*(*chunkCols); ++i)
    {
    ostr << tScalar[i] << " ";
    }
  ostr << endl;
  std::cerr << ostr.str();
}

void addToGrid(vtkUnstructuredGrid* grid2d, vtkUnstructuredGrid* grid3d, 
               double lonRad, double latRad, GridType gridType)
{
  // 2d grid
  double lonDeg = toDegrees(lonRad);
  int lonIndex = round (lonDeg / g_lonStep);          // interval [0, 360]
  double latDeg = toDegrees(latRad);
  int latIndex = round ((90 + latDeg) / g_latStep);   // interval [-90, 90]
  double
    lonMinus = lonDeg - g_lonStep / 2,
    lonPlus = lonDeg + g_lonStep / 2,
    latMinus = latDeg - g_latStep / 2,
    latPlus = latDeg + g_latStep / 2;
  double xMinus, xPlus, yMinus, yPlus, z;
  if (gridType == RECTILINEAR)
    {
    xMinus = lonMinus;
    xPlus = lonPlus;
    yMinus = latMinus;
    yPlus = latPlus;
    z = 0;
    }
  else
    {
    
    }
  int cornerIndex[4][2] = {{lonIndex, latIndex},
                           {lonIndex, latIndex + 1},
                           {lonIndex + 1, latIndex + 1},
                           {lonIndex + 1, latIndex}};
  double corner[4][2] = {{lonMinus, latMinus},
                         {lonMinus, latPlus},
                         {lonPlus, latPlus},
                         {lonPlus, latMinus}};
  vtkIdType cornerId[4];
  for (int p = 0; p < 4; ++p)
    {
    int pointIdIndex = 
      cornerIndex[p][0] + cornerIndex[p][1] * (g_dim[0] + 1);
    vtkIdType id = g_pointId2d[pointIdIndex];
    if (id < 0)
      {
      id = g_pointId2d[pointIdIndex] = 
        grid2d->GetPoints()->InsertNextPoint(corner[p][0], corner[p][1], 0);
      }
    cornerId[p] = id;
    }
  vtkIdType currentCell = grid2d->InsertNextCell(VTK_QUAD, 4, cornerId);
  // std::cerr << "cellId(" << lonIndex << ", " << latIndex << ") = " 
  //           << currentCell << endl;
  g_cellId2d[lonIndex + latIndex * g_dim[0]] = currentCell;
  // 3d grid
  for (int j = 0; j < g_dim[2]; ++j)
    {
    double levMinus, levPlus;
    if (j == 0)
      {
      double step = (g_level[1] - g_level[0]);
      levMinus = g_level[j] - step / 2;
      levPlus = g_level[j] + step / 2;
      }
    else if (j == g_dim[2] - 1)
      {
      double step = (g_level[g_dim[2] - 1] - g_level[g_dim[2] - 2]);
      levMinus = g_level[j] - step / 2;
      levPlus = g_level[j] + step / 2;
      }
    else
      {
      levMinus = g_level[j] - (g_level[j] - g_level[j-1]) / 2;
      levPlus = g_level[j] + (g_level[j+1] - g_level[j]) / 2;
      }
    int cornerIndex[8][3] = {{lonIndex,   latIndex,   j},
                             {lonIndex+1, latIndex,   j},
                             {lonIndex,   latIndex+1, j},
                             {lonIndex+1, latIndex+1, j},
                             {lonIndex,   latIndex,   j+1},
                             {lonIndex+1, latIndex,   j+1},
                             {lonIndex,   latIndex+1, j+1},
                             {lonIndex+1, latIndex+1, j+1}};

    double corner[8][3] = {{lonMinus, latMinus, levMinus},
                           {lonPlus, latMinus, levMinus},
                           {lonMinus, latPlus, levMinus},
                           {lonPlus, latPlus, levMinus},
                           {lonMinus, latMinus, levPlus},
                           {lonPlus, latMinus, levPlus},
                           {lonMinus, latPlus, levPlus},
                           {lonPlus, latPlus, levPlus}};

    vtkIdType cornerId[8];
    for (int p = 0; p < 8; ++p)
      {
      int pointIdIndex = 
        cornerIndex[p][0] + cornerIndex[p][1] * (g_dim[0] + 1) + 
        cornerIndex[p][2] * (g_dim[0] + 1) * (g_dim[1] + 1);
      vtkIdType id = g_pointId3d[pointIdIndex];
      if (id < 0)
        {
        id = g_pointId3d[pointIdIndex] =
          grid3d->GetPoints()->InsertNextPoint(
            corner[p][0], corner[p][1], corner[p][2]);
        // std::cerr << "InsertNextPoint(" 
        //           << corner[p][0] << ", " << corner[p][1] << ", "
        //           << corner[p][2] << ")" 
        //           << "\npointIdIndex=" << pointIdIndex << endl;
        }
      cornerId[p] = id;
      }
    currentCell = grid3d->InsertNextCell(VTK_VOXEL, 8, cornerId);
    g_cellId3d[lonIndex + latIndex * g_dim[0] + j * g_dim[0] * g_dim[1]] = 
      currentCell;
    }
}

};

extern "C" void cxx_create_grid_(
  int* dim, double* lonCoord, double* latCoord, double* levCoord,
  int* nCells2d, int* maxNcols,
  int* myRank)
{
  printCreateGrid(dim, lonCoord, latCoord, levCoord, nCells2d, maxNcols, myRank);

  g_rank = *myRank;
  g_chunkCapacity = *maxNcols;
  g_nCells2d = *nCells2d;
  g_lonStep = lonCoord[1] - lonCoord[0];
  g_latStep = latCoord[1] - latCoord[0];
  std::copy(dim, dim + 3, g_dim);
  if (!g_coprocessorData)
    {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
    }

  // 2d grid
  vtkSmartPointer<vtkUnstructuredGrid> grid2d =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points2d = vtkSmartPointer<vtkPoints>::New();
  grid2d->SetPoints(points2d);
  grid2d->GetCells()->Initialize();
  grid2d->Allocate(g_nCells2d);
  g_rankArrayIndex2d = addAttribute<vtkIntArray>(grid2d, "Rank", g_nCells2d);
  g_psArrayIndex = addAttribute<vtkDoubleArray>(grid2d, "PS", g_nCells2d);

  vtkCPInputDataDescription* idd = g_coprocessorData->
    GetInputDescriptionByName("input");
  if (!idd)
    {
    std::ostringstream ostr;
    ostr << "No vtkCPInputDataDescription for 'input'";
    vtkGenericWarningMacro(<< ostr.str());
    return;
    }
  idd->SetGrid(grid2d);
  g_cellId2d = new vtkIdType[dim[0] * dim[1]];
  for (int lat = 0; lat < dim[1]; ++lat)
    for(int lon = 0; lon < dim[0]; ++lon)
      {
      g_cellId2d[lon + lat * g_dim[0]] = -1;
      }
  g_pointId2d = new vtkIdType[(dim[0] + 1) * (dim[1] + 1)];
  for (int lat = 0; lat < dim[1] + 1; ++lat)
    for(int lon = 0; lon < dim[0] + 1; ++lon)
      {
      g_pointId2d[lon + lat * (g_dim[0] + 1)] = -1;
      }
  

  // 3D grid
  g_level = new double[g_dim[2]];
  std::copy(levCoord, levCoord + g_dim[2], g_level);
  vtkSmartPointer<vtkUnstructuredGrid> grid3d =
    vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> points3d = vtkSmartPointer<vtkPoints>::New();
  grid3d->SetPoints(points3d);
  grid3d->GetCells()->Initialize();
  grid3d->Allocate(g_nCells2d * g_dim[2]);
  g_rankArrayIndex3d = addAttribute<vtkIntArray>(
    grid3d, "Rank", g_nCells2d * g_dim[2]);
  g_tArrayIndex = addAttribute<vtkDoubleArray>(
    grid3d, "T", g_nCells2d * g_dim[2]);
  g_uArrayIndex = addAttribute<vtkDoubleArray>(
    grid3d, "U", g_nCells2d * g_dim[2]);
  g_vArrayIndex = addAttribute<vtkDoubleArray>(
    grid3d, "V", g_nCells2d * g_dim[2]);
  idd = g_coprocessorData->GetInputDescriptionByName("input3D");
  if (!idd)
    {
    std::ostringstream ostr;
    ostr << "No vtkCPInputDataDescription for 'input3D'";
    vtkGenericWarningMacro(<< ostr.str());
    return;
    }
  idd->SetGrid(grid3d);
  g_cellId3d = new vtkIdType[g_dim[0] * g_dim[1] * g_dim[2]];
  for (int lev = 0; lev < g_dim[2]; ++lev)
    for (int lat = 0; lat < g_dim[1]; ++lat)
      for (int lon= 0;lon < g_dim[0]; ++lon)
        {
        g_cellId3d[lon + lat*g_dim[0] + lev*g_dim[0]*g_dim[1]] = -1;
        }
  g_pointId3d = new vtkIdType[(g_dim[0] + 1) * (g_dim[1] + 1) * (g_dim[2] + 1)];
  for (int lev = 0; lev < g_dim[2] + 1; ++lev)
    for (int lat = 0; lat < g_dim[1] + 1; ++lat)
      for (int lon= 0; lon < g_dim[0] + 1; ++lon)
        {
        g_pointId3d[lon + lat*(g_dim[0]+1) + lev*(g_dim[0]+1)*(g_dim[1]+1)] = -1;
        }
}

extern "C" void cxx_add_chunk_(
  int* nstep, int* chunkSize,
  double* lonRad, double* latRad,
  double* psScalar, double *tScalar, double* uScalar, double* vScalar)
{
  // get the grids from Catalyst
  vtkCPInputDataDescription* idd = g_coprocessorData->
    GetInputDescriptionByName("input");
  if (!idd)
    {
    return;
    }
  vtkUnstructuredGrid* grid2d =
    vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());
  idd = g_coprocessorData->GetInputDescriptionByName("input3D");
  if (!idd)
    {
    return;
    }
  vtkUnstructuredGrid* grid3d =
    vtkUnstructuredGrid::SafeDownCast(idd->GetGrid());

  if (*nstep == 0)
    {
    // create points and cells
    for (int i = 0; i < *chunkSize; ++i)
      {
      addToGrid(grid2d, grid3d, lonRad[i], latRad[i], RECTILINEAR);
      }
    }

  printAddChunk(
    nstep, grid2d->GetNumberOfCells(), grid2d->GetPoints()->GetNumberOfPoints(), 
    grid2d->GetCellData()->GetNumberOfArrays(),
    chunkSize, lonRad, latRad, psScalar, tScalar);


  for (int i = 0; i < *chunkSize; ++i)
    {
    // 2d attributes
    int lonIndex = round (toDegrees(lonRad[i]) / g_lonStep);
    int latIndex = round((90 + toDegrees(latRad[i])) / g_latStep);
    vtkIdType cellId = g_cellId2d[lonIndex + latIndex * g_dim[0]];
    if (cellId == -1)
      {
      vtkGenericWarningMacro(<< "Invalid 2D cell at: "
                           << lonIndex << ", " << latIndex << endl);
      exit(13);
      }
    vtkDoubleArray::SafeDownCast(
      grid2d->GetCellData()->GetArray(g_psArrayIndex))
      ->SetValue(cellId, psScalar[i]);
    vtkIntArray::SafeDownCast(
      grid2d->GetCellData()->GetArray(g_rankArrayIndex2d))
      ->SetValue(cellId, g_rank);

    // 3d attributes
    for (int j = 0; j < g_dim[2]; ++j)
      {
      cellId = g_cellId3d[
        lonIndex + latIndex * g_dim[0] + j * g_dim[0] * g_dim[1]];
      if (cellId == -1)
        {
        vtkGenericWarningMacro(<< "Invalid 3D cell at: "
                               << lonIndex << ", "
                               << latIndex << ", "
                               << j << endl);
        exit(13);
        }
      vtkIntArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_rankArrayIndex3d))
        ->SetValue(cellId, g_rank);
      vtkDoubleArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_tArrayIndex))
        ->SetValue(cellId, tScalar[i + j*g_chunkCapacity]);
      vtkDoubleArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_uArrayIndex))
        ->SetValue(cellId, uScalar[i + j*g_chunkCapacity]);
      vtkDoubleArray::SafeDownCast(
        grid3d->GetCellData()->GetArray(g_vArrayIndex))
        ->SetValue(cellId, vScalar[i + j*g_chunkCapacity]);
      }
    }
}

extern "C" void cxx_finalize_()
{
  if (g_level)
    {
    delete[] g_level;
    }
  if (g_cellId2d)
    {
    delete[] g_cellId2d;
    }
  if (g_cellId3d)
    {
    delete[] g_cellId3d;
    }
  if (g_pointId2d)
    {
    delete[] g_pointId2d;
    }
  if (g_pointId3d)
    {
    delete[] g_pointId3d;
    }
}


// make sure you pass a zero terminated string
extern "C" void cxx_coprocessorinitializewithpython_(const char* pythonScriptName)
{
  if (!g_coprocessor)
    {
    g_coprocessor = vtkCPProcessor::New();
    g_coprocessor->Initialize();
    // python pipeline
    vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline =
      vtkSmartPointer<vtkCPPythonScriptPipeline>::New();
    pipeline->Initialize(pythonScriptName);
    g_coprocessor->AddPipeline(pipeline);
    }
  if (!g_coprocessorData)
    {
    g_coprocessorData = vtkCPDataDescription::New();
    g_coprocessorData->AddInput("input");
    g_coprocessorData->AddInput("input3D");
    }
}

extern "C" void cxx_coprocessorfinalize_()
{
  if (g_coprocessor)
    {
    g_coprocessor->Delete();
    g_coprocessor = NULL;
    }
  if (g_coprocessorData)
    {
    g_coprocessorData->Delete();
    g_coprocessorData = NULL;
    }
}

extern "C" void cxx_requestdatadescription_(int* timeStep, double* time,
                                            int* coprocessThisTimeStep)
{
  if(!g_coprocessorData || !g_coprocessor)
    {
    vtkGenericWarningMacro("Problem in needtocoprocessthistimestep."
                           << "Probably need to initialize.");
    *coprocessThisTimeStep = 0;
    return;
    }
  vtkIdType tStep = *timeStep;
  g_coprocessorData->SetTimeData(*time, tStep);
  if(g_coprocessor->RequestDataDescription(g_coprocessorData))
    {
    *coprocessThisTimeStep = 1;
    g_isTimeDataSet = true;
    }
  else
    {
    *coprocessThisTimeStep = 0;
    g_isTimeDataSet = false;
    }
}

extern "C" void cxx_needtocreategrid_(int* needGrid)
{
  if(!g_isTimeDataSet)
    {
    vtkGenericWarningMacro("Time data not set.");
    *needGrid = 0;
    return;
    }

  // assume that the grid is not changing so that we only build it
  // the first time, otherwise we clear out the field data
  if(g_coprocessorData->GetInputDescriptionByName("input")->GetGrid())
    {
    *needGrid = 0;
    }
  else
    {
    *needGrid = 1;
    }
}

//-----------------------------------------------------------------------------
extern "C" void cxx_coprocess_()
{
  if(!g_isTimeDataSet)
    {
    vtkGenericWarningMacro("Time data not set.");
    }
  else
    {
    g_coprocessor->CoProcess(g_coprocessorData);
    }
  // Reset time data.
  g_isTimeDataSet = false;
}
