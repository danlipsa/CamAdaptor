// Adaptor for getting Fortran simulation code into ParaView CoProcessor.

// CoProcessor specific headers
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
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

// coprocessor data
vtkCPProcessor* g_coprocessor;             // catalyst coprocessor
vtkCPDataDescription* g_coprocessorData;   // input, sinput, input3D, sinput3D
bool g_isTimeDataSet;                      // is time data set?

// shared between grids
int g_rank = -1;
int g_chunkCapacity;  // maximum number of (vertical) columns in a chunk
int g_nCells2d;       // total number of 2D cells on a MPI processor
int g_dim[3];         // lon x lat x lev
double g_lonStep;     // longitude step in degrees
double g_latStep;     // latitude step in degrees
double* g_level;      // level values

// rectilinear (2D, 3D) and structured (spherical) (2D, 3Da) grids
class Grid;
Grid* g_grid;
Grid* g_sgrid;
};

namespace
{
double toDegrees(double rad)
{
  return rad / PI * 180.0;
}

double toRadians(double deg)
{
  return deg * PI / 180.0;
}

void rotateAroundYDeg(double degQ, double v[3], double r[3])
{
  double q = toRadians(degQ);
  r[0] = v[2]*sin(q) + v[0]*cos(q);
  r[1] = v[1];
  r[2] = v[2]*cos(q) - v[0]*sin(q);
}

void rotateAroundZDeg(double degQ, double v[3], double r[3])
{
  double q = toRadians(degQ);
  r[0] = v[0]*cos(q) - v[1]*sin(q);
  r[1] = v[0]*sin(q) + v[1]*cos(q);
  r[2] = v[2];
}

void sphericalToCartesian(double corner[3])
{
  double rotLon = corner[0];
  double rotLat = corner[1];
  double src[3] = {corner[2], 0, 0};
  double r[3];
  rotateAroundYDeg(-rotLat, src, r);
  rotateAroundZDeg(rotLon, r, corner);
  // std::cerr << "rotLon=" << rotLon
  //           << " rotLat=" << rotLat
  //           << " x=" << src[0]
  //           << " xp=" << corner[0]
  //           << " yp=" << corner[1]
  //           << " zp=" << corner[2] << endl;
}

// Creates and accumulates data for a 2D and 3D grid. It can generate 
// either a rectilinear or a spherical grid.
class Grid
{
public:
  enum Type
  {
    RECTILINEAR,
    SPHERE
  };

  // Creates a 2D and a 3D grid of the specified 'type'
  Grid(Type type, const char* name2d, const char* name3d) : 
    RankArrayIndex2d(-1), CoordArrayIndex2d(-1),
    PSArrayIndex(-1), CellId2d(NULL), PointId2d(NULL),
    
    RankArrayIndex3d(-1), CoordArrayIndex3d(-1),
    TArrayIndex(-1), UArrayIndex(-1), VArrayIndex(-1),
    CellId3d(NULL), PointId3d(NULL), GridType(type)
  {
    this->Name2d = name2d;
    this->Name3d = name3d;
    this->Grid2d = CreateGrid2d(name2d);
    this->Grid3d = CreateGrid3d(name3d);
  }

  // Deletes data used to build the grids. Note that the grid memory is managed by
  // the Catalyst Coprocessor.
  ~Grid()
  {
  if (this->CellId2d)
    {
    delete[] this->CellId2d;
    }
  if (this->CellId3d)
    {
    delete[] this->CellId3d;
    }
  if (this->PointId2d)
    {
    delete[] this->PointId2d;
    }
  if (this->PointId3d)
    {
    delete[] this->PointId3d;
    }
  }

  // Creates the 2D grid and the data used to add attributes to the grid
  vtkUnstructuredGrid* CreateGrid2d(const char* name)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid2d =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points2d = vtkSmartPointer<vtkPoints>::New();
    grid2d->SetPoints(points2d);
    grid2d->GetCells()->Initialize();
    grid2d->Allocate(g_nCells2d);
    this->RankArrayIndex2d = 
      addAttribute<vtkIntArray>(grid2d, "Rank", g_nCells2d, 1);
    this->CoordArrayIndex2d = addAttribute<vtkDoubleArray>(
      grid2d, "Coord", g_nCells2d, 2);
    this->PSArrayIndex = 
      addAttribute<vtkDoubleArray>(grid2d, "PS", g_nCells2d, 1);  

    vtkCPInputDataDescription* idd = g_coprocessorData->
      GetInputDescriptionByName(name);
    if (!idd)
      {
      std::ostringstream ostr;
      ostr << "No vtkCPInputDataDescription for " << name;
      throw std::runtime_error(ostr.str());
      }
    idd->SetGrid(grid2d);
    this->CellId2d = new vtkIdType[g_dim[0] * g_dim[1]];
    for (int lat = 0; lat < g_dim[1]; ++lat)
      for(int lon = 0; lon < g_dim[0]; ++lon)
        {
        this->CellId2d[lon + lat * g_dim[0]] = -1;
        }
    this->PointId2d = new vtkIdType[(g_dim[0] + 1) * (g_dim[1] + 1)];
    for (int lat = 0; lat < g_dim[1] + 1; ++lat)
      for(int lon = 0; lon < g_dim[0] + 1; ++lon)
        {
        this->PointId2d[lon + lat * (g_dim[0] + 1)] = -1;
        }
    // the pointer is managed by vtkCPInputDataDescription
    return grid2d;
  }

  // Creates the 3D grid and the data used to add attributes to the grid
  vtkUnstructuredGrid* CreateGrid3d(const char* name)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid3d =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points3d = vtkSmartPointer<vtkPoints>::New();
    grid3d->SetPoints(points3d);
    grid3d->GetCells()->Initialize();
    grid3d->Allocate(g_nCells2d * g_dim[2]);
    this->RankArrayIndex3d = addAttribute<vtkIntArray>(
      grid3d, "Rank", g_nCells2d * g_dim[2], 1);
    this->CoordArrayIndex3d = addAttribute<vtkDoubleArray>(
      grid3d, "Coord", g_nCells2d * g_dim[2], 3);
    this->TArrayIndex = addAttribute<vtkDoubleArray>(
      grid3d, "T", g_nCells2d * g_dim[2], 1);
    this->UArrayIndex = addAttribute<vtkDoubleArray>(
      grid3d, "U", g_nCells2d * g_dim[2], 1);
    this->VArrayIndex = addAttribute<vtkDoubleArray>(
      grid3d, "V", g_nCells2d * g_dim[2], 1);
    vtkCPInputDataDescription* idd = 
      g_coprocessorData->GetInputDescriptionByName(name);
    if (!idd)
      {
      std::ostringstream ostr;
      ostr << "No vtkCPInputDataDescription for " << name;
      throw std::runtime_error(ostr.str());
      }
    idd->SetGrid(grid3d);
    this->CellId3d = new vtkIdType[g_dim[0] * g_dim[1] * g_dim[2]];
    for (int lev = 0; lev < g_dim[2]; ++lev)
      for (int lat = 0; lat < g_dim[1]; ++lat)
        for (int lon= 0;lon < g_dim[0]; ++lon)
          {
          this->CellId3d[lon + lat*g_dim[0] + lev*g_dim[0]*g_dim[1]] = -1;
          }
    this->PointId3d = new vtkIdType[
      (g_dim[0] + 1) * (g_dim[1] + 1) * (g_dim[2] + 1)];
    for (int lev = 0; lev < g_dim[2] + 1; ++lev)
      for (int lat = 0; lat < g_dim[1] + 1; ++lat)
        for (int lon= 0; lon < g_dim[0] + 1; ++lon)
          {
          this->PointId3d[lon + lat*(g_dim[0]+1) + 
                          lev*(g_dim[0]+1)*(g_dim[1]+1)] = -1;
          }
    return grid3d;
  }

  // Adds the points and the cells for a vertical column to the grid 
  void AddPointsAndCells(double lonRad, double latRad)
  {
    if (! this)
      return;
    // 2d grid
    double lonDeg = toDegrees(lonRad);
    int lonIndex = round (lonDeg / g_lonStep);          // interval [0, 360]
    double latDeg = toDegrees(latRad);
    int latIndex = round ((90 + latDeg) / g_latStep);   // interval [-90, 90]
    double
      lonMinus = lonDeg - g_lonStep / 2,
      lonPlus = lonDeg + g_lonStep / 2,
      latMinus,
      latPlus;
    latitudeMinusPlus(latIndex, &latMinus, &latPlus);
    int cornerIndex[4][2] = {{lonIndex, latIndex},
                             {lonIndex, latIndex + 1},
                             {lonIndex + 1, latIndex + 1},
                             {lonIndex + 1, latIndex}};
    double corner[4][3] = {{lonMinus, latMinus, 1},
                           {lonMinus, latPlus, 1},
                           {lonPlus, latPlus, 1},
                           {lonPlus, latMinus, 1}};
    if (this->GridType == SPHERE)
      {
      for (int i = 0; i < 4; ++i)
        {
        sphericalToCartesian(corner[i]);
        }
      }
    vtkIdType cornerId[4];
    for (int p = 0; p < 4; ++p)
      {
      int pointIdIndex = 
        cornerIndex[p][0] + cornerIndex[p][1] * (g_dim[0] + 1);
      vtkIdType id = this->PointId2d[pointIdIndex];
      if (id < 0)
        {
        id = this->PointId2d[pointIdIndex] = 
          this->Grid2d->GetPoints()->InsertNextPoint(corner[p][0], corner[p][1], 
                                                     corner[p][2]);
        }
      cornerId[p] = id;
      }
    vtkIdType currentCell = this->Grid2d->InsertNextCell(VTK_QUAD, 4, cornerId);
    // std::cerr << "cellId(" << lonIndex << ", " << latIndex << ") = " 
    //           << currentCell << endl;
    this->CellId2d[lonIndex + latIndex * g_dim[0]] = currentCell;
    // 3d grid
    for (int j = 0; j < g_dim[2]; ++j)
      {
      double levMinus, levPlus;
      levelMinusPlus(j, &levMinus, &levPlus);
      int cornerIndex[8][3] = {{lonIndex,   latIndex,   j},
                               {lonIndex+1, latIndex,   j},
                               {lonIndex+1, latIndex+1, j},
                               {lonIndex,   latIndex+1, j},
                               {lonIndex,   latIndex,   j+1},
                               {lonIndex+1, latIndex,   j+1},
                               {lonIndex+1, latIndex+1, j+1},
                               {lonIndex,   latIndex+1, j+1}};
      double corner[8][3] = {{lonMinus, latMinus, this->levelToRadius(levMinus)},
                             {lonPlus, latMinus, this->levelToRadius(levMinus)},
                             {lonPlus, latPlus, this->levelToRadius(levMinus)},
                             {lonMinus, latPlus, this->levelToRadius(levMinus)},
                             {lonMinus, latMinus, this->levelToRadius(levPlus)},
                             {lonPlus, latMinus, this->levelToRadius(levPlus)},
                             {lonPlus, latPlus, this->levelToRadius(levPlus)},
                             {lonMinus, latPlus, this->levelToRadius(levPlus)}};
      if (this->GridType == SPHERE)
        {
        for (int i = 0; i < 8; ++i)
          {
          sphericalToCartesian(corner[i]);
          }
        }
      vtkIdType cornerId[8];
      for (int p = 0; p < 8; ++p)
        {
        int pointIdIndex = 
          cornerIndex[p][0] + cornerIndex[p][1] * (g_dim[0] + 1) + 
          cornerIndex[p][2] * (g_dim[0] + 1) * (g_dim[1] + 1);
        vtkIdType id = this->PointId3d[pointIdIndex];
        if (id < 0)
          {
          id = this->PointId3d[pointIdIndex] =
            this->Grid3d->GetPoints()->InsertNextPoint(
              corner[p][0], corner[p][1], corner[p][2]);
          // std::cerr << "InsertNextPoint(" 
          //           << corner[p][0] << ", " << corner[p][1] << ", "
          //           << corner[p][2] << ")" 
          //           << "\npointIdIndex=" << pointIdIndex << endl;
          }
        cornerId[p] = id;
        }
      currentCell = this->Grid3d->InsertNextCell(VTK_HEXAHEDRON, 8, cornerId);
      this->CellId3d[lonIndex + latIndex * g_dim[0] + j * g_dim[0] * g_dim[1]] = 
        currentCell;
      }
  }

  // Sets attributes for a chunk (a list of vertical columns) to 
  // the 2D and 3D grids
  void SetAttributeValue(int chunkSize,
                         double* lonRad, double* latRad,
                         double* psScalar, double *tScalar, double* uScalar, 
                         double* vScalar)
  {
    for (int i = 0; i < chunkSize; ++i)
      {
      // 2d attributes
      double lonDeg = toDegrees(lonRad[i]);
      int lonIndex = round (lonDeg / g_lonStep);          // interval [0, 360]
      double latDeg = toDegrees(latRad[i]);
      int latIndex = round ((90 + latDeg) / g_latStep);   // interval [-90, 90]
      vtkIdType cellId = this->CellId2d[lonIndex + latIndex * g_dim[0]];
      if (cellId == -1)
        {
        vtkGenericWarningMacro(<< "Invalid 2D cell at: "
                               << lonIndex << ", " << latIndex << endl);
        exit(13);
        }
      vtkDoubleArray::SafeDownCast(
        this->Grid2d->GetCellData()->GetArray(this->PSArrayIndex))
        ->SetValue(cellId, psScalar[i]);
      this->Grid2d->GetCellData()->GetArray(this->CoordArrayIndex2d)
        ->SetTuple2(cellId, lonDeg, latDeg);
      vtkIntArray::SafeDownCast(
        this->Grid2d->GetCellData()->GetArray(this->RankArrayIndex2d))
        ->SetValue(cellId, g_rank);

      // 3d attributes
      for (int j = 0; j < g_dim[2]; ++j)
        {
        cellId = this->CellId3d[
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
          this->Grid3d->GetCellData()->GetArray(this->RankArrayIndex3d))
          ->SetValue(cellId, g_rank);
        vtkDoubleArray::SafeDownCast(
          this->Grid3d->GetCellData()->GetArray(this->TArrayIndex))
          ->SetValue(cellId, tScalar[i + j*g_chunkCapacity]);
        vtkDoubleArray::SafeDownCast(
          this->Grid3d->GetCellData()->GetArray(this->UArrayIndex))
          ->SetValue(cellId, uScalar[i + j*g_chunkCapacity]);
        vtkDoubleArray::SafeDownCast(
          this->Grid3d->GetCellData()->GetArray(this->VArrayIndex))
          ->SetValue(cellId, vScalar[i + j*g_chunkCapacity]);
        this->Grid3d->GetCellData()->GetArray(this->CoordArrayIndex3d)
          ->SetTuple3(cellId, lonDeg, latDeg, g_level[j]);
        }
      }
  }

  // Print routine for debugging
  void PrintAddChunk(
    int* nstep, int* chunkCols,
    double* lonRad, double* latRad, double* psScalar, double *tScalar)
  {
    std::ostringstream ostr;
    ostr << "add_chunk: " << *nstep
         << "\nchunkCols: " << *chunkCols
         << "\nnCells: " << this->Grid2d->GetNumberOfCells()
         << "\nnPoints: " << this->Grid2d->GetPoints()->GetNumberOfPoints()
         << "\nnCellArrays: " << this->Grid2d->GetCellData()->GetNumberOfArrays()
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

  // Returns the name for the 2D grid
  const std::string& GetName2d() const
  {
    return this->Name2d;
  }

  // Returns the name for the 3D grid
  const std::string& GetName3d() const
  {
    return this->Name3d;
  }

  // Returns the 2D grid
  vtkUnstructuredGrid* GetGrid2d() const
  {
    return this->Grid2d;
  }

  // Returns the 3D grid
  vtkUnstructuredGrid* GetGrid3d() const
  {
    return this->Grid3d;
  }

private:
  // Creates an array for attribute 'name' and adds it to 'grid'.
  template<typename T> static
  int addAttribute(vtkUnstructuredGrid* grid, const char* name,
                   vtkIdType size, int nComponents)
  {
    vtkSmartPointer<T> a = vtkSmartPointer<T>::New();
    a->SetNumberOfComponents(nComponents);
    a->SetNumberOfTuples(size);
    a->SetName(name);
    return grid->GetCellData()->AddArray(a);
  }
  // The level for a grid measures preasure which decreases with height.
  // This transformation sets 
  // - max pressure to 0 (surface level) and min pressure to 1 (SPHERE)
  // - negates the pressure so that globe surface (high pressure) is closer to 
  //   origin than high in the atmosphere (low pressure) (RECTILINEAR)
  // These transformations match the transformation in the history file.
  double levelToRadius(double level)
  {
    double maxLevel = g_level[g_dim[2] - 1];
    if (this->GridType == SPHERE)
      {
      return  (maxLevel - level) / maxLevel;
      }
    else
      {
      return - level;
      }
  }
  // Compute the location of the points surounding a cell at index 'j' for level
  static double levelMinusPlus(int j, double* levMinus, double* levPlus)
  {
    if (j == 0)
      {
      double step = (g_level[1] - g_level[0]);
      *levMinus = g_level[j];
      *levPlus = g_level[j] + step / 2;
      }
    else if (j == g_dim[2] - 1)
      {
      double step = (g_level[g_dim[2] - 1] - g_level[g_dim[2] - 2]);
      *levMinus = g_level[j] - step / 2;
      *levPlus = g_level[j];
      }
    else
      {
      *levMinus = g_level[j] - (g_level[j] - g_level[j-1]) / 2;
      *levPlus = g_level[j] + (g_level[j+1] - g_level[j]) / 2;
      }
  }
  // Compute the location of the points surouding a cell at index 'latIndex'
  // for latitude
  static double latitudeMinusPlus(int latIndex, double* latMinus, double* latPlus)
  {
    double latDeg = -90 + g_latStep * latIndex;
    if (latIndex == 0)
      {
      *latMinus = latDeg;
      *latPlus = latDeg + g_latStep / 2;
      }
    else if (latIndex == g_dim[1] - 1)
      {
      *latMinus = latDeg - g_latStep / 2;
      *latPlus = latDeg;
      }
    else
      {
      *latMinus = latDeg - g_latStep / 2;
      *latPlus = latDeg + g_latStep / 2;
      }
  }

private:
  Type GridType;

  // 2d grid
  std::string Name2d;         //grid name
  vtkUnstructuredGrid* Grid2d;//the grid
  int RankArrayIndex2d;       //index of the rank array
  int CoordArrayIndex2d;      //index of the coordinates array (lon x lat)
  int PSArrayIndex;           // index of the PS array
  vtkIdType* CellId2d;        // 2d array with cellId at lon x lat
  vtkIdType* PointId2d;       // 2d array with pointId at lon x lat

  // 3d grid
  std::string Name3d;         // grid name
  vtkUnstructuredGrid* Grid3d;// the grid
  int RankArrayIndex3d;       // index of the rank array
  int CoordArrayIndex3d;      // index of the coordinates array (lon x lat x lev)
  int TArrayIndex;            // index of the attribute array for T
  int UArrayIndex;            // index of the attribute array for U
  int VArrayIndex;            // index of the attribute array for V
  vtkIdType* CellId3d;        // 3d array with the cellId at lon x lat x lev
  vtkIdType* PointId3d;       // 3d array with the pointId at lon x lat x lev.
};

// Debugging printouts
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
};

// Creates the Grids for 2D, 3D rectilinear and 2D, 3D spherical
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
  g_level = new double[g_dim[2]];
  std::copy(levCoord, levCoord + g_dim[2], g_level);

  if (!g_coprocessorData)
    {
    vtkGenericWarningMacro("Unable to access CoProcessorData.");
    return;
    }
    try
      {
      g_grid = new Grid(Grid::RECTILINEAR, "input", "input3D");
      }
    catch(std::exception& e)
      {
      vtkGenericWarningMacro(<< e.what());
      delete g_grid;
      g_grid = NULL;
      }
    try
      {
      g_sgrid = new Grid(Grid::SPHERE, "sinput", "sinput3D");
      }
    catch(std::exception& e)
      {
      vtkGenericWarningMacro(<< e.what());
      delete g_sgrid;
      g_sgrid = NULL;
      }
}

// for timestep 0: creates the points and cells for the grids.
// for all timesteps: copies data from the simulation to Catalyst.
extern "C" void cxx_add_chunk_(
  int* nstep, int* chunkSize,
  double* lonRad, double* latRad,
  double* psScalar, double *tScalar, double* uScalar, double* vScalar)
{
  if (*nstep == 0)
    {
    for (int i = 0; i < *chunkSize; ++i)
      {
      if (g_grid)
        {
        g_grid->AddPointsAndCells(lonRad[i], latRad[i]);
        }
      if (g_sgrid)
        {
        g_sgrid->AddPointsAndCells(lonRad[i], latRad[i]);
        }
      }
    }
  if (g_grid)
    {
    g_grid->PrintAddChunk(
      nstep, chunkSize, lonRad, latRad, psScalar, tScalar);
    g_grid->SetAttributeValue(*chunkSize, lonRad, latRad,
                              psScalar, tScalar, uScalar, vScalar);
    }
  if (g_sgrid)
    {
    g_sgrid->SetAttributeValue(*chunkSize, lonRad, latRad,
                               psScalar, tScalar, uScalar, vScalar);
    }
}

// Deletes global data
extern "C" void cxx_finalize_()
{
  if (g_level)
    {
    delete[] g_level;
    }
  if (g_grid)
    {
    delete g_grid;
    }
  if (g_sgrid)
    {
    delete g_sgrid;
    }
}


// Initializes the Catalyst Coprocessor
// WARNING: Make sure you pass a zero terminated string
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
    g_coprocessorData->AddInput("sinput");
    g_coprocessorData->AddInput("sinput3D");
    }
}

// Deletes the Catalyt Coprocessor and data
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

// Checks if Catalyst needs to coprocess data
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

// Checks if the grids need to be created
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
  vtkCPInputDataDescription* idd = 
    g_coprocessorData->GetInputDescriptionByName("input");
  if(idd == NULL || idd->GetGrid() == NULL)
    {
    *needGrid = 1;
    }
  else
    {
    *needGrid = 0;
    }
}

// calls the coprocessor
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
