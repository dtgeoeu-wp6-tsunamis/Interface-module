import vtk
from vtk.util import numpy_support
import numpy as np
import seissolxdmf
import time
import argparse
from pyproj import Transformer

# TODO: potentially set up inputCRS correctly

"""
Module for the SeisSol donor functionalities.

Contains the following functionalities:

* seissolxdmf                                 class definition
* project_coordinates                     project coordinates to WGS84 coordinates
* setUp_grid_interpolation             set up grid and interpolation structures 
* get_interpolation                         perform the interpolation on the data
* interpolate_Seissol2structured    main routine
* get_seissol                                   parent routine that is called from the main donorModel routine

"""

class seissolxdmfExtended(seissolxdmf.seissolxdmf):
  def generateVtkObject(self):
    """Filling in vtk arrays with data from hdf5 file."""

    connect = self.ReadConnect()
    nElements, ndim2 = connect.shape

    xyz = self.ReadGeometry()
    points = vtk.vtkPoints()
    if ndim2 == 3:
        print("surface output, assuming the grid is at z=0")
        xyz[:, 2] = 0.0
    points.SetData(numpy_support.numpy_to_vtk(xyz))

    vtkCells = vtk.vtkCellArray()
    connect2 = np.zeros((nElements, ndim2 + 1), dtype=np.int64)
    # number of points in the cell
    connect2[:, 0] = ndim2
    connect2[:, 1:] = connect
    vtkCells.SetCells(nElements, numpy_support.numpy_to_vtkIdTypeArray(connect2))

    if ndim2 == 4:
        unstrGrid3d = vtk.vtkUnstructuredGrid()
        unstrGrid3d.SetPoints(points)
        unstrGrid3d.SetCells(vtk.VTK_TETRA, vtkCells)
        return unstrGrid3d
    elif ndim2 == 3:
        myPolydata = vtk.vtkPolyData()
        myPolydata.SetPoints(points)
        myPolydata.SetPolys(vtkCells)
        return myPolydata
    else:
        raise NotImplementedError
  
  
  
def project_coordinates(x, y, inputcrs):
  """
  This function takes the input x and y coordinates and returns the transformed WGS84 coordinates

  :param x: x-coordinates
  :param y: y-coordinates
  :param inputcrs: input CRS
  """
  transformer = Transformer.from_crs(inputcrs, "epsg:4326", always_xy=True)

    # transform x-coordinates
  y_lowerRow = np.repeat(y[0], len(x))
  xnew, ynew = transformer.transform(x, y_lowerRow)
  x_proj = xnew  
    
    # transform y-coordinates
  x_leftColumn = np.repeat(x[0], len(y))
  xnew, ynew = transformer.transform(x_leftColumn, y)
  y_proj = ynew

  return x_proj, y_proj
  
  
  
def setUp_grid_interpolation(coord_min, coord_max, dx, inputCRS):
  """
  Sets up the grid (image) for interpolation using VTK and probe filter. Returns the probe filter and shape for reshaping (needed within the interpolation)
  
  :param coord_min: minimum coordinates (bottom left corner)
  :param coord_min: maximum coordinates (top right corner)
  :param dx:  spatial resolution 
  :param inputCRS:  CRS of the input 2d mesh 
  """
  
  # set up x and y coordinates
  x = np.arange(coord_min[0], coord_max[0] + dx, dx)
  y = np.arange(coord_min[1], coord_max[1] + dx, dx)

  z = np.array([0])   # ensure that the mesh is 2D
  xx, yy = np.meshgrid(x, y)

  # project the x and y coordinates to lat/lon
  transformer = Transformer.from_crs(inputCRS, "epsg:4326", always_xy=True)
  x_proj, y_proj = project_coordinates(x, y, inputCRS)

  # Create grid image volume
  imageSize = [x.shape[0], y.shape[0], z.shape[0]]
  imageOrigin = [coord_min[0], coord_min[1], coord_min[2]]
  imageSpacing = [dx, dx, dx]
  
  imageData = vtk.vtkImageData()
  imageData.SetDimensions(imageSize)
  imageData.SetOrigin(imageOrigin)
  imageData.SetSpacing(imageSpacing)

  # Create the interpolation filter
  probeFilter = vtk.vtkProbeFilter()

  probeFilter.SetInputData(imageData)
  probeFilter.SpatialMatchOn()
  
  return probeFilter, xx.shape, x_proj, y_proj



def get_interpolation(sx, unstrGrid3d, probeFilter, projDataShape, timestep, varName):
  """
  Routine that calculates the interpolation via probe filter.

  :param sx: seissolxdmf file   
  :param unstrGrid3d: unstructed grid object (derived from XDMF file)
  :param probeFilter:  probe filter object
  ;param projDataShape: shape that the data has to be reshaped to after interpolation
  :param time: time index
  :param varName: name of variable to be interpolated
  """
  
  # Set up data structure
  scalars = vtk.vtkFloatArray()

  W = sx.ReadData(varName, timestep)
  scalars = numpy_support.numpy_to_vtk(
      num_array=W, deep=True, array_type=vtk.VTK_FLOAT
  )
  unstrGrid3d.GetCellData().SetScalars(scalars)
  
  # Create the CellDataToPointData filter
  cellToPointFilter = vtk.vtkCellDataToPointData()
  cellToPointFilter.SetInputData(unstrGrid3d)
  cellToPointFilter.Update()

  # Get the output grid with point data
  outputGrid = cellToPointFilter.GetOutput()

  # Perform the interpolation
  probeFilter.SetSourceData(outputGrid)
  start = time.time()
  probeFilter.Update()
  stop = time.time()
  print(f"{varName} {timestep}: done probe filter in {stop - start} s")

  polyout = probeFilter.GetOutput()
  projData = polyout.GetPointData().GetScalars()
  projDataNp = numpy_support.vtk_to_numpy(projData).reshape(projDataShape)
  
  return projDataNp
  
  
  
def interpolate_Seissol2structured(sx, dx, coord_min, coord_max, include_horizontal, instants=[]):
  """
  Interpolate the SeisSol XDMF to VTK
  :param sx: seissolxdmf file 
  :param dx: spatial resolution
  :param coord_min: minimum coordinates for box
  :param coord_max: maximum coordinates for box
  :param include_horizontal: handle whether to interpolate only the vertical component or not
  :param instants: time steps to include in the out netCDF. If not provided, all time steps will be included
  
  returns deformation data and coordinates
  """
  unstrGrid3d = sx.generateVtkObject()

  nTime = sx.ReadNdt()  # number of time steps in the Seissol file
  if not instants:
    # if no specific instant is provided, use all in the seissol output
    instants = list(range(0, nTime))
  else:
    # check if instants provided are within the possible ones
    if not all(x in range(0, nTime) for x in instants):
      print("Instants provided are outside the range of timesteps in SeisSol's output")
  
  # Choose which data interpolate (only vertical or all components) based on the input handle  
  if (include_horizontal):
    data = ['u1', 'u2', 'u3']
  else:
    data = ['u3']    
  
  # Set up the inputCRS   
  inputCRS = "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=26.25 +lat_0=37.75"  # CRS of the input 2d mesh

  # Create probe filter and get projected coordinates
  probeFilter, projDataShape, x_proj, y_proj = setUp_grid_interpolation(coord_min, coord_max, dx, inputCRS)

  # List for interpolated data
  probedData = []

  # Perform interpolation (over each timestep and each variable)
  print("Interpolation is performed using VTK probe filter")
  for timestep in instants:
      for varName in data:
          projDataNp = get_interpolation(sx, unstrGrid3d, probeFilter, projDataShape, timestep, varName)
          probedData.append(projDataNp)
 
  return probedData, x_proj, y_proj



def get_seissol(filename, spatial_resolution, include_horizontal):
  """
  Actual part that will be used by the main donor model functionality to get the deformation data
  :param filename: filename for the SeisSol data (has to be an XDMF file)
  :param spatial_resolution: spatial_resolution in meters
  :param include_horizontal: boolean handle whether the output will only contain the vertical deformation (deformation) or also include the horizontal deformation
  """

  sx = seissolxdmfExtended(filename) # get seissolxdmf from provided XDMF file
# get x, y and z interval min/max and round to nearest 1000 (if not done, the interpolation does nothing)
  geom = sx.ReadGeometry()
  coordinate_min = np.round(geom.min(0) +  spatial_resolution, -4) 
  coordinate_max = np.round(geom.max(0) -  spatial_resolution, -4) 

  donor_deformation, donor_x, donor_y = interpolate_Seissol2structured(sx, spatial_resolution, coordinate_min, coordinate_max, include_horizontal)

  return donor_deformation, donor_x, donor_y
