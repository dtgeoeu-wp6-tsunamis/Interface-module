import numpy as np
import os
import time
import seissolxdmf
import vtk

from pyproj import Transformer
from vtk.util import numpy_support

#TODO Lon/Lat for the bottom left coordinates of the domain should be provided by some way within the SeisSol output. Will be used for inputCRS. 

"""
Module for the SeisSol donor functionalities.

Contains the following functionalities:

* seissolxdmf                                 class definition
* get_seissol_time                          function to get time values of SeisSol output
* project_coordinates                     project coordinates to WGS84 coordinates
* setUp_grid_interpolation             set up grid and interpolation structures 
* get_interpolation                         perform the interpolation on the data
* interpolate_seissol2structured    main routine
* get_seissol                                  parent routine that is called from the main donorModel routine
"""

# Global parameters
basicCRS = 'epsg:4326' # basic lat-lon coordinate system

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


class seissolxdmfExtended(seissolxdmf.seissolxdmf):
  def generateVtkObject(self):
    """Filling in vtk arrays with data from hdf5 file."""

    connect = self.ReadConnect()
    nElements, ndim2 = connect.shape

    xyz = self.ReadGeometry()
    points = vtk.vtkPoints()
    if ndim2 == 3:
        print("Surface output, assuming the grid is at z=0".center(column_size))
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
  
  
  
def get_seissol_time(sx):
  """
  Function to read the time data from the seissolxdmf file
  
  :param sx: seissolxdmf file
  """
  
  time_array = np.zeros(sx.ReadNdt())
  root = sx.tree.getroot()
  i = 0
  for Property in root.findall("Domain/Grid/Grid/Time"):
    if i == 0:
      time_array[i] = float(Property.get("Value"))
      i = 1
    else:
      time_array[i] = float(Property.get("Value")) 
      i += 1
  return time_array



def project_coordinates(x, y, inputCRS):
  """
  This function takes the input x and y coordinates and returns the transformed WGS84 coordinates.

  :param x: x-coordinates
  :param y: y-coordinates
  :param inputCRS: input CRS
  """
  transformer = Transformer.from_crs(inputCRS, basicCRS, always_xy=True)

  # move x- and y-coordinates so that the lower left corner lies at [0,0]
  x_tmp = x + np.abs(x[0])
  y_tmp = y + np.abs(y[0])
  
  # transform x-coordinates
  y_lowerRow = np.repeat(y_tmp[0], len(x))
  xnew, ynew = transformer.transform(x_tmp, y_lowerRow)
  x_proj = xnew  
    
    # transform y-coordinates
  x_leftColumn = np.repeat(x_tmp[0], len(y))
  xnew, ynew = transformer.transform(x_leftColumn, y_tmp)
  y_proj = ynew

  return x_proj, y_proj
  
  
  
def setUp_grid_interpolation(coord_min, coord_max, dx, inputCRS):
  """
  Sets up the grid (image) for interpolation using VTK and probe filter. Returns the probe filter and shape for reshaping (needed within the interpolation).
  
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
  print(f"{varName} {timestep}: done probe filter in {stop - start} s.".center(column_size))

  polyout = probeFilter.GetOutput()
  projData = polyout.GetPointData().GetScalars()
  projDataNp = numpy_support.vtk_to_numpy(projData).reshape(projDataShape)
  
  return projDataNp
  
  
  
def interpolate_seissol2structured(sx, dx, coord_min, coord_max, inputCRS, include_horizontal):
  """
  Interpolate the SeisSol XDMF to VTK.
  
  :param sx: seissolxdmf file 
  :param dx: spatial resolution
  :param coord_min: minimum coordinates for box
  :param coord_max: maximum coordinates for box
  :param inputCRS:  CRS of the input 2d mesh 
  :param include_horizontal: handle whether to interpolate only the vertical component or not
  
  returns deformation data and coordinates
  """
  unstrGrid3d = sx.generateVtkObject()

  nTime = sx.ReadNdt()  # number of time steps in the Seissol file
  
  # Choose which data interpolate (only vertical or all components) based on the input handle  
  is_new_format = 'u3' in sx.ReadAvailableDataFields()
  if (include_horizontal):
    data = ['u1', 'u2', 'u3'] if is_new_format else ['U', 'V', 'W']
  else:
    data = ['u3'] if is_new_format else ['W']   
  
  # Create probe filter and get projected coordinates
  probeFilter, projDataShape, x_proj, y_proj = setUp_grid_interpolation(coord_min, coord_max, dx, inputCRS)

  # Read time data from seissolxdmf file
  seissol_time = get_seissol_time(sx)
 
  # List for interpolated data
  probedData = []
  
  # Perform interpolation (over each timestep and each variable)
  print("Interpolation is performed using VTK probe filter.".center(column_size))
  for timestep in range(nTime):
    for varName in data:
      projDataNp = get_interpolation(sx, unstrGrid3d, probeFilter, projDataShape, timestep, varName)
      probedData.append(projDataNp)
 
  return probedData, x_proj, y_proj, seissol_time



def get_seissol(filename, spatial_resolution, CRS_reference_coordinates, include_horizontal):
  """
  Actual part that will be used by the main donor model functionality to get the deformation data.
  
  :param filename: filename for the SeisSol data (has to be an XDMF file)
  :param spatial_resolution: spatial_resolution in meters
  :param CRS_reference_coordinates: CRS reference coordinates (list of longitude and latitude of lower left corner of the domain)  
  :param include_horizontal: boolean handle whether the output will only contain the vertical deformation (deformation) or also include the horizontal deformation
  """
  
  # CRS of the 2d mesh
  inputCRS = "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=" + CRS_reference_coordinates[0] + \
                      " +lat_0=" + CRS_reference_coordinates[1]  
  
  print("Getting output data from SeisSol.\n".center(column_size))
  
  sx = seissolxdmfExtended(filename) # get seissolxdmf from provided XDMF file
  
  # Get x, y and z interval min/max 
  geom = sx.ReadGeometry()
  coordinate_min = geom.min(0)
  coordinate_max = geom.max(0) 
  
  donor_deformation, donor_x, donor_y, donor_time = interpolate_seissol2structured(sx, spatial_resolution, coordinate_min, coordinate_max, inputCRS, include_horizontal)

  return donor_deformation, donor_x, donor_y, donor_time
