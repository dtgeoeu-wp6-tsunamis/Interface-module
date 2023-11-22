import vtk
from vtk.util import numpy_support
import numpy as np
import seissolxdmf
import time
import argparse
from pyproj import Transformer

# TODO: separate returns for horizontal deformation
# TODO: set up inputCRS correctly

"""
Module for the SeisSol donor functionalities.
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
  
  
  
def projectCoordinates(x, y, inputcrs):
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
  
  
  
def interpolateSeissol2structured(sx, dx, coord_min, coord_max, include_horizontal, instants=[]):
  """
  Interpolate the SeisSol XDMF to VTK
  :param sx: seissolxdmf file 
  :param dx: spatial resolution
  :param coord_min: minimum coordinates for box
  :param coord_max: maximum coordinates for box
  :param include_horizontal: handle whether to interpolate only the vertical component or not
  :param instants: time steps to include in the out netCDF. If not provided, all time steps will be included
  
  returns uplift data and coordinates
  """
  unstrGrid3d = sx.generateVtkObject()

  ndt = sx.ReadNdt()  # number of time steps in the Seissol file
  if not instants:
    # if not specific instant provided, use all in the seissol output
    instants = list(range(0, ndt))
  else:
    # check if instants provided are within the possible ones
    if not all(x in range(0, ndt) for x in instants):
      print("Instants provided are outside the range of timesteps in SeisSol's output")
  
  # Choose which data interpolate (only vertical or all components) based on the input handle  
  if (include_horizontal):
    data = ['u1', 'u2', 'u3']
  else:
    data = ['u3']    
  
## Set up the mesh    
  # set up x and y coordinates
  x = np.arange(coord_min[0], coord_max[0] + dx, dx)
  y = np.arange(coord_min[1], coord_max[1] + dx, dx)

  z = np.array([0])   # ensure that the mesh is 2D
  xx, yy = np.meshgrid(x, y)

  inputCRS = "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=26.25 +lat_0=37.75"  # CRS of the input 2d mesh

  # project the x and y coordinates to lat/lon
  transformer = Transformer.from_crs(inputCRS, "epsg:4326", always_xy=True)
  x_proj, y_proj = projectCoordinates(x, y, inputCRS)

  # Create grid image volume
  image1Size = [x.shape[0], y.shape[0], z.shape[0]]
  image1Origin = [coord_min[0], coord_min[1], coord_min[2]]
  image1Spacing = [dx, dx, dx]
  
  imageData1 = vtk.vtkImageData()
  imageData1.SetDimensions(image1Size)
  imageData1.SetOrigin(image1Origin)
  imageData1.SetSpacing(image1Spacing)

  # Perform the interpolation
  probeFilter = vtk.vtkProbeFilter()

  probeFilter.SetInputData(imageData1)
  probeFilter.SpatialMatchOn()

  probedData = []

  for idt in instants:

      for var in data:
          print("update scalars")
          scalars = vtk.vtkFloatArray()

          W = sx.ReadData(var, idt)
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
          print("start probe filter")
          probeFilter.Update()
          stop = time.time()
          print(f"{var} {idt}: done probe filter in {stop - start} s")

          polyout = probeFilter.GetOutput()
          projData = polyout.GetPointData().GetScalars()
          projDataNp = numpy_support.vtk_to_numpy(projData).reshape(xx.shape)
          probedData.append(projDataNp)
         
  return probedData, x_proj, y_proj



def get_seissol(filename, spatial_resolution, include_horizontal):
  """
  Actual part that will be used by the main donor model functionality to get the uplift data
  :param filename: filename for the SeisSol data (has to be an XDMF file)
  :param spatial_resolution: spatial_resolution in meters
  :param include_horizontal: boolean handle whether the output will only contain the vertical deformation (uplift) or also include the horizontal deformation
  """

  sx = seissolxdmfExtended(filename) # get seissolxdmf from provided XDMF file
# get x, y and z interval min/max and round to nearest 1000 (if not done, the interpolation does nothing)
  geom = sx.ReadGeometry()
  coordinate_min = np.round(geom.min(0) +  spatial_resolution, -4) 
  coordinate_max = np.round(geom.max(0) -  spatial_resolution, -4) 

  donor_uplift, donor_x, donor_y = interpolateSeissol2structured(sx, spatial_resolution, coordinate_min, coordinate_max, include_horizontal)

  return donor_uplift, donor_x, donor_y
