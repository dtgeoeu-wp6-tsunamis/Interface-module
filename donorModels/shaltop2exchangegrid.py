import numpy as np
import os
import pyvista as pv
import time

from pathlib import Path
from pyproj import Transformer
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import KDTree


#TODO Lon/Lat for the bottom left coordinates of the domain should be provided by some way within the SHALTOP output. Will be used for inputCRS. 

"""
Module for the SHALTOP donor functionalities.

This module reads the water colum height from a provided SHALTOP data directory. The correct have to be computed from the mesh normal vectors due to SHALTOP's local coordinates.
Assumes a nearest neighbor approach for the meshing for efficiency reasons. The SHALTOP HySEA coupling written by Alexis Marboeuf allows for more options (see GitHub).

Contains the following functionalities:

* get_currentnormal            outputs current normal vector and ensures that it points upwards
* get_mesh                          calculates the mesh from given vertices
* get_newvertices                calculates new vertices based on the normal vectors
* get_normaldistances         calculates the distance between old and new vertices
* project_coordinates           project coordinates to WGS84 coordinates
* calculate_correct_height    calculates the correct water column height (for a lon/lat grid)
* obtain_shaltop_data           routine to obtain the data from SHALTOP  
* interpolate_shaltop_data    interpolate SHALTOP data to exchange grid
* get_shaltop                         parent routine that is called from the main donorModel routine
"""

# Global parameters
basicCRS = 'epsg:4326' # basic lat-lon coordinate system

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


def get_currentnormal(mesh, idx):
  """
  Get normal vector for a number of normal vectors. Also ensures that the normal vector points updwards.
  
  :param mesh: input mesh
  :param idx:     current index
  """
  
  normal = mesh["Normals"][idx]  
  if (normal[2] < 0.0): normal = -normal
  return normal



def get_mesh(vertices):
  """
  Gets the mesh data from input vertices.
  
  :param vertices: vertices 
  """
  
  mesh = pv.PolyData(vertices)
  meshD = mesh.delaunay_2d()
  meshT = meshD.extract_surface()
  meshFaces = np.insert(meshT.faces.reshape(-1,4)[:, 1:],0,3,axis=1)
  mesh = pv.PolyData(vertices, meshFaces)
  return mesh



def get_newvertices(mesh_normals, height):
  """
  Calculates new vertices with the provided normal vectors.
  
  :param mesh_normals: normal vectors
  :param height:  water column height
  """

  new_vertices = np.zeros((mesh_normals.n_points,3))
  for i in range(mesh_normals.n_points):
    normal = get_currentnormal(mesh_normals, i)
    new_vertices[i,:] = mesh_normals.points[i] + normal * height[i]
    
  return new_vertices


def get_normaldistances(vertices, new_vertices, mesh_normals):
  """
  Calculate the distances between the old and new vertices.
  
  :param vertices:            old vertices
  :param new_vertices:    new vertices
  :param mesh_normals:  normal vectors
  """
  
  # Create distance storage
  mesh_normals["distances"] = np.empty(mesh_normals.n_points)

  # Calculate distances (from tree structure)
  vertices2D = np.delete(vertices, 2, axis=1)
  vertices22D = np.delete(new_vertices, 2, axis=1)
  tree = KDTree(vertices22D)
  distances, idx = tree.query(vertices2D)
  mesh_normals["distances"] = np.abs(vertices[:,2] - new_vertices[idx,2])

  # Remove NaNs
  mask = np.isnan(mesh_normals["distances"])
  mesh_normals["distances"][mask] = 0.0  
  
  return mesh_normals["distances"]


  
def project_coordinates(x, y, inputCRS):
  """
  This function takes the input x and y coordinates and returns the transformed WGS84 coordinates.

  :param x: x-coordinates
  :param y: y-coordinates
  :param inputCRS: input CRS
  """
  transformer = Transformer.from_crs(inputCRS, basicCRS, always_xy=True)

    # transform x-coordinates
  y_lowerRow = np.repeat(y[0], len(x))
  xnew, ynew = transformer.transform(x, y_lowerRow)
  x_proj = xnew  
    
    # transform y-coordinates
  x_leftColumn = np.repeat(x[0], len(y))
  xnew, ynew = transformer.transform(x_leftColumn, y)
  y_proj = ynew

  return x_proj, y_proj



def calculate_correct_height(shaltop_local_deformation, bathymetry, Ntime, Nx, Ny, xmin, ymin, dx, dy):
  """
  This function calculates the correct water column height for each node.
  
  :param shaltop_local_deformation:  Local deformation
  :param bathymetry:   Bathymetry data
  :param Ntime:   Number of timesteps
  :param Nx:        Number of points in x-direction
  :param Ny:        Number of points in y-direction
  :param xmin:    Smallest x values
  :param ymin:    Smallest y values
  :param dx:        Resolution in x-direction
  :param dy:        Resolution in y-direction
  """

  shaltop_deformation = np.zeros((Ntime, Ny, Nx))
  previous_deformation = np.zeros((Ny, Nx))
  frmt = len(str(Ntime)) # formatter for terminal output  
  
  for time in range(Ntime):
    if (time % 10 == 0): print(f"Reading timestep {time+1:{frmt}d} out of {Ntime}.".center(column_size))
    
    current_deformation = shaltop_local_deformation[time]
    
    # Indices where deformation is larger than zero
    posy, posx = np.where(current_deformation > 0.0)
    wsize = np.shape(posy)[0]
    
    # Define coordinates for found indices
    x = [xmin+dx*i for i in posx]
    y = [ymin+dy*(Ny-1-j) for j in posy]
    current_bathymetry = bathymetry[posy, posx]
    current_height = current_deformation[posy, posx]
    
    # Build vertices for all found indices
    vertices = np.array([[x[i], y[i], current_bathymetry[i]] for i in range(0, wsize)])
    
    # Get mesh from vertices
    mesh = get_mesh(vertices)
    
    # Compute normal vectors
    mesh_normals = mesh.compute_normals(point_normals=True, cell_normals=False, auto_orient_normals=False, flip_normals=False)
    
    # Compute new vertices (with normal vectors)
    new_vertices = get_newvertices(mesh_normals, current_height)
    
    # Get new mesh from updated vertices
    mesh2 = get_mesh(new_vertices)
    
    # Calculate distances between normals of both vertices types
    mesh_normals["distances"] = get_normaldistances(vertices, new_vertices, mesh_normals)
    
    # Calculate deformation relative to previous timestep
    tmp_deformation = np.zeros((Ny, Nx))
    tmp_deformation[posy,posx] = mesh_normals["distances"]
    if time > 0:
      shaltop_deformation[time,:,:] = - (previous_deformation - tmp_deformation)
    
    # Save deformation for next timestep
    previous_deformation = tmp_deformation
    
  return shaltop_deformation



def obtain_shaltop_data(shaltop_data_path):
  """
  Reads the data from SHALTOP. We assume that the data is stored in the directory given by shaltop_data_path and that the scenario input files are stored in its parent directory.
  Part of the routine calculates the correct water column height.
  
  :param shaltop_data_path:  path to the directory where the SHALTOP output is stored
  """

  start = time.time()
  shaltop_input_path = Path(shaltop_data_path).parent
  
  shaltop_parameters = np.loadtxt(os.path.join(shaltop_data_path, 'plotmat.dat'))
  Nx = int(shaltop_parameters[2,0])
  Ny = int(shaltop_parameters[2,1])
  
  shaltop_time = np.loadtxt(os.path.join(shaltop_data_path, 'time_im.d'))
  Ntime = np.shape(shaltop_time)[0]
  
  # Load local x- and y-coordinates
  local_x = np.loadtxt(os.path.join(shaltop_data_path, 'grillexw.d'))
  local_y = np.loadtxt(os.path.join(shaltop_data_path, 'grilleyw.d'))
  
  # Computing error for initial time
  bathymetry = np.reshape(np.loadtxt(os.path.join(shaltop_input_path, 'z.d')), (Ny, Nx))
  
  # Domain [xmin, xmax]
  xmin = np.min(local_x)
  xmax = np.max(local_x)
  ymin = np.min(local_y)
  ymax = np.max(local_y)
  
  # Computation of local height:
  dx = (xmax - xmin) / Nx
  dy = (ymax - ymin) / Ny
  
  # Get local height data from SHALTOP output files
  shaltop_deformation_file = open(os.path.join(shaltop_data_path, 'rho.bin'),'r')
  shaltop_local_deformation = np.fromfile(shaltop_deformation_file, dtype='float32', count=Ntime*Nx*Ny).reshape((Ntime, Ny, Nx))
  
  # Calculate correct height for each timestep
  shaltop_deformation = calculate_correct_height(shaltop_local_deformation, bathymetry, Ntime, Nx, Ny, xmin, ymin, dx, dy)
  stop = time.time()    
  print(f"Data has been obtained. It took {stop - start} s.".center(column_size))
  
  return shaltop_deformation, local_x, local_y, shaltop_time



def interpolate_shaltop_data(shaltop_deformation, spatial_resolution, shaltop_x, shaltop_y, Ntime):
  """
  Interpolation routine to interpolate the SHALTOP data to the new grid. Makes use of RectBivariateSpline.
  
  :param shaltop_deformation: correct deformation from SHALTOP
  :param spatial_resolution: spatial_resolution in meters
  :param shaltop_x:  x-coordinates used within SHALTOP
  :param shaltop_y:  y-coordinates used within SHALTOP 
  :param Ntime:   Number of timesteps
  """
  
  interpolated_x = np.arange(np.min(shaltop_x), np.max(shaltop_x) + spatial_resolution, spatial_resolution)
  interpolated_y = np.arange(np.min(shaltop_y), np.max(shaltop_y) + spatial_resolution, spatial_resolution)

  interpolated_deformation = []
  for time in range(Ntime):
    # Create interpolation class (new for each timestep)
    # Note that x and y-coordinates are flipped due to the storage in the netCDF file
    interpolator = RectBivariateSpline(shaltop_y, shaltop_x, shaltop_deformation[time])
    
    interpolated_deformation.append(interpolator(interpolated_y, interpolated_x))

  return interpolated_deformation, interpolated_x, interpolated_y
  
  
  
def get_shaltop(donor_output_path, spatial_resolution, CRS_reference_coordinates):
  """
  Main donor model functionality to get the deformation data.

  :param donor_output_path: path to the directory where the SHALTOP output is stored
  :param spatial_resolution: spatial_resolution in meters
  :param CRS_reference_coordinates: CRS reference coordinates (list of longitude and latitude of lower left corner of the domain)
  """
  
  # CRS of the 2d mesh
  inputCRS = "+proj=tmerc +datum=WGS84 +k=0.9996 +lon_0=" + CRS_reference_coordinates[0] + \
                      " +lat_0=" + CRS_reference_coordinates[1]  
  
  print("Getting output data from SHALTOP.\n".center(column_size))
  
  shaltop_deformation, shaltop_x, shaltop_y, donor_time = obtain_shaltop_data(donor_output_path)
  
  # Interpolate Bingclaw data to new grid
  print("Starting the interpolation (SHALTOP).".center(column_size))
  start = time.time()

  donor_deformation, interpolated_x, interpolated_y = interpolate_shaltop_data(shaltop_deformation, spatial_resolution, shaltop_x, shaltop_y, len(donor_time))

  donor_x, donor_y = project_coordinates(interpolated_x, interpolated_y, inputCRS)
  stop = time.time()    
  print(f"The interpolation took {stop - start} s.\n".center(column_size))
  
  return donor_deformation, donor_x, donor_y, donor_time



