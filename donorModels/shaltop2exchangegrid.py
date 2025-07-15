import numpy as np
import os
import time

from pathlib import Path
from pyproj import Transformer
from scipy.interpolate import RectBivariateSpline
from scipy.spatial import KDTree
from scipy import interpolate


#TODO Lon/Lat for the bottom left coordinates of the domain should be provided by some way within the SHALTOP output. Will be used for inputCRS. 

"""
Module for the SHALTOP donor functionalities.

This module reads the water colum height from a provided SHALTOP data directory. The correct height has to be computed from the mesh normal vectors due to SHALTOP's local coordinates.

The SHALTOP data directory produces several output files and the following are required:
* plotmat.dat (to read parameter data)
* time_im.d (times for which data is written),
* grillexw.d and grilleyw.d (provide the x- and y-coordinates)
* z.bin (for bathymetry) and finally 
* rho.bin (which contains the actual height data, but in local coordinates)

Contains the following functionalities:

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
  
def project_coordinates(x, y, projection):
  """
  This function takes the input x and y coordinates and returns the transformed WGS84 coordinates.

  :param x: x-coordinates
  :param y: y-coordinates
  :param projection: grid information for converting cartesian coordinates
  """
  gridinfo = projection[1:-1].split(",")

  inputCRS = gridinfo[2]

  # Perform the transform from the resolution in m to degree
  transformer = Transformer.from_crs(inputCRS, basicCRS, always_xy=True)
  ytmp = y[0]*np.ones(np.shape(x))
  x_proj,_ = transformer.transform(x,ytmp)
  xtmp = x[0]*np.ones(np.shape(y))
  _,y_proj = transformer.transform(xtmp,y)

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
  #previous_deformation = np.zeros((Ny, Nx))
  original_mass = np.zeros((Ny, Nx))
  frmt = len(str(Ntime)) # formatter for terminal output  

  # gradients: upward normal at each bathy point is (-sinx,-siny,cost)
  dzdy,dzdx = np.gradient(bathymetry,dy,dx)
  aux=np.sqrt(1.0+dzdy*dzdy+dzdx*dzdx)
  cost=1./aux
  sinx=dzdx/aux
  siny=dzdy/aux
  
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

    # Define outward normals for found indices: (-sinx,-siny,cost)
    costr = cost[posy,posx]
    sinxr = sinx[posy,posx]
    sinyr = siny[posy,posx]
    normal = np.transpose(np.stack([-sinxr,-sinyr,costr]))
    
    # Build vertices for all found indices and create arrays for vertices representing the surface from Shaltop
    vertices = np.array([[x[i], y[i], current_bathymetry[i]] for i in range(0, wsize)])
    vertices2 = np.zeros((wsize,3))
    H = np.zeros(wsize)

    for i in range(wsize):
      vertices2[i,:] = vertices[i,:] + normal[i,:] * current_height[i]
      H[i] = vertices2[i,2]

    tmp_deformation = np.zeros((Ny, Nx))
    # Nearest neighboor: faster but less accurate than the cubic interpolation
    # vertices2D = np.delete(vertices,2,axis=1)
    # vertices22D = np.delete(vertices2,2,axis=1)
    # tree = KDTree(vertices22D)
    # distances, idx = tree.query(vertices2D)
    # tmp_deformation[posy,posx] = np.abs(vertices[:,2] - vertices2[idx,2])
    f = interpolate.griddata(np.delete(vertices2,2,axis=1),H,np.delete(vertices,2,axis=1),method='cubic')
    tmp_deformation[posy,posx] = np.maximum(f.reshape(-1)-current_bathymetry,np.zeros(np.shape(current_bathymetry)))
    # Need to remove NaN from the cubic interpolation
    mask = np.isnan(tmp_deformation)
    tmp_deformation[mask] = 0.0
    if time == 0:
      original_mass[posy,posx] = tmp_deformation[posy,posx]
    else:
      shaltop_deformation[time,:,:] = original_mass - tmp_deformation
    # # Calculate deformation relative to previous timestep
    # if time > 0:
    #   shaltop_deformation[time,:,:] = tmp_deformation - previous_deformation
    # # Save deformation for next timestep
    # previous_deformation = tmp_deformation
    
  return shaltop_deformation



def obtain_shaltop_data(shaltop_data_path,projection):
  """
  Reads the data from SHALTOP. We assume that the data is stored in the directory given by shaltop_data_path and that the scenario input files are stored in its parent directory.
  Part of the routine calculates the correct water column height.
  
  :param shaltop_data_path:  path to the directory where the SHALTOP output is stored
  """

  start = time.time()
  shaltop_input_path = Path(shaltop_data_path).parent

  gridinfo = projection[1:-1].split(",")
  
  # data2 is the default directory of Shaltop outputs in the Shaltop directory
  shaltop_parameters = np.loadtxt(os.path.join(shaltop_data_path, 'data2/plotmat.dat'))
  Nx = int(shaltop_parameters[2,0])
  Ny = int(shaltop_parameters[2,1])
  
  shaltop_time = np.loadtxt(os.path.join(shaltop_data_path, 'data2/time_im.d'))
  Ntime = np.shape(shaltop_time)[0]
  
  # Load local x- and y-coordinates
  local_x = np.loadtxt(os.path.join(shaltop_data_path, 'data2/grillexw.d'))
  local_y = np.loadtxt(os.path.join(shaltop_data_path, 'data2/grilleyw.d'))

  local_x = local_x + float(gridinfo[3])
  local_y = local_y + float(gridinfo[4])

  # Computing error for initial time
  bathymetry_file = open(os.path.join(shaltop_data_path, 'data2/z.bin'),'r')
  bathymetry = np.fromfile(bathymetry_file, dtype='float32', count=Nx*Ny).reshape((Ny, Nx))
  ####### Flip due to error while initializing Shaltop. This will be fixed for the next scenarios. ######
  bathymetry = np.flip(bathymetry,axis=0)
  
  # Domain [xmin, xmax]
  xmin = np.min(local_x)
  xmax = np.max(local_x)
  ymin = np.min(local_y)
  ymax = np.max(local_y)
  
  # Computation of local height:
  dx = (xmax - xmin) / Nx
  dy = (ymax - ymin) / Ny
  
  # Get local height data from SHALTOP output files
  shaltop_deformation_file = open(os.path.join(shaltop_data_path, 'data2/rho.bin'),'r')
  shaltop_local_deformation = np.fromfile(shaltop_deformation_file, dtype='float32', count=Ntime*Nx*Ny).reshape((Ntime, Ny, Nx))
  ###### Flip due to error while initializing Shaltop. This will be fixed for the next scenarios. ######
  shaltop_local_deformation = np.flip(shaltop_local_deformation,axis=1)
  
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
  
  
  
def get_shaltop(donor_output_path, spatial_resolution, projection):
  """
  Main donor model functionality to get the deformation data.

  :param donor_output_path: path to the directory where the SHALTOP output is stored
  :param spatial_resolution: spatial_resolution in meters
  *param projection: projection parameters (in Proj4 format) for converting from Cartesian to geographic (lon, lat) coordinates
  """
  
  print("Getting output data from SHALTOP.\n".center(column_size))
  
  shaltop_deformation, shaltop_x, shaltop_y, donor_time = obtain_shaltop_data(donor_output_path,projection)
  
  # Interpolate Bingclaw data to new grid
  print("Starting the interpolation (SHALTOP).".center(column_size))
  start = time.time()

  donor_deformation, interpolated_x, interpolated_y = interpolate_shaltop_data(shaltop_deformation, spatial_resolution, shaltop_x, shaltop_y, len(donor_time))

  donor_x, donor_y = project_coordinates(interpolated_x, interpolated_y, projection)
  stop = time.time()    
  print(f"The interpolation took {stop - start} s.\n".center(column_size))

  ####################

  # netCDFFile = netCDF4.Dataset('/home/marboeuf/shalbing-to-hysea/outputs/mscen_v0.141_x0_15.471_y0_38.004/shaltop_out/shaltopverif.nc', "w", format="NETCDF4")

  # netCDFFile.Conventions = "CF-1.5"
  # netCDFFile.GDAL = "GDAL 3.4.1, released 2021/12/27"
  # netCDFFile.history = "created by script"
  # netCDFFile.NCO = "4.7.2"
  # netCDFFile.nco_openmp_thread_number = 1.0

  # lat = netCDFFile.createDimension("lat", len(donor_y))
  # lon = netCDFFile.createDimension("lon", len(donor_x))
  # timev = netCDFFile.createDimension("time", len(donor_time))

  # latitudes = netCDFFile.createVariable("lat","f4",("lat",))
  # latitudes.long_name = "latitude"
  # latitudes.standard_name = "latitude"
  # latitudes.units = "degrees_north"
  # longitudes = netCDFFile.createVariable("lon","f4",("lon",))
  # longitudes.long_name = "longitude"
  # longitudes.standard_name = "longitude"
  # longitudes.units = "degrees_east"
  # times = netCDFFile.createVariable("times","f4",("time",))
  # times.long_name = "time"
  # times.standard_name = "time"
  # times.units = "seconds"
  # z = netCDFFile.createVariable("z","f4",("time","lat","lon",),fill_value=-9999.0)
  # z.units = "GDAL Band Number 1"

  # latitudes[:] = donor_y
  # longitudes[:] = donor_x

  # # For an unknown reason, ASCII GRD file store data with a reversed order for latitudes compared to netCDF
  # times[:] = donor_time
  # z[:,:,:] = donor_deformation

  # netCDFFile.close()

  ####################

  return donor_deformation, donor_x, donor_y, donor_time



