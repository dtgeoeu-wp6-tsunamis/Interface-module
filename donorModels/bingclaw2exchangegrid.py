import numpy as np
from pyproj import Transformer
import os
import time
from scipy.interpolate import RectBivariateSpline

"""
Module for the Bingclaw donor functionalities.

The Bingclaw files are formatted with an ESRI ASCII format. They are comprised of fort.a, fort.q and fort.t files that represent auxiliary variables, main output and time output. Only fort.q and fort.t files are needed. Although AMR is mentioned with Bingclaw, is is not utilized. Consequently, a fixed grid size is assumed.

Contains the following functionalities:

* read_bingclaw_time               read time from all fort.t files     
* read_bingclaw_header           read header of first fort.q files to set up Bingclaw data structure
* read_bingclaw_variables        read variables from all fort.q files       
* read_bingclaw_ESRIASCII       main routine to read bingclaw files
* interpolate_bingclawdata       perform interpolation to new grid
* get_bingclaw                          parent routine that is called from the main donorModel routine
"""

# Global parameters
basicCRS = 'epsg:4326' # basic lat-lon coordinate system


def read_bingclaw_time(datapath, filenames):
  """
  Read the time data for a set of Bingclaw files. The time files have the following structure:
  
  --------------------
  time
  meqn
  ngrids
  naux
  ndim
  nghost
  --------------------
  
  Returns number of timesteps and time data for each file.
  
  :param datapath:  path to the time data  
  :param filenames: list of filenames for the time data
  """
  
  # Get number of timesteps
  Ntime = len(filenames)

  time_data = np.zeros(Ntime)      
  for time in range(Ntime):
    
    # Read only the first line of a new file
    with open(os.path.join(datapath, filenames[time])) as f:
      first_line = f.readline().split()
      time_data[time] = float(first_line[0])
      
  return Ntime, time_data



def read_bingclaw_header(filename):
  """
  Read the header for a Bingclaw file which always has the following structure:
  
  --------------------
  grid_number           can be ignored
  AMR_level               can be ignored
  mx                          number of points in x/longitude-direction
  my                          number of points in y/latitude-direction
  xlow                        lowest longitude value
  ylow                        lowest latitude value
  dx                           resolution for longitude
  dy                           resolution for latitude
  --------------------
  
  Returns the last six variables.
  
  :param filename: name of the file for which the header should be read
  """

  # Open file and read the the header
  bingclaw_file = open(filename, 'r')
  header = []
  for i in range(8):
    header.append(bingclaw_file.readline().split())

  # Get respective outputs from header
  Nx = int(header[2][0])
  Ny = int(header[3][0])
  x_min = float(header[4][0])
  y_min = float(header[5][0])
  dx = float(header[6][0])
  dy = float(header[7][0])
  
  bingclaw_file.close()
  
  return Nx, Ny, x_min, y_min, dx, dy
  
  

def read_bingclaw_variables(Ntime, Nx, Ny, datapath, filenames):
  """
  Read the bingclaw data. The data is formatted in the following way:
  (h, u_p , v_p , hu, hv, γ, η), where γ is the accumulated shear strain at the bottom and η = h + b is the (sea) surface elevation. 
  
  The output starts from xlow and ylow in x-direction and needs to be transposed, as a result to match the netCDF file structure.
  
  :param Ntime: Number of timesteps
  :param Nx:      Number of points in x/longitude-direction
  :param Ny:      Number of points in y/latitude-direction
  :param datapath:  path to the time data  
  :param filenames: list of filenames for the time data
  """
  
  bingclaw_deformation = np.zeros((Ntime, Ny, Nx))
  
  for time in range(Ntime):
    if (time % 10 == 0): print(f"Reading timestep {time} out of {Ntime}.")
    bingclaw_file = open(os.path.join(datapath, filenames[time]), 'r')
    bingclaw_data = np.loadtxt(bingclaw_file, skiprows=8)
    
    #raster_h = np.reshape(bingclaw_data[:,0], (Nx, Ny))
    raster_η = np.reshape(bingclaw_data[:,-1], (Nx, Ny))
    bingclaw_deformation[time] = np.transpose(raster_η) 

  return bingclaw_deformation


        
def read_bingclaw_ESRIASCII(donor_output_path):
  """     
  Main functionality to read the data from the ESRI ASCII Bingclaw files.
  
  :param donor_output_path: path to the direction where the ESRI ASCII files are stored
  """
  
  # Get all output files for Bingclaw (variable and time data)
  script_path = os.getcwd()
  binglcaw_path = os.path.join(script_path, donor_output_path)
  data_files = []
  data_files += [each for each in sorted(os.listdir(binglcaw_path)) if each.startswith('fort.q')]
  time_files = []
  time_files += [each for each in sorted(os.listdir(binglcaw_path)) if each.startswith('fort.t')]
  
  # Read time data
  Ntime, bingclaw_time = read_bingclaw_time(donor_output_path, time_files)

  # Read grid structure from the first data file
  Nx, Ny, x_min, y_min, dx, dy = read_bingclaw_header(os.path.join(donor_output_path, data_files[0]))
  
  # Create grid from header values  
  x_max = x_min + (Nx - 1) * dx
  y_max = y_min + (Ny - 1) * dy
  bingclaw_x = np.linspace(x_min, x_max, Nx)
  bingclaw_y = np.linspace(y_min, y_max, Ny)

  start = time.time()
  bingclaw_deformation = read_bingclaw_variables(Ntime, Nx, Ny, donor_output_path, data_files)
  stop = time.time()    
  print(f"Data has been read. It took {stop - start} s.")
    
  return Ntime, bingclaw_deformation, bingclaw_x, bingclaw_y



def interpolate_bingclawdata(Ntime, bingclaw_x, bingclaw_y, bingclaw_deformation, donor_x, donor_y):
  """
  Interpolation routine to interpolate the Bingclaw data to the new grid. Makes use of RectBivariateSpline.
  
  :param Ntime: Number of timesteps
  :param bingclaw_x:  x/longitude from Bingclaw
  :param bingclaw_y:  y/latitude from Bingclaw
  :param bingclaw_deformation: deformation from Bingclaw
  :param donor_x:  x/longitude for new grid
  :param donor_y:  y/latitude for new grid
  """
  
  interpolated_deformation = []
  
  for time in range(Ntime):
    # Create interpolation class (new for each timestep)
    # Note that x and y-coordinates are flipped due to the storage in the netCDF file
    interpolator = RectBivariateSpline(bingclaw_y, bingclaw_x, bingclaw_deformation[time])

    interpolated_deformation.append(interpolator(donor_y, donor_x))
    
  return interpolated_deformation



def get_bingclaw(donor_output_path, spatial_resolution):
  """
  Main donor model functionality to get the deformation data.
  
  :param donor_output_path: path to the direction where the ESRI ASCII files are stored
  :param spatial_resolution: spatial_resolution in meters
  """

  print("Getting output data from Bingclaw.\n")

  # Get projected spatial resolution
  inputCRS = "+proj=tmerc +datum=WGS84"  # CRS
  transformer = Transformer.from_crs(inputCRS, basicCRS, always_xy=True)
  projection_resolution = transformer.transform(spatial_resolution, 1.)[0]

  # Get Bingclaw data
  Ntime, bingclaw_deformation, bingclaw_x, bingclaw_y = read_bingclaw_ESRIASCII(donor_output_path)

  # Use Bingclaw coordinates to create new grid coordinates (with projection_resolution)
  donor_x = np.arange(np.min(bingclaw_x), np.max(bingclaw_x) + projection_resolution, projection_resolution)
  donor_y = np.arange(np.min(bingclaw_y), np.max(bingclaw_y) + projection_resolution, projection_resolution)
  
  # Interpolate Bingclaw data to new grid
  print("Starting the interpolation (Bingclaw).")
  start = time.time()
  donor_deformation = interpolate_bingclawdata(Ntime, bingclaw_x, bingclaw_y, bingclaw_deformation, donor_x, donor_y)
  stop = time.time()    
  print(f"The interpolation took {stop - start} s\n")

  return donor_deformation, donor_x, donor_y
