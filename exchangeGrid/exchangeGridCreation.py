import numpy as np
import os
import time

from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline

from .interpolateBathy import get_interpolatedBathy

"""
Functionality to either choose to interpolate the 

Contains the following functionalities:

* donor2bathyDomain    write donor deformation to a larger grid similar to the bathymetry grid
* createExchangeGrid    functionality from which the exchange grid operation are performed
"""

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


def donor2bathyDomain(bathy_file, donor_x, donor_y, donor_deformation):
  """
  Write the given donor deformation to a structured mesh similar to the bathymetry grid but defined with the spatial resolution provided via the terminal input.
  
  :param bathy_file: file for the bathymetry (has to be larger than the mesh provided by the donor mesh)
  :param donor_x: x/lon-coordinates of the donor mesh 
  :param donor_y: y/lat-coordinates of the donor mesh 
  :param donor_deformation: deformation to calculate the actual (new) bathymetry 
  """

  # Load bathy and donor data
  bathy_data = Dataset(bathy_file, 'r', format='NETCDF4')

  # Get coordinate values from bathymetry file
  try:
    bathy_x = bathy_data.variables['x']
    bathy_y = bathy_data.variables['y']
  except:
    bathy_x = bathy_data.variables['lon']
    bathy_y = bathy_data.variables['lat']

  x_min = np.min(bathy_x)
  x_max = np.max(bathy_x)
  y_min = np.min(bathy_y)
  y_max = np.max(bathy_y)

  start = time.time()
  print("Starting the exchange grid creation.".center(column_size))
  
  # Get spatial resolution from the difference of some donor coordinates (donor grid is uniform!)
  dx = donor_x[1] - donor_x[0]

  # Calculate new number of points for both x- and y-coordinates
  new_Nx = int((bathy_x[-1] - bathy_x[0]) / dx) + 1
  new_Ny = int((bathy_y[-1] - bathy_y[0]) / dx) + 1

  # Calculate indices at which the "old"/donor data will start and end on the "new"/exchange grid
  donor_x_index_start = int((donor_x[0] - x_min)/dx)
  donor_x_index_end = donor_x_index_start + len(donor_x)
  donor_y_index_start = int((donor_y[0] - y_min)/dx)
  donor_y_index_end = donor_y_index_start + len(donor_y)
      
  # Create new x-coordinate array
  exchange_grid_x = np.zeros(new_Nx)
  exchange_grid_x[donor_x_index_start:donor_x_index_end] = donor_x
  
  # Calculate the new coordinate points recursively from the donor coordinates
  for idx in range(donor_x_index_start):
    exchange_grid_x[donor_x_index_start-idx-1] = exchange_grid_x[donor_x_index_start-idx] - dx
  for idx in range(new_Nx-donor_x_index_end):
    exchange_grid_x[donor_x_index_end+idx] = exchange_grid_x[donor_x_index_end+idx-1] + dx
    
  # Create new y-coordinate array  
  exchange_grid_y = np.zeros(new_Ny)
  exchange_grid_y[donor_y_index_start:donor_y_index_end] = donor_y
  
  # Calculate the new coordinate points recursively from the donor coordinates
  for idx in range(donor_y_index_start):
    exchange_grid_y[donor_y_index_start-idx-1] = exchange_grid_y[donor_y_index_start-idx] - dx
  for idx in range(new_Ny-donor_y_index_end):
    exchange_grid_y[donor_y_index_end+idx] = exchange_grid_y[donor_y_index_end+idx-1] + dx
    
  # Get number of deformation timesteps
  Ntime = np.shape(donor_deformation)[0] 
  
  # Create exchange grid deformation array and fill it with the "old"/donor data (with the indices as before)
  exchange_grid_deformation = np.zeros((Ntime, new_Ny, new_Nx))
  exchange_grid_deformation[:, donor_y_index_start:donor_y_index_end, donor_x_index_start:donor_x_index_end] = donor_deformation

  stop = time.time()
  print(f"The interpolation took {stop - start} s\n".center(column_size))
  
  return exchange_grid_deformation, exchange_grid_x, exchange_grid_y



#**********************************************************************************
# Generic function to create the exchange grid
def createExchangeGrid(bathy_file, donor_x, donor_y, donor_deformation, only_donor_domain):
  """
  Write the given donor deformation to a structured mesh similar to the bathymetry grid.
  
  :param bathy_file: file for the bathymetry (has to be larger than the mesh provided by the donor mesh)
  :param donor_x: x/lon-coordinates of the donor mesh 
  :param donor_y: y/lat-coordinates of the donor mesh 
  :param donor_deformation: deformation to calculate the actual (new) bathymetry 
  :param only_donor_domain: handle to only use the domain given by the donor model
  """
  
  # Check if only donor domain should be used
  if (only_donor_domain):
    exchange_grid_deformation = donor_deformation
    exchange_grid_x = donor_x
    exchange_grid_y = donor_y
  else:
    exchange_grid_deformation, exchange_grid_x, exchange_grid_y = donor2bathyDomain(bathy_file, donor_x, donor_y, donor_deformation)
  interpolated_bathymetry = get_interpolatedBathy(bathy_file, exchange_grid_x, exchange_grid_y, exchange_grid_deformation)
  return exchange_grid_deformation, exchange_grid_x, exchange_grid_y, interpolated_bathymetry
