import numpy as np
import os
import time

from netCDF4 import Dataset
from scipy.interpolate import RectBivariateSpline

"""
Functionality to get the interpolated bathymetry (deformation taken into account)

Contains the following functionalities:

* get_interpolatedBathy   interpolate the bathymetry data to exchange grid and add the deformation for each timestep
"""

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


def get_interpolatedBathy(bathy_file, donor_x, donor_y, donor_deformation):
  """
  Interpolate the given bathymetry to the structured mesh given by donor_x and donor_y.
  
  :param bathy_file: file for the bathymetry (has to be larger than the mesh provided by the donor mesh)
  :param donor_x: x/lon-coordinates of the donor mesh 
  :param donor_y: y/lat-coordinates of the donor mesh 
  :param donor_deformation: deformation to calculate the actual (new) bathymetry 
  """

  # load bathy and SeisSol data
  bathy_data = Dataset(bathy_file, 'r', format='NETCDF4')

  # get values from bathymetry file
  try:
    bathy_x = bathy_data.variables['x']
    bathy_y = bathy_data.variables['y']
    bathy = bathy_data.variables['z']
  except:
    bathy_x = bathy_data.variables['lon']
    bathy_y = bathy_data.variables['lat']
    bathy = bathy_data.variables['elevation']

  # Create interpolation class 
  # Note that x and y-coordinates are flipped due to the storage in the netCDF file
  interpolator = RectBivariateSpline(bathy_y, bathy_x, bathy)

  # Get number of deformation timesteps
  Ntime = np.shape(donor_deformation)[0] 

  # Interpolate bathymetry to SeisSol grid
  start = time.time()
  print("Starting the interpolation (bathymetry).".center(column_size))
  
  # Calculate updated bathymetry (add the deformation for every timestep)
  interpolated_bathy = []
  for t in range(Ntime):
    interpolated_bathy.append(interpolator(donor_y, donor_x) - donor_deformation[t])
  stop = time.time()
  print(f"The interpolation took {stop - start} s\n".center(column_size))
  
  return interpolated_bathy 
