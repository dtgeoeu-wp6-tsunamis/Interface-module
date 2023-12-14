import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline
import time

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
  bathy_x = bathy_data.variables['x']
  bathy_y = bathy_data.variables['y']
  bathy = bathy_data.variables['z']

  xg, yg = np.meshgrid(bathy_x, bathy_y)

  # Create interpolation class 
  # Note that x and y-coordinates are flipped due to the storage in the netCDF file
  interpolator = RectBivariateSpline(bathy_y, bathy_x, bathy)

  # Get number of deformation timesteps
  Ntime = np.shape(donor_deformation)[0] 

  # Interpolate bathymetry to SeisSol grid
  start = time.time()
  print("Starting the interpolation (bathymetry).")
  
  # Calculate updated bathymetry (add the deformation for every timestep)
  interpolated_bathy = []
  for t in range(Ntime):
    interpolated_bathy.append(interpolator(donor_y, donor_x) - donor_deformation[t])
  stop = time.time()
  print(f"The interpolation took {stop - start} s\n")
  
  return interpolated_bathy 
