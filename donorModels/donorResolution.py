import numpy as np

from netCDF4 import Dataset
from pyproj import Transformer

"""
Functionality to check the spatial resolution and output the bathymetry.

Contains the following functionalities:

* donor_chooseResolution  checks resolution (if 0) and outputs bathymetry
"""


def donor_chooseResolution(spatial_resolution, bathy_file):
  """
  Checks whether the given spatial resolution (im m) is 0. If yes, uses the spatial resolution of the bathymetry

  :param spatial_resolution: spatial resolution for which the data will be given as a structured mesh (resolution is the same in both horizontal directions)
  :param bathy_file: file for the bathymetry (has to be larger than the mesh provided by the donor mesh)
  """

  # Load bathy and donor data
  bathy_data = Dataset(bathy_file, 'r', format='NETCDF4')  
  
  # Get values from bathymetry file
  try:
    bathy_x = bathy_data.variables['x']
    bathy_y = bathy_data.variables['y']
    bathy = bathy_data.variables['z']
  except:
    bathy_x = bathy_data.variables['lon']
    bathy_y = bathy_data.variables['lat']
    try:
      bathy = bathy_data.variables['elevation']
    except:
      bathy = bathy_data.variables['z']
  
  # Check spatial resolution and set it equal to the bathymetry resolution if 0
  if (spatial_resolution == 0.0):
    
    # Calculate bathymetry resolution (uniform grid is assumed!); ensure that dx = dy
    bathy_x_res = bathy_x[1]-bathy_x[0]
    bathy_y_res = bathy_y[1]-bathy_y[0]

    # check that the grid is uniform
    dx = np.diff(bathy_x[:])
    dy = np.diff(bathy_y[:])
    dx, dy = np.abs(dx), np.abs(dy)

    if not np.allclose(dx, dx[0]) or not np.allclose(dy, dy[0]):
        raise ValueError("Bathymetry grid is not uniform!")

    if (np.abs(bathy_x_res - bathy_y_res) <= 1E-8):
      bathy_resolution = bathy_x_res
    else: 
      raise ValueError("The provided bathymetry file has different resolutions for x- and y-coordinates. Please provide a bathymetry file that has the same resolution in both directions.")

    return bathy_resolution, bathy
  
  else:
    return spatial_resolution, bathy
