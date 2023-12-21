import numpy as np
from netCDF4 import Dataset
from datetime import datetime

"""
Functionality to write bathymetry data to a netCDF file usable by HySEA (T-HySEA 1.1.1)

Contains the following functionalities:

* write_interpolatedBathy     write interpolated bathymetry output (netCDF file)
"""


def write_interpolatedBathy(interpolated_bathy, eg_x, eg_y, casename, Ntime=1):
  """
  Write the interpolated bathymetry to a netCDF file
  
  :param interpolated_bathy: interpolated bathymetry data
  :param eg_x: x/lon-coordinates of the exchange grid 
  :param eg_y: y/lat-coordinates of the donor mesh 
  :param casename: prefix (character) for the name of the netCDF file
  :param Ntime: timesteps for bathymetry
  """
  
  new_bathyfile = casename + "_bathymetry_interp.nc"
  ds = Dataset(new_bathyfile, 'w', format='NETCDF4')
  ds.Conventions = "CF-1.7"
  ds.title = "Interpolated bathymetry data"
  ds.history = "based on bathymetry for HySEA"
  today = datetime.today()
  ds.description = "Created " + today.strftime("%d/%m/%y")
  lon_dim = ds.createDimension('x', len(eg_x))
  lat_dim = ds.createDimension('y', len(eg_y))
  time_dim = ds.createDimension('time', None)

  time = ds.createVariable('time', 'f4', ('time',))
  time.units = 'time step'
  longitude = ds.createVariable('x', 'f8', ('x',))
  longitude.long_name = 'x'
  longitude.actual_range = [eg_x[0], eg_x[-1]]
  latitude = ds.createVariable('y', 'f8', ('y',))
  latitude.long_name = 'y'
  latitude.actual_range = [eg_y[0], eg_y[-1]]
  longitude[:] = eg_x
  latitude[:] = eg_y
  
  bathy_nc = ds.createVariable('z', 'f4', ('time', 'y', 'x',))
  bathy_nc.long_name = 'z'
  bathy_nc.fill_values = np.nan
  bathy_nc.actual_range = [np.min(interpolated_bathy), np.max(interpolated_bathy)]
  
  for t in range(Ntime):
    time[t] = t
    bathy_nc[t,:,:] = interpolated_bathy[t]
  ds.close()
