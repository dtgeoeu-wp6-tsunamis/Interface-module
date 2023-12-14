import numpy as np
from netCDF4 import Dataset
from datetime import datetime

"""
Functionality to write deformation data to a netCDF file usable by HySEA (T-HySEA 1.1.1)

Contains the following functionalities:

* write2hysea   write HySEA output (netCDF file)

"""

def write2hysea(eg_deformation, eg_x, eg_y, donor, filtername, casename, spatial_resolution):
  """
  Write deformation to a netCDF file viable for HySEA
  
  :param eg_deformation: deformation data
  :param eg_x: x/lon-coordinates of the exchange grid
  :param eg_y: y/lat-coordinates of the exchange grid
  :param donor: donor model that was used 
  :param filtername: prefix (character) for the name of the netCDF file  
  :param casename: prefix (character) for the name of the netCDF file
  :param spatial_resolution: spatial resolution. Will be used as prefix (character) for the name of the netCDF file  
  """
  if (filtername == 'none'):
    filtername = 'unfiltered'
  deformationfile = f"{casename}_{donor}_{filtername}_dx{int(spatial_resolution)}_deformation.nc"
  Nrow = len(eg_y)
  Ncolumn = len(eg_x)
  Ntime = np.shape(eg_deformation)[0]

  # Create the netCDF file and fill it
  ds = Dataset(deformationfile, 'w', format='NETCDF4')
  ds.title = f"{donor} model outputs converted to a structured mesh by interpolation"
  ds.history = "File written using netCDF4 Python module"
  today = datetime.today()
  ds.description = "Created " + today.strftime("%d/%m/%y")
  lon_dim = ds.createDimension('x', Ncolumn)
  lat_dim = ds.createDimension('y', Nrow)
  time_dim = ds.createDimension('time', None)
  
  time = ds.createVariable('time', 'f4', ('time',))
  time.units = 'time step'
  latitude = ds.createVariable('y', 'f8', ('y',))
  latitude.units = 'degrees north (WGS84)'
  latitude.long_name = 'latitude'
  longitude = ds.createVariable('x', 'f8', ('x',))
  longitude.units = 'degrees east (WGS84)'
  longitude.long_name = 'longitude'
  longitude[:] = eg_x
  latitude[:] = eg_y
  
  z = ds.createVariable('z', 'f4', ('time', 'y', 'x'))
  
  for t in range(Ntime):
    time[t] = t
    z[t,:,:] = eg_deformation[t]
  
  ds.close()
