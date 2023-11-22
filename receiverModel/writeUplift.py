import numpy as np
from .uplift2hysea import write2hysea

"""
Functionality to write uplift data. As of now, only HySEA is a viable receiver.
"""

def write_uplift(eg_uplift, eg_x, eg_y, receiver, donor, filtername, casename, spatial_resolution):
  """
  Choose receiver model to write the uplift for
  
  :param eg_uplift: uplift data
  :param eg_x: x/lon-coordinates of the exchange grid
  :param eg_y: y/lat-coordinates of the exchange grid
  :param receiver: receiver model for which the data is written 
  :param donor: donor model that was used 
  :param filtername: prefix (character) for the name of the netCDF file  
  :param casename: prefix (character) for the name of the netCDF file  
  :param spatial_resolution: spatial resolution. Will be used as prefix (character) for the name of the netCDF file  
  """
  if (receiver == 'hysea'):
    write2hysea(eg_uplift, eg_x, eg_y, donor, filtername, casename, spatial_resolution)
  else:
    raise NotImplementedError("The provided receiver model is unknown. Possible receiver models include: hysea")

