import numpy as np

from .deformation2hysea import write2hysea

"""
Functionality to write deformation data. As of now, only HySEA is a viable receiver.

Contains the following functionalities:

* write_deformation      parent routine; choose receiver model for which the output is written.
"""


def write_deformation(eg_deformation, eg_x, eg_y, eg_time, receiver, donor, filtername, casename, spatial_resolution):
  """
  Choose receiver model to write the deformation for
  
  :param eg_deformation: deformation data
  :param eg_x: x/lon-coordinates of the exchange grid
  :param eg_y: y/lat-coordinates of the exchange grid
  :param eg_time: time for each timestep on the exchange grid
  :param receiver: receiver model for which the data is written 
  :param donor: donor model that was used 
  :param filtername: prefix (character) for the name of the netCDF file  
  :param casename: prefix (character) for the name of the netCDF file  
  :param spatial_resolution: spatial resolution. Will be used as prefix (character) for the name of the netCDF file  
  """
  
  if (receiver == 'hysea'):
    write2hysea(eg_deformation, eg_x, eg_y, eg_time, donor, filtername, casename, spatial_resolution)
  else:
    raise NotImplementedError("The provided receiver model is unknown. Possible receiver models include: hysea")

