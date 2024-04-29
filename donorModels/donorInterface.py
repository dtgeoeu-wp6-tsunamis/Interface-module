import numpy as np
from .seissol2exchangegrid import get_seissol
from .shaltop2exchangegrid import get_shaltop
from .bingclaw2exchangegrid import get_bingclaw
from .donorResolution import donor_chooseResolution

"""
Functionality to get deformation data from the respective donor model. As of now, BingClaw, SeisSol and Shaltop are viable donors.

Contains the following functionalities:

* get_donorModel    choose donor model and get the output from said model
"""


def get_donorModel(choose_donormodel, donor_output, spatial_resolution, CRS_reference, bathy_file, include_horizontal=False):
  """
  This function provides the data for the respective donor model.
  
  :param choose_donormodel: donor model 
  :param donor_output: name of the output file or path from the donor model
  :param spatial_resolution: spatial resolution for which the data will be given as a structured mesh (resolution is the same in both horizontal directions)
  :param CRS_reference: CRS reference coordinates (list of longitude and latitude of lower left corner of the domain)  
  :param bathy_file: file for the bathymetry (has to be larger than the mesh provided by the donor mesh)
  :param include_horizontal: handle whether to include horizontal deformations (False by default)
  """
  
  # Perform resolution check
  resolution, donor_bathy =  donor_chooseResolution(spatial_resolution, bathy_file)
  
  # Choose corresponding donor model and get the deformation data
  if (choose_donormodel == 'bingclaw'):
    donor_deformation, donor_x, donor_y, donor_time = get_bingclaw(donor_output, resolution)
    return donor_deformation, donor_x, donor_y, donor_time, donor_bathy
  if (choose_donormodel == 'seissol'):
    donor_deformation, donor_x, donor_y, donor_time = get_seissol(donor_output, resolution, CRS_reference, include_horizontal)
    return donor_deformation, donor_x, donor_y, donor_time, donor_bathy
  if (choose_donormodel == 'shaltop'):
    donor_deformation, donor_x, donor_y, donor_time = get_shaltop(donor_output, resolution, CRS_reference)
    return donor_deformation, donor_x, donor_y, donor_time, donor_bathy
  else:
    raise NotImplementedError("The provided donor model is unknown. Possible donor models include: bingclaw, seissol or shaltop")

