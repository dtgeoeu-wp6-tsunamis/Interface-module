import numpy as np
from .seissol2exchangegrid import get_seissol

"""
Functionality to get deformation data from the respective donor model. As of now, only SeisSol is a viable donor.
"""

def get_donorModel(choose_donormodel, filename, spatial_resolution, include_horizontal=False):
  """
  This function provides the data for the respective donor model.
  
  :param choose_donormodel: donor model 
  :param filename: name of the output file from the donor model
  :param spatial_resolution: spatial resolution for which the data will be given as a structured mesh (resolution is the same in both horizontal directions)
  :param include_horizontal: handle whether to include horizontal deformations (False by default)
  """
  # choose corresponding donor model and get the deformation data
  if (choose_donormodel == 'seissol'):
    donor_deformation, donor_x, donor_y = get_seissol(filename, spatial_resolution, include_horizontal)
    return donor_deformation, donor_x, donor_y
  else:
    raise NotImplementedError("The provided donor model is unknown. Possible donor models include: seissol")
