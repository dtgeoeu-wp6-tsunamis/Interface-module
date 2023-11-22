import numpy as np
from scipy import interpolate
from math import sqrt, ceil
from scipy import fft
import time

"""
Functionalities to filter the uplift data. Consists of no filter as of now.
"""


# Generic filtering function. Needs type of filter (character: none or kajiura), as well as uplift and bathymetry as inputs
def filter_uplift(choose_filter, uplift, bathymetry, spatial_resolution):
  if (choose_filter == 'none'):
      print('Uplift data will not be filtered')
      return uplift
  else:
      print('No known filter was used. Will opt for no filtering as a default.')
      return uplift

