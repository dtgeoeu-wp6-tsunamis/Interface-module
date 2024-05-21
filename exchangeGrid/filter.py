import numpy as np
import os
import time

from math import sqrt, ceil
from scipy import interpolate, fft

"""
Functionalities to filter the deformation data. Consists of a Kajiura or no filter.

Contains the following functionalities:

* precompute_R          precompute sea surface response as in Kajiura 1963
* precompute_σ          precompute σ for apply_kajiura_fft
* compute_filter          compute filter matrix for the use in apply_kajiura_fft
* apply_kajiura_fft       calculate FFT to get sea surface response
* use_kajiura_filter      main Kajiura filter routine
* filter_deformation    parent routine where the filter is chosen (currently: kajiura or none)
"""

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


#**********************************************************************************
# The following routines are functionalities for the Kajiura filter which are based on a julia version from LMU

def precompute_R():
  """
  Function to (pre)calculate the initial surface displacement R (follows the formula given in Kajiura 1963)
  """
  
  Δr = .01
  l_r = np.arange(0.,10 + Δr, Δr)
  
  # Auxiliary function for approximation of the infinite series
  def R(r):
      # R(r) ≈ 1/π * Σ [(-1)^n * (2n + 1) / ((2n + 1)^2 + r^2)^(3/2)] for n = 0 to 10000
      result = 0.0
      for n in range(10001):
          frac = (-1) ** n * (2 * n + 1) / ((2 * n + 1) ** 2 + r ** 2) ** 1.5
          result += frac
      return result / np.pi
    
  l_R = [R(r) for r in l_r]  
  itp = interpolate.interp1d(l_r, l_R, kind='nearest', fill_value='extrapolate')
  
  return itp



def precompute_σ(h_min, h_max, Δx, Δy, precalc_R, n_h=20.0):
  """
  Function to precompute σ, which will be used later.
                         h^2
  σ = ———————————————————————————————————————
        Σ   Σ  R(√[(nΔx)^2 + (mΔy)^2]/h)ΔxΔy
       n∈N m∈M
  
  N = {n ∈ Z| |nΔx| ≤ n_h * h}
  M = {m ∈ Z| |mΔy| ≤ n_h * h}
  
  Here, R(r) is defined as above.
  As the domain of n and m are point symmetric around 0, only ~1/4th of the sum element have to be evaluated:
  * One quadrant (n > 0, m > 0), multiplied by 4
  * One half of the x-axis (n > 0), multiplied by 2
  * One hald of the y-axis (m > 0), multiplied by 2
  * The point at the origin (n, m) = (0, 0)  
  
  :param h_min: minimum water column height (if negative, 0 will be used)
  :param h_max: maximum water column height
  :param Δx: spatial resolution (in longitude direction)
  :param Δy: spatial resolution (in latitude direction)
  :param precalc_R: precomputed values for initial surface displacement
  :param n_h: size of box in which σ will be calculated (default: 20)
  """
  
  Δh = (h_max - max(0.0, h_min)) / 10
  if Δh == 0.0:
      Δh = 1.0
  l_h = np.arange(max(0.0001, h_min), h_max + Δh, Δh)

  # Auxiliary function to calculate σ
  def σ(h):
    σ_inv = 0.0
    n_max = ceil(n_h * h / Δx)
    m_max = ceil(n_h * h / Δy)

    N = np.arange(1, n_max + 1)
    M = np.arange(1, m_max + 1)

    for n in N:
        for m in M:
            σ_inv += 4 * precalc_R(sqrt((n * Δx) ** 2 + (m * Δy) ** 2) / h)

    for n in N:
        σ_inv += 2 * precalc_R(n * Δx / h)

    for m in M:
        σ_inv += 2 * precalc_R(m * Δy / h)

    σ_inv += precalc_R(0)

    return h ** 2 / (σ_inv * Δx * Δy)

  l_σ = [σ(h) for h in l_h]

  itp = interpolate.interp1d(l_h, l_σ, kind='cubic', fill_value='extrapolate')
  return itp



def compute_filter(filter_depth, h_max, Δx, Δy, nx, ny, n_h, precalc_R):
  """
  Routine to compute the filter matrix for the FFT
  
  :param filter_depth: Kajiura filter depth
  :param h_max: maximum water column height
  :param Δx: spatial resolution (in longitude direction)
  :param Δy: spatial resolution (in latitude direction)
  :param nx: number of points in longitude direction
  :param ny: number of points in latitude direction
  :param n_h: size of box in which σ will be calculated (default: 20)  
  :param precalc_R: precomputed values for initial surface displacement
  """
  
  filter_nx_half = int(np.ceil(n_h * h_max / Δx / 2))
  filter_ny_half = int(np.ceil(n_h * h_max / Δy / 2))

  filter_matrix = np.zeros((ny, nx), dtype=float)

  assert filter_nx_half < nx
  assert filter_ny_half < ny

  for x in range(-filter_nx_half, filter_nx_half + 1):
      for y in range(-filter_ny_half, filter_ny_half + 1):
          #if (h[y,x] > 0.0):
          filter_matrix[(ny + y) % ny, (nx + x) % nx] = precalc_R(sqrt((x * Δx) ** 2 + (y * Δy) ** 2) / filter_depth)  
  
  return filter_matrix        



def apply_kajiura_fft(bathymetry, deformation, η, h_max, Δx, Δy, precalc_σ, precalc_R, filter_depth, water_level=0.0, n_h=20.0):
  """
  Kajiura filter functionality with FFT
  
  :param bathymetry: bathymetry data (interpolated and updated with deformation)
  :param deformation: deformation data from the donor model
  :param η: sea surface response (output)
  :param h_max: maximum water column height
  :param Δx: spatial resolution (in longitude direction)
  :param Δy: spatial resolution (in latitude direction)
  :param precalc_σ: precomputed values for σ
  :param precalc_R: precomputed values for initial surface displacement
  :param filter_depth: Kajiura filter depth
  :param water_level: water level displacement (default: zero)
  :param n_h: size of box in which σ will be calculated (default: 20)  
  """
    
  ny, nx = η.shape
  η_aux = np.zeros((ny,nx))
  σ = precalc_σ(filter_depth)
    
  for x in range(nx):
      for y in range(ny):
          h_yx = max(0.0, water_level - bathymetry[y, x]) # set height to 0 on land
          
          if ((h_yx > 0.1) and (abs(deformation[y, x]) > 1e-8)): # only calculate sea surface response if in ocean
            η_aux[y, x] = σ * Δx * Δy / filter_depth ** 2 * deformation[y, x]
            
  # compute filter matrix
  filter_matrix = compute_filter(filter_depth, h_max, Δx, Δy, nx, ny, n_h, precalc_R)
  
  # FFT 
  Η = fft.fft2(η_aux)
  Filter = fft.fft2(filter_matrix)
  Η *= Filter
  η_complex = fft.ifft2(Η)
  η[:] = np.real(η_complex)
  
  

def use_kajiura_filter(deformation, bathymetry, spatial_resolution):
  """
  Functionality for the Kajiura filter
  
  :param deformation: deformation data from the donor model
  :param bathymetry: bathymetry data (interpolated and updated with deformation)
  :param spatial_resolution: spatial resolution for both longitude and latitude direction    
  """
  
  Ntime = np.shape(deformation)[0] # number of timesteps
  Ny = np.shape(deformation)[1]
  Nx = np.shape(deformation)[2]

  current_disp = np.zeros((Ny, Nx))
  current_η = np.zeros((Ny, Nx))
  current_deformation_diff = np.zeros((Ny, Nx))
  
  print(f"Precomputing parts for the filtering.".center(column_size))
  start = time.time()
  precalc_R = precompute_R()
  precalc_σ = precompute_σ(-np.max(bathymetry), -np.min(bathymetry), 
                                                      spatial_resolution, spatial_resolution, precalc_R)
  stop = time.time()
  print(f"Precomputations performed. It took {stop - start} seconds.".center(column_size))

  filtered_deformation = np.zeros_like(deformation)
  frmt = len(str(Ntime)) # formatter for terminal output

  print("Starting filtering of deformation data.".center(column_size))
  
  for t in range(Ntime):
    start = time.time()  
    
    current_bathymetry = bathymetry[t]
    current_deformation = deformation[t]
    
    # Get indices for location of largest deformation and use those to define Kajiura depth
    maxdeform_indices = np.unravel_index(np.argmax(np.abs(current_deformation), axis=None), (Ny, Nx))
    # Set up bathymetry at largest deformation as Kajiura depth
    kajiura_depth = np.abs(current_bathymetry[maxdeform_indices[0], maxdeform_indices[1]])
    
    print(f"Filter is in timestep {t+1:{frmt}d} of {Ntime} with a filtering depth of {kajiura_depth} m.".center(column_size))

    # Start filterting
    current_η[:] = 0.0
    current_deformation_diff[:] = current_deformation - current_disp  
    apply_kajiura_fft(current_bathymetry, current_deformation_diff, current_η, 
                               -np.min(bathymetry), spatial_resolution, spatial_resolution, 
                               precalc_σ, precalc_R, kajiura_depth)
    current_disp = deformation[t]

    filtered_deformation[t] = current_η

    stop = time.time()
    print(f"Timestep {t+1} took {stop - start} seconds.".center(column_size))

  return filtered_deformation



#**********************************************************************************
# Generic filtering function
def filter_deformation(choose_filter, deformation, bathymetry, spatial_resolution):
  """
  Generic filtering function. Needs type of filter (character: none or kajiura), as well as deformation and bathymetry as inputs.
   
  :param choose_filter: character for choosing the filter
  :param deformation: deformation data from the donor model
  :param bathymetry: bathymetry data (interpolated and updated with deformation)
  :param spatial_resolution: spatial resolution for both longitude and latitude direction
  """
  
  if (choose_filter == "kajiura"):
      print("Using a Kajiura filter for smoothing the deformation data.".center(column_size))
      start = time.time()
      filtered_deformation = use_kajiura_filter(deformation, bathymetry, spatial_resolution)
      stop = time.time()
      print(f"The Kajiura filter took {stop - start} seconds.\n".center(column_size))
      return filtered_deformation
  elif (choose_filter == "none"):
      print("Uplift data will not be filtered.\n".center(column_size))
      return deformation
  else:
      print("No known filter was used. Will opt for no filtering as a default.\n".center(column_size))
      return deformation

