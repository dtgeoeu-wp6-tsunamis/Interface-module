"""
Interface module for the WP6 tsunami workflow

This module consists of three main parts:
  * Reading/getting input from one of the donor models (SeisSol, SHALTOP, Bingclaw)
      and transferring the data onto a structured exchange grid
  * Data handling (e.g. filtering) on the exchange grid level
  * Writing the exchange grid data to a netCDF file readable by the receiver model (HySEA)

Many parts of the routines are based on work and codes provided by LMU, UMA and IPGP. As a result, this filtering module is a joint effort.
Input was provided by:
  * Thomas Ulrich (LMU)
  * Alex González (UMA)
  * Alexis Marboeuf (IPGP)

Latest changes made in 04/24 by M. Bänsch (UHAM)

*** Instructions for this module ***

This script can be run from the command line and needs some additional arguments.

How to run the script: 
  python interface_module.py --donor donor_model --projection proj4string donor_output bathy_file (--resolution resolution --only_donor_domain --receiver receiver_model --filter filter --filtering_depth filtering_depth --casename casename --include_horizontal_deformation)

Arguments that need/can to be provided:
  * --donor donor_model           where donor_model = seissol, shaltop, bingclaw (all lower case!) 
  * --projection proj4_string     parameters (in Proj4 format) for converting from Cartesian to geographic (lon, lat) coordinates
  * donor_output                       name of the output file or path from the donor model
  * bathy_file                             name of the bathymetry file. Domain has to be larger compared to the domain from the donor model
  * --resolution resolution:          (optional) spatial resolution the donor output will be interpolated to (will be used for both x- and y-coordinates; has to be provided in meters)
  * --only_donor_domain          (optional) handle to only use the domain given by the donor model (False by default)
  * --receiver receiver_model     (optional) receiver model (as of now, only hysea is available)
  * --filter filter                          (optional) filter for the deformation data where filter = none, kajiura; default: none
  * --filtering_depth filtering_depth (optional) filtering/Kajiura depth (in m) that is used within the filter.
  * --casename casename        (optional) string to append the filename with 
  * --include_horizontal_deformation     (optional) handle whether to include horizontal deformations (False by default; only for SeisSol)

Regarding the donor_output, note that SeisSol requires a filename, while Bingclaw requires the path to the directory where the (ESRI ASCII) output files are located
"""

# Generic modules that are needed 
import argparse
import time
import donorModels.donorInterface as donorInterface
import exchangeGrid.exchangeGridCreation as exchangeGridCreation
import exchangeGrid.filter as filtering
import receiverModel.writeInterpolatedBathy as writeBathy
import receiverModel.writeDeformation as writeUplift
import numpy as np
import matplotlib.pyplot as plt
import os, sys

#TODO: include functionality for parameter file ?

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


# Define arguments the script has to be called with
parser = argparse.ArgumentParser(description="Interface module for the DT-GEO WP6 tsunami workflow (source to wave-filter)")
parser.add_argument("-d", "--donor",
    help="donor model; seissol, shaltop, bingclaw (all lower case)",
    required=True,)
parser.add_argument("-p", "--projection",
    type=str,
    help="projection parameters (in Proj4 format) for converting from Cartesian to geographic (lon, lat) coordinates",
    default = "epsg:4326")
parser.add_argument("donor_output", help="name of the output file(s) from the donor model")
parser.add_argument("bathy_file", help="name of the bathymetry file")
parser.add_argument(
    "-r", "--receiver",
    help="receiver model; hysea (all lower case)",
    default="hysea",)
parser.add_argument("--resolution", help="spatial resolution for both horizontal directions (in m)", 
    default=0.0)
parser.add_argument("--only_donor_domain", 
    help="handle to only use the domain given by the donor model; default: False", 
    default=False)
parser.add_argument("-f", "--filter", 
    help="filter for the deformation data where filter = none, kajiura; default: none",
    default='none')
parser.add_argument("--filtering_depth", 
    help="filtering/Kajiura depth (in m) that is used within the filter; default: 0",
    default=0.0)
parser.add_argument("-c", "--casename", 
    help="string to append the filename with",
    default='src2waveOut')
parser.add_argument("--include_horizontal_deformation", 
    help="horizontal deformation handle (for SeisSol)", 
    default=False)

args = parser.parse_args()
spatial_resolution = float(args.resolution)
only_donor_domain = args.only_donor_domain
incl_horizontal = args.include_horizontal_deformation
filtername = args.filter
filtering_depth = float(args.filtering_depth)
casename = args.casename
projection = args.projection

print(asterisk_fill + "\n")
print("WP6 Interface module\n".center(column_size))

"""
Stage 1: Get data from donor model 
"""
print(asterisk_fill + "\n")

print("Entering Stage 1: getting the data from the donor model.\n".center(column_size))
print(asterisk_fill + "\n")


start = time.time()

# ------------- CHANGES FOR KAJIURA TESTING -------------
# Commented out line below as it is not needed in simple geometry example for testing
### donor_deformation, donor_x, donor_y, donor_time, donor_bathy, eg_resolution = donorInterface.get_donorModel(args.donor, args.donor_output, spatial_resolution, projection, args.bathy_file, incl_horizontal)
#------------- END of CHANGES FOR KAJIURA TESTING -------------

stop = time.time()
print((f"Stage 1 completed. It took {stop - start} seconds.\n").center(column_size))
print(asterisk_fill + "\n")


"""
Stage 2: Create exchange grid from donor and bathymetry coordinates, interpolate bathymetry data to same grid as the deformation (if necessary). Filter the deformation if desired.
"""

print("Entering Stage 2: processing the data.\n".center(column_size))
print(asterisk_fill + "\n")


start = time.time()

# ------------- CHANGES FOR KAJIURA TESTING -------------
# Example of command to run script for testing Kajiura with simple 2D geometry:
# > python interface_module.py --donor seissol fake/tmp bathyfake.grd --casename kajiuratesting --filter kajiura --filtering_depth 2000

# Commented out line below as all these output values will be imposed
### eg_tmp_deformation, eg_x, eg_y, eg_bathymetry =  exchangeGridCreation.createExchangeGrid(args.bathy_file, donor_x, donor_y, donor_deformation, only_donor_domain)

# Create computational mesh 
eg_res = 100    # Resolution of the mesh (m)
eg_x = np.arange(-25000,25100,eg_res)   # Crete array of x coordinates (m)
eg_y = np.arange(-25000,25100,eg_res)   # Crete array of y coordinates (m)

# Create bathymetry matrix
# and impose a different bathymetry value at the side of the box to allow for the fit in function precompute_σ to work
# If bathymetry has only one value, then the fit function in precompute_σ fails as it has only one point and cannot compute the fit
eg_bathy = np.zeros((1, eg_x.shape[0], eg_y.shape[0]))
eg_bathy[0, :, :] = -2000       # Bathymetry (m)           
eg_bathy[0, :, 0:10] = eg_bathy[0,0,0]+200    # Impose a different bathymetry value at the side of the box (m)

# Set up initial sea floor deformation at the center of the domain
eg_tmp_def = np.zeros((1, eg_x.shape[0], eg_y.shape[0]))
a = 10000       # x extent of the initial sea floor deformation (m)
a_el_half = int(a/2/eg_res)
x_center = int((eg_x.shape[0]-1)/2)
x1 = x_center-a_el_half     
x2 = x_center+a_el_half+1   

b = 10000       # y extent of the initial sea floor deformation (m)
b_el_half = int(b/2/eg_res)
y_center = int((eg_y.shape[0]-1)/2)
y1 = y_center-b_el_half
y2 = y_center+b_el_half+1

print("x and y coordinates of initial deformation")
print(eg_x[x1:x2])
print(eg_y[y1:y2])

eg_tmp_def[0, y1:y2, x1:x2] = 2 # Assign initial sea floor deformation (m) 

"""
# PLOT BATHYMETRY
plt.figure(figsize=(6, 6))
plt.imshow(eg_bathymetry[0], cmap='viridis', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('Deformation')
plt.show()
"""
# Commented out line below to use imposed values of bathymetry and deformation
### eg_deformation = filtering.filter_deformation(filtername, eg_tmp_deformation, eg_bathymetry, eg_resolution, filtering_depth)
eg_deformation = filtering.filter_deformation(filtername, eg_tmp_def, eg_bathy, eg_res, filtering_depth)

print(f"Maximum amplitude of filtered deformation is {eg_deformation.max()} m")

# Plot initial seafloor deformation and filtered deformation
plt.figure(figsize=(20, 6))
ax1 = plt.subplot(131)
plt.imshow(eg_tmp_def[0], cmap='viridis', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('Seafloor deformation (m)')
ax1.set_xlabel("x elements")
ax1.set_ylabel("y elements")

ax2 = plt.subplot(132)
plt.imshow(eg_deformation[0], cmap='viridis', interpolation='nearest')
plt.colorbar(label='Value')
plt.title('Filtered deformation (m)')
ax2.set_xlabel("x elements")
ax2.set_ylabel("y elements")

ax3 = plt.subplot(133)
#plt.plot(eg_y, eg_deformation[0,y_center,:])
plt.plot(eg_x, eg_tmp_def[0,y_center,:], color = "blue", label='Seafloor def x')
#plt.plot(eg_y, eg_tmp_def[0,:,x_center], color="orange", label='Seafloor def y')
plt.plot(eg_x, eg_deformation[0,y_center,:], color = "blue", linestyle='dashed',label='Filtered def along x')
#plt.plot(eg_y, eg_deformation[0,:,x_center], color="orange", linestyle='dashed',label='Filtered def along y')
plt.title("Seafloor vs Filtered deformation along x profile")
plt.legend()
ax3.set_xlabel("x (m)")
ax3.set_ylabel("Deformation (m)")
#plt.ylim((0,1))


plt.show()
sys.exit()
#------------- END of CHANGES FOR KAJIURA TESTING -------------

stop = time.time()
print(f"Stage 2 completed. It took {stop - start} seconds.\n".center(column_size))
print(asterisk_fill + "\n")



"""
Stage 3: Write interpolated bathymetry and deformation to corresponding netCDF files
"""
print("Entering Stage 3: writing output.\n".center(column_size))
print(asterisk_fill + "\n")

start = time.time()

# if (spatial_resolution > 0.0):
writeBathy.write_interpolatedBathy(args.receiver, eg_bathymetry, eg_x, eg_y, casename,  Ntime=np.shape(eg_deformation)[0])
writeUplift.write_deformation(eg_deformation, eg_x, eg_y, donor_time, args.receiver, args.donor, filtername, casename, spatial_resolution)

stop = time.time()
print(f"Stage 3 completed. It took {stop - start} seconds.\n".center(column_size))
print(asterisk_fill + "\n")


