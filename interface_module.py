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

Versions: 
0.1 (10/23) M. Bänsch (UHAM)    Initial version (SeisSol to HySEA)
0.2 (12/23) M. Bänsch (UHAM)    Added Kajiura Filter; some renaming

*** Instructions for this module ***

This script can be run from the command line and needs some additional arguments.

How to run the script: 
  python interface_module.py --donor donor_model donor_output bathy_file resolution (--receiver receiver_model --filter filter --casename casename)

Arguments that need/can to be provided:
  * --donor donor_model           where donor_model = seissol, shaltop, bingclaw (all lower case!) 
  * donor_output                       name of the output file or path from the donor model
  * bathy_file                             name of the bathymetry file. Domain has to be larger compared to the domain from the donor model
  * resolution                             spatial resolution for the interpolation (will be used for both x- and y-coordinates)
  * --receiver receiver_model     (optional) receiver model (as of now, only hysea is available)
  * --filter filter                          (optional) filter for the deformation data where filter = none, kajiura; default: none
  * --casename casename        (optional) string to append the filename with 

Regarding the donor_output, note that SeisSol requires a filename, while Bingclaw requires the path to the directory where the (ESRI ASCII) output files are located
"""

# Generic modules that are needed 
import argparse
import time
import donorModels.donorInterface as donorInterface
import exchangeGrid.interpolateBathy as interpolateBathy
import exchangeGrid.filter as filtering
import receiverModel.writeInterpolatedBathy as writeBathy
import receiverModel.writeDeformation as writeUplift
import numpy as np

#TODO: include functionality for parameter file ?
#TODO: include time for output files

# Define arguments the script has to be called with
parser = argparse.ArgumentParser(
    description="Interface module for the DT-GEO WP6 tsunami workflow (source to wave-filter)"
)
parser.add_argument(
    "--donor",
    nargs=1,
    help="donor model; seissol, shaltop, bingclaw (all lower case)",
    required=True,
)
parser.add_argument("donor_output", help="name of the output file(s) from the donor model")
parser.add_argument("bathy_file", help="name of the bathymetry file")
parser.add_argument(
    "--receiver",
    nargs=1,
    help="receiver model; hysea (all lower case)",
    default="hysea",
)
parser.add_argument("resolution", help="spatial resolution for both horizontal directions (in m)" )
parser.add_argument("--include_horizontal_deformation", 
    help="spatial resolution for the interpolation (will be used for both x- and y-coordinates)", 
    default=False)
parser.add_argument("--filter", 
    help="filter for the deformation data where filter = none, kajiura; default: none",
    default='none')
parser.add_argument("--casename", 
    help="string to append the filename with",
    default='src2waveOut')

args = parser.parse_args()
spatial_resolution = float(args.resolution)
incl_horizontal = args.include_horizontal_deformation
filtername = args.filter
casename = args.casename

print("************************************\n")
print("       WP6 Interface module       \n")

"""
Stage 1: Get data from donor model 
"""
print("************************************\n")
print("Entering Stage 1: getting the data from the donor model.\n")
print("************************************\n")

start = time.time()

donor_deformation, donor_x, donor_y = donorInterface.get_donorModel(args.donor[0], args.donor_output, spatial_resolution, incl_horizontal)

stop = time.time()
print(f"\nStage 1 completed (it took {stop - start} seconds).\n")
print("************************************\n")

"""
Stage 2: interpolate bathymetry data to same grid as the deformation. Filter the deformation if desired.
"""

print("Entering Stage 2: processing the data.\n")
print("************************************\n")

start = time.time()

interpolated_bathymetry = interpolateBathy.get_interpolatedBathy(args.bathy_file, donor_x, donor_y, donor_deformation)
eg_deformation = filtering.filter_deformation(filtername, donor_deformation, interpolated_bathymetry, spatial_resolution)

stop = time.time()
print(f"\nStage 2 completed (it took {stop - start} seconds).\n")
print("************************************\n")

"""
Stage 3: Write interpolated bathymetry and deformation to corresponding netCDF files
"""
print("Entering Stage 3: writing output.\n")
print("************************************")

start = time.time()

writeBathy.write_interpolatedBathy(interpolated_bathymetry, donor_x, donor_y, casename,  Ntime=np.shape(eg_deformation)[0])
writeUplift.write_deformation(eg_deformation, donor_x, donor_y, args.receiver, args.donor[0], filtername, casename, spatial_resolution)

stop = time.time()
print(f"\nStage 3 completed (it took {stop - start} seconds).\n")
print("************************************\n")
