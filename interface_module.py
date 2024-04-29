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
  python interface_module.py --donor donor_model --CRS_reference lon, lat donor_output bathy_file resolution (--receiver receiver_model --filter filter --casename casename)

Arguments that need/can to be provided:
  * --donor donor_model           where donor_model = seissol, shaltop, bingclaw (all lower case!) 
  * --CRS_reference                    CRS coordinates reference (lon, lat of lower left corner of domain)
  * donor_output                       name of the output file or path from the donor model
  * bathy_file                             name of the bathymetry file. Domain has to be larger compared to the domain from the donor model
  * resolution                             (optional) spatial resolution the donor output will be interpolated to (will be used for both x- and y-coordinates; has to be provided in meters)
  * --receiver receiver_model     (optional) receiver model (as of now, only hysea is available)
  * --filter filter                          (optional) filter for the deformation data where filter = none, kajiura; default: none
  * --casename casename        (optional) string to append the filename with 

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
import os

#TODO: include functionality for parameter file ?
#TODO: Check if CRS_reference can be removed after SeisSol and SHALTOP have a CRS reference; In this case, the CRS_reference in each donor module needs to be changed

# Some definitions for a nice print on the terminal
column_size = os.get_terminal_size().columns
asterisk_fill = "*" * column_size


# Define arguments the script has to be called with
parser = argparse.ArgumentParser(description="Interface module for the DT-GEO WP6 tsunami workflow (source to wave-filter)")
parser.add_argument("-d", "--donor",
    help="donor model; seissol, shaltop, bingclaw (all lower case)",
    required=True,)
parser.add_argument("-crs", "--CRS_reference", 
    metavar=("[lon, lat]"),
    nargs=2,
    type=str,
    help="CRS reference coordinates (list of longitude and latitude of lower left corner of the domain); required for SHALTOP and SeisSol.",
    default = ["0", "0"])
parser.add_argument("donor_output", help="name of the output file(s) from the donor model")
parser.add_argument("bathy_file", help="name of the bathymetry file")
parser.add_argument(
    "-r", "--receiver",
    help="receiver model; hysea (all lower case)",
    default="hysea",)
parser.add_argument("--resolution", help="spatial resolution for both horizontal directions (in m)", 
    default=0.0)
parser.add_argument("--include_horizontal_deformation", 
    help="horizontal deformation handle (for SeisSol)", 
    default=False)
parser.add_argument("-f", "--filter", 
    help="filter for the deformation data where filter = none, kajiura; default: none",
    default='none')
parser.add_argument("-c", "--casename", 
    help="string to append the filename with",
    default='src2waveOut')

args = parser.parse_args()
spatial_resolution = float(args.resolution)
incl_horizontal = args.include_horizontal_deformation
filtername = args.filter
casename = args.casename
CRS_reference = args.CRS_reference

print(asterisk_fill + "\n")
print("WP6 Interface module\n".center(column_size))

"""
Stage 1: Get data from donor model 
"""
print(asterisk_fill + "\n")

print("Entering Stage 1: getting the data from the donor model.\n".center(column_size))
print(asterisk_fill + "\n")


start = time.time()

donor_deformation, donor_x, donor_y, donor_time, donor_bathy = donorInterface.get_donorModel(args.donor, args.donor_output, spatial_resolution, CRS_reference, args.bathy_file, incl_horizontal)

stop = time.time()
print((f"Stage 1 completed. It took {stop - start} seconds.\n").center(column_size))
print(asterisk_fill + "\n")


"""
Stage 2: Create exchange grid from donor and bathymetry coordinates, interpolate bathymetry data to same grid as the deformation (if necessary). Filter the deformation if desired.
"""

print("Entering Stage 2: processing the data.\n".center(column_size))
print(asterisk_fill + "\n")


start = time.time()


eg_tmp_deformation, eg_x, eg_y, eg_bathymetry =  exchangeGridInterface.createExchangeGrid(args.bathy_file, donor_x, donor_y, donor_deformation)
eg_deformation = filtering.filter_deformation(filtername, eg_tmp_deformation, eg_bathymetry, spatial_resolution)

stop = time.time()
print(f"Stage 2 completed. It took {stop - start} seconds.\n".center(column_size))
print(asterisk_fill + "\n")



"""
Stage 3: Write interpolated bathymetry and deformation to corresponding netCDF files
"""
print("Entering Stage 3: writing output.\n".center(column_size))
print(asterisk_fill + "\n")

start = time.time()

if (spatial_resolution > 0.0):
  writeBathy.write_interpolatedBathy(args.receiver, eg_bathymetry, eg_x, eg_y, casename,  Ntime=np.shape(eg_deformation)[0])
writeUplift.write_deformation(eg_deformation, eg_x, eg_y, donor_time, args.receiver, args.donor, filtername, casename, spatial_resolution)

stop = time.time()
print(f"Stage 3 completed. It took {stop - start} seconds.\n".center(column_size))
print(asterisk_fill + "\n")


