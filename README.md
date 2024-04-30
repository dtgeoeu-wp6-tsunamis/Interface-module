# Interface module for the WP6 tsunami workflow

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

## Instructions for this module

This script can be run from the command line and needs some additional arguments.

How to run the script: 
  python interface_module.py --donor donor_model --CRS_reference lon, lat donor_output bathy_file (--resolution resolution --only_donor_domain --receiver receiver_model --filter filter --casename casename --include_horizontal_deformation)

Arguments that need/can to be provided:
  * --donor donor_model           where donor_model = seissol, shaltop, bingclaw (all lower case!) 
  * --CRS_reference                    CRS coordinates reference (lon, lat of lower left corner of domain)
  * donor_output                       name of the output file or path from the donor model
  * bathy_file                             name of the bathymetry file. Domain has to be larger compared to the domain from the donor model
  * --resolution resolution          (optional) spatial resolution the donor output will be interpolated to (will be used for both x- and y-coordinates; has to be provided in meters)
  * --only_donor_domain          (optional) handle to only use the domain given by the donor model (False by default)
  * --receiver receiver_model     (optional) receiver model (as of now, only hysea is available)
  * --filter filter                          (optional) filter for the deformation data where filter = none, kajiura; default: none
  * --casename casename        (optional) string to append the filename with 
  * --include_horizontal_deformation     (optional) handle whether to include horizontal deformations (False by default; only for SeisSol)

# Required Python packages
To run the Interface module, the following Python packages have to be installed on the system. 

Standard Packages:
  * argparse
  * datetime
  * math
  * numpy
  * os
  * pathlib
  * pyproj
  * scipy
  * time 

Non-standard packages:
  * netCDF4
  * pyvista
  * seissolxdmf
  * vtk

# To Dos:
  * CRS references for SeisSol and SHALTOP need to be implemented and then included in the interface module
  * (potentially) include parameter file support
