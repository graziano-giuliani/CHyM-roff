#!/usr/bin/env python3

import sys
import os
import numpy as np
from netCDF4 import Dataset

# https://github.com/dengwirda/reach
from reach.filter_hydrosheds import filter_hydrosheds

# This MUST be configured.

resolution = 12000.0

try:
    ds = Dataset(sys.argv[1],'r')
    xlon = ds.variables["lon"][:]
    ylat = ds.variables["lat"][:]
    ds.close( )
except:
    xlon = np.arange(-7.0,48.6,0.1)
    ylat = np.arange(27.0,63.6,0.1)

nlon = xlon.size
nlat = ylat.size

class base: pass

tag_ = "filtered"

args = base()
args.shp_file = os.path.join( # location of Hy-SHEDS
        "dat", "HydroRIVERS_v10_shp", "HydroRIVERS_v10.shp")

args.out_file = os.path.join( # location of filtered
        "out", tag_, "filtered.shp")

args.spc_file = os.path.join( # location for spacing
        "out", tag_, "resolution.nc")

args.msh_tags = os.path.join( # location for meshes
        "msh", "HydroRIVERS_v10_msh", tag_) 

dir_ = os.path.dirname(os.path.abspath(args.spc_file))

if (not os.path.exists(dir_)): os.makedirs(dir_)

args.sph_size = 6371220.  # earth radius
args.flt_endo = True  # strip no-outlet rivers??   
args.box_xmin = np.amin(xlon)
args.box_ymin = np.amin(ylat)
args.box_xmax = np.amax(xlon)
args.box_ymax = np.amax(ylat)

data = Dataset(args.spc_file, "w")
data.createDimension("nlon", nlon)
data.createDimension("nlat", nlat)
data.createVariable("xlon", "f8", ("nlon"))
data["xlon"][:] = xlon
data.createVariable("ylat", "f8", ("nlat"))
data["ylat"][:] = ylat
data.createVariable("vals", "f4", ("nlat", "nlon"))
data["vals"][:] = resolution  # in [m]
data.close()

filter_hydrosheds(args)
