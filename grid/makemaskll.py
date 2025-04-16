#!/usr/bin/env python3

import sys
import datetime
import numpy as np
from netCDF4 import Dataset

try:
    src = Dataset(sys.argv[1],'r')
except:
    try:
        src = Dataset('landfile.nc','r')
    except:
        print('Cannot open land file.')
        sys.exit(1)

try:
    nlon = src.dimensions['lon'].size
    nlat = src.dimensions['lat'].size
except:
    nlon = src.dimensions['x'].size
    nlat = src.dimensions['y'].size

dst = Dataset('maskfile.nc', 'w')

dst.title = 'CHyM maskfile'
dst.datetime = datetime.datetime.utcnow( ).isoformat( )

# Define dimension
dst.createDimension("lon",nlon)
dst.createDimension("lat",nlat)
# Define variables
lon = dst.createVariable("lon", 'f4', ["lat", "lon"])
lon.standard_name = "longitude"
lon.long_name = "Longitude"
lon.units = "degrees_east"
lat = dst.createVariable("lat", 'f4', ["lat", "lon"])
lat.standard_name = "latitude"
lat.long_name = "Latitude"
lat.units = "degrees_north"
mask = dst.createVariable("mask", 'u1', ["lat", "lon"])
mask.standard_name = "sea_binary_mask"
mask.long_name = "Sea binary mask"
mask.units = "1"
mask.coordinates = "lat lon"

# Write variables in file
lat[:] = src.variables['lat'][:]
lon[:] = src.variables['lon'][:]
ocean = np.logical_or(src.variables['luc'][:] == 15,
                      src.variables['luc'][:] == 14)
mask[:] = np.where(ocean, 1, 0)
dst.close( )
