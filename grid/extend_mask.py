#!/usr/bin/env python3

import os
import sys
from netCDF4 import Dataset
import numpy as np

mit = Dataset('mitgcm_mask.nc', 'r')
lat = mit.variables["lat"][:]
lon = mit.variables["lon"][:]
mask = mit.variables["mask"][:]

nx,ny = mask.shape

src = Dataset('gridfile.nc','r')
glon = src.variables["grid_center_lon"][:]
glat = src.variables["grid_center_lat"][:]

dst = Dataset('maskfile.nc','w')
dst.setncatts(mit.__dict__)
dst.createDimension('lon', len(src.dimensions['grid_xsize']))
dst.createDimension('lat', len(src.dimensions['grid_ysize']))

replace_names = { 'x' : 'lon',
		  'y' : 'lat' }
for name, variable in mit.variables.items():
    x = dst.createVariable(name, variable.datatype,
		    (replace_names[x] for x in variable.dimensions))
    x.setncatts(mit[name].__dict__)

dst.variables["lon"][:] = glon
dst.variables["lat"][:] = glat

dmask = np.zeros_like(glon, dtype='u1')

print(np.shape(dmask))
print(np.shape(mask))

dmask[35:368+35,:611] = mask[:,49:]

dst.variables["mask"][:] = dmask

mit.close( )
src.close( )
dst.close( )
