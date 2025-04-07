#!/usr/bin/env python3

import os
import sys
from netCDF4 import Dataset
import numpy as np
import regionmask as rg
import geopandas as gp

src = Dataset(sys.argv[1], 'r+')
luc = src.variables["luc"][:].data
lat = src.variables["lat"][:]
lon = src.variables["lon"][:]
mask = np.where(luc == 15,1,0)

dst = Dataset('maskfile.nc','w')
dst.setncatts(src.__dict__)

for name, dimension in src.dimensions.items():
    dst.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else None))

for name, variable in src.variables.items():
   if name == 'luc':
       x = dst.createVariable('mask', 'u1', variable.dimensions)
       x[:] = mask
       x.standard_name = "sea_binary_mask"
       x.units = "1"
   else:
       x = dst.createVariable(name, variable.datatype, variable.dimensions)
       x.setncatts(src[name].__dict__)
       x[:] = src[name][:]

dst.close( )
src.close( )
