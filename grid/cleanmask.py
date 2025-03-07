#!/usr/bin/env python3

import os
import sys
from netCDF4 import Dataset
import numpy as np
import regionmask as rg
import geopandas as gp

basins = '/home/esp-shared-b/HydroSheds/hydrosheds/basins'
med = 'hybas_eu_lev02_v1c.shp'
afr = 'hybas_af_lev03_v1c.shp'

ds = Dataset(sys.argv[1], 'r+')

maskdata = ds.variables["mask"][:].data
lat = ds.variables["lat"][:]
lon = ds.variables["lon"][:]

cleanmask = np.where(maskdata > 1,0,maskdata)

def roundres(op,v,r):
    if op == min:
        return (int(op(v)/r)-1)*r
    else:
        return (int(op(v)/r)+1)*r

def orthemall(mm):
    res = mm[0].notnull( )
    for m in mm[1:]:
        res = np.logical_or(res,m.notnull( ))
    return res

medbas = gp.read_file(os.path.join(basins,med))
afrbas = gp.read_file(os.path.join(basins,afr))

memask = medbas[medbas.HYBAS_ID==2020000010]
bsmask = medbas[medbas.HYBAS_ID==2020003440]
namask1 = afrbas[afrbas.HYBAS_ID==1030029810]
namask2 = afrbas[afrbas.HYBAS_ID==1030031860]
namask3 = afrbas[afrbas.HYBAS_ID==1030034170]
namask4 = afrbas[afrbas.HYBAS_ID==1030034270]

maskocean = (cleanmask == 1)
landmask = (cleanmask == 0)
maskmed = rg.mask_geopandas(memask,lon,lat)
maskbls = rg.mask_geopandas(bsmask,lon,lat)
masknaf1 = rg.mask_geopandas(namask1,lon,lat)
masknaf2 = rg.mask_geopandas(namask2,lon,lat)
masknaf3 = rg.mask_geopandas(namask3,lon,lat)
masknaf4 = rg.mask_geopandas(namask4,lon,lat)

totalmask = orthemall((maskmed,maskbls,masknaf1,masknaf2,masknaf3,masknaf4))
bothmask = np.logical_and(landmask,totalmask)
fullmask = np.logical_or(maskocean,bothmask)
cleanmask = np.where(fullmask, cleanmask, 2)

ds.variables["mask"][:] = cleanmask
ds.close( )
