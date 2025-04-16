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

cleanmask = np.where(maskdata > 0, 1, maskdata)

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
cleanmask = np.where(bothmask, 2, cleanmask)

cleanmask[146:149,70] = 2
cleanmask[284:289,219] = 2
cleanmask[290,221:223] = 2
cleanmask[286,247] = 2
cleanmask[268,259] = 2
cleanmask[103,209] = 2
cleanmask[233,198] = 2
cleanmask[229,205] = 2
cleanmask[236,199] = 2
cleanmask[235,200] = 2
cleanmask[206,231] = 2
cleanmask[142,197:199] = 2
cleanmask[178,331] = 2
cleanmask[176,331] = 2
cleanmask[177,353] = 2
cleanmask[88,210] = 2
cleanmask[55:57,446] = 2
cleanmask[57,450] = 2
cleanmask[57,452] = 2
cleanmask[58,451:453] = 2
cleanmask[204,364] = 2
cleanmask[166,336] = 2
cleanmask[209,398] = 2
cleanmask[211,398] = 2
cleanmask[153,408] = 2
cleanmask[153,408] = 2
cleanmask[147,419] = 2
cleanmask[349,438] = 2
cleanmask[349,438] = 2
cleanmask[355,489] = 2
cleanmask[356,489] = 2
cleanmask[358,449] = 2
cleanmask[355,448] = 2
cleanmask[352,449] = 2
cleanmask[353,449] = 2
cleanmask[355:357,448] = 2
cleanmask[387,529] = 2
cleanmask[205,364] = 2
cleanmask[177,331] = 2
cleanmask[177,331] = 2
cleanmask[247,284] = 2
cleanmask[56,448] = 2
cleanmask[55,447:449] = 2
cleanmask[54,447:449] = 2
cleanmask[53,446:449] = 2
cleanmask[52,447:450] = 2
cleanmask[51,447:452] = 2
cleanmask[51,447:452] = 2
cleanmask[50,449] = 2
cleanmask[50:45,450] = 2
cleanmask[184,392:395] = 2


ds.variables["mask"][:] = cleanmask
ds.close( )
