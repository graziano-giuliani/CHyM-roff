#!/usr/bin/env python3

from osgeo import gdal
import numpy as np
from netCDF4 import Dataset
from cdo import *
import sys
import os

gdal.DontUseExceptions()

global_dem = 'hyd_glo_dem_15s.tif'
src = gdal.Open(global_dem)
_, dem_xres, _, _, _, dem_yres  = src.GetGeoTransform()
del src

global_luc = 'gbogegeo20.tif'
src = gdal.Open(global_luc)
_, luc_xres, _, _, _, luc_yres  = src.GetGeoTransform()
del src

grid = Dataset('gridfile.nc')

lat = grid.variables['grid_corner_lat'][:]
lon = grid.variables['grid_corner_lon'][:]

dem_xres = np.abs(dem_xres)
dem_yres = np.abs(dem_yres)
minlat = float((np.floor(np.min(lat)/dem_yres) - 0.5)*dem_yres)
maxlat = float((np.ceil(np.max(lat)/dem_yres) + 0.5) *dem_yres)
minlon = float((np.floor(np.min(lon)/dem_xres) - 0.5)*dem_xres)
maxlon = float((np.ceil(np.max(lon)/dem_xres) + 0.5) *dem_xres)
opts = {'projWin' : [ minlon, maxlat, maxlon, minlat ],
        'format' : 'NetCDF',
        'outputType': gdal.GDT_Float32 }
dst = gdal.Translate('HydroSheds_15s.nc',global_dem,**opts)
del dst

luc_xres = np.abs(luc_xres)
luc_yres = np.abs(luc_yres)
minlat = float((np.floor(np.min(lat)/luc_yres) - 0.5)*luc_yres)
maxlat = float((np.ceil(np.max(lat)/luc_yres) + 0.5) *luc_yres)
minlon = float((np.floor(np.min(lon)/luc_xres) - 0.5)*luc_xres)
maxlon = float((np.ceil(np.max(lon)/luc_xres) + 0.5) *luc_xres)
opts = {'projWin' : [ minlon, maxlat, maxlon, minlat ],
        'format' : 'NetCDF',
        'outputType': gdal.GDT_Int16 }
dst = gdal.Translate('GLCC_gboggeo20.nc',global_luc,**opts)
del dst

ds = Dataset('HydroSheds_15s.nc',mode='a')
ds.renameVariable('Band1','dem')
ds.variables['dem'].standard_name = 'elevation'
ds.variables['dem'].long_name = 'Digital Elevation Model'
ds.variables['dem'].units = 'm'
ds.close( )

ds = Dataset('GLCC_gboggeo20.nc',mode='a')
ds.renameVariable('Band1','luc')
ds.variables['luc'].standard_name = 'class_type'
ds.variables['luc'].long_name = 'Class Type'
ds.variables['luc'].units = '1'
ds.close( )

cdo = cdo.Cdo( )
cdo.remaplaf('gridfile.nc',input='GLCC_gboggeo20.nc',output='landfile.nc',
             option='-f nc4 -z zip_9')
