#!/usr/bin/env python3

import sys
import os
import numpy as np
import pandas as pd
import geopandas as gp
from shapely.geometry import Point, Polygon
from netCDF4 import Dataset
from tqdm import tqdm

def coord_to_point(alon, alat):
    lon = alon.data.flatten( )
    lat = alat.data.flatten( )
    return gp.points_from_xy(lon, lat, crs = 'epsg:4326')

argerror = False
try:
    gridfile = Dataset(sys.argv[1], 'r')
except:
    try:
        gridfile = Dataset('gridfile.nc', 'r')
    except:
        print('No gridfile provided.')
        argerror = True
try:
    demfile = Dataset(sys.argv[2], 'r')
except:
    try:
        demfile = Dataset('HydroSheds_15s_Mediterraneo.nc', 'r')
    except:
        print('No demfile provided')
        argerror = True
try:
    maskfile = Dataset(sys.argv[3], 'r')
except:
    try:
        maskfile = Dataset('maskfile.nc', 'r')
    except:
        print('No maskfile provided')
        argerror = True

try:
    rivers = gp.read_file(sys.argv[3])
except:
    try:
        rivers = gp.read_file('filtered.shp')
    except:
        print('No river shapes provided')
        argerror = True

if argerror:
    print('Usage: ', sys.argv[0], ' grid.nc dem.nc mask.nc rivers.shp')
    sys.exit(-1)

print('Reading input DEM data:')
lon, lat = np.meshgrid(demfile.variables['lon'][:],
                       demfile.variables['lat'][:])

print('   Creating input GeoDataFrame...')
in_coords = coord_to_point(lon, lat)
print('   ...coords...')
in_values = demfile.variables['dem'][:].data.flatten( )
in_values = np.where(in_values < 0.0, 0.0, in_values)
in_values = np.where(in_values > 9000.0, 0.0, in_values)
in_dem = gp.GeoDataFrame(data = {'dem' : in_values,
                                 'geometry': in_coords},
                         crs = 'epsg:4326')
print(    'Done.')
demfile.close( )

print('Reading mask.')
mask = maskfile.variables['mask'][:]
maskfile.close( )

print('Reading output grid.')
gridlat = gridfile.variables['grid_center_lat'][:]
gridlon = gridfile.variables['grid_center_lon'][:]

bl_lon = gridfile.variables['grid_corner_lon'][:, :, 0].flatten( )
bl_lat = gridfile.variables['grid_corner_lat'][:, :, 0].flatten( )
br_lon = gridfile.variables['grid_corner_lon'][:, :, 1].flatten( )
br_lat = gridfile.variables['grid_corner_lat'][:, :, 1].flatten( )
tr_lon = gridfile.variables['grid_corner_lon'][:, :, 2].flatten( )
tr_lat = gridfile.variables['grid_corner_lat'][:, :, 2].flatten( )
tl_lon = gridfile.variables['grid_corner_lon'][:, :, 3].flatten( )
tl_lat = gridfile.variables['grid_corner_lat'][:, :, 3].flatten( )

bl = coord_to_point(gridfile.variables['grid_corner_lon'][:, :, 0],
                    gridfile.variables['grid_corner_lat'][:, :, 0])
br = coord_to_point(gridfile.variables['grid_corner_lon'][:, :, 1],
                    gridfile.variables['grid_corner_lat'][:, :, 1])
tr = coord_to_point(gridfile.variables['grid_corner_lon'][:, :, 2],
                    gridfile.variables['grid_corner_lat'][:, :, 2])
tl = coord_to_point(gridfile.variables['grid_corner_lon'][:, :, 3],
                    gridfile.variables['grid_corner_lat'][:, :, 3])
pointlist = [x for x in zip(bl, br, tr, tl)]
polygons = [x for x in map(Polygon, [p for p in pointlist])]

gridfile.close( )

print('Creating output DataFrame...')
nlat,nlon = mask.shape
mask1d = mask.data.flatten( )
out_dem = pd.DataFrame(
        data = {'dem' : np.zeros(nlat*nlon),
                'low' : np.zeros(nlat*nlon),
                'mask' : mask1d,
                'bl_lon' : bl_lon,
                'bl_lat' : bl_lat,
                'br_lon' : br_lon,
                'br_lat' : br_lat,
                'tr_lon' : tr_lon,
                'tr_lat' : tr_lat,
                'tl_lon' : tl_lon,
                'tl_lat' : tl_lat,
                })
print('Done')

print('Created output dem data, size ', nlat, ',', nlon)
print('    Filling it...')

def process_point(row):
    if row[2] == 0:
        a = in_dem
        poly = Polygon([Point(row[3],row[4]),
                        Point(row[5],row[6]),
                        Point(row[7],row[8]),
                        Point(row[9],row[10])])
        b = gp.GeoDataFrame(geometry = [poly], crs = 'epsg:4326')
        pos = pd.Index(gp.sjoin(a, b, predicate='within').index)
        if pos.size > 0:
            within = in_dem.loc[pos,'dem'].to_numpy( )
            within = within[within > 0]
            if within.size > 0:
                row[0] = np.maximum(np.quantile(within,0.5), 5.0)
                row[1] = np.maximum(np.quantile(within,0.1), 2.0)
            else:
                row[0] = 5.0
                row[1] = 2.0
        else:
            row[0] = 5.0
            row[1] = 2.0
    return row

tqdm.pandas()
out_dem = out_dem.progress_apply(process_point,
        raw = True, axis = 1, result_type = 'expand')

print('Data filling complete.')

print('Simplify rivers...')
simpler = rivers.simplify(0.1)
print('Done')

print('Etching rivers....')
geo_dem = gp.GeoDataFrame(data = out_dem[['dem','low']],
        geometry = coord_to_point(gridlon, gridlat), crs = 'epsg:4326')
geo_dem["x"] = geo_dem["geometry"].x
geo_dem["y"] = geo_dem["geometry"].y
a = gp.GeoDataFrame(geometry = polygons, crs = 'epsg:4326')
b = gp.GeoDataFrame(geometry = simpler, crs = 'epsg:4326')
pos = pd.Index(pd.unique(gp.sjoin(a, b, predicate='intersects').index))
geo_dem.loc[pos,'dem'] = geo_dem.loc[pos,'low']

dem_out = geo_dem.set_index(["y","x"])["dem"].to_numpy( )
dem_out = np.where(dem_out < 0.0, 0.0, dem_out)
print('Done!')

print('Writing output results....')
ncout = Dataset('demfile.nc', 'w')
ncout.createDimension('lon', nlon);
ncout.createDimension('lat', nlat);
lat = ncout.createVariable('lat', 'f4', ('lat', 'lon'))
lat.units = 'degrees_north'
lat.standard_name = 'latitude'
lat.long_name = 'Latitude'
lon = ncout.createVariable('lon', 'f4', ('lat', 'lon'))
lon.units = 'degrees_east'
lon.standard_name = 'longitude'
lon.long_name = 'Longitude'
dem = ncout.createVariable('dem', 'f4', ('lat', 'lon'))
dem.units = 'm'
dem.long_name = 'Digital Elevation model'
dem.standard_name = 'elevation'
dem.coordinates = 'lat lon'
lat[:] = gridlat
lon[:] = gridlon
dem[:] = dem_out.reshape((nlat, nlon))
ncout.close();

print('Done!')
