# Procedure

This file describes the procedure to create input files for the CHyM preproc.

# NOTE!

IT IS USER TASK TO EXAMINE THE python SCRIPTS AND VERIFY ALL THE PACKAGES
IMPORTED ARE AVAILABLE ON HER SYSTEM.

## Create the Gridfile

Create grid file using the python script, either onto regular
latitude-longitude grid or an ORCA NEMO tripolar grid

### Regular lat-lon grid

Edit the **makegrid.yaml** input file for extremes and resolution.

Run the python script **makegrid.py** to create the **gridfile.nc**

### ORCA NEMO Mediterranenan Grid

I have an input file, **namelist_R12**, for the CHyM Mediterranean
domain in the **MED12-ocean-mit** repository [here](https://github.com/graziano-giuliani/MED12-ocean-mit).

In the **MED12-ocean-mit/grid** directory, run:

    ./create_coordinates namelist_R12

and copy the file **1_coordinates_ORCA_R12.nc** in this directory as the file
**orca_coordinates.nc**.

Run the python script **makegrid.py** to create the **gridfile.nc** using the
*orca_grid* settings.

## Input datasets for the Mediterranean basin

### Source global datasets:

The HydroSheds, HydroRivers and HydroBasins data can be obtained from [Hydrosheds](https://www.hydrosheds.org)

The GLCC data can be otianed from the USGS [here](https://doi.org/10.5066/F7GB230D). The classification required is the Geo-Biosphere Observatories GEO20.

### Interpolate on regional domain

The global files are high resolution and BIG, so the first step is to cut a
smaller area and interpolate the land cover area onto the chosen grid.

The user can find a python script doing these steps. It needs the user to
have installed [cdo](https://code.mpimet.mpg.de/projects/cdo) and
[GDAL](https://gdal.org/en/stable).

Once the two input global dataset files
are copied or linked in this directory (*hyd_glo_dem_15s.tif* and
*gbogegeo20.tif*) and the *gridfile.nc* is present, the user can run the
interpolation script:

    python3 interpolate.py

The result should be the the two files:

1. *HydroSheds_15s.nc* : contains "windowed" data on the original 15s grid covering the area defined in the *gridfile.nc*
2. *landfile.nc* : GLCC GBO data interpolated on the target grid from the GLCC *30s* AVHRR daraset using the Largest Area Fraction cdo algorithm.

#### Optional:

The user may want to create a clipped down version of the global HydroEivers shape file (same source as Hydrosheds). For the Mediterranean (NOT USED!) this would ampunt to:

     ogr2ogr -clipsrc -7.0 27.0 48.5 63.5 filtered.shp HydroRIVERS_v10.shp 

## Create land/ocean mask

### ORCA Mediterranenan Grid

For the ORCA Med GRID, we can use the MITgcm mask file to create a common
land/ocean mask. After running the bathymetry step of the MITgcm processing,
the created mask has been extended and "corrected" to remove differences in
the two coastlines. The user may want to look in the scripts I have used to
replicate my setting OR adapt to her own needs.

1. Extend the MITgcm Mediterranean mask.

    python3 extend_mask.py

2. Fine tuning the mask. It is both to match the coastlines between CHyM and MITgcm, and to reduce computational cost by removing basins not feeding into the Mediterranean-Black Sea basins. It needs the user to download also Basin shapefiles from the HydroSheds. The resulting maskfile will contain the following mask:

* 0 : Points that should be masked out.
* 1 : Ocean points for BOTH MITgcm and CHyM.
* 2 : land points the CHyM should work on.

    python3 cleanmask.py

### Regular lat/lon grid

Run the script **createmask_ll.py**, which creates the mask based on the
landuse category in the file **landfile.nc**

    python3 createmask_ll.py

## Create River Network enforcing (optional)

It is based on reach python package, here my fork:

  * https://github.com/graziano-giuliani/reach-streamfilter

I have a script located in this directory **create_rivernet.py**

It must be run in the **reach-streamfilter** directory. Download the input
data:

     https://data.hydrosheds.org/file/HydroRIVERS/HydroRIVERS_v10_shp.zip

and uncompress it in the dat directory. Edit the python script to reflect
your domain, and run it. Copy all the data in the **out/filtered** data in
this directory. Create a raster data out of that:
   
     gdal_rasterize -a DIST_DN_KM -tr 0.025 0.025 -of NetCDF \
       -ot Float32 -a_nodata -1 filtered.shp rivernet_ll.nc
     ncrename -v Band1,network rivernet_ll.nc
     cdo remaplaf,gridfile.nc rivernet_ll.nc rivernet.nc

## Create the conditioned DEM for the model

Once all data are in the directory, the script **makedem.py** should be
executed to create the input **demfile.nc**

The script is using settings for the Mediterranean, but the user can modify
it to her liking.

The algorithm used for the coarsening of the input high resolution DEM file
is the following:

1. Read all input datasets:
  1. The *gridfile.nc* created above
  2. The *HydroSheds_15s.nc* created above
  3. The *maskfile.nc* created above
  4. The *HydroRIVERS_v10.shp* or the clipped down file created above.
2. Loop over the gridcell of the coarse destination grid and for all the selected ($mask == 2$) points falling in the polygon enclosed by the gridcell compute the "high" ($50$ percentile) and "low" ($10$ percentile) values.
3. Select one value or the other if the gridcell intersects a river basin as defined by the shapefile vectors simplified to a target resolution.
4. Write the final (hopefully still) conditioned coarse DEM file.

Once this is done, the input files for the CHyM preproc are ready:

  * gridfile.nc : File containing the grid geolocation informations
  * maskfile.nc : CHyM land sea mask
  * landfile.nc : Land use categories
  * demfile.nc : Conditioned Digital elevation model topography

Proceed!
