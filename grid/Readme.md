# Procedure

This file describes the procedure to create input files for the CHyM preproc.

## Create the Gridfile

Create grid file using the python script, either onto regular
latitude-longitude grid or an ORCA NEMO tripolar grid

### ORCA NEMO Mediterranenan Grid

I have an input file, **namelist_R12_chym**, for the CHyM Mediterranean
domain in the **MED12-ocean-mit** repository here:

    https://github.com/graziano-giuliani/MED12-ocean-mit

In the **MED12-ocean-mit/grid** directory, run:

    ./create_coordinates namelist_R12_chym

and copy the file **1_coordinates_ORCA_R12.nc** in this directory as the file
**orca_coordinates.nc**.

### Regular lat-lon grid

Edit the **makegrid.yaml** input file for extremes and resolution.

### Create the grid

Run the python script **makegrid.py** to create the **gridfile.nc**

## Input datasets for the Mediterranean basin

### Source global datasets:

The HydroSheds, HydroRivers and HydroBasins data can be obtained from:

    https://www.hydrosheds.org/

The GLCC data can be otianed from the USGS:

    https://doi.org/10.5066/F7GB230D

### Interpolate on regional domain

From Global files, this are the commands I have used for the Mediterranean
grid [GDAL 3.5.2]:

    gdal_translate -ot Float32 -of netCDF \
       -projwin -7.0 63.5 48.5 27.0 \
       hyd_glo_dem_15s.tif HydroSheds_15s_Mediterraneo.nc

     gdal_translate -ot Int16 -of netCDF \
       -projwin -7.0 63.5 48.5 27.0 \
       gbogegeo20.tif GLCC_gboggeo20_Mediterraneo.nc

Some renaming and interpolation [nco 5.1.3, cdo 2.2.0]:

     ncrename -h -v Band1,luc GLCC_gboggeo20_Mediterraneo.nc
     ncatted -h -a standard_name,luc,c,c,class_type \
                -a long_name,luc,m,c,"Class Type" \
                -a units,luc,c,c,1 GLCC_gboggeo20_Mediterraneo.nc
     cdo remaplaf,gridfile.nc GLCC_gboggeo20_Mediterraneo.nc landfile.nc

     ncrename -h -v Band1,dem HydroSheds_15s_Mediterraneo.nc
     ncatted -h -a standard_name,dem,c,c,elevation \
                -a long_name,dem,m,c,"Digital Elevation Model" \
                -a units,dem,c,c,m HydroSheds_15s_Mediterraneo.nc

Create clipped river network shape file:

     ogr2ogr -clipsrc -7.0 27.0 48.5 63.5 filtered.shp HydroRIVERS_v10.shp 

## Create land/ocean mask

### ORCA Mediterranenan Grid

For the ORCA Med GRID, we can use the MITgcm mask file to create a common
land/ocean mask. After running the bathymetry step of the MITgcm processing,
remap the ocean mask onto the CHyM grid:

    cdo remapnn,gridfile.nc mitgcm_mask.nc maskfile.nc

### Regular lat/lon grid

Run the script **makemaskll.py**, which creates the mask based on the
landuse category in the file **landfile.nc**

### Cleanup the mask and select only relevant basins

To reduce computational cost, we can create a digital mask to select only
the relevant basins. The script in this directory, **cleanmask.py**, is tailored over the Mediterranean, but the user can look in hydrobasins file for the code
of the river basins in her area of interest. This script sets:

 * 0 : land points the CHyM should work on
 * 1 : Ocean points
 * 2 : Points that should be masked out.

To run it:

     python3 cleanmask.py maskfile.nc

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

Once all data are in the directory, the script **makedem.py** shoudl be
executed to create the input **demfile.nc**

Once this is done, the input files for the CHyM preproc are ready:

  * gridfile.nc : File containing the grid geolocation informations
  * demfile.nc : Conditioned Digital elevation model topography
  * landfile.nc : Land use categories
  * maskfile.nc : CHyM land sea mask
