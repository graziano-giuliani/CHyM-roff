# Procedure

This file describes the procedure to create input files for the CHyM preproc.

## Gridfile

Create grid file using the python script, either onto regular
latitude-longitude or an ORCA NEMO tripolar grid (needs coordinates file).

Run the python script **makegrid.py** to create the **gridfile.nc**

### ORCA Mediterranenan Grid

I have an input file, **namelist_R12_chym**, for the CHyM Mediterranean domain.
In the **MED12-ocean-mit/grid** directory, run:

    ./create_coordinates namelist_R12_chym

and copy the file **1_coordinates_ORCA_R12.nc** in this directory as the file
**orca_coordinates.nc**.

### Regulat latitude-longitude grid

Edit the **makegrid.yaml** input file for extremes and resolution.

## Create input datasets for the Mediterranean basin

### Source global datasets:

The Hydorsheds and HydroBasins data can be obtained from:

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

Some renaming [nco 5.1.3]:

     ncrename -h -v Band1,luc GLCC_gboggeo20_Mediterraneo.nc
     ncatted -h -a standard_name,luc,c,c,class_type \
                -a long_name,luc,m,c,"Class Type" \
                -a units,luc,c,c,1 GLCC_gboggeo20_Mediterraneo.nc

     ncrename -h -v Band1,dem HydroSheds_15s_Mediterraneo.nc
     ncatted -h -a standard_name,dem,c,c,elevation \
                -a long_name,dem,m,c,"Digital Elevation Model" \
                -a units,dem,c,c,m HydroSheds_15s_Mediterraneo.nc

Interpolate data on CHyM grid [cdo 2.2.0]:

     cdo remapdis,gridfile.nc,25 HydroSheds_15s_Mediterraneo.nc demfile.nc
     cdo remaplaf,gridfile.nc GLCC_gboggeo20_Mediterraneo.nc landfile.nc

## Create land/ocean mask

### ORCA Mediterranenan Grid

For the ORCA Med GRID, we can use the MITgcm mask file to create a common
land/ocean mask. The procedure will mark:

 * 0 : land points the CHyM should work on
 * 1 : Ocean points
 * 2 : Points that should be masked out.

     cdo remapnn,gridfile.nc mitgcm_mask.nc maskfile.nc
     python3 cleanmask.py maskfile.nc

### Regular Latitude-Longitude grid

The maskfile is created from the **demfile.nc**

     cdo 

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

