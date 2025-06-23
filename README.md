# CHyM runoff

The river routing model CHyM-roff is derived from hydrological model CHyM.

The code in this repository is a fork of the original work in ICTP by
Fabio Di Sante modified by Graziano Giuliani to be used in the evaluation
of the RegCM-ES1-1 ICTP coupled Regional Model over the Mediterranean
basin [in this repository](https://github.com/graziano-giuliani/MED12-ocean-mit)

The code here has its own pre-processing (extracted from the original CETEMPS
CHyM code [here](https://github.com/graziano-giuliani/CHyM), which I have
forked for simple bug fixing.

Major modification with respect to Fabio code is the introduction of:

    1. Scripts and examples on creating input data for the CHym-roff *preproc* program, usig the data available from [Hydrosheds](https://www.hydrosheds.org/products/hydrosheds) and the GBOGEO classification of the [GLCC AVHRR data](https://doi.org/10.5066/F7GB230D)
    2. The *preproc* program itself, which creates the static data for the drainage network used to compute the discarge.
    3. The capability to configure the programs using standard Fortran namelist
    4. The Manning coefficients are liberally derived from the values in [HEC HERAS Hydraulic Reference Manual](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/basic-data-requirements/geometric-data/energy-loss-coefficients)
    5. A simple seasonal scheme to account for irrigation withdrawal with monthly coefficients and a small configurable water loss.
    6. Revised I/O layer to be used in offline coupling with MITgcm
    7. General code revision and optimization: you can diff with Fabio [repository](https://github.com/fdisante/CHyM-roff) for the full list of modifications.

To compile the code, the user needs NetCDF and MPI libraries. The configure script should be able to locate the required bit. The procedure to create executables is:

    1. *autoreconf -f-i*
    2. *./configure*
    3. *make install*

The two binary files will be created in the *bin* directory. The *run* directory contains example namelist files.
