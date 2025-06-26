# CHyM runoff

The river routing model CHyM-roff is derived from hydrological model CHyM.

The code in this repository is a fork of the original work in ICTP by Fabio Di Sante modified by [Graziano Giuliani](https://www.ictp.it/member/graziano-giuliani) to be used in the evaluation of the RegCM-ES1-1 ICTP coupled Regional Model over the Mediterranean basin [in this repository](https://github.com/graziano-giuliani/MED12-ocean-mit)

The code in this repository has its own pre-processing (extracted from the original CETEMPS CHyM code [here](https://github.com/graziano-giuliani/CHyM)), which I have forked for simple bug fixing.

Major modification with respect to Fabio code is the introduction of:

1. Scripts and examples on creating input data for the CHym-roff *preproc* program, usig the data available from [Hydrosheds](https://www.hydrosheds.org/products/hydrosheds) and the GBOGEO classification of the [GLCC AVHRR data](https://doi.org/10.5066/F7GB230D)
2. The *preproc* program itself, which creates the static data for the drainage network used to compute the discarge.
3. The capability to configure the programs using standard Fortran namelist
4. The Manning coefficients are liberally derived from the values in [HEC HERAS Hydraulic Reference Manual](https://www.hec.usace.army.mil/confluence/rasdocs/ras1dtechref/6.6/basic-data-requirements/geometric-data/energy-loss-coefficients)
5. A simple seasonal scheme to account for irrigation withdrawal with monthly coefficients and a small configurable water loss.
6. Revised I/O layer to be used in offline coupling with MITgcm
7. General code revision and optimization: you can diff with Fabio [repository](https://github.com/fdisante/CHyM-roff) for the full list of modifications.

For the description of the Cellular Automata method creating the river network, the user can look at [Coppola et al. paper](https://www.tandfonline.com/doi/abs/10.1623/hysj.52.3.579).

The time integration in CHyM-roff requires input data for the *runoff* physical variable, which is a common product of any Land model component of most NWP or GCM models, and is generally defined as the portion of precipitation that doesn't infiltrate into the soil or evaporate and is "expected" to flow towards streams, rivers, and other bodies of water.

Because the CHyM model does not have a description of the ground but only takes care of the water transmission phase, the infiltration part of the original CHyM model is not present in the CHyM-roff model, which simulates only the momentum equation:

$$q = \frac{\sqrt{S} R^{\frac{2}{3}}}{n}$$

where $n$ is a function of the Manning coefiicient $M$:

$$n(M) = \frac{M}{\delta}$$

with $\delta$ configurable parameter, currently set to $\delta = 5.5$ ($cpar8$ in the preproc namelist).

$q$ is the flow of water discharge, $S$ is the slope, $R$ is the hydraulic radius, linear function of the drained area $D$ as in:

$$R = \alpha + \beta \max(D,D_{min})^{\gamma}$$

with $\alpha$, $\beta$, $\gamma$ and $D_{min}$ calibration coefficients. Values used are:

$\alpha = 0.0015$ ($cpar2$ in the preproc namelist), $\beta = 0.05$ ($cpar3$ in the preproc namelist), $\gamma = \frac{1}{3}$ ($cpar4$ in the preproc namelist) and $D_{min} = 100 km^2$ ($cpar6$ in the preproc namelist).

If the drain area $D$ is less than the threshold value $D_{min}$, the $n$ function is modified to be:

$$n(M) = \frac{M}{1+(\delta-1) \frac{1+(D-D_{min})}{D_{min}}}$$

Once the flow $q$ is computed, using the continuity equation for the total water $A$, stating that the water per unit element changes at a rate equal to the difference between outflow and inflow plus any source term:

$$\frac{\partial A}{\partial t} + \frac{\partial q}{\partial x} = q_c$$

can be used to compute the river discharge $Q$ by time integration with a configurable timestep as a fraction ($step$ in chymroff namelist) of the input runoff data time step ($dstep$ in the chymroff namelist). $q_c$ is the water per length unit sourced from the input runoff.

$$Q = q A(x,t)$$

The optional irrigation loss mod introduced acts by reducing the water per unit lenght available for flow by a factor changing on a monthly basis for the gridcells where a "crop type" category is present (classes $30,31,35,36,37,38,39,76,92,93,94,95,96$):

$$q_c = \left( 1 - \frac{r_i}{dt} \right) q_c$$

and allowing for a loss of the water reservoir:

$$\left.A right|_{i} = left( 1 - \frac{r_l r_i}{dt} \right) A$$

where $r_i$ is a factor changing with calendar month, and $r_l$ is a constant.
The values can be tuned in the chymroff namelist ($irloss$ and $irmonfac$). Setting all to $0$ disables the mod.

## Steps to run

The *preproc* program pre-computes the flow rate for every point of the domain, and the *chymroff* model is taking care of the time integration with configurable timestep, and can be run in parallel. The *chymroff* program needs the output of the *preproc* program as its input, together with the interpolated runoff dataset.

The sequence of action to obtain a successful simulation are the following.

### Download the source code from GitHub

    git clone https://github.com/graziano-giuliani/CHyM-roff.git

or fork the repository if you plan to modify the code, or download a *zip* snapshot of the code.

### Configure and compile the binary programs

To compile the code, the user needs NetCDF and MPI libraries. The configure script should be able to locate the required bit. The procedure to create executables is:

    1. autoreconf -f-i
    2. ./configure
    3. make install

The two binary files will be created in the *bin* directory. The *run* directory contains example namelist files.

### Configure the grid and create input data

Refer to the README file in the *grid* directory for this steps.

In principle, given the global coverage of the input data, the model can run on any region of the world. The user needs to interpolate global data and create the following netCDF format input files for the *preproc* program:

1. *gridfile.nc* Geolocation and geometric information
2. *landfile.nc* Landuse file with per gridcell information
3. *maskfile.nc* Land/ocean mask, it is important if coupling with an ocean model to have matching nad masks.
4. *demfile.nc* Digital Elevation Model data
5. Optional : *riverfile.nc* Rasterized river mask

### Run the *preproc* program

Change into the run directory and copy over or link the above described input file from the *grid* directory and the model binary files from the *bin* directory. The file names are not important, the use can set their values in the namelist files.

The model parameters, above described, can be changed in the input namelist file. The name of the file is not relevant, it is important the value specified on the command line argument when running the program:

    ./preproc preproc.namelist

If the program works, the file specified in the namelist as *outfile* should have been created. It is recommended the user checks the variable inside the file with a visual inspection to verify the consistency of the river drainage network. More sofisticated GIS analysis tools can be used to compare against a vectorial representation of the expected "real" river network, but usually a visual comparison is enough to identify macroscopic differences.

The parameters used by the CA algorithms can be modified in according to what is presented in the Coppola paper above linked.

### Download and interpolate runoff dataset

The user is responsible to prepare the runoff data for the CHyM-roff. Refer to the README file in the run/input directory for example guideline on how to do this for the ERA5 ECMWF reanalysis.

In the *run/input/global* the user can find an example python script to download the ECMWF ERA5 daily runoff data. The model expect input data to be on the same grid as defined on the *gridfile.nc*.

Common is to use the program [cdo](https://code.mpimet.mpg.de/projects/cdo) for its flexibility and ease of use in performing the required interpolation to create the input file.

An example is:

      cdo remapcon,gridfile.nc global/runoff_1996.nc runoff_1996.nc

The data are expected to be divided in yearly files, containg the runoff data.

### Running the *chymroff* program

The user needs to edit the namelist input file to setup an eventual conversion factor from the input data units to the expected runoff units (mm/s or kg m-2 s-1). The user can find the example in the *run* directory for the ERA5 model.

The yearly input file names and the data frequency herein is configurable in the namelist file, together with the name of the static file produced by the *preproc* program.

The user can restart the run from a previous checkpoint file (*rst*) by changing to the value $1$ the value of the $isread$ entry.
The restart file production and frequency is also configurable, and the output file basic naming prefixed part. The user is responsible to create the output file path if it does not already exist.

The user must finally define the time integration window and calendar. Once the input is ready and the model is configured, the user can run the model in parallel in the *run* directory:

    mpirun ./chymroff chymroff.namelist

Similar consideration as above for the name of the namelist file. The output files will be created in the requested output path with the requested naming convention. A simple dataset which can be used to verify results is the RivDIS v1.1 available [here](https://doi.org/10.3334/ORNLDAAC/199). Example plotting script is in the run directory as *plot.py*

# Questions?

Please drop me an email if you need help: [ggiulian\@ictp.it](mailto:ggiulian\@ictp.it)

