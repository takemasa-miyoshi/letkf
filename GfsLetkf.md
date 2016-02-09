

# Getting started #

## Download our LETKF package ##

Download the latest release version from the [Source tab](http://code.google.com/p/miyoshi/source/checkout).

## Prepare required executables and libraries ##

Our LETKF package provides the LETKF assimilation code coupled to the National Centers for Environmental Prediction (NCEP) Global Forecasting System (GFS) model. It requires several executables and libraries from the NCEP, which mainly belong to two NCEP models: GFS and Gridpoint Statistical Interpolation (GSI). Compiling these libraries and models on your computer may not be easy and may take long time.

  * **GFS**: The GFS model is an operational global NWP model developed by the Environmental Modeling Center (EMC) at the NCEP. It is one of the major state-of-the-art operational NWP models over the world and provides main model guidance for the weather forecast in the United States. It is not a community model, and some of the code may not be open to public, so users interested to run the GFS-LETKF system need to contact NCEP to obtain the GFS model code and the associated libraries.
    * NCEP/EMC GFS model page: http://www.emc.ncep.noaa.gov/index.php?branch=GFS
    * When you compile the GFS model, you may need to compile the [Earth System Modeling Framework (ESMF)](http://www.earthsystemmodeling.org/) first.

  * **GSI**: The GSI is the primary data assimilation system for the GFS model, based on 3-dimensional variational method (3DVar). In our GFS-LETKF implementation, there is an option to use the GSI as an observation operator. Note that we only use the GSI as the observation operator, and the GSI 3DVar solution is not computed in the GFS-LETKF. Besides, if you do not plan to assimilate the satellite radiance data, in the GFS-LETKF package there is also a set of built-in observation operators that does not rely on the GSI. In this case, you do not need to install the GSI.
    * The GSI has been a community model, available at http://www.dtcenter.org/com-GSI/users/

### List of required executables and libraries ###

  * **NCEP executables**:
> After successfully compiling the GFS and GSI, you are asked to put all executable files together in the '`$EXECGLOBAL`' directory configured in '`gfs/run/configure.sh`'.
| **Name** | **Purpose** | **Called from** |
|:---------|:------------|:----------------|
| `global_fcst` | GFS main program | `gfs/run/run_gfs.sh`<br><code>gfs/run/cycle.sh</code><br><code>gfs/run/fcst.sh</code> <br>
<tr><td> <code>global_gsi</code> </td><td> GSI main program </td><td> <code>gfs/run/run_gsi.sh</code><br><code>gfs/run/cycle.sh</code> </td></tr>
<tr><td> <code>global_sighdr</code> </td><td> Read the header of the GFS sigma-level input/output files </td><td> <code>gfs/run/run_gfs.sh</code><br><code>gfs/run/run_gsi.sh</code> </td></tr>
<tr><td> <code>global_sfchdr</code> </td><td> Read the header of the GFS surface input/output files </td><td> <code>gfs/run/run_gfs.sh</code><br><code>gfs/run/run_gsi.sh</code> </td></tr>
<tr><td> <code>global_chgres</code> </td><td> Change the resolution of the GFS input/output files </td><td> <code>gfs/run/run_chgres.sh</code> </td></tr></li></ul></tbody></table>

<ul><li><b>NCEP libraries required for compiling the GFS-LETKF code</b>:<br>
</li></ul><blockquote>These libraries are required for compiling several programs in the GFS-LETKF package. The path of these libraries are configured in '<code>gfs/configure.user</code>' :<br>
<table><thead><th> <b>Name</b> </th><th> <b>Purpose</b> </th><th> <b>List of files required</b> </th><th> <b>Link</b> </th></thead><tbody>
<tr><td> <code>bacio</code> </td><td> Basic I/O library for the GFS </td><td> <code>$(BACIO_LIB)/libbacio_4.a</code> </td><td>             </td></tr>
<tr><td> <code>sigio</code> </td><td> I/O library for the GFS sigma-level files </td><td> <code>$(SIGIO_LIB)/libsigio_4.a</code><br><code>$(SIGIO_INC)/sigio_module.mod</code><br><code>$(SIGIO_INC)/sigio_r_module.mod</code> </td><td>             </td></tr>
<tr><td> <code>sfcio</code> </td><td> I/O library for the GFS surface-level files </td><td> <code>$(SFCIO_LIB)/libsfcio_4.a</code><br><code>$(SFCIO_INC)/sfcio_module.mod</code> </td><td>             </td></tr>
<tr><td> <code>sp</code> </td><td> Spectral transform library </td><td> <code>$(SP_LIB)/libsp_4.a</code> </td><td> <a href='http://www.nco.ncep.noaa.gov/pmb/docs/libs/splib/ncep_splib.shtml'>http://www.nco.ncep.noaa.gov/pmb/docs/libs/splib/ncep_splib.shtml</a> </td></tr>
<tr><td> <code>bufrlib</code> </td><td> I/O library for the BUFR file format </td><td> <code>$(BUFR_LIB)/libbufrlib.a</code> </td><td> <a href='http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/'>http://www.nco.ncep.noaa.gov/sib/decoders/BUFRLIB/</a> </td></tr></blockquote></tbody></table>

<ul><li><b>Other NCEP libraries that may be required for compiling the GFS and GSI</b>:<br>
</li></ul><blockquote>These libraries may be required when you compile the GFS and GSI, but they are NOT directly required in the GFS-LETKF.<br>
<table><thead><th> <b>Name</b> </th><th> <b>Purpose</b> </th></thead><tbody>
<tr><td> <code>ip</code> </td><td> General interpolation library </td></tr>
<tr><td> <code>landsfcutil</code> </td><td> Land surface library </td></tr>
<tr><td> <code>gfsio</code> </td><td> GFS I/O library </td></tr>
<tr><td> <code>nemsio</code> </td><td> NEMS I/O library </td></tr>
<tr><td> <code>w3lib</code> </td><td> I/O library for the GRIB file format </td></tr></blockquote></tbody></table>

<ul><li><b>Other required libraries</b>:<br>
</li></ul><blockquote>LAPACK may be required when it is not provided by the Fortran compiler. For example, if you use the Intel Fortran compiler with its Math Kernel Library (MKL), the LAPACK is included and you do not need to install LAPACK by yourself.<br>
<table><thead><th> <b>Name</b> </th><th> <b>Purpose</b> </th><th> <b>Link</b> </th></thead><tbody>
<tr><td> <code>LAPACK</code> </td><td> Linear algebra package </td><td> <a href='http://www.netlib.org/lapack/'>http://www.netlib.org/lapack/</a> </td></tr></blockquote></tbody></table>

<h2>Prepare fix files for the GFS and GSI</h2>

To run the GFS and GSI, they require several "fix" data files (e.g., boundary conditions) and CRTM coefficients. These files need to be put in three directories: '<code>$FIXGLOBAL</code>', '<code>$FIXGSI</code>', and '<code>$FIXCRTM</code>', respectively, configured in '<code>gfs/run/configure.sh</code>'.<br>
<br>
The GFS fix files can be found at <a href='http://www.nco.ncep.noaa.gov/pmb/codes/nwprod/fix/'>http://www.nco.ncep.noaa.gov/pmb/codes/nwprod/fix/</a> , but they are much more than enough to drive a GFS model. It is recommended to select only the necessary files when putting into '<code>$FIXGLOBAL</code>'.<br>
<br>
<h2>Test compiling the GFS-LETKF programs</h2>

When the required libraries are ready, we can do a test compilation of all GFS-LETKF programs.<br>
<br>
<ul><li>Go to the GFS main directory.<br>
<pre><code>$ cd gfs<br>
</code></pre>
</li><li>Edit the '<code>configure.user</code>' file, specify the compiler flags and paths of the required libraries in your system. In the '<code>arch</code>' directory, there are a few examples of the '<code>configure.user</code>' file using the PGI or Intel compilers.<br>
</li><li>Compile all GFS-LETKF programs.<br>
<pre><code>$ make<br>
</code></pre></li></ul>

<h3>List of the GFS-LETKF programs</h3>

If the compilation is successful, you will obtain a number of executable files listed below:<br>
<table><thead><th> <b>File path and name</b> </th><th> <b>Purpose</b> </th></thead><tbody>
<tr><td> <code>common/datetime</code> </td><td> A utility for date and time computation </td></tr>
<tr><td> <code>common/enssize</code> </td><td> To simply print the ensemble size </td></tr>
<tr><td> <code>letkf/efsoXXX</code> </td><td> EFSO mean program </td></tr>
<tr><td> <code>letkf/letkfXXX</code> </td><td> LETKF main program </td></tr>
<tr><td> <code>letkf/meanXXX</code> </td><td> To compute the ensemble mean from grid files </td></tr>
<tr><td> <code>letkf/obsope</code> </td><td> Built-in observation operator for non-radiance data </td></tr>
<tr><td> <code>obs/dec_prcp</code> </td><td> To convert the gridded precipitation data to the LETKF observation format </td></tr>
<tr><td> <code>obs/dec_prepbufr</code> </td><td> To convert the NCEP PREPBUFR data to the LETKF observation format </td></tr>
<tr><td> <code>readdiag_conv</code> </td><td> To convert the GSI diagnostic files to the LETKF observation format with model background values </td></tr>
<tr><td> <code>obs/superob</code>  </td><td> A superobing/thinning utility </td></tr>
<tr><td> <code>ssio/grd2ss</code>  </td><td> To convert a grid file to GFS sig/sfc (spectral) files, and also cycle forecasts into analyses </td></tr>
<tr><td> <code>ssio/grdctl</code>  </td><td> To print GrADS CTL files according to the current grid settings </td></tr>
<tr><td> <code>ssio/ss2grd</code>  </td><td> To convert GFS sig/sfc (spectral) files to a grid file </td></tr>
<tr><td> <code>ssio/ss2grdp</code> </td><td> To convert GFS sig/sfc (spectral) files to a grid file in pressure coordinates </td></tr>
<tr><td> <code>ssio/sscycle</code> </td><td> To cycle forecasts into analyses </td></tr>
<tr><td> <code>util/gfsmeanXXX</code> </td><td> To compute the ensemble mean from GFS sig/sfc (spectral) files </td></tr>
<tr><td> <code>verify/verify</code> </td><td> A utility to perform verification against observations or other model analyses </td></tr></tbody></table>

Note that <code>XXX</code> is the ensemble size that the program is used for.<br>
<br>
If some executable files are missing, check the errors and try to fix them.<br>
<br>
<h1>Code overview</h1>

<h2>List of scripts</h2>
All run scripts are located in the '<code>gfs/run</code>' directory, so you will do every task in this directory. The purposes of the scripts are described as follows. Note that for most of the scripts, you can show the usage help by directly executing them with no arguments.<br>
<br>
<table><thead><th> <b>Script</b> </th><th> <b>Purpose</b> </th></thead><tbody>
<tr><td> <code>configure.sh</code> </td><td> Main configurations for all GFS-LETKF scripts. </td></tr>
<tr><td> <code>datetime.sh</code> </td><td> Date and time functions. It uses the Unix 'date' command.<br>If the 'date' command does not work in your system, replace '<code>datetime.sh</code>' by '<code>datetime.sh.no_unix_date</code>' that does not require the 'date' command. </td></tr>
<tr><td> <code>distribution.sh</code> </td><td> Functions to adaptively distribute members on nodes. </td></tr>
<tr><td> <code>stageinout.sh</code> </td><td> Functions to copy files between the server node and local disks on computing nodes. </td></tr>
<tr><td> <code>get_ncepobs.sh</code> </td><td> Download NCEP conventional observation data from the CISL Research Data Archive. </td></tr>
<tr><td> <code>get_cfsr.sh</code> </td><td> Download NCEP CFSR data. </td></tr>
<tr><td> <code>run_chgres.sh</code> </td><td> Change the resolution of a series of GFS initial conditions. </td></tr>
<tr><td> <code>mss2grd.sh</code> </td><td> Convert GFS sig/sfc files to GrADS grd/grdp files. </td></tr>
<tr><td> <code>outdir.sh</code> </td><td> Create necessary subdirectories and files in the output (<code>$OUTDIR</code>) directory. </td></tr>
<tr><td> <code>init.sh</code> </td><td> Prepare a initial ensemble from a series of analyses at different times. </td></tr>
<tr><td> <code>init2.sh</code> </td><td> Prepare a initial ensemble from outputs of another experiments. </td></tr>
<tr><td> <code>init3.sh</code> </td><td> Prepare a series of initial mean analyses from another source, which is useful to run forecast experiments. </td></tr>
<tr><td> <code>run_gfs.sh</code> </td><td> Prepare a temporary directory for a GFS run, and may run the model. </td></tr>
<tr><td> <code>run_gsi.sh</code> </td><td> Prepare a temporary directory for a GSI run, and may run the model. </td></tr>
<tr><td> <code>cycle.sh</code> </td><td> Run a GFS-LETKF forecast/analysis cycle. </td></tr>
<tr><td> <code>mcycle.sh</code> </td><td> Run multiple cycles. </td></tr>
<tr><td> <code>fcst.sh</code> </td><td> Run (ensemble) forecasts and may also perform forecast verification. </td></tr>
<tr><td> <code>mfcst.sh</code> </td><td> Run multiple-cycle ensemble mean forecasts. </td></tr>
<tr><td> <code>efsofcst.sh</code> </td><td> Run multiple-cycle ensemble forecasts for the EFSO computation. </td></tr>
<tr><td> <code>verify.sh</code> </td><td> Compute verification of (ensemble) forecasts. </td></tr>
<tr><td> <code>mverify.sh</code> </td><td> Run multiple-cycle verification. </td></tr>
<tr><td> <code>efso.sh</code> </td><td> Compute the ensemble forecast sensitivity to observations (EFSO). </td></tr>
<tr><td> <code>mefso.sh</code> </td><td> Run multiple-cycle EFSO. </td></tr>
<tr><td> <code>pbs_example.sh</code> </td><td> An example script to submit a parallel job with the PBS jobs scheduling system. </td></tr></tbody></table>

<h2>Main configuration script</h2>

The '<code>configure.sh</code>' is the main configuration file. The other scripts are run based on this file. The explanation of all the configurable variables are provided along with the file (using code comments). Selected variables are explained here:<br>
<br>
<table><thead><th> <b>Variable</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> <code>OUTDIR</code> </td><td> GFS-LETKF input and output directory. </td></tr>
<tr><td> <code>SYSNAME</code> </td><td> A unique name in the machine, which is used to identify multiple GFS-LETKF runs. </td></tr>
<tr><td> <code>LTMP1</code><br><code>LTMP2</code> </td><td> Level 1 and 2 local temporary directories on computing nodes. They are only used when <code>$SHAREDISK = 0</code>. They must exist on all computing nodes. The capacity requirement of the level 1 temporary directory is not too big but the I/O speed is very important; therefore, a RAMdisk (<code>/dev/shm</code> in Linux machines) could be assigned to this variable if the memory on your machine is sufficient. The capacity requirement of the level 2 temporary directory is bigger than <code>$LTMP1</code> and the I/O speed is less important than <code>$LTMP1</code>. The <code>$LTMP1</code> and the <code>$LTMP2</code> can be assigned to the same directory if you do not want to separate the location of different temporary files. </td></tr>
<tr><td> <code>TMP1</code> </td><td> Temporary directory on the server machine. Accessing to this directory from other computing nodes is not required. </td></tr>
<tr><td> <code>TMPMPI</code> </td><td> Temporary directory on the server machine with shared access from all computing nodes. </td></tr>
<tr><td> <code>SHAREDISK</code> </td><td> Do we use local disks on computing nodes to store runtime temporary file?<br>0: Yes, use local disks (<code>$LTMP1</code>, <code>$LTMP2</code>) to store temporary files.<br>1: No, the <code>$OUTDIR</code> is a shared disk that can be accessed from all computing nodes. Use this shared disk to store temporary files.</td></tr>
<tr><td> <code>EXECGLOBAL</code> </td><td> Directory of the <a href='#List_of_required_executables_and_libraries.md'>NCEP executable files</a> (GFS, GSI... etc.). </td></tr>
<tr><td> <code>FIXGLOBAL</code> </td><td> Directory of GFS fix files. </td></tr>
<tr><td> <code>FIXGSI</code> </td><td> Directory of GSI fix files, only required when using GSI as the observation operator (<code>$OBSOPE_OPT = 2</code>). </td></tr>
<tr><td> <code>FIXCRTM</code> </td><td> Directory of CRTM fix files, only required when using GSI as the observation operator (<code>$OBSOPE_OPT = 2</code>). </td></tr>
<tr><td> <code>ANLGFS</code> </td><td> Directory of reference model files in the GFS sig/sfc formats. </td></tr>
<tr><td> <code>ANLGRD</code> </td><td> Directory of reference model files in the sigma-level grid format (grd). </td></tr>
<tr><td> <code>ANLGRDP</code><br><code>ANLGRDP2</code> </td><td> Directory of reference model files in the pressure-level grid format (grdp). They are used for verification. The verification results against the model data in <code>$ANLGRDP</code> will be stored in '<code>$OUTDIR/verfa1</code>'; the verification results against the model data in <code>$ANLGRDP2</code> will be stored in '<code>$OUTDIR/verfa2</code>'. </td></tr>
<tr><td> <code>INITGFS</code> </td><td> Directory of arbitrary initial condition files in the GFS sig/sfc formats. It can be the same as <code>$ANLGFS</code> if the quality of the reference model data in that directory is good. See explanation in the <a href='#Prepare_the_GFS_initial_conditions.md'>latter section</a>. </td></tr>
<tr><td> <code>OBS</code> </td><td> Directory of observation data in the LETKF observation format, required when using the built-in observation operator (<code>$OBSOPE_OPT = 1</code>). </td></tr>
<tr><td> <code>OBSNCEP</code> </td><td> Directory of observation data in the NCEP BUFR format, required when using GSI as the observation operator (<code>$OBSOPE_OPT = 2</code>). </td></tr>
<tr><td> <code>MEMBER</code> </td><td> Ensemble size. It should be the same as the '<code>nbv</code>' variable in '<code>common/common_letkf.f90</code>'. </td></tr>
<tr><td> <code>MIN_NP_GFS</code><br><code>MIN_NP_GSI</code> </td><td> Minimum numbers of CPU cores required to run the GFS and GSI. This limits are to avoid using up all available memory per node/core. </td></tr>
<tr><td> <code>MAX_NP_GFS</code><br><code>MAX_NP_GSI</code> </td><td> Maximum numbers of cores suggested to run the GFS and GSI. This limits are to avoid poor parallelization efficiency when using too many nodes/cores. </td></tr>
<tr><td> <code>OBSOPE_OPT</code> </td><td> Observation operator options:<br>1: use the LETKF built-in observation operators.<br>2: use the GSI as the observation operator. </td></tr>
<tr><td> <code>THIN_OPT</code> </td><td> Superobing/thinning options:<br>-- Options below (1-2) is for <code>$OBSOPE_OPT = 1</code><br>1: No superobing/thinning.<br>2: Use superobed/thinned observations processed by the '<code>gfs/obs/superob</code>' program before the LETKF assimilation.<br>-- Options below (3-4) are for <code>$OBSOPE_OPT = 2</code><br>3: use thinned observations for satellite radiance observations only, processed by the GSI.<br>4: use thinned observations for both conventional and satellite radiance observations, processed by the GSI. </td></tr>
<tr><td> <code>ADAPTINFL</code> </td><td> Adaptive inflation options:<br>0: OFF, no adaptive inflation.<br>1: ON, using inflation parameter 1 cycle ago as the prior.<br>2: ON, using inflation parameter 2 cycles ago as the prior (leap-frog the adaptive inflation fields). </td></tr>
<tr><td> <code>FCSTLEN</code> </td><td> GFS forecast length when running the (ensemble) forecasts ('<code>gfs/run/fcst.sh</code>'). </td></tr>
<tr><td> <code>OUT_OPT</code><br><code>FOUT_OPT</code><br><code>OBSOUT_OPT</code><br><code>LOG_OPT</code> </td><td> How detail do you want to keep for the regular and diagnostic output files in '<code>$OUTDIR</code>'? See details in the '<code>configure.sh</code>' script. </td></tr>
<tr><td> <code>MPIBIN</code> </td><td> The path of the '<code>mpiexec</code>' command to run a parallel program. </td></tr>
<tr><td> <code>BUFRBIN</code> </td><td> The path of the '<code>grabbufr</code>' command. </td></tr></tbody></table>

<h2>Variables in the source code</h2>

There are other variables that can not be configured in the '<code>configure.sh</code>', but need to be set in the Fortran source code before compiling the GFS-LETKF programs. After changing the values of these variables, the GFS-LETKF programs need to be re-compiled. Some important variables in the source code are listed and explained below:<br>
<br>
<table><thead><th> <b>Source file</b> </th><th> <b>Variable</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> <code>common/common_letkf.f90</code> </td><td> <code>nbv</code> </td><td> Ensemble size      </td></tr>
<tr><td> <code>gfs/common/common_gfs.f90</code> </td><td> <code>nlon</code><br><code>nlat</code> </td><td> Longitude and latitude grid numbers of the GFS model. They correspond to the spectral truncation wavenumbers:<br>For T62, <code>$nlon</code> = 192, <code>$nlat</code> = 94<br>For T126, <code>$nlon</code> = 384, <code>$nlat</code> = 190 </td></tr>
<tr><td> <code>gfs/common/common_gfs.f90</code> </td><td> <code>nlev</code> </td><td> Number of the vertical levels of the GFS model </td></tr>
<tr><td> <code>gfs/common/common_gfs.f90</code> </td><td> <code>gfs_jcap</code> </td><td> Spectral truncation wavenumber of the GFS model </td></tr>
<tr><td> <code>gfs/common/common_gfs_pres.f90</code> </td><td> <code>nlevp</code> </td><td> Number of the vertical levels in the GFS-LETKF pressure-level outputs </td></tr>
<tr><td> <code>gfs/common/common_gfs_pres.f90</code> </td><td> <code>levp</code> </td><td> List of the vertical levels in the GFS-LETKF pressure-level outputs </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>omb_output</code> </td><td> Whether output the (observation - background) statistics? </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>oma_output</code> </td><td> Whether output the (observation - analysis) statistics? </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>obsgues_output</code> </td><td> Whether output the observation values in the ensemble model background [H(X<sup>b</sup>)]? </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>obsanal_output</code> </td><td> Whether output the observation values in the ensemble model analyses [H(X<sup>a</sup>)]? This is used for the EFSO computation, and it requires considerable additional computational time of the LETKF main program. </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>sigma_obs</code> </td><td> Horizontal localization length scale </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>sigma_obsv</code> </td><td> Vertical localization length scale </td></tr>
<tr><td> <code>gfs/letkf/letkf_obs.f90</code> </td><td> <code>gross_error</code> </td><td> Quality control with gross errors </td></tr>
<tr><td> <code>gfs/letkf/letkf_tools.f90</code> </td><td> <code>cov_infl_mul</code> </td><td> Multiplicative  covariance inflation parameter:<br>> 0: globally constant inflation<br>< 0: 3D inflation values input from the '<code>infl_mul.grd</code>' file </td></tr>
<tr><td> <code>gfs/letkf/letkf_tools.f90</code> </td><td> <code>sp_infl_mul</code> </td><td> Additive covariance inflation parameter </td></tr>
<tr><td> <code>gfs/letkf/letkf_tools.f90</code> </td><td> <code>var_local</code> </td><td> Variable localization matrix </td></tr>
<tr><td> <code>gfs/obs/superob.f90</code> </td><td> <code>obmethod_g</code><br><code>obmethod_v</code><br><code>obmethod_t</code><br><code>obmethod_h</code> </td><td> Superobing and thinning settings. See details in the code comments. </td></tr>
<tr><td> <code>gfs/verify/verify.f90</code> </td><td> <code>nvrf_obs</code><br><code>nvrf_ana</code> </td><td> Numbers of the observation datasets and the model analysis datasets used for the verification </td></tr>
<tr><td> <code>gfs/verify/verify.f90</code> </td><td> <code>narea</code><br><code>vlon1</code><br><code>vlon2</code><br><code>vlat1</code><br><code>vlat2</code> </td><td> Settings of the verification regions </td></tr></tbody></table>

<h2>Data formats</h2>

Below are several data formats appeared in the GFS-LETKF system:<br>
<br>
<ul><li><b>GFS sigma/surface files</b>
</li></ul><blockquote>GFS model inputs and outputs.<br>Abbreviation: <b>sig/sfc</b></blockquote>

<ul><li><b>GrADS grid files in sigma coordinate (same as the GFS model levels)</b>
</li></ul><blockquote>Gridded files in model levels that can be read by the LETKF main program and plotted with the GrADS software.<br>Abbreviation: <b>grd</b></blockquote>

<ul><li><b>GrADS grid files in pressure coordinate</b>
</li></ul><blockquote>Gridded files in pressure levels, which is defined by '<code>nlevp</code>' and '<code>levp</code>' in '<code>gfs/common/common_gfs_pres.f90</code>'. It can be plotted with the GrADS software.<br>Abbreviation: <b>grdp</b></blockquote>

<ul><li><b>NCEP PREPBUFR files</b>
</li></ul><blockquote>NCEP PREPBUFR observation date format.<br>Abbreviation: <b>prepbufr</b></blockquote>

<ul><li><b>LETKF observation format</b>
</li></ul><blockquote>A special observation data format used by the LETKF code.<br>Abbreviation: <b>letkfobs</b></blockquote>

<blockquote>The observation data are saved in a simple binary structure. For a single observation, it consists of these columns, all of which are 4-bytes floating variables:<br>
<table><thead><th> <b>Column</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> 1             </td><td> Variable type; see '<code>id_*_obs</code>' variables in '<code>gfs/common/common_obs_gfs.f90</code>' </td></tr>
<tr><td> 2             </td><td> Longitude (degree) </td></tr>
<tr><td> 3             </td><td> Latitude (degree)  </td></tr>
<tr><td> 4             </td><td> u,v,t,tv,q,rh: level (hPa)<br>ps: station elevation (m) </td></tr>
<tr><td> 5             </td><td> Observation value:<br>wind (m/s)<br>temperature (K)<br>specific humidity (kg/kg)<br>relative humidity (%)<br>surface pressure (hPa) </td></tr>
<tr><td> 6             </td><td> Observation error (unit same as the observation value) </td></tr>
<tr><td> 7             </td><td> Observation platform type; see the '<code>obtypelist</code>' array in '<code>gfs/common/common_obs_gfs.f90</code>' </td></tr></blockquote></tbody></table>

<ul><li><b>LETKF observation format with additional columns</b>
</li></ul><blockquote>It is an extension of the "letkfobs" format with more columns to save the observation value in the model.<br>Abbreviation: <b>letkfobs2</b></blockquote>

<blockquote>However, the meaning of the additional columns varies in different cases. In the observation data processed by '<code>obsope</code>' or '<code>readdiag_*</code>' programs, or the input data of the LETKF main program:<br>
<table><thead><th> <b>Column</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> 8             </td><td> Observation time relative to analysis time (hour; can be positive or negative) </td></tr>
<tr><td> 9             </td><td> h(x) observation in model background (unit same as observation value except surface pressure in Pa) </td></tr>
<tr><td> 10            </td><td> quality control mark (1 = pass; others = do not pass) </td></tr></blockquote></tbody></table>

<blockquote>In the observation diagnostics output by the LETKF program (<code>$OUTDIR/obs/obsdiag</code>):<br>
<table><thead><th> <b>Column</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> 8             </td><td> Observation time relative to analysis time (hour; can be positive or negative) </td></tr>
<tr><td> 9             </td><td> h(x) observation in model, either background or analysis (unit same as observation value except surface pressure in Pa) </td></tr>
<tr><td> 10            </td><td> 0 (no meaning)     </td></tr></blockquote></tbody></table>

<blockquote>In the observation departure statistics output by the LETKF program (<code>$OUTDIR/obs/obsdep</code>):<br>
<table><thead><th> <b>Column</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> 8             </td><td> Observation time relative to analysis time (hour; can be positive or negative) </td></tr>
<tr><td> 9             </td><td> (Observation - background) (O-B) or (Observation - analysis) (O-A) </td></tr>
<tr><td> 10            </td><td> 0 (no meaning)     </td></tr></blockquote></tbody></table>

<blockquote>In the observation sensitivity estimates output by the EFSO program (<code>$OUTDIR/efso</code>):<br>
<table><thead><th> <b>Column</b> </th><th> <b>Description</b> </th></thead><tbody>
<tr><td> 8             </td><td> Observation time relative to analysis time (hour; can be positive or negative) </td></tr>
<tr><td> 9             </td><td> Observation sensitivity estimated by the EFSO (J/kg) </td></tr>
<tr><td> 10            </td><td> 0 (no meaning)     </td></tr></blockquote></tbody></table>

<ul><li><b>GSI diagnostic files</b>
</li></ul><blockquote>The format of the GSI diagnostic outputs<br>Abbreviation: <b>gsidiag</b></blockquote>

<h2>Flow chart</h2>

<img src='http://miyoshi.googlecode.com/svn/wiki/images/gfs-letkf-flowchart.png' />

In the flow chart, the rectangles represent any kind of files with their formats shown in square brackets (<a href='#Data_formats.md'>explanation of the data formats and their abbreviations</a>). Those rectangles are connected by arrow lines that represent program execution, with corresponding program file names shown in the bold italic font next to the arrow lines. Three main components of the system – the data assimilation cycle, the observation processing module, and the model forecast and verification modules – are boxed by the red dashed rectangles. The data assimilation cycle is illustrated in the lower part of the figure: the 9-hour ensemble GFS model integration is executed based on the GFS sigma/surface file formats (sig/sfc), and the LETKF analysis is executed based on the grid file format (grd). The purpose of conducting 9-hour forecasts is to perform a 4-dimensional LETKF (4D-LETKF) which assimilates asynchronous observation data at their right time within a window from hour 3 to hour 9.<br>
<br>
The source of the observation data is from the NCEP PREPBUFR dataset that not only provides the observed values but also the observation errors associated with each observation. These observation errors will be used in our system. There are two routes of the observation processing:<br>
<ol><li>The route 1 shown in green arrows uses the built-in observation operators that can only process conventional (non-radiance) observation data.<br>
</li><li>The route 2 shown in blue arrows uses the GSI as the observation operator.</li></ol>

A set of reference model analysis data (gray rectangle) is needed in order to provide updated values of some prognostic variables that are not able to be analyzed by the atmospheric data assimilation system, such as ozone concentration and sea surface temperature (SST). It can also be used for the verification.<br>
<br>
<h1>Prepare model and observation data</h1>

We need to prepare a set of reference model data in order to update some "boundary conditions" such as SST and ozone concentration during the data assimilation cycle. Therefore, this dataset should cover the entire experiment period. It can also be used for verification. The NCEP has put their Climate Forecast System version 2 Reanalysis (CFSR) <a href='http://nomads.ncdc.noaa.gov/data.php#cfs'>online</a> over the period of 1979–2011. They provide the model initial conditions in the GFS sig/sfc format at the <a href='http://nomads.ncdc.noaa.gov/modeldata/cmd_LIC/'>T126</a> and the <a href='http://nomads.ncdc.noaa.gov/modeldata/cmd_HIC/'>T382</a> resolutions. It is convenient to use this dataset as the reference model data for the GFS-LETKF.<br>
<br>
There are scripts to download the CFSR (sig/sfc) data, to change the resolution, and to convert them to the gridded data in both model levels (grd) and pressure levels (grdp). In the following demonstration, everything will be run at the T62 resolution, so we are going to download the T126 data and convert them to T62.<br>
<br>
<h2>Download the reference model data</h2>

<ul><li>Go to the GFS run directory.<br>
<pre><code>$ cd gfs/run<br>
</code></pre>
</li><li>Edit the '<code>configure.sh</code>' file. Set '<code>$OUTDIR</code>' to an existing empty directory, and set other variables according to the <a href='#Main_configuration_script.md'>variable description</a>. Pay attention to the '<code>$TMP1</code>', '<code>$EXECGLOBAL</code>', '<code>$FIXGLOBAL</code>', '<code>$ANLGFS</code>', '<code>$ANLGRD</code>', and '<code>$ANLGRDP</code>' variables that will be used in this section. Create empty directories for '<code>$ANLGFS</code>', '<code>$ANLGRD</code>', and '<code>$ANLGRDP</code>'.<br>
<pre><code>$ mkdir $OUTDIR $ANLGFS $ANLGRD $ANLGRDP<br>
</code></pre>
</li><li>Create a directory to store the T126 CFSR data.<br>
<pre><code>$ mkdir $CFSR_T126<br>
</code></pre>
</li><li>Use the '<code>get_cfsr.sh</code>' script to download the CFSR data. Note that you can run it without any argument to see the usage help.<br>
<pre><code>$ ./get_cfsr.sh<br>
<br>
[get_cfsr.sh] Download NCEP CFSR.<br>
                 *use settings in 'configure.sh'<br>
<br>
Usage: ./get_cfsr.sh STIME [ETIME] [GFS_DIR]<br>
<br>
  STIME    Start time (format: YYYYMMDD)<br>
  ETIME    End   time (format: YYYYMMDD)<br>
           (default: same as STIME)<br>
  GFS_DIR  Directory to save the CFSR sig/sfc files<br>
           (default: $ANLGFS in 'configure.sh')<br>
<br>
$ ./get_cfsr.sh 20080101 20080110 $CFSR_T126<br>
</code></pre>
</li><li>Then you will see the CFSR data from January 1 to 10, 2008 are saved to the '<code>$CFSR_T126</code>' directory following the '<code>YYYYMMDDHH.sig</code>' and '<code>YYYYMMDDHH.sfc</code>' convention.</li></ul>

<h2>Change the resolution of the reference model data</h2>

<ul><li>Use the '<code>run_chgres.sh</code>' script to change the resolution of the CFSR data from T126 to T62. The '<code>$ANLGFS</code>' below should be the same directory as the variable in '<code>configure.sh</code>'.<br>
<pre><code>$ ./run_chgres.sh<br>
<br>
[run_chgres.sh] Change the resolution of GFS initial conditions.<br>
                *use settings in 'configure.sh'<br>
<br>
Usage: ./run_chgres.sh STIME ETIME JCAP_NEW LONB_NEW LATB_NEW GFS_INDIR GFS_OUTDIR [FORMAT]<br>
<br>
  STIME       Start time (format: YYYYMMDDHH)<br>
  ETIME       End   time (format: YYYYMMDDHH)<br>
  JCAP_NEW    New JCAP (spectral truncation)<br>
  LONB_NEW    New LONB (number of longitudes)<br>
  LATB_NEW    New LATB (number of latitudes)<br>
  GFS_INDIR   Input  directory of GFS sig/sfc files<br>
  GFS_OUTDIR  Output directory of GFS sig/sfc files<br>
  FORMAT      Input data format<br>
              0: Do not specify the input data format<br>
              1: GFS/GDAS data<br>
              2: CFSR data [- 2010-12-31 18]<br>
              3: CFSR data [2011-01-01 00 -]<br>
              (default: 1)<br>
<br>
$ ./run_chgres.sh 2008010100 2008011018 62 192 94 $CFSR_T126 $ANLGFS 2<br>
</code></pre>
</li><li>Then you will see the CFSR data from January 1 to 10, 2008 are converted to the T62 resolution, saved to the '<code>$ANLGFS</code>' directory.</li></ul>

<h2>Convert the reference model data to gridded data</h2>

<ul><li>Use the '<code>mss2grd.sh</code>' script to convert the T62 CFSR sig/sfc data (in '<code>$ANLGFS</code>' directory) to the gridded data in both model levels and pressure levels (grd/grdp; in '<code>$ANLGRD</code>' and '<code>$ANLGRDP</code>' directories). Pay attention to the '<code>$ANLGFS</code>', '<code>$ANLGRD</code>', and '<code>$ANLGRDP</code>' variables in the '<code>configure.sh</code>' file.<br>
<pre><code>$ ./mss2grd.sh<br>
<br>
[mss2grd.sh] Convert NCEP sig/sfc files to GrADS grd/grdp files.<br>
             *use settings in 'configure.sh'<br>
<br>
Usage: ./mss2grd.sh STIME [ETIME]<br>
<br>
  STIME  Start time (format: YYYYMMDDHH)<br>
  ETIME  End   time (format: YYYYMMDDHH)<br>
         (default: same as STIME)<br>
<br>
$ ./mss2grd.sh 2008010100 2008011018<br>
</code></pre>
</li><li>Then you will see the CFSR data from January 1 to 10, 2008 are converted to the gridded data following the '<code>YYYYMMDDHH.grd</code>' convention, saved to the '<code>$ANLGRD</code>' and '<code>$ANLGRDP</code>' directories.<br>
</li><li>You can create GrADS CTL files by the '<code>grdctl</code>' program and use the GrADS software to visualize those data.<br>
<pre><code>$ ../ssio/grdctl<br>
 <br>
Usage: ./grdctl DSET OPTIONS T_START T_INT T_NUM GRD_TYPE<br>
 <br>
  GRD_TYPE: 's': Sigma-level output<br>
            'x': Extended sigma-level output<br>
            'p': Pressure-level output<br>
 <br>
FORTRAN STOP<br>
$ ../ssio/grdctl "%y4%m2%d2%h2.grd" "template byteswapped" 00Z01Jan2008 6hr 10000 x &gt; $ANLGRD/yyyymmddhhx.ctl<br>
$ ../ssio/grdctl "%y4%m2%d2%h2.grd" "template byteswapped" 00Z01Jan2008 6hr 10000 p &gt; $ANLGRDP/yyyymmddhhp.ctl<br>
</code></pre></li></ul>

<h2>Prepare the GFS initial conditions</h2>

One of the simplest ways to create an initial ensemble is using a combination of initial conditions at different times. In this way, a series of reference model data, such as the CFSR, at an unrelated time period could be used to form an initial ensemble. However, when we choose the CFSR to do so, a consistent temperature bias is observed near the tropopause which can be as large as -8 K, especially near the polar region. It is an unacceptable huge bias that can significantly degrade the LETKF data assimilation performance. For this reason, it is recommended to use the operational GFS model analyses to form the initial ensemble.<br>
<br>
The latest one month operational GFS model analyses can be download from <a href='http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/'>http://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/</a> . Earlier data may not be found via the public Internet. Please put these data into the '<code>$INITGFS</code>' directory configured in '<code>configure.sh</code>'. Note that these data also need to be converted to the T62 resolution before using them.<br>
<br>
<h2>Download the NCEP observation data</h2>

<ul><li>Get and compile the '<code>grabbufr</code>' utility that is used to pre-process the NCEP BUFR data. It can be downloaded from <a href='http://ftp.emc.ncep.noaa.gov/gc_wmb/jwoollen/grabbufr/'>here</a>. Please compile it and name the executable file as '<code>grabbufr</code>', put into the '<code>$BUFRBIN</code>' directory set in '<code>configure.sh</code>'.</li></ul>

<ul><li>Edit the '<code>configure.sh</code>' file. Pay attention to the '<code>$OBS</code>', '<code>$OBSNCEP</code>', and '<code>$BUFRBIN</code>' variables. Create empty directories for '<code>$OBS</code>' and '<code>$OBSNCEP</code>'.<br>
<pre><code>$ mkdir $OBS $OBSNCEP<br>
</code></pre>
</li><li>We are going to download the NCEP PREPBUFR data from the <a href='http://rda.ucar.edu/datasets/ds337.0/'>CISL Research Data Archive</a>. You need to register (free) on this website to get an account.<br>
</li><li>Use the '<code>get_ncepobs.sh</code>' script to download the NCEP PREPBUFR data and convert them to the LETKF observation format.<br>
<pre><code>$ ./get_ncepobs.sh<br>
<br>
[get_ncepobs.sh] Download NCEP conventional observation data.<br>
                 *use settings in 'configure.sh'<br>
<br>
Usage: ./get_ncepobs.sh EMAIL PASSWD STIME [ETIME] [IF_DECODE]<br>
<br>
  EMAIL      UCAR/DSS account (register at http://rda.ucar.edu )<br>
  PASSWD     UCAR/DSS account password<br>
  STIME      Start time (format: YYYYMMDD)<br>
  ETIME      End   time (format: YYYYMMDD)<br>
             (default: same as STIME)<br>
  IF_DECODE  Convert PREPBUFR to LETKF obs format or not<br>
             0: No,  store only PREPBUFR format<br>
             1: Yes, decode to LETKF obs format<br>
             (default: Yes)<br>
<br>
$ ./get_ncepobs.sh EMAIL PASSWD 20080101 20080110 1<br>
</code></pre>
</li><li>Then you will see the NCEP PREPBUFR data from January 1 to 10, 2008 are saved in the '<code>$OBSNCEP</code>' directory, and they are also converted to the LETKF observation format saved in the '<code>$OBS</code>' directory.</li></ul>

<h1>Run data assimilation cycles</h1>

Here we demonstrate how to run the data assimilation cycle with 32 members.<br>
<br>
<h2>Create the initial ensemble</h2>

<ul><li>Check the '<code>configure.sh</code>' file, set the '<code>$MEMBER</code>' variable to 32 if it is not.<br>
</li><li>Edit the '<code>common/common_letkf.f90</code>' file. Set the '<code>nbv</code>' variable to 32 (default is 20).<br>
</li><li>Recompile all GFS-LETKF programs.<br>
<pre><code>$ cd ..   (go to the 'gfs' top directory)<br>
$ make<br>
$ cd run<br>
</code></pre>
</li><li>The '<code>init.sh</code>' script takes a series of the model analysis fields at an unrelated time period, saved in the '<code>$INITGFS</code>' directory in '<code>configure.sh</code>', to form the initial ensemble. Below we use the operational GFS model analyses from 00Z December 15, 2011 to 00Z January 15, 2012 with a 24-hour interval to form a 32-member initial ensemble. Note that the 24-hour interval is recommended in order to avoid the issue with the diurnal cycle.<br>
<pre><code>$ ./init.sh<br>
<br>
[init.sh] Prepare a initial ensemble from analyses at different time.<br>
          *use settings in 'configure.sh'<br>
<br>
Usage: ./init.sh STIME RTIME [RINT_H]<br>
<br>
  STIME   Initial time of the ensemble (format: YYYYMMDDHH)<br>
  RTIME   An arbitrary time that the first member is taken from (format: YYYYMMDDHH)<br>
  RINT_H  Interval between each members (hour) (default: 24)<br>
<br>
  For example, if RTIME = 1991010100, RINT_H = 24,<br>
  then the initial ensemble are constructed by analyses at<br>
  1991010100, 1991010200, 1991010300, ...<br>
<br>
$ ./init.sh 2008010100 2011121500 24<br>
</code></pre>
</li><li>Then you will see the output directory structure is created in the '<code>$OUTDIR</code>' directory and the initial ensemble is prepared in the '<code>$OUTDIR/anal/mmm</code>' directories.</li></ul>

<h2>Run forecast/analysis cycles</h2>

<ul><li>Now we are going to run parallel programs. Before running parallel programs, we need to create the "machinefile" listing the computing nodes available for use. In the '<code>gfs/run/machine</code>' directory, there is a example of the machinefile. One CPU core corresponds to one line. Note: if your system has a job scheduling system, skip this two steps and see the following explanation.<br>
<pre><code>$ cp machine/machinefile_4nodes machinefile<br>
</code></pre>
</li><li>Use the '<code>mcycle.sh</code>' script to run multiple forecast/analysis cycles. Here we run 4 cycles (1 day) for demonstration.<br>
<pre><code>$ ./mcycle.sh<br>
<br>
[mcycle.sh] Run multiple cycles.<br>
            *use settings in 'configure.sh'<br>
<br>
Usage: ./mcycle.sh STIME [ETIME]<br>
<br>
  STIME   Time of the first cycle (format: YYYYMMDDHH)<br>
  ETIME   Time of the last  cycle (format: YYYYMMDDHH)<br>
          (default: same as STIME)<br>
<br>
$ ./mcycle.sh 2008010100 2008010118 &amp;<br>
</code></pre>
</li><li>If you use a cluster with a job scheduling system, you need to submit jobs via the system. Instead of doing the above two steps, you need to create a script to submit the job according to your job scheduling system. Copy the machinefile the scheduling system gives to '<code>gfs/run/machinefile</code>' before running the GFS-LETKF scripts. An example of the script for the PBS system is provided at '<code>pbs_example.sh</code>'.<br>
</li><li>The data assimilation cycles take a while to finish. You can monitor the progress and (possible) error messages in the '<code>gfs/run/log</code>' directory while the execution. Using 32 cores of the AMD Opteron(tm) 2378 CPU at 2.4 GHz (8 X quad-core processors), it takes roughly 15 minutes for one 6-hour forecast/analysis cycle, and 1 hour for 4 cycles.<br>
</li><li>Finally you will see the results in the '<code>$OUTDIR</code>' directory. You can use the GrADS software to visualize the results in the '<code>$OUTDIR/guesg</code>', '<code>$OUTDIR/analg</code>' (model levels) and the '<code>$OUTDIR/guesgp</code>', '<code>$OUTDIR/analgp</code>' (pressure levels) directories.</li></ul>

<h1>Run (ensemble) forecasts and verification</h1>

Here we demonstrate how to run the 5-day forecasts initialized from the ensemble mean analysis every cycle.<br>
<br>
<ul><li>Use the '<code>mfcst.sh</code>' script to run multiple 5-day forecasts initialized from the 4-cycle ensemble means we previously got with the data assimilation cycles. The four cycles can be run parallelly in order to accelerate the computation. The last argument '<code>1</code>' means performing verification as well after the forecasts finish. Type the '<code>mfcst.sh</code>' command without any argument to see the usage help.<br>
<pre><code>$ ./mfcst.sh<br>
<br>
[mfcst.sh] Run multiple ensemble mean forecasts.<br>
           *use settings in 'configure.sh'<br>
<br>
Usage: ./mfcst.sh STIME [ETIME] [CYCLES] [IF_VERF]<br>
<br>
  STIME    Time of the first cycle (format: YYYYMMDDHH)<br>
  ETIME    Time of the last  cycle (format: YYYYMMDDHH)<br>
           (default: same as STIME)<br>
  CYCLES   Number of forecast cycles run in parallel<br>
           (default: 1)<br>
  IF_VERF  Run verification or not<br>
           0: Do not run verification<br>
           1: Run verification<br>
           (default: 0)<br>
<br>
$ ./mfcst.sh 2008010106 2008010200 4 1 &amp;<br>
</code></pre>
</li><li>You can monitor the progress and (possible) error messages in the '<code>gfs/run/log</code>' directory while the execution. This ensemble mean forecasts and verification takes roughly 15 minutes using the same configuration of the 32 cores.<br>
</li><li>You will see the results in the '<code>$OUTDIR/fcst*</code>' directory. You can use the GrADS software to visualize the results in the '<code>$OUTDIR/fcstg</code>', '<code>$OUTDIR/fcstgp</code>' (sorted according to the initial times) and the '<code>$OUTDIR/fcstv</code>', '<code>$OUTDIR/fcstvp</code>' (sorted according to the forecast hours) directories. The results of verification against (rawinsonde) observations are saved in the '<code>$OUTDIR/verfo1</code>' directory and verification against reference model analyses are saved in the '<code>$OUTDIR/verfa1</code>' directory.</li></ul>

<h1>Compute the EFSO</h1>

Here we demonstrate how to compute the ensemble forecast sensitivity to observations (EFSO) with the 6 hour evaluation time. To compute the EFSO at time t = 0, evaluated at t = 6 h, we will use the forecasts from the ensemble mean at t = 0 and t = -6 h, both valid at t = 6 h. This should be already done by the '<code>mfcst.sh</code>' step. Besides, it also needs the <b>ensemble</b> forecasts from all members at t = 0, valid at t = 6 h, which can be done by '<code>efsofcst.sh</code>' separately. In addition, we further need the observation values in the analysis at t = 0 [H(X<sup>a</sup>)], which need to be output from the LETKF program with a flag turned on.<br>
<br>
<ul><li>We are going to compute the EFSO at 12Z January 1, 2008. We need to compute the LETKF analysis step again at this time if the [H(X<sup>a</sup>)] output flag was not turned on. Edit the '<code>gfs/letkf/letkf_obs.f90</code>' file, set both the '<code>obsgues_output</code>' and '<code>obsanal_output</code>' variables to '<code>.TRUE.</code>'.<br>
</li><li>Recompile all GFS-LETKF programs.<br>
<pre><code>$ cd ..   (go to the 'gfs' top directory)<br>
$ make<br>
$ cd run<br>
</code></pre>
</li><li>Edit the '<code>configure.sh</code>' file, set '<code>EFSOT</code>' to '<code>6</code>' and set '<code>$OBSOUT_DIAG</code>' to '<code>1</code>'.<br>
</li><li>Run the LETKF data assimilation cycle again with the [H(X<sup>a</sup>)] output flag turned on.<br>
<pre><code>$ ./mcycle.sh 2008010106 &amp;<br>
</code></pre>
</li><li>Then the [H(X<sup>a</sup>)] data will be saved to '<code>$OUTDIR/obs/obsdiag/2008010112</code>'. They are in the <a href='#Data_formats.md'>"letkfobs2" format</a>.<br>
</li><li>Assuming the ensemble mean forecasts ('<code>mfcst.sh</code>') have been conducted, use the '<code>efsofcst.sh</code>' script to conduct the 6-hour ensemble forecasts.<br>
<pre><code>$ ./efsofcst.sh<br>
<br>
[efsofcst.sh] Run ensemble forecasts for EFSO computation.<br>
              *use settings in 'configure.sh'<br>
<br>
Usage: ./efsofcst.sh STIME [ETIME] [IF_MEAN] [CYCLES]<br>
<br>
  STIME    Time of the first cycle (format: YYYYMMDDHH)<br>
  ETIME    Time of the last  cycle (format: YYYYMMDDHH)<br>
           (default: same as STIME)<br>
  IF_MEAN  Also run the ensemble mean forecast?<br>
           0: No<br>
           1: Yes<br>
           (default: 0)<br>
  CYCLES   Number of forecast cycles run in parallel<br>
           (default: 1)<br>
<br>
$ ./efsofcst.sh 2008010112 2008010112 &amp;<br>
</code></pre>
</li><li>You can monitor the progress and (possible) error messages in the '<code>gfs/run/log</code>' directory.<br>
</li><li>Compute the EFSO at 12Z January 1, 2008 by the '<code>mefso.sh</code>' script.<br>
<pre><code>$ ./mefso.sh<br>
<br>
[mefso.sh] Run EFSO for multiple cycles.<br>
            *use settings in 'configure.sh'<br>
<br>
Usage: ./mefso.sh STIME [ETIME] [EFT] [LOCADV_RATE] [WMOIST]<br>
<br>
  STIME        Time of the first cycle (format: YYYYMMDDHH)<br>
  ETIME        Time of the last  cycle (format: YYYYMMDDHH)<br>
               (default: same as STIME)<br>
  EFT          Evaluation forecast time (hours)<br>
               (default: 24)<br>
  LOCADV_RATE  Localization advection rate relative to the phase velocity (winds)<br>
               0: No advection<br>
               (default: 0)<br>
  WMOIST       Wight for the moist term in the energy norm (dimensionless)<br>
               (default: 1)<br>
<br>
$ ./mefso.sh 2008010112 2008010112 6 &amp;<br>
</code></pre>
</li><li>You can monitor the progress and (possible) error messages in the '<code>gfs/run/log</code>' directory, and then the EFSO results will be stored in the '<code>$OUTDIR/efso</code>' directory, also in the <a href='#Data_formats.md'>"letkfobs2" format</a>.