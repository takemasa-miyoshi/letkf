

# SPEEDY model #
  * **SPEEDY** is a simplified AGCM (Atmospheric General Circulation Model) developed by [Molteni (2003)](http://dx.doi.org/10.1007/s00382-002-0268-2). It has only 96x48x7 grid points for Earth's global atmosphere. It is feasible for single processor PC (like your desktop or laptop PC) with a Fortran 90 compiler and shell scripting environment. I have tested with a Linux PC with Intel Celeron D 2.16 GHz processor (about 5 years old), as well as very recent Intel Xeon quad-core 2.5 GHz processors. It takes about 2 seconds for 6-hr forecast with the Celeron D processor, and about 0.4 seconds with the Xeon processor (single process).
  * [Franco Molteni](http://users.ictp.it/~moltenif/index.html) and [Fred Kucharski](http://users.ictp.it/~kucharsk/) should be acknowledged for the use of SPEEDY model.

# Getting started #
  1. Download the latest release version from the [Source tab](http://code.google.com/p/miyoshi/source/checkout).
  1. Run the SPEEDY model
    * Dig into the model directory
```
$cd speedy/model/run
```
    * Execute run\_first.sh
```
$sh run_first.sh
```
> > This executes 1-year integration. It may take time according to your computer. With an Intel Xeon 2.5GHz Quad-core processor, it takes about 6 min with a single process (without making use of multi-core function).
    * Edit and execute run\_cycle.sh<br>Edit the final date setting for your convenience. Just for testing, one day may be sufficient.<br>
<pre><code>$sh run_cycle.sh<br>
</code></pre>
<ul><li>At this point, you already got the "nature run" of the SPEEDY model. You can visualize the data using <a href='http://www.iges.org/grads/'>GrADS software</a>.<br>
<pre><code>$cd ../DATA<br>
$grads<br>
ga-&gt; open nature/yyyymmddhh.ctl<br>
</code></pre>
</li></ul><ol><li>Generate simulated observation data<br>
<ul><li>Change directory to the observation tools directory<br>
<pre><code>$cd ../obs<br>
</code></pre>
</li><li>Compile the program<br>
<pre><code>$sh make_obsmake.sh<br>
</code></pre>
</li><li>Edit and execute obsmake.sh<br>
</li></ul><blockquote>Edit the final date setting shorter than the one that you generated the nature run.<br>
<pre><code>$sh obsmake.sh<br>
</code></pre>
</blockquote></li><li>Run the SPEEDY-LETKF forecast/assimilation cycle experiment<br>
<ul><li>Change directory to the LETKF directory<br>
<pre><code>$cd ../letkf<br>
</code></pre>
</li><li>Compile the program<br>
<pre><code>$sh make_letkf.sh<br>
</code></pre>
</li><li>Change directory to the run directory<br>
<pre><code>$cd run<br>
</code></pre>
</li><li>Initialize the experiment<br>
<pre><code>$sh init.sh<br>
</code></pre>
</li><li>Edit and execute run_cycle.sh<br>
</li></ul><blockquote>Edit the final date setting shorter than the one that you generated the observation data.<br>
<pre><code>$sh run_cycle.sh<br>
</code></pre>
It may take time according to your computing environment. LETKF data assimilation is executed first, then, 20-member ensemble prediction follows. With a computer with 2 quad-core processors, it takes about 2-3 sec for LETKF assimilation and 4 sec for 20-member ensemble prediction with 4 parallel processes (making use of multi-core processors).<br>
</blockquote></li><li>Monitor the standard output<br>
</li></ol><blockquote>You will see in your console as follows:<br>
<pre><code>== PARTIAL OBSERVATIONAL DEPARTURE (gues) ==============================<br>
           U           V           T           Q          PS          RH<br>
  -1.071E-01   8.660E-03   2.719E-03  -6.378E-05  -3.453E+00  -9.990E+33<br>
   1.950E+00   1.877E+00   1.291E+00   1.042E-03   1.629E+02  -9.990E+33<br>
== NUMBER OF OBSERVATIONS ==============================================<br>
           U           V           T           Q          PS          RH<br>
        5024        5035        5040        5040        1005           0<br>
========================================================================<br>
== PARTIAL OBSERVATIONAL DEPARTURE (anal) ==============================<br>
           U           V           T           Q          PS          RH<br>
  -6.608E-02   2.802E-03  -2.012E-04  -6.380E-05  -1.053E+00  -9.990E+33<br>
   1.249E+00   1.251E+00   1.162E+00   1.024E-03   1.165E+02  -9.990E+33<br>
== NUMBER OF OBSERVATIONS ==============================================<br>
           U           V           T           Q          PS          RH<br>
        5024        5035        5040        5040        1005           0<br>
========================================================================<br>
</code></pre>
The upper section indicates first guess (gues), i.e., before assimilation, and the bottom section indicates analysis (anal), i.e., after assimilation. The floating numbers indicate bias (first row) and root mean square errors (rmse, second row) relative to observations (not nature run!). For example, The first guess of U (zonal wind) has -0.1071 m/s bias and 1.950 m/s rmse, that are reduced to be -0.06608 and 1.249 respectively after the LETKF analysis. If you see decrease of numbers in the analysis, that means LETKF works properly.<br>
<pre><code>### TIMER(MONIT_MEAN):      2.32      0.05<br>
</code></pre>
This line indicates how many seconds it took for the LETKF computation. In this case, 2.32 seconds.</blockquote>

<h1>LETKF parameter settings</h1>
<table><thead><th> <b>parameter</b> </th><th> <b>source file</b> </th><th> <b>variable name</b> </th><th> <b>unit</b> </th></thead><tbody>
<tr><td> ensemble size    </td><td> common/common_letkf.f90 </td><td> nbv                  </td><td> members     </td></tr>
<tr><td> horizontal localization </td><td> speedy/letkf/letkf_obs.f90 </td><td> sigma_obs            </td><td> meters      </td></tr>
<tr><td> vertical localization </td><td> speedy/letkf/letkf_obs.f90 </td><td> sigma_obsv           </td><td> ln p        </td></tr>
<tr><td> temporal localization </td><td> speedy/letkf/letkf_obs.f90 </td><td> sigma_obst           </td><td> time slots  </td></tr>
<tr><td> multiplicative inflation </td><td> speedy/letkf/letkf_tools.f90 </td><td> cov_infl_mul         </td><td> positive: fixed inflation (1.00 for no inflation), negative: adaptive inflation </td></tr></tbody></table>

<h1>Script settings</h1>
<table><thead><th> <b>source file</b> </th><th> <b>variable name</b> </th><th> <b>setting</b> </th></thead><tbody>
<tr><td> speedy/obs/obsmake.sh </td><td> EXP                  </td><td> name of obs network </td></tr>
<tr><td> speedy/letkf/run/init.sh </td><td> OBS                  </td><td> name of obs network </td></tr>
<tr><td> speedy/letkf/run/init.sh </td><td> EXP                  </td><td> name of experiment </td></tr>
<tr><td> speedy/letkf/run/run_cycle.sh </td><td> OBS                  </td><td> name of obs network </td></tr>
<tr><td> speedy/letkf/run/run_cycle.sh </td><td> EXP                  </td><td> name of experiment </td></tr>