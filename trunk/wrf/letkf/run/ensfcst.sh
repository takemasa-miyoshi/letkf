#!/bin/sh
# ===================================================================
# ensfcst.sh
#   This script runs the WRF model
# ===================================================================
### input for this shell
TMPDIR=$1
IYMDH=$2
FYMDH=$3
MEM=$4
MPI_PER_MEMBER=$5
FNLDIR=$6
INIT_DATA=$7
FLAG=$8
# ==============================
#    PREAMBLE 
# ==============================
syy=`echo $IYMDH | cut -c1-4`
smm=`echo $IYMDH | cut -c5-6`
sdd=`echo $IYMDH | cut -c7-8`
shh=`echo $IYMDH | cut -c9-10`
smn=00

eyy=`echo $FYMDH | cut -c1-4`
emm=`echo $FYMDH | cut -c5-6`
edd=`echo $FYMDH | cut -c7-8`
ehh=`echo $FYMDH | cut -c9-10`
emn=00

wkdir_root=${TMPDIR}/gues/${MEM}
io_data_fname=fnl_${syy}${smm}
io_data_type=GFS

wk_common_geog=${TMPDIR}/model/GEOG
wk_common_wps=${TMPDIR}/model/${MEM}/WPS
wk_common_wrf=${TMPDIR}/model/${MEM}/WRF
wk_common_post=${TMPDIR}/model/${MEM}/ARWpost

mpibin=/usr/local/mpich2/bin

start_step=0
end_step=410
num_frames=1	# set 1 or 1000
num_mpi=${MPI_PER_MEMBER}
window_length=6
ft_lag=3

ft=`expr ${window_length} + ${ft_lag}`

source ${TMPDIR}/model/util.sh
# ==============================
#    STEPS 
# ==============================
# ------------------------------
#    STEP 100: Init
# ------------------------------  
step=100_Init
# --- set model directory ---
if test ! -d ${TMPDIR}/model/${MEM}/WPS -o ! -d ${TMPDIR}/model/${MEM}/WRF
then		
test -d ${TMPDIR}/model/${MEM} || mkdir -p ${TMPDIR}/model/${MEM}
cp -RL ${TMPDIR}/model/WPS     ${TMPDIR}/model/${MEM}
cp -RL ${TMPDIR}/model/WRF     ${TMPDIR}/model/${MEM}		
cp -RL ${TMPDIR}/model/ARWpost     ${TMPDIR}/model/${MEM}		
fi
# --- set work directory ---
anal_hour=`expr ${ft_lag} + ${window_length} / 2`
anal_min=`expr ${anal_hour} \* 60`
date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${anal_min} > anal.$$
read ayy amm add ahh amn < anal.$$
rm -f anal.$$

wkdir=${wkdir_root}/${ayy}${amm}${add}${ahh}${amn}
test -d ${wkdir} || mkdir -p ${wkdir}
cd ${wk_common_wps}
# --- set dummy time
case ${ehh} in
03|09|15|21)
min=180
date_edit ${eyy} ${emm} ${edd} ${ehh} ${emn} ${min} > dummy.txt
;;
*)
echo "${eyy} ${emm} ${edd} ${ehh} ${emn}" > dummy.txt
;;
esac
read dyy dmm ddd dhh dmn < dummy.txt
rm -f dummy.txt	
# --- make namelist ---
cat << EOF >  namelist.wps	
&share
 wrf_core = 'ARW',
 max_dom = 1,
 start_date = "${syy}-${smm}-${sdd}_${shh}:${smn}:00"
 end_date   = "${dyy}-${dmm}-${ddd}_${dhh}:${dmn}:00"
 interval_seconds = 10800
 io_form_geogrid = 2,
/

&geogrid
 parent_id         =   1,
 parent_grid_ratio =   1,
 i_parent_start    =   1,
 j_parent_start    =   1,
 e_we              = 137,
 e_sn              = 109,
 geog_data_res     = '10m',
 dx = 60000,
 dy = 60000,
 map_proj = 'mercator',
 ref_lat   =  30.0,
 ref_lon   = 140.0,
 truelat1  =  22.5,
 truelat2  =  60.0,
 stand_lon = 140.0,
 geog_data_path = "${wk_common_geog}"
/

&ungrib
 out_format = 'WPS',
 prefix = "${io_data_type}",
/

&metgrid
 fg_name = "${io_data_type}",
 io_form_metgrid = 2,
/
EOF
# ------------------------------
#    STEP 200: GEOGRID
# ------------------------------ 
step=200
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
# --- run the pgm ---
./geogrid.exe > /dev/null 2> /dev/null
fi
# ------------------------------
#    STEP 210: UNGRIB
# ------------------------------ 
step=210
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
# --- link Vtable ---
vtable=ungrib/Variable_Tables/Vtable.${io_data_type}
if test ! -f ${vtable}
then
echo "xxx   not found ${vtable}   xxx"
exit
fi	
rm -f Vtable
ln -s ungrib/Variable_Tables/Vtable.${io_data_type} Vtable
# --- link iodata ---
date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} 1440 > next.$$
read nyy nmm ndd nhh nmn < next.$$
rm -f next.$$	
csh link_grib.csh $FNLDIR/fnl_${syy}${smm}${sdd}*  $FNLDIR/fnl_${nyy}${nmm}${ndd}* 
# --- run the pgm ---
rm -f ${io_data_type}:*	
./ungrib.exe > /dev/null 2> /dev/null
# --- end trunsaction ---
rm -f GRIBFILE.???
fi
# ------------------------------
#    STEP 220: METGRID
# ------------------------------ 
step=220
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
# --- run the pgm ---
rm -f met_em.d01.*
./metgrid.exe > /dev/null 2> /dev/null
init_file=met_em.d01.${syy}-${smm}-${sdd}_${shh}:${smn}:00.nc
# --- overwrite initial data for first run ---
if test -f ${INIT_DATA} -a $FLAG -eq 0
then
echo "@@@ Overwrite Initial Data after METGRID @@@"
cp ${INIT_DATA} ./${init_file}
fi		
# --- end transaction ---
cp ${init_file} ${wkdir}
fi
# ------------------------------
#    STEP 300: Init
# ------------------------------
step=300
cd ${wk_common_wrf}/test/em_real
# --- make namelist ---
cat << EOF >  namelist.input	
 &time_control
 run_days                            = 0,
 run_hours                           = ${ft},
 run_minutes                         = 0,
 run_seconds                         = 0,
 start_year                          = ${syy}, 2000, 2000,
 start_month                         = ${smm},   01,   01,
 start_day                           = ${sdd},   24,   24,
 start_hour                          = ${shh},   12,   12,
 start_minute                        = ${smn},   00,   00,
 start_second                        = 00,   00,   00,
 end_year                            = ${eyy}, 2000, 2000,
 end_month                           = ${emm},   01,   01,
 end_day                             = ${edd},   25,   25,
 end_hour                            = ${ehh},   12,   12,
 end_minute                          = ${emn},   00,   00,
 end_second                          = 00,   00,   00,
 interval_seconds                    = 10800
 input_from_file                     = .true.,.true.,.true.,
 history_interval                    = 60,  60,   60,
 frames_per_outfile                  = ${num_frames}, 1000, 1000,
 restart                             = .false.,
 restart_interval                    = 5000,
 io_form_history                     = 2
 io_form_restart                     = 2
 io_form_input                       = 2
 io_form_boundary                    = 2
 debug_level                         = 0
 /

 &domains
 time_step                           = 180,
 time_step_fract_num                 = 0,
 time_step_fract_den                 = 1,
 max_dom                             = 1,
 e_we                                = 137,    112,   94,
 e_sn                                = 109,     97,    91,
 e_vert                              = 40,    28,    28,
 p_top_requested                     = 1000,
 num_metgrid_levels                  = 27,
 num_metgrid_soil_levels             = 4,
 dx                                  = 60000, 10000,  3333.33,
 dy                                  = 60000, 10000,  3333.33,
 grid_id                             = 1,     2,     3,
 parent_id                           = 0,     1,     2,
 i_parent_start                      = 1,     31,    30,
 j_parent_start                      = 1,     17,    30,
 parent_grid_ratio                   = 1,     3,     3,
 parent_time_step_ratio              = 1,     3,     3,
 feedback                            = 1,
 smooth_option                       = 0
 /

 &physics
 mp_physics                          = 3,     3,     3,
 ra_lw_physics                       = 1,     1,     1,
 ra_sw_physics                       = 1,     1,     1,
 radt                                = 60,    30,    30,
 sf_sfclay_physics                   = 1,     1,     1,
 sf_surface_physics                  = 2,     2,     2,
 bl_pbl_physics                      = 1,     1,     1,
 bldt                                = 0,     0,     0,
 cu_physics                          = 1,     1,     0,
 cudt                                = 5,     5,     5,
 isfflx                              = 1,
 ifsnow                              = 0,
 icloud                              = 1,
 surface_input_source                = 1,
 num_soil_layers                     = 4,
 sf_urban_physics                    = 0,     0,     0,
 maxiens                             = 1,
 maxens                              = 3,
 maxens2                             = 3,
 maxens3                             = 16,
 ensdim                              = 144,
 /

 &fdda
 /

 &dynamics
 w_damping                           = 0,
 diff_opt                            = 1,
 km_opt                              = 4,
 diff_6th_opt                        = 0,      0,      0,
 diff_6th_factor                     = 0.12,   0.12,   0.12,
 base_temp                           = 290.
 damp_opt                            = 0,
 zdamp                               = 5000.,  5000.,  5000.,
 dampcoef                            = 0.2,    0.2,    0.2
 khdif                               = 0,      0,      0,
 kvdif                               = 0,      0,      0,
 non_hydrostatic                     = .true., .true., .true.,
 moist_adv_opt                       = 1,      1,      1,     
 scalar_adv_opt                      = 1,      1,      1,     
 /

 &bdy_control
 spec_bdy_width                      = 5,
 spec_zone                           = 1,
 relax_zone                          = 4,
 specified                           = .true., .false.,.false.,
 nested                              = .false., .true., .true.,
 /

 &grib2
 /

 &namelist_quilt
 nio_tasks_per_group = 0,
 nio_groups = 1,
 /
EOF
# ------------------------------
#    STEP 310: REAL
# ------------------------------ 
step=310
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
# -- remove logs ---
rm -f  rsl_error.* rsl_out.*
# --- input files ---
rm -f met_em*.nc
ln -s ${wk_common_wps}/met_em*.nc .
# --- run the pgm ---
${mpibin}/mpiexec -n ${num_mpi} ./real.exe < /dev/null > /dev/null 2> /dev/null
# --- use analysis ---
if test -e ${INIT_DATA} -a $FLAG -eq 2
then
echo "@@@ Overwrite Initial Data after REAL @@@"
rm -f wrfinput_d01
cp ${INIT_DATA} ./wrfinput_d01	
elif test -e ${INIT_DATA} -a $FLAG -eq 3
then
echo "@@@ Merge Initial Data after REAL (based on other analysis) @@@"
cp wrfinput_d01 wrfinput_d01.org 
ln -s wrfinput_d01 input.nc
ln -s ${INIT_DATA} anal.grd
$TMPDIR/ncio/init_merge > log.init_merge${MEM} 2> /dev/null
rm -f input.nc anal.grd
elif test -e ${INIT_DATA} -a $FLAG -eq 4
then
echo "@@@ Merge Initial Data after REAL (based on previous forecast) @@@"
mv wrfinput_d01 wrfinput_d01.org 
fcst_result=$TMPDIR/gues/${MEM}/${syy}${smm}${sdd}${shh}${smn}/wrfout_d01_${syy}-${smm}-${sdd}_${shh}:${smn}:00
cp ${fcst_result} wrfinput_d01
ln -s wrfinput_d01 input.nc
ln -s ${INIT_DATA} anal.grd
$TMPDIR/ncio/init_merge > log.init_merge${MEM} 2> /dev/null
rm -f input.nc anal.grd	
fi
fi
# ------------------------------
#    STEP 320: WRF
# ------------------------------ 
step=320
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
# -- remove logs ---
rm -f rsl_error.* rsl_out.*
rm -f wrfout_d01*
# --- run the pgm ---
${mpibin}/mpiexec -n ${num_mpi} ./wrf.exe < /dev/null > /dev/null 2> /dev/null
# --- end transaction ---
rm -f met_em*.nc
rm -f log_error log_out
i=0
while test ${i} -lt ${num_mpi}
do
cat rsl.error.000${i} >> log_error
cat rsl.out.000${i}   >> log_out
i=`expr ${i} + 1`
done
mv wrfinput_d01 ${wkdir}
mv wrfout_d01* ${wkdir}
mv log_error ${wkdir}
mv log_out   ${wkdir}
mv namelist.input  ${wkdir}
mv namelist.output ${wkdir}
fi
# ------------------------------
#    STEP 400: POST
# ------------------------------ 
step=400
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
cd ${wk_common_post}
it=0
# --- output ft = 3 - 9
#	while [ ${it} -le ${window_length} ]
#	do
#		diff_hour=` expr ${it} + ${ft_lag}`
#		min=`expr ${diff_hour} \* 60`
#		date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${min} > vdate.txt
#		read vyy vmm vdd vhh vmn < vdate.txt
#		rm -f vdate.txt
# --- output ft = 0 - 9
while test ${it} -le ${ft}
do
diff_hour=${it}
min=`expr ${diff_hour} \* 60`
date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${min} > vdate.txt
read vyy vmm vdd vhh vmn < vdate.txt
rm -f vdate.txt
# --- set file date ---
if test ${num_frames} -eq 1
then
fyy=${vyy}
fmm=${vmm}
fdd=${vdd}
fhh=${vhh}
fmn=${vmn}									
else
fyy=${syy}
fmm=${smm}
fdd=${sdd}
fhh=${shh}
fmn=${smn}		
fi
# --- make namelist (model level) ---
cat << EOF >  namelist.ARWpost
&datetime
 start_date = "${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00",
 end_date   = "${vyy}-${vmm}-${vdd}_${vhh}:${vmn}:00",
 interval_seconds = 3600,
 tacc = 0,
 debug_level = 0,
/

&io
 io_form_input  = 2,
 input_root_name = "./wrfout_d01_${fyy}-${fmm}-${fdd}_${fhh}:${fmn}:00"
 output_root_name = "./${vyy}${vmm}${vdd}${vhh}${vmn}"
 plot = 'list'
 fields = 'U, V, W, QVAPOR, PSFC, RAINC, RAINNC, HGT, height, geopt, tk, rh, rh2, u10m, v10m, slp'
 output_type = 'grads'
 mercator_defs = .true.
/
 split_output = .true.
 frames_per_outfile = 2

 output_type = 'grads'
 output_type = 'v5d'

 plot = 'all'
 plot = 'list' 
 plot = 'all_list'
! Below is a list of all available diagnostics
 fields = 'height,geopt,theta,tc,tk,td,td2,rh,rh2,umet,vmet,pressure,u10m,v10m,wdir,wspd,wd10,ws10,slp,mcape,mcin,lcl,lfc,cape,cin,dbz,max_dbz,clfr'
 

&interp
 interp_method = 1,
 interp_levels = 1000.,975.,950.,925.,900.,850.,800.,700.,600.,500.,400.,300.,250.,200.,150.,100.,
/
extrapolate = .true.

 interp_method = 0,     ! 0 is model levels, -1 is nice height levels, 1 is user specified pressure/height

 interp_levels = 1000.,950.,900.,850.,800.,750.,700.,650.,600.,550.,500.,450.,400.,350.,300.,250.,200.,150.,100.,
 interp_levels = 0.25, 0.50, 0.75, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,	
EOF
# --- input files ---
rm -f wrfout_d01_${fyy}-${fmm}-${fdd}_${fhh}:${fmn}:00
ln -s ${wkdir}/wrfout_d01_${fyy}-${fmm}-${fdd}_${fhh}:${fmn}:00 .
# --- run the pgm ---
./ARWpost.exe > /dev/null 2> /dev/null
# --- end transaction
test -d ${wkdir}/plev || mkdir -p ${wkdir}/plev
rm -f wrfout_d01_${fyy}-${fmm}-${fdd}_${fhh}:${fmn}:00
mv ${vyy}${vmm}${vdd}${vhh}${vmn}.ctl ${wkdir}/plev
mv ${vyy}${vmm}${vdd}${vhh}${vmn}.dat ${wkdir}/plev
it=`expr ${it} + 1`
done
fi
# ------------------------------
#    STEP 410: CONVERT
# ------------------------------ 
step=410
if test ${step} -ge ${start_step} -a ${step} -le ${end_step}
then
cd ${wkdir}
it=3
while test ${it} -le ${ft}
do
diff_hour=${it}
min=`expr ${diff_hour} \* 60`
date_edit ${syy} ${smm} ${sdd} ${shh} ${smn} ${min} > vdate.txt
read vyy vmm vdd vhh vmn < vdate.txt
rm -f vdate.txt
ip=`expr ${it} - 2`
if test ${ip} -lt 10
then
ip="0${ip}"
fi
# --- set file date ---
if [ ${num_frames} -eq 1 ]; then
fyy=${vyy}
fmm=${vmm}
fdd=${vdd}
fhh=${vhh}
fmn=${vmn}
else
echo "xxx  JOB:${step} : Not work when num_frames = ${num_frames}  xxx"
continue
fi
# --- remove files first ---
rm -f input.nc const.grd guess.grd
# --- input file ---
ln -s wrfout_d01_${fyy}-${fmm}-${fdd}_${fhh}:${fmn}:00	input.nc
# --- output files ---		
#ln -s const.grd const.grd
ln -s gs${ip}${MEM}.grd guess.grd
# --- run the pgm ---
$TMPDIR/ncio/nc2grd > @log.convert_${fyy}${fmm}${fdd}${fhh} 2> /dev/null
# --- end transaction ---
rm -f input.nc guess.grd #const.grd 
it=`expr ${it} + 1`
done
fi
