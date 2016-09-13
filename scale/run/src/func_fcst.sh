#!/bin/bash
#===============================================================================
#
#  Steps of 'fcst.sh'
#  October 2014, created   Guo-Yuan Lien
#
#===============================================================================

setting () {
#-------------------------------------------------------------------------------
# define steps

nsteps=3
stepname[1]='Run SCALE pp'
stepexecdir[1]="$TMPRUN/scale_pp"
stepexecname[1]="scale-rm_pp_ens"
stepname[2]='Run SCALE init'
stepexecdir[2]="$TMPRUN/scale_init"
stepexecname[2]="scale-rm_init_ens"
stepname[3]='Run ensemble forecasts'
stepexecdir[3]="$TMPRUN/scale"
stepexecname[3]="scale-rm_ens"
#stepname[4]='Run verification'
#stepexecdir[4]="$TMPRUN/verify"
#stepexecname[4]="verify"

#-------------------------------------------------------------------------------
# usage help string

USAGE="
[$myname] Run ensemble forecasts and (optional) verifications.

Configuration files:
  config.main
  config.cycle

Steps:
$(for i in $(seq $nsteps); do echo "  ${i}. ${stepname[$i]}"; done)

Usage: $myname [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP TIME_LIMIT]

  STIME       Time of the first cycle (format: YYYY[MMDDHHMMSS])
  ETIME       Time of the last  cycle (format: YYYY[MMDDHHMMSS])
               (default: same as STIME)
  MEMBERS     List of forecast members ('mean' for ensemble mean)
               all:     Run all members including ensemble mean (default)
               mems:    Run all members but not including ensemble mean
               '2 4 6': Run members 2, 4, 6
  CYCLE       Number of forecast cycles run in parallel
               (default: 1)
  CYCLE_SKIP  Run forecasts every ? cycles
               (default: 1)
  IF_VERF     Run verification? [Not finished!]
               0: No (default)
               1: Yes
              * to run the verification, a shared disk storing observations
                and reference model analyses needs to be used
  IF_EFSO     Use EFSO forecast length and output interval? [Not finished!]
               0: No (default)
               1: Yes
  ISTEP       The initial step in the first cycle from which this script starts
               (default: the first step)
  FSTEP       The final step in the last cycle by which this script ends
               (default: the last step)
  TIME_LIMIT  Requested time limit (only used when using a job scheduler)
               (default: 30 minutes)
"

#if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
#  echo "$USAGE"
#  exit 0
#fi

#-------------------------------------------------------------------------------
# set parameters from command line

STIME=${1:-$STIME}; shift
ETIME=${1:-$ETIME}; shift
MEMBERS=${1:-$MEMBERS}; shift
CYCLE=${1:-$CYCLE}; shift
CYCLE_SKIP=${1:-$CYCLE_SKIP}; shift
IF_VERF=${1:-$IF_VERF}; shift
IF_EFSO=${1:-$IF_EFSO}; shift
ISTEP=${1:-$ISTEP}; shift
FSTEP=${1:-$FSTEP}; shift
TIME_LIMIT="${1:-$TIME_LIMIT}"

#-------------------------------------------------------------------------------
# if some necessary parameters are not given, print the usage help and exit

#if [ -z "$STIME" ]; then
#  echo "$USAGE" >&2
#  exit 1
#fi

#-------------------------------------------------------------------------------
# error detection

#if ((MACHINE_TYPE == 10 && ONLINE_STGOUT != 0)); then
#  echo "[Error] $myname: When \$MACHINE_TYPE = 10, \$ONLINE_STGOUT needs to be 0." >&2
#  exit 1
#fi

#... more detections...

#-------------------------------------------------------------------------------
# assign default values to and standardize the parameters

STIME=$(datetime $STIME)
ETIME=$(datetime ${ETIME:-$STIME})
if [ -z "$MEMBERS" ] || [ "$MEMBERS" = 'all' ]; then
  MEMBERS="mean $(printf "$MEMBER_FMT " $(seq $MEMBER))"
elif [ "$MEMBERS" = 'mems' ]; then
  MEMBERS=$(printf "$MEMBER_FMT " $(seq $MEMBER))
else
  tmpstr=''
  for m in $MEMBERS; do
    if [ "$m" = 'mean' ] || [ "$m" = 'sprd' ]; then
      tmpstr="$tmpstr$m "
    else
      tmpstr="$tmpstr$(printf $MEMBER_FMT $((10#$m))) "
      (($? != 0)) && exit 1
    fi
  done
  MEMBERS="$tmpstr"
fi
CYCLE=${CYCLE:-0}
CYCLE_SKIP=${CYCLE_SKIP:-1}
IF_VERF=${IF_VERF:-0}
IF_EFSO=${IF_EFSO:-0}
ISTEP=${ISTEP:-1}
FSTEP=${FSTEP:-$nsteps}
TIME_LIMIT=${TIME_LIMIT:-"0:30:00"}

#-------------------------------------------------------------------------------
# common variables

if ((BDY_FORMAT >= 1)); then
  if ((BDYCYCLE_INT % BDYINT != 0)); then
    echo "[Error] \$BDYCYCLE_INT needs to be an exact multiple of \$BDYINT" >&2
    exit 1
  fi
  BDY_STARTFRAME_MAX=$((BDYCYCLE_INT / BDYINT))
  if [ -z "$PARENT_REF_TIME" ]; then
    PARENT_REF_TIME=$STIME
    for bdy_startframe in $(seq $BDY_STARTFRAME_MAX); do
      if ((BDY_FORMAT == 1)) && [ -s "$DATA_BDY_SCALE/${PARENT_REF_TIME}/hist/meanf/history.pe000000.nc" ]; then
        break
      elif ((BDY_FORMAT == 2 && BDY_ROTATING == 1)) && [ -s "$DATA_BDY_WRF/${PARENT_REF_TIME}/mean/wrfout_${PARENT_REF_TIME}" ]; then
        break
      elif ((BDY_FORMAT == 2 && BDY_ROTATING != 1)) && [ -s "$DATA_BDY_WRF/mean/wrfout_${PARENT_REF_TIME}" ]; then
        break
      elif ((bdy_startframe == BDY_STARTFRAME_MAX)); then
        echo "[Error] Cannot find boundary files." >&2
        exit 1
      fi
      PARENT_REF_TIME=$(datetime $PARENT_REF_TIME -${BDYINT} s)
    done
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

print_setting () {
#-------------------------------------------------------------------------------

for vname in DIR OUTDIR DATA_TOPO DATA_TOPO_BDY_SCALE DATA_LANDUSE DATA_BDY_SCALE \
             DATA_BDY_SCALE_PREP DATA_BDY_WRF DATA_BDY_NICAM OBS OBSNCEP TOPO_FORMAT \
             LANDUSE_FORMAT LANDUSE_UPDATE BDY_FORMAT BDY_ENS BDYINT BDYCYCLE_INT \
             PARENT_REF_TIME OCEAN_INPUT OCEAN_FORMAT OBSNUM WINDOW_S WINDOW_E \
             LCYCLE LTIMESLOT MEMBER NNODES PPN THREADS SCALE_NP \
             STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP \
             FCSTLEN FCSTOUT MAKEINIT OUT_OPT TOPOOUT_OPT LANDUSEOUT_OPT BDYOUT_OPT \
             LOG_OPT LOG_TYPE; do
  printf '  %-20s = %s\n' $vname "${!vname}"
done

#-------------------------------------------------------------------------------
}

#===============================================================================

staging_list () {
#-------------------------------------------------------------------------------
# TMPDAT

if ((TMPDAT_MODE == 1)); then
#-------------------
  echo "[Error] \$TMPDAT_MODE == 1 not available in this version!" >&2
  exit 1
#  safe_init_tmpdir $TMPDAT
#  safe_init_tmpdir $TMPDAT/exec
##  ln -fs $MODELDIR/scale-rm_pp $TMPDAT/exec
##  ln -fs $MODELDIR/scale-rm_init $TMPDAT/exec
##  ln -fs $MODELDIR/scale-rm $TMPDAT/exec
#  ln -fs $ENSMODEL_DIR/scale-rm_pp_ens $TMPDAT/exec
#  ln -fs $ENSMODEL_DIR/scale-rm_init_ens $TMPDAT/exec
#  ln -fs $ENSMODEL_DIR/scale-rm_ens $TMPDAT/exec
#  ln -fs $COMMON_DIR/pdbash $TMPDAT/exec
#  ln -fs $DATADIR/rad $TMPDAT/rad
#  ln -fs $DATADIR/land $TMPDAT/land
#  ln -fs $DATADIR/topo $TMPDAT
#  ln -fs $DATADIR/landuse $TMPDAT

#  if ((DATA_BDY_TMPLOC == 1)); then
#    if ((BDY_FORMAT == 2)); then
#      ln -fs $DATA_BDY_WRF $TMPDAT/bdyorg
#    fi
#  fi

#  safe_init_tmpdir $TMPDAT/conf
#  ln -fs $SCRP_DIR/config.* $TMPDAT/conf
#-------------------
else
#-------------------
  cat >> $STAGING_DIR/stagein.dat << EOF
${ENSMODEL_DIR}/scale-rm_pp_ens|exec/scale-rm_pp_ens
${ENSMODEL_DIR}/scale-rm_init_ens|exec/scale-rm_init_ens
${ENSMODEL_DIR}/scale-rm_ens|exec/scale-rm_ens
${COMMON_DIR}/pdbash|exec/pdbash
${SCRP_DIR}/config.nml.scale_pp|conf/config.nml.scale_pp
${SCRP_DIR}/config.nml.scale_init|conf/config.nml.scale_init
${SCRP_DIR}/config.nml.scale|conf/config.nml.scale
${SCRP_DIR}/config.nml.ensmodel|conf/config.nml.ensmodel
${DATADIR}/rad|rad
${DATADIR}/land|land
EOF
#${MODELDIR}/scale-rm_pp|exec/scale-rm_pp
#${MODELDIR}/scale-rm_init|exec/scale-rm_init
#${MODELDIR}/scale-rm|exec/scale-rm

  if [ "$TOPO_FORMAT" != 'prep' ]; then
    if ((DISK_MODE_TOPO_LANDUSE_DB == 2)); then
      echo "${DATADIR}/topo/${TOPO_FORMAT}/Products|topo/${TOPO_FORMAT}/Products|s" >> $STAGING_DIR/stagein.dat
    else
      echo "${DATADIR}/topo/${TOPO_FORMAT}/Products|topo/${TOPO_FORMAT}/Products" >> $STAGING_DIR/stagein.dat
    fi
  fi
  if [ "$LANDUSE_FORMAT" != 'prep' ]; then
    if ((DISK_MODE_TOPO_LANDUSE_DB == 2)); then
      echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products|landuse/${LANDUSE_FORMAT}/Products|s" >> $STAGING_DIR/stagein.dat
    else
      echo "${DATADIR}/landuse/${LANDUSE_FORMAT}/Products|landuse/${LANDUSE_FORMAT}/Products" >> $STAGING_DIR/stagein.dat
    fi
  fi

  if [ "$STG_TYPE" = 'K' ] || [ "$STG_TYPE" = 'K_rankdir' ]; then
    echo "${COMMON_DIR}/datetime|exec/datetime" >> $STAGING_DIR/stagein.dat
  fi
#-------------------
fi

#-------------------------------------------------------------------------------
# TMPOUT

if ((TMPOUT_MODE == 1)); then
#-------------------
  echo "[Error] \$TMPOUT_MODE == 1 not available in this version!" >&2
  exit 1
#  mkdir -p $(dirname $TMPOUT)
#  ln -fs $OUTDIR $TMPOUT

#  lcycles=$((LCYCLE * CYCLE_SKIP))
#  time=$STIME
#  while ((time <= ETIME)); do
#    for c in $(seq $CYCLE); do
#      time2=$(datetime $time $((lcycles * (c-1))) s)
#      if ((time2 <= ETIME)); then
#        #-------------------
#        if [ "$TOPO_FORMAT" = 'prep' ]; then
#          ln -fs ${DATA_TOPO} $TMPOUT/${time2}/topo
#        fi
#        #-------------------
#        if [ "$LANDUSE_FORMAT" = 'prep' ]; then
#          if ((LANDUSE_UPDATE == 1)); then
#            ln -fs ${DATA_LANDUSE}/${time2} $TMPOUT/${time2}/landuse
#          else
#            ln -fs ${DATA_LANDUSE} $TMPOUT/${time2}/landuse
#          fi
#        fi
#        #-------------------
#        if ((BDY_FORMAT == 0)); then
#          ln -fs ${DATA_BDY_SCALE_PREP}/${time2} $TMPOUT/${time2}/bdy
#        fi
#        #-------------------
#      fi
#    done
#    time=$(datetime $time $((lcycles * CYCLE)) s)
#  done

#  if ((DATA_BDY_TMPLOC == 2)); then
#    if ((BDY_FORMAT == 2)); then
#      ln -fs $DATA_BDY_WRF $TMPOUT/bdyorg
#    fi
#  fi

#  if ((BDY_FORMAT == 1)) || ((BDY_FORMAT == -1)); then
#    if ((DATA_BDY_TMPLOC == 1)); then
#      bdyorgf="$TMPDAT/bdyorg"
#    elif ((DATA_BDY_TMPLOC == 2)); then
#      bdyorgf="$TMPOUT/bdyorg"
#    fi
#    mkdir -p $bdyorgf

#    find_catalogue=0
#    for ibdy in $(seq $nfiles_all); do
#      time_bdy=${history_times_all[$ibdy]}

#      if ((find_catalogue == 0)); then
#        time_catalogue=$(datetime $time_bdy -$BDYCYCLE_INT s)
#        if [ -s "$DATA_BDY_SCALE/${time_catalogue}/log/scale/latlon_domain_catalogue.txt" ]; then
#          pathin="$DATA_BDY_SCALE/${time_catalogue}/log/scale/latlon_domain_catalogue.txt"
#          ln -fs ${pathin} ${bdyorgf}/latlon_domain_catalogue.txt
#          find_catalogue=1
#        fi
#      fi

#      if ((BDY_ENS == 1)); then
#        for m in $(seq $fmember); do
#          mem=${name_m[$m]}
#          [ "$mem" = 'mean' ] && mem='meanf'
#          mkdir -p ${bdyorgf}/${time_bdy}/${name_m[$m]}
#          for ifile in $(ls $DATA_BDY_SCALE/${time_bdy}/gues/${mem}/history.*.nc 2> /dev/null); do
#            pathin="$ifile"
#            ln -fs ${pathin} ${bdyorgf}/${time_bdy}/${name_m[$m]}/$(basename $ifile)
#          done
#        done
#      else
#        mkdir -p ${bdyorgf}/${time_bdy}/mean
#        for ifile in $(ls $DATA_BDY_SCALE/${time_bdy}/gues/meanf/history.*.nc 2> /dev/null); do
#          pathin="$ifile"
#          ln -fs ${pathin} ${bdyorgf}/${time_bdy}/mean/$(basename $ifile)
#        done
#      fi
#    done

#    if ((find_catalogue == 0)); then
#      echo "[Error] Cannot find a lat/lon domain catalogue file." >&2
#      exit 1
#    fi
#  fi
#-------------------
else
#-------------------
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  loop=0
  while ((time <= ETIME)); do
    loop=$((loop+1))
    if ((ONLINE_STGOUT == 1)); then
      stgoutstep="stageout.loop.${loop}"
    else
      stgoutstep='stageout.out'
    fi

    for c in $(seq $CYCLE); do
      time2=$(datetime $time $((lcycles * (c-1))) s)
      if ((time2 <= ETIME)); then
        #-------------------
        # stage-in
        #-------------------

        # anal
        #-------------------
        if ((MAKEINIT != 1)); then
          for m in $(seq $fmember); do
            mm=$(((c-1) * fmember + m))
            for q in $(seq $mem_np); do
              path="${time2}/anal/${name_m[$mm]}/init$(printf $SCALE_SFX $((q-1)))"
              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((mm-1)*mem_np+q))]}
            done
          done
        fi

        # anal_ocean
        #-------------------
#        if ((OCEAN_INPUT == 1)) && ((OCEAN_FORMAT == 0)); then
#          for m in $(seq $fmember); do
#            mm=$(((c-1) * fmember + m))
#            for q in $(seq $mem_np); do
#              path="${time2}/anal/${name_m[$mm]}/init_ocean$(printf $SCALE_SFX $((q-1)))"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((mm-1)*mem_np+q))]}
#            done
#          done
#        fi

        # topo
        #-------------------
        if ((loop == 1)) && [ "$TOPO_FORMAT" = 'prep' ]; then
          for m in $(seq $fmember); do
            mm=$(((c-1) * fmember + m))
            for q in $(seq $mem_np); do
              pathin="${DATA_TOPO}/const/topo/topo$(printf $SCALE_SFX $((q-1)))"
              path="const/topo/topo$(printf $SCALE_SFX $((q-1)))"
              echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((mm-1)*mem_np+q))]}
            done
          done
        fi

        # topo (bdy_scale)
        #-------------------
        if ((loop == 1 && BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
#          for ifile in $(ls ${DATA_TOPO_BDY_SCALE}/topo.*.nc 2> /dev/null); do
#            pathin="$ifile"
#            path="bdytopo/const/$(basename $ifile)"
#            if ((DISK_MODE_DATA_BDY == 2)); then
#              echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
#            else
#              echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
#            fi
#          done
          pathin="${DATA_TOPO_BDY_SCALE}"
          path="bdytopo/const"
          if ((DISK_MODE_DATA_BDY == 2)); then
            echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
          else
            echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
          fi
        fi

        # landuse
        #-------------------
        if ((loop == 1 || LANDUSE_UPDATE == 1)) && [ "$LANDUSE_FORMAT" = 'prep' ]; then
          for m in $(seq $fmember); do
            mm=$(((c-1) * fmember + m))
            for q in $(seq $mem_np); do
              if ((LANDUSE_UPDATE == 1)); then
                pathin="${DATA_LANDUSE}/${time2}/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
                path="${time2}/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
              else
                pathin="${DATA_LANDUSE}/const/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
                path="const/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
              fi
              echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((mm-1)*mem_np+q))]}
            done
          done
        fi

        # bdy (prepared)
        #-------------------
        if ((BDY_FORMAT == 0)); then
          if ((BDY_ENS == 0)); then
            for m in $(seq $fmember); do
              mm=$(((c-1) * fmember + m))
              for q in $(seq $mem_np); do
                pathin="${DATA_BDY_SCALE_PREP}/${time2}/bdy/mean/boundary$(printf $SCALE_SFX $((q-1)))"
                path="${time2}/bdy/mean/boundary$(printf $SCALE_SFX $((q-1)))"
                echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((mm-1)*mem_np+q))]}
              done
            done
          elif ((BDY_ENS == 1)); then
            for m in $(seq $fmember); do
              mm=$(((c-1) * fmember + m))
              for q in $(seq $mem_np); do
                pathin="${DATA_BDY_SCALE_PREP}/${time2}/bdy/${name_m[$m]}/boundary$(printf $SCALE_SFX $((q-1)))"
                path="${time2}/bdy/${name_m[$m]}/boundary$(printf $SCALE_SFX $((q-1)))"
                echo "${pathin}|${path}" >> $STAGING_DIR/stagein.out.${mem2node[$(((mm-1)*mem_np+q))]}
              done
            done
          fi
        fi

        #-------------------
        # stage-out
        #-------------------

#        #++++++
#        if ((SIMPLE_STGOUT == 1)); then
#        #++++++

          # anal
          #-------------------
          if ((MAKEINIT == 1)); then
            path="${time2}/anal"
            echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
          fi

          # topo
          #-------------------
          if ((loop == 1 && TOPOOUT_OPT <= 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
            path="const/topo"
            echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
          fi

          # landuse
          #-------------------
          if ((loop == 1 || LANDUSE_UPDATE == 1)) && ((LANDUSEOUT_OPT <= 1)) && [ "$LANDUSE_FORMAT" != 'prep' ]; then
            if ((LANDUSE_UPDATE == 1)); then
              path="${time2}/landuse"
            else
              path="const/landuse"
            fi
            echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
          fi

          # bdy
          #-------------------
          if ((BDY_FORMAT != 0)); then
            if ((BDY_ENS == 1 && BDYOUT_OPT <= 1)); then
#              for m in $(seq $fmember); do
##                mm=$(((c-1) * fmember + m))
#                path="${time2}/bdy/${name_m[$m]}"
#                echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
#              done
              path="${time2}/bdy"
              echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
            elif ((BDYOUT_OPT <= 2)); then
              path="${time2}/bdy/mean"
              echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
            fi
          fi

          # fcst
          #-------------------
#          for m in $(seq $fmember); do
##            mm=$(((c-1) * fmember + m))
#            path="${time2}/fcst/${name_m[$m]}"
#            echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
#          done
          path="${time2}/fcst"
          echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}

          ### anal_ocean [mean]

          # log
          #-------------------
          if [ "$MPI_TYPE" = 'K' ]; then
            log_zeros='0'
          else
            log_zeros='000000'
          fi

          if ((LOG_OPT <= 2)); then
            if ((LOG_TYPE == 1)); then
              if ((c == 1)); then
                path="${time2}/log/fcst_scale_pp/${name_m[1]}_pp.conf"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale_pp/${name_m[1]}_LOG.pe000000"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale_pp/NOUT.${log_zeros}"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale_init/${name_m[1]}_init.conf"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale_init/${name_m[1]}_LOG.pe000000"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                if ((BDY_ENS == 1)); then
                  path="${time2}/log/fcst_scale_init/NOUT-1.${log_zeros}"
                else
                  path="${time2}/log/fcst_scale_init/NOUT.${log_zeros}"
                fi
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
              fi
            else
              path="${time2}/log/fcst_scale_pp"
              echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
              path="${time2}/log/fcst_scale_init"
              echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
            fi
          fi
          if ((LOG_OPT <= 3)); then
            if ((LOG_TYPE == 1)); then
              if ((c == 1)); then
                path="${time2}/log/fcst_scale/${name_m[1]}_run.conf"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale/${name_m[1]}_LOG.pe000000"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale/NOUT-1.${log_zeros}"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
                path="${time2}/log/fcst_scale/latlon_domain_catalogue.txt"
                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.1
              fi
            else
              path="${time2}/log/fcst_scale"
              echo "${OUTDIR}/${path}|${path}|d" >> $STAGING_DIR/${stgoutstep}
            fi
          fi

#        #++++++
#        else
#        #++++++
#          for m in $(seq $fmember); do
#            mm=$(((c-1) * fmember + m))
#            #-------------------

#            for q in $(seq $mem_np); do
#              #-------------------

#              # bdy [members]
#              #-------------------
#              if ((BDYOUT_OPT <= 1)) && ((BDY_ENS == 1)); then
#                path="${time2}/bdy/${name_m[$mm]}/boundary$(printf $SCALE_SFX $((q-1)))"
#                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+q))]}
#              fi

#              # anal
#              #-------------------
#              if ((MAKEINIT == 1)); then
#                path="${time2}/anal/${name_m[$mm]}/init$(printf $SCALE_SFX $((q-1)))"
#                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+q))]}
#              fi

#              # anal_ocean
#              #-------------------
#  #            if ((OCEAN_INPUT == 1)) && ((MAKEINIT != 1)); then
#  #              path="${time2}/anal/${name_m[$mm]}/init_ocean$(printf $SCALE_SFX $((q-1)))"
#  #              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+q))]}
#  #            fi

#              # fcst [history]
#              #-------------------
#              if ((OUT_OPT <= 2)); then
#                path="${time2}/fcst/${name_m[$mm]}/history$(printf $SCALE_SFX $((q-1)))"
#                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+q))]}
#              fi

#              # fcst [restart]
#              #-------------------
#              if ((OUT_OPT <= 1)); then
#                path="${time2}/fcst/${name_m[$mm]}/init_$(datetime ${time2} $FCSTLEN s)$(printf $SCALE_SFX $((q-1)))"
#                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+q))]}
#              fi

#              #-------------------
#            done

#            # log [scale_init: members]
#            #-------------------
#            if ((BDY_FORMAT > 0)) && ((LOG_OPT <= 2)) && ((BDY_ENS == 1)); then
#              path="${time2}/log/scale_init/${name_m[$mm]}_fcst_LOG${SCALE_LOG_SFX}"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+1))]}
#            fi

#            # log [scale]
#            #-------------------
#            if ((LOG_OPT <= 3)); then
#              path="${time2}/log/scale/${name_m[$mm]}_fcst_LOG${SCALE_LOG_SFX}"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+1))]}
#              path="${time2}/log/scale/${name_m[$mm]}_latlon_domain_catalogue.txt"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$(((mm-1)*mem_np+1))]}
#            fi

#  #          if ((LOG_OPT <= 1)); then
#  #            # perturb bdy log
#  #          fi

#            #-------------------
#          done
#          #-------------------

#          if ((repeat_mems <= fmember)); then
#            tmpidx=0                            # mm=1
#          else
#            tmpidx=$((((c-1)*fmember)*mem_np))  # mm=$(((c-1) * fmember + 1))
#          fi

#          for q in $(seq $mem_np); do
#            #-------------------

#            # topo
#            #-------------------
#            if ((TOPOOUT_OPT <= 1)); then
#              path="${time2}/topo/topo$(printf $SCALE_SFX $((q-1)))"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#            fi

#            # landuse
#            #-------------------
#            if ((LANDUSEOUT_OPT <= 1)); then
#              path="${time2}/landuse/landuse$(printf $SCALE_SFX $((q-1)))"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#            fi

#            # bdy [mean]
#            #-------------------
#            if ((BDYOUT_OPT <= 2)) && ((BDY_ENS != 1)); then
#              path="${time2}/bdy/mean/boundary$(printf $SCALE_SFX $((q-1)))"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#            fi

#            # anal_ocean [mean]
#            #-------------------
#            if ((OCEAN_INPUT == 1)) && ((MAKEINIT != 1)); then
#              path="${time2}/anal/mean/init_ocean$(printf $SCALE_SFX $((q-1)))"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#            fi

#            # log [scale_pp/scale_init/scale]
#            #-------------------
#            if ((LOG_OPT <= 4)); then
#              if [ "$TOPO_FORMAT" != 'prep' ] || [ "$LANDUSE_FORMAT" != 'prep' ] && ((BDY_FORMAT != 0)); then
#                path="${time2}/log/scale_pp/NOUT-$(printf $PROCESS_FMT $((tmpidx+q-1)))"
#                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#              fi
#              if ((BDY_FORMAT != 0)); then
#                path="${time2}/log/scale_init/NOUT-$(printf $PROCESS_FMT $((tmpidx+q-1)))"
#                echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#              fi
#              path="${time2}/log/scale/NOUT-$(printf $PROCESS_FMT $((tmpidx+q-1)))"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+q))]}
#            fi

#            #-------------------
#          done

#          # log [scale_pp]
#          #-------------------
#          if [ "$TOPO_FORMAT" != 'prep' ] || [ "$LANDUSE_FORMAT" != 'prep' ]; then
#            if ((LOG_OPT <= 2)); then
#              path="${time2}/log/scale_pp/fcst_LOG${SCALE_LOG_SFX}"
#              echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+1))]}
#            fi
#          fi

#          # log [scale_init: mean]
#          #-------------------
#          if ((BDY_FORMAT > 0)) && ((LOG_OPT <= 2)) && ((BDY_ENS != 1)); then
#            path="${time2}/log/scale_init/mean_fcst_LOG${SCALE_LOG_SFX}"
#            echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+1))]}
#          fi

#          # log [scale: catalogue]
#          #-------------------
#          path="${time2}/log/scale/latlon_domain_catalogue.txt"
#          echo "${OUTDIR}/${path}|${path}" >> $STAGING_DIR/${stgoutstep}.${mem2node[$((tmpidx+1))]}

#        #++++++
#        fi # ((SIMPLE_STGOUT == 1))
#        #++++++

        #-------------------
      fi
    done

    time=$(datetime $time $((lcycles * CYCLE)) s)
  done

  #-------------------
  # stage-in
  #-------------------

  # bdy
  #-------------------
  if ((BDY_FORMAT >= 1)); then
    if ((BDY_FORMAT == 1)); then
      if [ -s "$DATA_BDY_SCALE/${PARENT_REF_TIME}/log/scale/latlon_domain_catalogue.txt" ]; then
        pathin="$DATA_BDY_SCALE/${PARENT_REF_TIME}/log/scale/latlon_domain_catalogue.txt"
        path="bdyorg/latlon_domain_catalogue.txt"
        if ((DISK_MODE_DATA_BDY == 2)); then
          echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
        else
          echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
        fi
      else
        echo "[Error] Cannot find a lat/lon domain catalogue file at" >&2
        echo "        '$DATA_BDY_SCALE/${PARENT_REF_TIME}/log/scale/latlon_domain_catalogue.txt'" >&2
        exit 1
      fi
    fi

    nbdy_all=0
    time=$STIME
    while ((time <= ETIME)); do
      for c in $(seq $CYCLE); do
        time2=$(datetime $time $((lcycles * (c-1))) s)
        if ((time2 <= ETIME)); then

          if ((BDY_FORMAT == 1 && BDY_ROTATING == 1)); then
            bdy_setting $time2 $FCSTLEN - $BDYINT $PARENT_REF_TIME
          else
            bdy_setting $time2 $FCSTLEN $BDYCYCLE_INT $BDYINT $PARENT_REF_TIME
          fi

          for ibdy in $(seq $nbdy); do
            time_bdy=${bdy_times[$ibdy]}

            bdy_processed=0
            for ibdy2 in $(seq $nbdy_all); do
              if ((${bdy_times_all[$ibdy2]} == $time_bdy)); then
                bdy_processed=1
                break
              fi
            done

            if ((bdy_processed == 0)); then
              nbdy_all=$((nbdy_all+1))
              bdy_times_all[${nbdy_all}]=$time_bdy
            fi

            if ((bdy_processed == 0 || BDY_ROTATING == 1)); then
              if ((BDY_FORMAT == 1)); then

                if ((BDY_ENS == 1)); then
                  for m in $(seq $fmember); do
                    mem=${name_m[$m]}
                    if [ "$BDY_SCALE_DIR" = 'hist' ] && [ "$mem" = 'mean' ]; then
                      mem='meanf'
                    fi
#                    for ifile in $(ls $DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/${mem}/history.*.nc 2> /dev/null); do
#                      pathin="$ifile"
#                      if ((BDY_ROTATING == 1)); then
#                        path="bdyorg/${time_bdy}/${name_m[$m]}/${time_bdy}/$(basename $ifile)"
#                      else
#                        path="bdyorg/const/${name_m[$m]}/${time_bdy}/$(basename $ifile)"
#                      fi
#                      if ((DISK_MODE_DATA_BDY == 2)); then
#                        echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
#                      else
#                        echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
#                      fi
#                    done
                    pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/${mem}"
                    if ((BDY_ROTATING == 1)); then
                      path="bdyorg/${time_bdy}/${name_m[$m]}/${time_bdy}"
                    else
                      path="bdyorg/const/${name_m[$m]}/${time_bdy}"
                    fi
                    if ((DISK_MODE_DATA_BDY == 2)); then
                      echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
                    else
                      echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
                    fi
                  done
                else
#                  for ifile in $(ls $DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/meanf/history.*.nc 2> /dev/null); do
#                    pathin="$ifile"
#                    if ((BDY_ROTATING == 1)); then
#                      path="bdyorg/${time_bdy}/mean/${time_bdy}/$(basename $ifile)"
#                    else
#                      path="bdyorg/const/mean/${time_bdy}/$(basename $ifile)"
#                    fi
#                    if ((DISK_MODE_DATA_BDY == 2)); then
#                      echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
#                    else
#                      echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
#                    fi
#                  done
                  if [ "$BDY_SCALE_DIR" = 'hist' ]; then
                    pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/meanf"
                  else
                    pathin="$DATA_BDY_SCALE/${time_bdy}/${BDY_SCALE_DIR}/mean"
                  fi
                  if ((BDY_ROTATING == 1)); then
                    path="bdyorg/${time_bdy}/mean/${time_bdy}"
                  else
                    path="bdyorg/const/mean/${time_bdy}"
                  fi
                  if ((DISK_MODE_DATA_BDY == 2)); then
                    echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
                  else
                    echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
                  fi
                fi

              elif ((BDY_FORMAT == 2)); then

                if ((BDY_ENS == 1)); then
                  for m in $(seq $fmember); do
                    if ((BDY_ROTATING == 1)); then
                      pathin="$DATA_BDY_WRF/${time2}/${name_m[$m]}/wrfout_${time_bdy}"
                      path="bdyorg/${time2}/${name_m[$m]}/wrfout_${time_bdy}"
                    else
                      pathin="$DATA_BDY_WRF/${name_m[$m]}/wrfout_${time_bdy}"
                      path="bdyorg/const/${name_m[$m]}/wrfout_${time_bdy}"
                    fi
                    if ((DISK_MODE_DATA_BDY == 2)); then
                      echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
                    else
                      echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
                    fi
                  done
                else
                  if ((BDY_ROTATING == 1)); then
                    pathin="$DATA_BDY_WRF/${time2}/mean/wrfout_${time_bdy}"
                    path="bdyorg/${time2}/mean/wrfout_${time_bdy}"
                  else
                    pathin="$DATA_BDY_WRF/mean/wrfout_${time_bdy}"
                    path="bdyorg/const/mean/wrfout_${time_bdy}"
                  fi
                  if ((DISK_MODE_DATA_BDY == 2)); then
                    echo "${pathin}|${path}|s" >> $STAGING_DIR/stagein.dat
                  else
                    echo "${pathin}|${path}" >> $STAGING_DIR/stagein.dat
                  fi
                fi

              fi
            fi # ((bdy_processed == 0 || BDY_ROTATING == 1))
          done # [ ibdy in $(seq $nbdy) ]

        fi # ((time2 <= ETIME))
      done
      time=$(datetime $time $((lcycles * CYCLE)) s)
    done
  fi # ((BDY_FORMAT >= 1))

  #-------------------

#-------------------
fi

### EFSO outputs...

#-------------------------------------------------------------------------------
}

#===============================================================================

enspp_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[1]}: Pre-processing script start" >&2
fi

if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
  echo "  ... skip this step (use prepared topo and landuse files)"
  exit 1
elif ((BDY_FORMAT == 0)); then
  echo "  ... skip this step (use prepared boundaries)"
  exit 1
elif ((LANDUSE_UPDATE != 1 && loop > 1)); then
  echo "  ... skip this step (already done in the first cycle)"
  exit 1
fi

if ((BDY_FORMAT == 1)); then
  if ((DISK_MODE_DATA_BDY == 2)); then
    bdycatalogue=${TMPDAT_S}/bdyorg/latlon_domain_catalogue.txt
    bdytopo=${TMPDAT_S}/bdytopo/const/topo
  else
    bdycatalogue=${TMPDAT_L}/bdyorg/latlon_domain_catalogue.txt
    bdytopo=${TMPDAT_L}/bdytopo/const/topo
  fi
fi

if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=$rcycle
#  MEMBER_RUN=1
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= fmember ? $((repeat_mems*rcycle)) : $((fmember*rcycle))))
fi

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_scale_pp_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale_pp $MEMBER_RUN $iter
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[1]}: Pre-processing script end" >&2
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  m=$(((it-1)*parallel_mems+g))
  if ((m >= 1 && m <= MEMBER_RUN)); then
    if ((TMPRUN_MODE <= 2)); then
      c=$m
    else
      c=$((repeat_mems <= fmember ? $(((m-1)/repeat_mems+1)) : $(((m-1)/fmember+1))))
    fi

    if (pdrun $g $PROC_OPT); then
      bash $SCRP_DIR/src/pre_scale_pp.sh $MYRANK ${stimes[$c]} ${name_m[$m]} \
           $TMPRUN/scale_pp/$(printf '%04d' $m) $TMPDAT \
           fcst ${bdytopo} ${bdycatalogue}
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

enspp_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
  return 1
elif ((BDY_FORMAT == 0)); then
  return 1
fi

if ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=$rcycle
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= fmember ? $((repeat_mems*rcycle)) : $((fmember*rcycle))))
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  m=$(((it-1)*parallel_mems+g))
  if ((m >= 1 && m <= MEMBER_RUN)); then
    if ((TMPRUN_MODE <= 2)); then
      c=$m
    else
      c=$((repeat_mems <= fmember ? $(((m-1)/repeat_mems+1)) : $(((m-1)/fmember+1))))
    fi

    if (pdrun $g $PROC_OPT); then
      bash $SCRP_DIR/src/post_scale_pp.sh $MYRANK ${stimes[$c]} \
           ${name_m[$m]} $TMPRUN/scale_pp/$(printf '%04d' $m) $LOG_OPT fcst
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[1]}: $it: Post-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensinit_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[2]}: Pre-processing script start" >&2
fi

if ((BDY_FORMAT == 0)); then
  echo "  ... skip this step (use prepared boundaries)"
  exit 1
fi

if ((DISK_MODE_DATA_BDY == 2)); then
  bdyorgf=${TMPDAT_S}/bdyorg
else
  bdyorgf=${TMPDAT_L}/bdyorg
fi

if ((BDY_ENS == 1)); then
  MEMBER_RUN=$((fmember*rcycle))
elif ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=$rcycle
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= fmember ? $((repeat_mems*rcycle)) : $((fmember*rcycle))))
fi

mkinit=0
if ((loop == 1)); then
  mkinit=$MAKEINIT
fi

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_scale_init_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale_init $MEMBER_RUN $iter
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[2]}: Pre-processing script end" >&2
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  m=$(((it-1)*parallel_mems+g))
  if ((m >= 1 && m <= MEMBER_RUN)); then
    if ((BDY_ENS == 1)); then
      c=$(((m-1)/fmember+1))
      mem_bdy=${name_m[$m]}
    elif ((TMPRUN_MODE <= 2)); then
      c=$m
      mem_bdy='mean'
    else
      c=$((repeat_mems <= fmember ? $(((m-1)/repeat_mems+1)) : $(((m-1)/fmember+1))))
      mem_bdy='mean'
    fi

    if ((BDY_FORMAT == 1 && BDY_ROTATING == 1)); then
      bdy_setting ${stimes[$c]} $FCSTLEN - $BDYINT $PARENT_REF_TIME
    else
      bdy_setting ${stimes[$c]} $FCSTLEN $BDYCYCLE_INT $BDYINT $PARENT_REF_TIME
    fi
    bdy_time_list=''
    for ibdy in $(seq $nbdy); do
      bdy_time_list="${bdy_time_list}${bdy_times[$ibdy]} "
    done

    if ((LANDUSE_UPDATE == 1)); then
      time_l=${stimes[$c]}
    else
      time_l='const'
    fi

    if (pdrun $g $PROC_OPT); then
      bash $SCRP_DIR/src/pre_scale_init.sh $MYRANK \
           $TMPOUT/const/topo/topo $TMPOUT/${time_l}/landuse/landuse \
           ${bdyorgf} ${stimes[$c]} $mkinit ${name_m[$m]} $mem_bdy \
           $TMPRUN/scale_init/$(printf '%04d' $m) \
           "$bdy_time_list" $ntsteps $ntsteps_skip fcst
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensinit_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

if ((BDY_FORMAT == 0)); then
  return 1
fi

if ((BDY_ENS == 1)); then
  MEMBER_RUN=$((fmember*rcycle))
elif ((TMPRUN_MODE <= 2)); then # shared run directory: only run one member per cycle
  MEMBER_RUN=$rcycle
else # local run directory: run multiple members as needed
  MEMBER_RUN=$((repeat_mems <= fmember ? $((repeat_mems*rcycle)) : $((fmember*rcycle))))
fi

mkinit=0
if ((loop == 1)); then
  mkinit=$MAKEINIT
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  m=$(((it-1)*parallel_mems+g))
  if ((m >= 1 && m <= MEMBER_RUN)); then
    if ((BDY_ENS == 1)); then
      c=$(((m-1)/fmember+1))
    elif ((TMPRUN_MODE <= 2)); then
      c=$m
    else
      c=$((repeat_mems <= fmember ? $(((m-1)/repeat_mems+1)) : $(((m-1)/fmember+1))))
    fi

    if (pdrun $g $PROC_OPT); then
      if ((BDY_ENS == 1)); then
        bash $SCRP_DIR/src/post_scale_init.sh $MYRANK ${stimes[$c]} \
             $mkinit ${name_m[$m]} $TMPRUN/scale_init/$(printf '%04d' $m) $LOG_OPT fcst
      else
        bash $SCRP_DIR/src/post_scale_init.sh $MYRANK ${stimes[$c]} \
             $mkinit mean $TMPRUN/scale_init/$(printf '%04d' $m) $LOG_OPT fcst
      fi
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[2]}: $it: Post-processing script (member) start" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensfcst_1 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Pre-processing scripts"
#echo

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[3]}: Pre-processing script start" >&2
fi
 
MEMBER_RUN=$((fmember*rcycle))

if (pdrun all $PROC_OPT); then
  bash $SCRP_DIR/src/pre_scale_node.sh $MYRANK \
       $mem_nodes $mem_np $TMPRUN/scale $MEMBER_RUN $iter
fi

mkinit=0
if ((loop == 1)); then
  mkinit=$MAKEINIT
fi

if ((MYRANK == 0)); then
  echo "[$(datetime_now)] ${time}: ${stepname[3]}: Pre-processing script end" >&2
fi

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Pre-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  m=$(((it-1)*parallel_mems+g))
  if ((m >= 1 && m <= fmembertot)); then
    c=$(((m-1)/fmember+1))

#    if ((PERTURB_BDY == 1)); then
#      ...
#    fi

    ocean_base='-'
    if ((OCEAN_INPUT == 1)); then
      if ((mkinit != 1 || OCEAN_FORMAT != 99)); then
        ocean_base="$TMPOUT/${stimes[$c]}/anal/mean/init_ocean"  ### always use mean???
      fi
    fi

    if ((BDY_ENS == 1)); then
      bdy_base="$TMPOUT/${stimes[$c]}/bdy/${name_m[$m]}/boundary"
    else
      bdy_base="$TMPOUT/${stimes[$c]}/bdy/mean/boundary"
    fi

    if ((BDY_FORMAT == 1 && BDY_ROTATING == 1)); then
      bdy_setting ${stimes[$c]} $FCSTLEN - $BDYINT $PARENT_REF_TIME
    else
      bdy_setting ${stimes[$c]} $FCSTLEN $BDYCYCLE_INT $BDYINT $PARENT_REF_TIME
    fi

    if ((LANDUSE_UPDATE == 1)); then
      time_l=${stimes[$c]}
    else
      time_l='const'
    fi

    if (pdrun $g $PROC_OPT); then
      bash $SCRP_DIR/src/pre_scale.sh $MYRANK ${name_m[$m]} \
           $TMPOUT/${stimes[$c]}/anal/${name_m[$m]}/init $ocean_base $bdy_base \
           $TMPOUT/const/topo/topo $TMPOUT/${time_l}/landuse/landuse \
           ${stimes[$c]} $FCSTLEN $FCSTLEN $FCSTOUT $TMPRUN/scale/$(printf '%04d' $m) \
           fcst $bdy_start_time
    fi
  fi

  if ((MYRANK == 0)); then
     echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Pre-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

ensfcst_2 () {
#-------------------------------------------------------------------------------

#echo
#echo "* Post-processing scripts"
#echo

for it in $(seq $its $ite); do
  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Post-processing script (member) start" >&2
  fi

  g=${proc2group[$((MYRANK+1))]}
  m=$(((it-1)*parallel_mems+g))
  if ((m >= 1 && m <= fmembertot)); then
    c=$(((m-1)/fmember+1))

#    if ((PERTURB_BDY == 1)); then
#      ...
#    fi

    if (pdrun $g $PROC_OPT); then
      bash $SCRP_DIR/src/post_scale.sh $MYRANK ${stimes[$c]} \
           ${name_m[$m]} $FCSTLEN $TMPRUN/scale/$(printf '%04d' $m) $LOG_OPT $OUT_OPT fcst
    fi
  fi

  if ((MYRANK == 0)); then
    echo "[$(datetime_now)] ${time}: ${stepname[3]}: $it: Post-processing script (member) end" >&2
  fi
done

#-------------------------------------------------------------------------------
}

#===============================================================================

finalization () {
#-------------------------------------------------------------------------------

if ((LOG_TYPE >= 3)); then
  lcycles=$((LCYCLE * CYCLE_SKIP))
  time=$STIME
  while ((time <= ETIME)); do
    for c in $(seq $CYCLE); do
      time2=$(datetime $time $((lcycles * (c-1))) s)

      if ((LOG_OPT <= 2)) && [ -d "$OUTDIR/${time2}/log/fcst_scale_pp" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C $OUTDIR/${time2}/log -cf $OUTDIR/${time2}/log/fcst_scale_pp.tar fcst_scale_pp && rm -fr $OUTDIR/${time2}/log/fcst_scale_pp ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C $OUTDIR/${time2}/log -czf $OUTDIR/${time2}/log/fcst_scale_pp.tar.gz fcst_scale_pp && rm -fr $OUTDIR/${time2}/log/fcst_scale_pp ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C $OUTDIR/${time2}/log -cf $OUTDIR/${time2}/log/fcst_scale_pp.tar fcst_scale_pp && rm -fr $OUTDIR/${time2}/log/fcst_scale_pp
          elif ((LOG_TYPE == 4)); then
            tar -C $OUTDIR/${time2}/log -czf $OUTDIR/${time2}/log/fcst_scale_pp.tar.gz fcst_scale_pp && rm -fr $OUTDIR/${time2}/log/fcst_scale_pp
          fi
        fi
      fi

      if ((LOG_OPT <= 2)) && [ -d "$OUTDIR/${time2}/log/fcst_scale_init" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C $OUTDIR/${time2}/log -cf $OUTDIR/${time2}/log/fcst_scale_init.tar fcst_scale_init && rm -fr $OUTDIR/${time2}/log/fcst_scale_init ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C $OUTDIR/${time2}/log -czf $OUTDIR/${time2}/log/fcst_scale_init.tar.gz fcst_scale_init && rm -fr $OUTDIR/${time2}/log/fcst_scale_init ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C $OUTDIR/${time2}/log -cf $OUTDIR/${time2}/log/fcst_scale_init.tar fcst_scale_init && rm -fr $OUTDIR/${time2}/log/fcst_scale_init
          elif ((LOG_TYPE == 4)); then
            tar -C $OUTDIR/${time2}/log -czf $OUTDIR/${time2}/log/fcst_scale_init.tar.gz fcst_scale_init && rm -fr $OUTDIR/${time2}/log/fcst_scale_init
          fi
        fi
      fi

      if ((LOG_OPT <= 3)) && [ -d "$OUTDIR/${time2}/log/fcst_scale" ]; then
        if ((TAR_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= TAR_THREAD)); do
            sleep 1s
          done
          if ((LOG_TYPE == 3)); then
            ( tar -C $OUTDIR/${time2}/log -cf $OUTDIR/${time2}/log/fcst_scale.tar fcst_scale && rm -fr $OUTDIR/${time2}/log/fcst_scale ) &
          elif ((LOG_TYPE == 4)); then
            ( tar -C $OUTDIR/${time2}/log -czf $OUTDIR/${time2}/log/fcst_scale.tar.gz fcst_scale && rm -fr $OUTDIR/${time2}/log/fcst_scale ) &
          fi
        else
          if ((LOG_TYPE == 3)); then
            tar -C $OUTDIR/${time2}/log -cf $OUTDIR/${time2}/log/fcst_scale.tar fcst_scale && rm -fr $OUTDIR/${time2}/log/fcst_scale
          elif ((LOG_TYPE == 4)); then
            tar -C $OUTDIR/${time2}/log -czf $OUTDIR/${time2}/log/fcst_scale.tar.gz fcst_scale && rm -fr $OUTDIR/${time2}/log/fcst_scale
          fi
        fi
      fi

    done
    time=$(datetime $time $((lcycles * CYCLE)) s)
  done
  if ((TAR_THREAD > 1)); then
    wait
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
