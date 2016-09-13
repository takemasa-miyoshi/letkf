#!/bin/bash
#===============================================================================
#
#  Run ensemble forecasts and (optional) verifications.
#
#  August  2014, modified from GFS-LETKF, Guo-Yuan Lien
#  October 2014, modified,                Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    fcst.sh [STIME ETIME MEMBERS CYCLE CYCLE_SKIP IF_VERF IF_EFSO ISTEP FSTEP]
#
#  Use settings:
#    config.main
#    config.fcst
#    config.nml.scale_pp_topo
#    config.nml.scale_pp_landuse
#    config.nml.scale_init
#    config.nml.scale
#
#===============================================================================

cd "$(dirname "$0")"
myname='fcst.sh'
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.main || exit $?
. config.$myname1 || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_$myname1.sh || exit $?

#-------------------------------------------------------------------------------

if [ "$STG_TYPE" = 'K_rankdir' ]; then
  SCRP_DIR="."
  if ((TMPDAT_MODE <= 2)); then
    TMPDAT="../dat"
  else
    TMPDAT="./dat"
  fi
  if ((TMPRUN_MODE <= 2)); then
    TMPRUN="../run"
  else
    TMPRUN="./run"
  fi
  if ((TMPOUT_MODE <= 2)); then
    TMPOUT="../out"
  else
    TMPOUT="./out"
  fi
fi

echo "[$(datetime_now)] Start $myname $@" >&2

setting "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" || exit $?

echo
print_setting || exit $?

#-------------------------------------------------------------------------------

if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
    safe_init_tmpdir $TMP || exit $?
  fi
  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
    safe_init_tmpdir $TMPL || exit $?
  fi
fi

#===============================================================================
# Determine the distibution schemes

declare -a node
declare -a node_m
declare -a name_m
declare -a mem2node
declare -a mem2proc
declare -a proc2node
declare -a proc2group
declare -a proc2grpproc

#if [ "$STG_TYPE" = 'builtin' ] && ((&& ISTEP == 1)); then
if [ "$STG_TYPE" = 'builtin' ]; then
  safe_init_tmpdir $NODEFILE_DIR || exit $?
  distribute_fcst "$MEMBERS" $CYCLE machinefile $NODEFILE_DIR || exit $?
else
  distribute_fcst "$MEMBERS" $CYCLE - - $NODEFILE_DIR/distr || exit $?
fi

if ((CYCLE == 0)); then
  CYCLE=$parallel_mems
fi

#===============================================================================
# Determine the staging list and then stage in

if [ "$STG_TYPE" = 'builtin' ]; then
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  safe_init_tmpdir $STAGING_DIR || exit $?
  staging_list || exit $?
  if ((TMPDAT_MODE >= 2 || TMPOUT_MODE >= 2)); then
    pdbash node all $SCRP_DIR/src/stage_in.sh || exit $?
  fi
fi

#===============================================================================
# Run initialization scripts on all nodes

if ((TMPRUN_MODE <= 2)); then
  pdbash node one $SCRP_DIR/src/init_all_node.sh $myname1 $CYCLE || exit $?
else
  pdbash node all $SCRP_DIR/src/init_all_node.sh $myname1 $CYCLE || exit $?
fi

#===============================================================================
# Run cycles of forecasts

declare -a stimes
declare -a stimesfmt
lcycles=$((LCYCLE * CYCLE_SKIP))
s_flag=1
e_flag=0
time=$STIME
loop=0

#-------------------------------------------------------------------------------
while ((time <= ETIME)); do
#-------------------------------------------------------------------------------

  loop=$((loop+1))

  for c in $(seq $CYCLE); do
    time2=$(datetime $time $((lcycles * (c-1))) s)
    if (($(datetime $time2 $lcycles s) > ETIME)); then
      e_flag=1
    fi

    if ((time2 <= ETIME)); then
      stimes[$c]=$time2
      stimesfmt[$c]="$(datetime_fmt ${stimes[$c]})"
      rcycle=$c  # The "real" number of cycles
    else
      stimes[$c]=
      stimesfmt[$c]=
    fi
  done

#-------------------------------------------------------------------------------
# Write the header of the log file

  echo
  echo " +----------------------------------------------------------------+"
  echo " |                        SCALE-Forecasts                         |"
  echo " +----------------------------------------------------------------+"
  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then
      printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
    fi
  done
  echo " +----------------------------------------------------------------+"
  echo
  echo "  Number of cycles:         $rcycle"
  echo "  Forecast start time:"
  for c in $(seq $rcycle); do
    printf "    Cycle %-5s %s\n" "$c:" "${stimesfmt[$c]}"
  done
  echo
  echo "  Forecast length:          $FCSTLEN s"
  echo "  Output interval:          $FCSTOUT s"
  echo
  echo "  Nodes used:               $NNODES"
#  if ((MTYPE == 1)); then
    for n in $(seq $NNODES); do
      echo "    ${node[$n]}"
    done
#  fi
  echo
  echo "  Processes per node:       $PPN"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
  echo "  Number of members:        $fmember"
  for c in $(seq $rcycle); do
    echo "    Cycle $c:"
    for m in $(seq $fmember); do
      mm=$(((c-1) * fmember + m))
      echo "      ${name_m[$m]}: ${node_m[$mm]}"
    done
  done
  echo

#-------------------------------------------------------------------------------
# Call functions to run the job

  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then

      ######
      if ((s == 1)); then
        if [ "$TOPO_FORMAT" == 'prep' ] && [ "$LANDUSE_FORMAT" == 'prep' ]; then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (use prepared topo and landuse files)" >&2
          continue
        elif ((BDY_FORMAT == 0)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (use prepared boundary files)" >&2
          continue
        elif ((LANDUSE_UPDATE != 1 && loop > 1)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (already done in the first cycle)" >&2
          continue
        fi
      fi
      if ((s == 2)); then
        if ((BDY_FORMAT == 0)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (use prepared boundary files)" >&2
          continue
        fi
      fi
      ######

      echo "[$(datetime_now)] ${time}: ${stepname[$s]}" >&2

      enable_iter=0
      if ((s == 2 && BDY_ENS == 1)); then
        enable_iter=1
      elif ((s == 3)); then
        enable_iter=1
      fi

      stdout_dir="$TMPOUT/${stimes[1]}/log/fcst_$(basename ${stepexecdir[$s]})"

      if ((enable_iter == 1)); then
        for it in $(seq $nitmax); do
          if [ "$STG_TYPE" = 'K_rankdir' ]; then
            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: start" >&2

            mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT-${it}" ${stepexecdir[$s]} \
                    "$(rev_path ${stepexecdir[$s]})/fcst_step.sh" $loop $it || exit $?

            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: end" >&2
          else
            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: start" >&2

            mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT-${it}" . \
                    "$SCRP_DIR/fcst_step.sh" $loop $it || exit $?

            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: end" >&2
          fi
        done
      else
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" ${stepexecdir[$s]} \
                  "$(rev_path ${stepexecdir[$s]})/fcst_step.sh" $loop || exit $?
        else
          mpirunf proc ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" . \
                  "$SCRP_DIR/fcst_step.sh" $loop || exit $?
        fi
      fi

    fi
  done

#-------------------------------------------------------------------------------
# Online stage out

#  if ((ONLINE_STGOUT == 1)); then
#    if ((MACHINE_TYPE == 11)); then
#      touch $TMP/loop.${loop}.done
#    fi
#    if ((BUILTIN_STAGING && $(datetime $time $((lcycles * CYCLE)) s) <= ETIME)); then
#      if ((MACHINE_TYPE == 12)); then
#        echo "[$(datetime_now)] ${stimes[1]}: Online stage out"
#        bash $SCRP_DIR/src/stage_out.sh s $loop || exit $?
#        pdbash node all $SCRP_DIR/src/stage_out.sh $loop || exit $?
#      else
#        echo "[$(datetime_now)] ${stimes[1]}: Online stage out (background job)"
#        ( bash $SCRP_DIR/src/stage_out.sh s $loop ;
#          pdbash node all $SCRP_DIR/src/stage_out.sh $loop ) &
#      fi
#    fi
#  fi

#-------------------------------------------------------------------------------
# Write the footer of the log file

  echo " +----------------------------------------------------------------+"
  echo " |             SCALE-Forecasts successfully completed             |"
  echo " +----------------------------------------------------------------+"
  echo

#-------------------------------------------------------------------------------

  time=$(datetime $time $((lcycles * CYCLE)) s)
  s_flag=0

#-------------------------------------------------------------------------------
done
#-------------------------------------------------------------------------------

#===============================================================================
# Stage out

if [ "$STG_TYPE" = 'builtin' ]; then
  echo "[$(datetime_now)] Finalization (stage out)" >&2

  if ((TMPOUT_MODE >= 2)); then
    if ((ONLINE_STGOUT == 1)); then
      wait
      bash $SCRP_DIR/src/stage_out.sh s $loop || exit $?
      pdbash node all $SCRP_DIR/src/stage_out.sh $loop || exit $?
    else
      bash $SCRP_DIR/src/stage_out.sh s || exit $?
      pdbash node all $SCRP_DIR/src/stage_out.sh || exit $?
    fi
  fi

#  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
#    safe_rm_tmpdir $TMP
#  fi
#  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
#    safe_rm_tmpdir $TMPL
#  fi
fi

#===============================================================================

echo "[$(datetime_now)] Finish $myname $@" >&2

exit 0
