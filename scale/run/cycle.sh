#!/bin/bash
#===============================================================================
#
#  Run data assimilation cycles.
#
#  November 2014, modified from GFS-LETKF, Guo-Yuan Lien
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle.sh [STIME ETIME MEMBERS ISTEP FSTEP]
#
#  Use settings:
#    config.main
#    config.cycle
#    config.nml.scale_pp_topo
#    config.nml.scale_pp_landuse
#    config.nml.scale_init
#    config.nml.scale
#    config.nml.obsope
#    config.nml.letkf
#
#===============================================================================

cd "$(dirname "$0")"
myname='cycle.sh'
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.main || exit $?
. config.$myname1 || exit $?

. src/func_distribute.sh || exit $?
. src/func_datetime.sh || exit $?
. src/func_util.sh || exit $?
. src/func_$myname1.sh || exit $?

echo "[$(datetime_now)] ### 1" >&2

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

setting "$1" "$2" "$3" "$4" "$5" || exit $?

echo
print_setting || exit $?

echo "[$(datetime_now)] ### 2" >&2

#-------------------------------------------------------------------------------

if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
  if ((TMPDAT_MODE <= 2 || TMPRUN_MODE <= 2 || TMPOUT_MODE <= 2)); then
    safe_init_tmpdir $TMP || exit $?
  fi
  if ((TMPDAT_MODE == 3 || TMPRUN_MODE == 3 || TMPOUT_MODE == 3)); then
    safe_init_tmpdir $TMPL || exit $?
  fi
fi

echo "[$(datetime_now)] ### 3" >&2

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

#if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
if [ "$STG_TYPE" = 'builtin' ]; then
  safe_init_tmpdir $NODEFILE_DIR || exit $?
  distribute_da_cycle machinefile $NODEFILE_DIR - "$MEMBERS" || exit $?
else
  distribute_da_cycle - - $NODEFILE_DIR/distr "$MEMBERS" || exit $?
fi

echo "[$(datetime_now)] ### 4" >&2

#===============================================================================
# Determine the staging list and then stage in

if [ "$STG_TYPE" = 'builtin' ] && ((ISTEP == 1)); then
  echo "[$(datetime_now)] Initialization (stage in)" >&2

  safe_init_tmpdir $STAGING_DIR || exit $?
  staging_list || exit $?
  if ((TMPDAT_MODE >= 2 || TMPOUT_MODE >= 2)); then
    pdbash node all $SCRP_DIR/src/stage_in_init.sh || exit $?
    pdbash node all $SCRP_DIR/src/stage_in.sh || exit $?
  fi
fi

echo "[$(datetime_now)] ### 5" >&2

#===============================================================================
# Run initialization scripts on all nodes

if ((TMPRUN_MODE <= 2)); then
  pdbash node one $SCRP_DIR/src/init_all_node.sh $myname1 || exit $?
else
  pdbash node all $SCRP_DIR/src/init_all_node.sh $myname1 || exit $?
fi

echo "[$(datetime_now)] ### 6" >&2

#===============================================================================
# Run data assimilation cycles

s_flag=1
e_flag=0
time=$STIME
atime=$(datetime $time $LCYCLE s)
loop=0

#-------------------------------------------------------------------------------
while ((time <= ETIME)); do
#-------------------------------------------------------------------------------

  timefmt="$(datetime_fmt ${time})"
  loop=$((loop+1))
  if (($(datetime $time $LCYCLE s) > ETIME)); then
    e_flag=1
  fi
  obstime $time || exit $?

#-------------------------------------------------------------------------------
# Write the header of the log file

  echo "[$(datetime_now)] ### 7" >&2

  echo
  echo " +----------------------------------------------------------------+"
  echo " |                          SCALE-LETKF                           |"
  echo " +----------------------------------------------------------------+"
  for s in $(seq $nsteps); do
    if (((s_flag == 0 || s >= ISTEP) && (e_flag == 0 || s <= FSTEP))); then
      printf " | %2d. %-58s |\n" ${s} "${stepname[$s]}"
    fi
  done
  echo " +----------------------------------------------------------------+"
  echo
  echo "  Start time:               ${timefmt}"
  echo "  Forecast length:          $CYCLEFLEN s"
  echo "  Assimilation window:      $WINDOW_S - $WINDOW_E s ($((WINDOW_E-WINDOW_S)) s)"
  echo
  echo "  Observation timeslots:"
  for is in $(seq $slot_s $slot_e); do
    if ((is == slot_b)); then
      printf "  %4d - %s [base]\n" ${is} "${timefmt_sl[$is]}"
    else
      printf "  %4d - %s\n" ${is} "${timefmt_sl[$is]}"
    fi
  done
  echo
  echo "  Nodes used:               $NNODES"
  for n in $(seq $NNODES); do
    echo "    ${node[$n]}"
  done
  echo
  echo "  Processes per node:       $PPN"
  echo "  Total processes:          $totalnp"
  echo
  echo "  Nodes per SCALE run:      $mem_nodes"
  echo "  Processes per SCALE run:  $mem_np"
  echo
  echo "  Ensemble size:            $MEMBER"
  for m in $(seq $msprd); do
    echo "      ${name_m[$m]}: ${node_m[$m]}"
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
      if ((s == 4)); then
        if ((OBSOPE_RUN == 0)); then
          echo "[$(datetime_now)] ${time}: ${stepname[$s]} ...skipped (only use integrated observation operators)" >&2
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

      nodestr=proc
      if ((ENABLE_SET == 1)); then                                    ##
        if ((s == 3)); then
          nodestr='set1.proc'
        elif ((s == 4)); then
          nodestr='set2.proc'
        elif ((s == 5)); then
          nodestr='set3.proc'
        fi
      fi

      if ((s <= 3)); then
        stdout_dir="$TMPOUT/${time}/log/$(basename ${stepexecdir[$s]})"
      else
        stdout_dir="$TMPOUT/${atime}/log/$(basename ${stepexecdir[$s]})"
      fi

#echo "$stdout_dir" >&2
#echo ${stepexecdir[$s]} >&2
#echo $(rev_path ${stepexecdir[$s]}) >&2

      if ((enable_iter == 1)); then
        for it in $(seq $nitmax); do
          if [ "$STG_TYPE" = 'K_rankdir' ]; then
            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: start" >&2

            mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT-${it}" ${stepexecdir[$s]} \
                    "$(rev_path ${stepexecdir[$s]})/${myname1}_step.sh" "$time" $loop $it || exit $?

            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: end" >&2
          else
            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: start" >&2

            mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT-${it}" . \
                    "$SCRP_DIR/${myname1}_step.sh" "$time" $loop $it || exit $?

            echo "[$(datetime_now)] ${time}: ${stepname[$s]}: $it: end" >&2
          fi
        done
      else
        if [ "$STG_TYPE" = 'K_rankdir' ]; then

          mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" ${stepexecdir[$s]} \
                  "$(rev_path ${stepexecdir[$s]})/${myname1}_step.sh" "$time" "$loop" || exit $?
        else

          mpirunf $nodestr ${stepexecdir[$s]}/${stepexecname[$s]} ${stepexecname[$s]}.conf "${stdout_dir}/NOUT" . \
                  "$SCRP_DIR/${myname1}_step.sh" "$time" "$loop" || exit $?
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
#    if ((BUILTIN_STAGING && $(datetime $time $LCYCLE s) <= ETIME)); then
#      if ((MACHINE_TYPE == 12)); then
#        echo "[$(datetime_now)] ${time}: Online stage out"
#        bash $SCRP_DIR/src/stage_out.sh s $loop || exit $?
#        pdbash node all $SCRP_DIR/src/stage_out.sh $loop || exit $?
#      else
#        echo "[$(datetime_now)] ${time}: Online stage out (background job)"
#        ( bash $SCRP_DIR/src/stage_out.sh s $loop ;
#          pdbash node all $SCRP_DIR/src/stage_out.sh $loop ) &
#      fi
#    fi
#  fi

#-------------------------------------------------------------------------------
# Write the footer of the log file

  echo " +----------------------------------------------------------------+"
  echo " |               SCALE-LETKF successfully completed               |"
  echo " +----------------------------------------------------------------+"
  echo

#-------------------------------------------------------------------------------

  time=$(datetime $time $LCYCLE s)
  atime=$(datetime $time $LCYCLE s)
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
