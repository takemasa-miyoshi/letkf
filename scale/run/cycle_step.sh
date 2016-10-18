#!/bin/bash
#===============================================================================
#
#  Run one step of data assimilation cycles.
#
#-------------------------------------------------------------------------------
#
#  Usage:
#    cycle_step.sh [STEPFUNC MYRANK TIME LOOP ITER]
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
#myname=$(basename "$0")
#myname1=${myname%.*}

#===============================================================================
# Configuration

if (($# < 3)); then
  echo "[Error] $FUNCNAME: Insufficient arguments." >&2
  exit 1
fi

#-------------------------------------------------------------------------------

. config.main
res=$? && ((res != 0)) && exit $res
. config.cycle
res=$? && ((res != 0)) && exit $res

. src/func_distribute.sh
. src/func_datetime.sh
. src/func_util.sh
. src/func_cycle.sh

######
echo "[$(datetime_now)] cycle_step: $@: start"
######

#-------------------------------------------------------------------------------

setting

STEPFUNC="${1}"; shift
MYRANK="${1}"; shift
TIME="${1}"; shift
LOOP="${1}"; shift
ITER="${1:-0}"

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

distribute_da_cycle machinefile - $NODEFILE_DIR/distr "$MEMBERS"

#===============================================================================
# Run one step

time=$TIME
loop=$LOOP
iter=$ITER
if ((ITER == 0)); then
  its=1
  ite=$nitmax
else
  its=$iter
  ite=$iter
fi

atime=$(datetime $time $LCYCLE s)
timefmt="$(datetime_fmt ${time})"
obstime $time

#-------------------------------------------------------------------------------

#echo $STEPFUNC $MYRANK $TIME $LOOP 1>&2

$STEPFUNC
res=$? && ((res != 0)) && exit $res

#echo $STEPFUNC $MYRANK $TIME $LOOP ... done 1>&2

#===============================================================================

######
echo "[$(datetime_now)] cycle_step: end"
######

exit 0
