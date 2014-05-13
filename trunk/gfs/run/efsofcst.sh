#!/bin/bash
#===============================================================================
#
#  Run multiple cycles.
#  Created  October  2012, Guo-Yuan Lien
#  Modified June     2012, Guo-Yuan Lien
#
#===============================================================================

RUNDIR="$( cd "$( dirname "$0" )" && pwd )"
cd "$RUNDIR"
if [ -f configure.sh ]; then
  . configure.sh
else
  echo "[Error] $0: 'configure.sh' does not exist." 1>&2
  exit 1
fi
. datetime.sh

if [ "$#" -lt 1 ]; then
  cat 1>&2 << EOF

[efsofcst.sh] Run ensemble forecasts for EFSO computation.
              *use settings in 'configure.sh'

Usage: $0 STIME [ETIME] [IF_MEAN] [CYCLES]

  STIME    Time of the first cycle (format: YYYYMMDDHH)
  ETIME    Time of the last  cycle (format: YYYYMMDDHH)
           (default: same as STIME)
  IF_MEAN  Also run the ensemble mean forecast?
           0: No
           1: Yes
           (default: 0)
  CYCLES   Number of forecast cycles run in parallel
           (default: 1)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})
IF_MEAN=${3:-0}
CYCLES=${4:-1}

#===============================================================================

mkdir -p log

time=$STIME
cyc=1
istime=$time
istimef=${time:0:10}
while [ "$time" -le "$ETIME" ]; do

  if [ "$cyc" -eq "$CYCLES" ]; then
    if [ "$IF_MEAN" = '1' ]; then
      ./fcst.sh "$istime" all $cyc 1 0 1 > log/efsofcst_${istimef}.log 2>> log/efsofcst.err &
    else
      ./fcst.sh "$istime" mems $cyc 1 0 1 > log/efsofcst_${istimef}.log 2>> log/efsofcst.err &
    fi
    wait
#    if [ "$?" -ne 0 ]; then
#      exit $?
#    fi

    time=$(datetime $time $LCYCLE h)
    cyc=1
    istime=$time
    istimef=${time:0:10}
  else
    time=$(datetime $time $LCYCLE h)
    cyc=$((cyc+1))
  fi

done

cyc=$((cyc-1))
if [ "$cyc" -ge 1 ]; then
  if [ "$IF_MEAN" = '1' ]; then
    ./fcst.sh "$istime" all $cyc 1 0 1 > log/efsofcst_${istimef}.log 2>> log/efsofcst.err &
  else
    ./fcst.sh "$istime" mems $cyc 1 0 1 > log/efsofcst_${istimef}.log 2>> log/efsofcst.err &
  fi
  wait
fi

#===============================================================================

exit 0
