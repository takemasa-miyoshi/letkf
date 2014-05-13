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

[mfcst.sh] Run multiple ensemble mean forecasts.
           *use settings in 'configure.sh'

Usage: $0 STIME [ETIME] [CYCLES] [IF_VERF]

  STIME    Time of the first cycle (format: YYYYMMDDHH)
  ETIME    Time of the last  cycle (format: YYYYMMDDHH)
           (default: same as STIME)
  CYCLES   Number of forecast cycles run in parallel
           (default: 1)
  IF_VERF  Run verification or not
           0: Do not run verification
           1: Run verification
           (default: 0)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})
CYCLES=${3:-1}
IF_VERF=${4:-0}

#===============================================================================

mkdir -p log

time=$STIME
cyc=1
istime=$time
istimef=${time:0:10}
while [ "$time" -le "$ETIME" ]; do

  if [ "$cyc" -eq "$CYCLES" ]; then
    ./fcst.sh "$istime" mean $cyc 1 $IF_VERF > log/fcst_${istimef}.log 2>> log/fcst.err &
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
  ./fcst.sh "$istime" mean $cyc 1 $IF_VERF > log/fcst_${istimef}.log 2>> log/fcst.err &
  wait
fi

#===============================================================================

exit 0
