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

[mcycle.sh] Run multiple cycles.
            *use settings in 'configure.sh'

Usage: $0 STIME [ETIME]

  STIME   Time of the first cycle (format: YYYYMMDDHH)
  ETIME   Time of the last  cycle (format: YYYYMMDDHH)
          (default: same as STIME)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})

#===============================================================================

mkdir -p log

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  timef=${time:0:10}

  if [ "$time" = "$STIME" ] && [ "$time" = "$ETIME" ]; then
    mcycle_opt=0
  else
    if [ "$time" = "$STIME" ]; then
      mcycle_opt=1
    elif [ "$time" = "$ETIME" ]; then
      mcycle_opt=3
    else
      mcycle_opt=2
    fi
  fi

  ./cycle.sh "$time" all $mcycle_opt > log/cycle_${timef}.log 2>> log/cycle.err &
  wait
#  if [ "$?" -ne 0 ]; then
#    exit $?
#  fi

time=$(datetime $time $LCYCLE h)
done

#===============================================================================

exit 0
