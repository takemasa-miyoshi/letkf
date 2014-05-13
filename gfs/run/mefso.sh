#!/bin/bash
#===============================================================================
#
#  Run EFSO for multiple cycles.
#  Created  January  2014, Guo-Yuan Lien
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

[mefso.sh] Run EFSO for multiple cycles.
            *use settings in 'configure.sh'

Usage: $0 STIME [ETIME] [EFT] [LOCADV_RATE] [WMOIST]

  STIME        Time of the first cycle (format: YYYYMMDDHH)
  ETIME        Time of the last  cycle (format: YYYYMMDDHH)
               (default: same as STIME)
  EFT          Evaluation forecast time (hours)
               (default: 24)
  LOCADV_RATE  Localization advection rate relative to the phase velocity (winds)
               0: No advection
               (default: 0)
  WMOIST       Wight for the moist term in the energy norm (dimensionless)
               (default: 1)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})
EFT=${3:-24}
LOCADV_RATE="${4:-0}"
WMOIST=${5:-1}

#===============================================================================

mkdir -p log

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  timef=${time:0:10}

  ./efso.sh "$time" $EFT $LOCADV_RATE $WMOIST > log/efso_${timef}.log 2>> log/efso.err &
  wait
#  if [ "$?" -ne 0 ]; then
#    exit $?
#  fi

time=$(datetime $time $LCYCLE h)
done

#===============================================================================

exit 0
