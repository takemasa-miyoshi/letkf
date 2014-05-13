#!/bin/bash
#===============================================================================
#
#  Run multiple ensemble verification.
#  December 2012, Guo-Yuan Lien
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

[mverify.sh] Run multiple ensemble verification."
             *use settings in 'configure.sh'"

Usage: $0 STIME [ETIME]"

  STIME   First verifivation time (format: YYYYMMDDHH)
  ETIME   Last  verification time (format: YYYYMMDDHH)
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

  ./verify.sh "$time" mean > log/verify_${timef}.log 2>> log/verify.err &
  wait
#  if [ "$?" -ne 0 ]; then
#    exit $?
#  fi

time=$(datetime $time $LCYCLE h)
done

#===============================================================================

exit 0
