#!/bin/bash
#===============================================================================
#
#  Prepare a series of initial mean analyses from another source,
#  useful to run forecast experiments.
#
#  August  2013           Guo-Yuan Lien
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

[init3.sh] Prepare a series of initial mean analyses from another source.
           *use settings in 'configure.sh'

Usage: $0 STIME [ETIME] [ANL_DIR]

  STIME    Start time of this series of analyses (format: YYYYMMDDHH)
  RTIME    End   time of this series of analyses (format: YYYYMMDDHH)
           (default: STIME)
  ANL_DIR  Directory of the other source of analyses
           (default: \$ANLGFS in 'configure.sh')

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})
ANL_DIR=${3:-$ANLGFS}

if [ ! -d "$ANL_DIR" ]; then
  echo "[Error] $0: '$ANL_DIR' does not exist." 1>&2
  exit 1
fi

#===============================================================================

echo
echo "Prepare output directory..."

./outdir.sh "$STIME"

#===============================================================================

echo
echo "Prepare initial mean analyses..."
echo

time=$STIME
while [ $(datetime $time) -le $(datetime $ETIME) ]; do
  Syyyymmddhh=${time:0:10}
  echo "  $Syyyymmddhh"

  cp $ANL_DIR/$Syyyymmddhh.sig $OUTDIR/anal/mean/$Syyyymmddhh.sig
  cp $ANL_DIR/$Syyyymmddhh.sfc $OUTDIR/anal/mean/$Syyyymmddhh.sfc
time=$(datetime $time $LCYCLE h)
done

echo

#===============================================================================

exit 0
