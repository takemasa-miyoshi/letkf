#!/bin/bash
#===============================================================================
#
#  Prepare a initial ensemble from outputs of another experiments.
#  May     2012           Guo-Yuan Lien
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

if [ "$#" -lt 2 ]; then
  cat 1>&2 << EOF

[init2.sh] Prepare a initial ensemble from outputs of another experiments.
           *use settings in 'configure.sh'

Usage: $0 STIME S_OUTDIR

  STIME     Initial time of the ensemble (format: YYYYMMDDHH)
  S_OUTDIR  The output directory of the "parent experiment"

EOF
  exit 1
fi

STIME=$(datetime $1)
S_OUTDIR=$2

if [ ! -d "$S_OUTDIR" ]; then
  echo "[Error] $0: '$S_OUTDIR' does not exist." 1>&2
  exit 1
fi
if [ "$OUTDIR" -ef "$S_OUTDIR" ]; then
  echo "[Error] $0: \$OUTDIR ('$OUTDIR') and \$S_OUTDIR are the same." 1>&2
  exit 1
fi

#===============================================================================

echo
echo "Prepare output directory..."

./outdir.sh "$STIME"

#===============================================================================

echo
echo "Prepare initial members..."

Syyyymmddhh=${STIME:0:10}
Eyyyymmddhh=$(datetime $STIME $LCYCLE h | cut -c 1-10)

for m in `seq $MEMBER`; do
  mem=`printf '%03d' $m`
  cp $S_OUTDIR/anal/${mem}/$Syyyymmddhh.sig $OUTDIR/anal/${mem}/$Syyyymmddhh.sig
  cp $S_OUTDIR/anal/${mem}/$Syyyymmddhh.sfc $OUTDIR/anal/${mem}/$Syyyymmddhh.sfc
done

if [ -s "$S_OUTDIR/infl/${Syyyymmddhh}.grd" ]; then
  cp $S_OUTDIR/infl/${Syyyymmddhh}.grd $OUTDIR/infl/${Syyyymmddhh}.grd
fi
if [ -s "$S_OUTDIR/infl/${Eyyyymmddhh}.grd" ]; then
  cp $S_OUTDIR/infl/${Eyyyymmddhh}.grd $OUTDIR/infl/${Eyyyymmddhh}.grd
fi

echo

#===============================================================================

exit 0
