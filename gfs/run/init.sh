#!/bin/bash
#===============================================================================
#
#  Prepare a initial ensemble from analyses at different time.
#  October 2012,          Guo-Yuan Lien
#  April   2012, modified Guo-Yuan Lien
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

[init.sh] Prepare a initial ensemble from analyses at different time.
          *use settings in 'configure.sh'

Usage: $0 STIME RTIME [RINT_H]

  STIME   Initial time of the ensemble (format: YYYYMMDDHH)
  RTIME   An arbitrary time that the first member is taken from (format: YYYYMMDDHH)
  RINT_H  Interval between each members (hour) (default: 24)

  For example, if RTIME = 1991010100, RINT_H = 24,
  then the initial ensemble are constructed by analyses at
  1991010100, 1991010200, 1991010300, ...

EOF
  exit 1
fi

STIME=$(datetime $1)
RTIME=$(datetime $2)
RINT_H=${3:-24}

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_init"
tmpssio="$TMP1/${tmpsubdir}/ssio"

#===============================================================================

echo
echo "Prepare output directory..."

./outdir.sh "$STIME"

#===============================================================================

echo
echo "Prepare initial members..."

mkdir -p $tmpssio
rm -fr $tmpssio/*
cd $tmpssio
cp $DIR/ssio/sscycle .

Syyyymmddhh=${STIME:0:10}
m=1
time=$RTIME
for m in `seq $MEMBER`; do
  yyyymmddhh=${time:0:10}
  mem=`printf '%03d' $m`
  echo "  $yyyymmddhh -> member $mem"

  rm -f fort.*
  ln -fs $INITGFS/$yyyymmddhh.sig fort.11
  ln -fs $INITGFS/$yyyymmddhh.sfc fort.12
  cp $ANLGFS/$Syyyymmddhh.sig fort.21
  cp $ANLGFS/$Syyyymmddhh.sfc fort.22
  ./sscycle
  mv -f fort.21 $OUTDIR/anal/$mem/$Syyyymmddhh.sig
  mv -f fort.22 $OUTDIR/anal/$mem/$Syyyymmddhh.sfc
time=$(datetime $time $RINT_H h)
done

cd $DIR
rm -fr $TMP1/${tmpsubdir}

#===============================================================================

exit 0
