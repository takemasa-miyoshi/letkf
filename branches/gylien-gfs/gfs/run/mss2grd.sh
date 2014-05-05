#!/bin/bash
#===============================================================================
#
#  Convert NCEP sig/sfc files to GrADS grd/grdp files
#   (in sigma/pressure coordinate)
#
#  Created  April  2013, Guo-Yuan Lien
#  Modified July   2013, Guo-Yuan Lien
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

[mss2grd.sh] Convert NCEP sig/sfc files to GrADS grd/grdp files.
             *use settings in 'configure.sh'

Usage: $0 STIME [ETIME]

  STIME  Start time (format: YYYYMMDDHH)
  ETIME  End   time (format: YYYYMMDDHH)
         (default: same as STIME)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_mss2grd"
tmpssio="$TMP1/${tmpsubdir}/ssio"

#===============================================================================

mkdir -p $ANLGRD
mkdir -p $ANLGRDP

mkdir -p $tmpssio
rm -fr $tmpssio/*
cd $tmpssio
cp $DIR/ssio/ss2grd .
cp $DIR/ssio/ss2grdp .
cp $DIR/ssio/grdctl .

if [ ! -s "$ANLGRD/yyyymmddhhx.ctl" ]; then
  STIMEgrads=$(datetimegrads $STIME)
  ./grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 x > \
           $ANLGRD/yyyymmddhhx.ctl
fi
if [ ! -s "$ANLGRDP/yyyymmddhhp.ctl" ]; then
  STIMEgrads=$(datetimegrads $STIME)
  ./grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
           $ANLGRDP/yyyymmddhhp.ctl
fi

time=$STIME
while [ "$time" -le "$ETIME" ]; do
  yyyymmddhh=${time:0:10}
  echo "  $yyyymmddhh"

  rm -f fort.*
  ln -fs $ANLGFS/$yyyymmddhh.sig fort.11
  ln -fs $ANLGFS/$yyyymmddhh.sfc fort.12
  ./ss2grd
  mv fort.31 $ANLGRD/$yyyymmddhh.grd
  ./ss2grdp
  mv fort.31 $ANLGRDP/$yyyymmddhh.grd
time=$(datetime $time $LCYCLE h)
done

cd $DIR
rm -fr $TMP1/${tmpsubdir}

#===============================================================================

exit 0
