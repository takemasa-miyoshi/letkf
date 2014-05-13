#!/bin/bash
#===============================================================================
#
#  Create necessary subdirectories and files in the output directory.
#  April 2013, Guo-Yuan Lien
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

if [ "$1" = '-h' ]; then
  cat 1>&2 << EOF

[outdir.sh] Create necessary subdirectories and files in the output directory.
            *use settings in 'configure.sh'

Usage: $0 [INIT_TIME]

  INIT_TIME  Initial time of the experiment (format: YYYYMMDDHH)
             Used for creating GrADS ctl files
             If not given, do not create GrADS ctl files

EOF
  exit 1
fi

if [ -z "$1" ]; then
  INIT_TIME=''
else
  INIT_TIME=$(datetime $1)
fi

#===============================================================================

mkdir -p $OUTDIR
mkdir -p $OUTDIR/gues
mkdir -p $OUTDIR/guesg
mkdir -p $OUTDIR/guesgp
mkdir -p $OUTDIR/anal
mkdir -p $OUTDIR/analg
mkdir -p $OUTDIR/analgp
mkdir -p $OUTDIR/infl
mkdir -p $OUTDIR/log
mkdir -p $OUTDIR/obs

mkdir -p $OUTDIR/gues/mean
mkdir -p $OUTDIR/guesg/mean
mkdir -p $OUTDIR/guesg/sprd
mkdir -p $OUTDIR/guesgp/mean

mkdir -p $OUTDIR/anal/mean
mkdir -p $OUTDIR/analg/mean
mkdir -p $OUTDIR/analg/sprd
mkdir -p $OUTDIR/analgp/mean

if [ "$INIT_TIME" != '' ]; then
  STIMEgrads=$(datetimegrads $(datetime $INIT_TIME $LCYCLE h))
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
                   $OUTDIR/guesg/mean/yyyymmddhh.ctl
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
                   $OUTDIR/guesg/sprd/yyyymmddhh.ctl
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
                   $OUTDIR/guesgp/mean/yyyymmddhhp.ctl
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
                   $OUTDIR/analg/mean/yyyymmddhh.ctl
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
                   $OUTDIR/analg/sprd/yyyymmddhh.ctl
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
                   $OUTDIR/analgp/mean/yyyymmddhhp.ctl
  $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
                   $OUTDIR/infl/yyyymmddhh.ctl
fi

#===============================================================================

for m in `seq $MEMBER`; do
  mem=`printf '%03d' $m`
  mkdir -p $OUTDIR/gues/$mem
  mkdir -p $OUTDIR/guesg/$mem
  mkdir -p $OUTDIR/guesgp/$mem
  mkdir -p $OUTDIR/anal/$mem
  mkdir -p $OUTDIR/analg/$mem
  mkdir -p $OUTDIR/analgp/$mem

  if [ "$INIT_TIME" != '' ]; then
    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 x > \
                     $OUTDIR/guesg/$mem/yyyymmddhhx.ctl
    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
                     $OUTDIR/guesgp/$mem/yyyymmddhhp.ctl
    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 s > \
                     $OUTDIR/analg/$mem/yyyymmddhh.ctl
    $DIR/ssio/grdctl '%y4%m2%d2%h2.grd' 'template byteswapped' $STIMEgrads ${LCYCLE}hr 10000 p > \
                     $OUTDIR/analgp/$mem/yyyymmddhhp.ctl
  fi
done

#===============================================================================

exit 0
