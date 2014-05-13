#!/bin/bash
#===============================================================================
#
#  Download NCEP CFSR
#      May 2014, Guo-Yuan Lien
#
#===============================================================================

if [ -f configure.sh ]; then
  . configure.sh
else
  echo "[Error] $0: 'configure.sh' does not exist." 1>&2
  exit 1
fi
. datetime.sh

if [ "$#" -lt 1 ]; then
  cat 1>&2 << EOF

[get_cfsr.sh] Download NCEP CFSR.
                 *use settings in 'configure.sh'

Usage: $0 STIME [ETIME] [GFS_DIR]

  STIME    Start time (format: YYYYMMDD)
  ETIME    End   time (format: YYYYMMDD)
           (default: same as STIME)
  GFS_DIR  Directory to save the CFSR sig/sfc files
           (default: \$ANLGFS in 'configure.sh')

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})
GFS_DIR="${3:-$ANLGFS}"

DATAURL="http://nomads.ncdc.noaa.gov/modeldata/cmd_LIC" # T126
#DATAURL="http://nomads.ncdc.noaa.gov/modeldata/cmd_HIC" # T382
WGET="wget"

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_get_cfsr"
tmprun="$TMP1/${tmpsubdir}"

#===============================================================================

mkdir -p $tmprun
rm -fr $tmprun/*
mkdir -p $tmprun/download
cd $tmprun/download

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  yyyymmdd=${time:0:8}
  yyyy=${time:0:4}
  mm=${time:4:2}
  dd=${time:6:2}

  $WGET "${DATAURL}/$yyyy/$yyyy$mm/$yyyymmdd/cfs_reanalysis_CFS_LIC_$yyyymmdd.tar"
  tar xvf cfs_reanalysis_CFS_LIC_$yyyymmdd.tar sfcanl.gdas2.${yyyymmdd}00 \
                                               sfcanl.gdas2.${yyyymmdd}06 \
                                               sfcanl.gdas2.${yyyymmdd}12 \
                                               sfcanl.gdas2.${yyyymmdd}18 \
                                               siganl.gdas2.${yyyymmdd}00 \
                                               siganl.gdas2.${yyyymmdd}06 \
                                               siganl.gdas2.${yyyymmdd}12 \
                                               siganl.gdas2.${yyyymmdd}18
  rm -f cfs_reanalysis_CFS_LIC_$yyyymmdd.tar

  mv -f sfcanl.gdas2.${yyyymmdd}00 ${GFS_DIR}/${yyyymmdd}00.sfc
  mv -f sfcanl.gdas2.${yyyymmdd}06 ${GFS_DIR}/${yyyymmdd}06.sfc
  mv -f sfcanl.gdas2.${yyyymmdd}12 ${GFS_DIR}/${yyyymmdd}12.sfc
  mv -f sfcanl.gdas2.${yyyymmdd}18 ${GFS_DIR}/${yyyymmdd}18.sfc
  mv -f siganl.gdas2.${yyyymmdd}00 ${GFS_DIR}/${yyyymmdd}00.sig
  mv -f siganl.gdas2.${yyyymmdd}06 ${GFS_DIR}/${yyyymmdd}06.sig
  mv -f siganl.gdas2.${yyyymmdd}12 ${GFS_DIR}/${yyyymmdd}12.sig
  mv -f siganl.gdas2.${yyyymmdd}18 ${GFS_DIR}/${yyyymmdd}18.sig

time=$(datetime $time 1 d)
done

rm -rf $tmprun

#===============================================================================

exit 0
