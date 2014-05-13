#!/bin/bash
#===============================================================================
#
#  Prepare precipitation observations from TMPA data
#   -- May 2013, Guo-Yuan Lien
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

[get_tmpa.sh] Prepare precipitation observations from TMPA data.
              *use settings in 'configure.sh'

Usage: $0 STIME [ETIME]

  STIME      Start time (format: YYYYMMDDHH)
  ETIME      End   time (format: YYYYMMDDHH)
             (default: same as STIME)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime ${2:-$STIME})

#TMPA_INT=3
#HDFdata_dir="/data/letkf04/gylien/sat/TMPA-3B42"
data_dir="/data/letkf04/gylien/sat/TMPA-3B42-dat-T62_accu"
prcp_obs_dir="/data/letkf04/gylien/obs/TMPA_obs_T62_accu"

mkdir -p $data_dir
mkdir -p $prcp_obs_dir

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_get_tmpa"
tmprun="$TMP1/${tmpsubdir}"

#===============================================================================

mkdir -p $tmprun
rm -fr $tmprun/*
cd $tmprun
cp $DIR/run/ncl/3B42_hdf2dat_regrid.ncl .
cp $DIR/obs/dec_prcp .

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  timef="${time:0:8}.${time:8:2}"
#  yyyymmddhh_m1=$(datetime $time -1 h | cut -c 1-10)
  yyyymmddhh=${time:0:10}
#  yyyymmddhh_p1=$(datetime $time 1 h | cut -c 1-10)

#  yyyy=${time:0:4}
#  mm=${time:4:2}
#  dd=${time:6:2}

#  echo "[3B42.${timef}.7.HDF]"

  rm -f *.hdf *.dat fort.*
#  gzip -c -d ${HDFdata_dir}/3B42.${timef}.7.HDF.Z > 3B42.in.hdf

#  ncl 3B42_hdf2dat_regrid.ncl >& /dev/null
#  mv -f 3B42.out.dat ${data_dir}/3B42.${timef}.7.dat

  ln -s ${data_dir}/${yyyymmddhh}.grd tmpa.dat
  ./dec_prcp
  mv fort.90 $prcp_obs_dir/tmpa_obs${yyyymmddhh}.dat

time=$(datetime $time $LCYCLE h)
done

rm -rf $tmprun

#===============================================================================

exit 0
