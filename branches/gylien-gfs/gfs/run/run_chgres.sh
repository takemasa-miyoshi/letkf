#!/bin/bash
#===============================================================================
#
#  Script to run 'global_chgres'
#  Change the resolution of GFS initial conditions
#  September 2013  created, Guo-Yuan Lien
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

if [ "$#" -lt 7 ]; then
  cat 1>&2 << EOF

[run_chgres.sh] Change the resolution of GFS initial conditions.
                *use settings in 'configure.sh'

Usage: $0 STIME ETIME JCAP_NEW LONB_NEW LATB_NEW GFS_INDIR GFS_OUTDIR [FORMAT]

  STIME       Start time (format: YYYYMMDDHH)
  ETIME       End   time (format: YYYYMMDDHH)
  JCAP_NEW    New JCAP (spectral truncation)
  LONB_NEW    New LONB (number of longitudes)
  LATB_NEW    New LATB (number of latitudes)
  GFS_INDIR   Input  directory of GFS sig/sfc files
  GFS_OUTDIR  Output directory of GFS sig/sfc files
  FORMAT      Input data format
              0: Do not specify the input data format
              1: GFS/GDAS data
              2: CFSR data [- 2010-12-31 18]
              3: CFSR data [2011-01-01 00 -]
              (default: 1)

EOF
  exit 1
fi

STIME=$(datetime $1)
ETIME=$(datetime $2)
JCAP_NEW=$3
LONB_NEW=$4
LATB_NEW=$5
GFS_INDIR="$6"
GFS_OUTDIR="$7"
FORMAT=${8:-1}

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_chgres"
tmprun="$TMP1/${tmpsubdir}"

chgres_script="/homes/metogra/gylien/work/gfs-letkf/branches/gylien-20140429/gfs/script/global_chgres.sh"

#===============================================================================

mkdir -p $tmprun
rm -fr $tmprun/*
mkdir -p $tmprun/tmp
cd $tmprun
cp -f $chgres_script .

export EXECGLOBAL=$EXECGLOBAL
export FIXGLOBAL=$FIXGLOBAL
export DATA=$tmprun/tmp

export JCAP=$JCAP_NEW
export LEVS=0
export LONB=$LONB_NEW
export LATB=$LATB_NEW

if [ "$FORMAT" -eq 1 ]; then
  export IDVC=2
  export IDSL=0
  export IDVM=0
elif [ "$FORMAT" -eq 2 ]; then
  export IDVC=3
  export IDSL=2
  export IDVM=11
elif [ "$FORMAT" -eq 3 ]; then
  export IDVC=2
  export IDSL=1
  export IDVM=1
fi

export LSOIL=4
export IVSSFC=200509
export LANDICE_OPT=2
export CLIMO_FIELDS_OPT=3

#===============================================================================

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  yyyymmddhh=${time:0:10}
  echo
  echo "[${yyyymmddhh}]"
  echo

  export SIGINP=${GFS_INDIR}/${yyyymmddhh}.sig
  export SFCINP=${GFS_INDIR}/${yyyymmddhh}.sfc
  export SIGOUT=${GFS_OUTDIR}/${yyyymmddhh}.sig
  export SFCOUT=${GFS_OUTDIR}/${yyyymmddhh}.sfc

  bash global_chgres.sh
  if [ "$?" -ne 0 ]; then
    exit $?
  fi

time=$(datetime $time $LCYCLE h)
done

rm -rf $tmprun

#===============================================================================

exit 0
