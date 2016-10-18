#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE topo creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 5)); then
  cat >&2 << EOF

[pre_scale_pp.sh] Prepare a temporary directory for scale-rm_pp run.

Usage: $0 MYRANK STIME MEM TMPDIR DATADIR [SCPCALL COPYTOPO CATALOGUE]

  MYRANK   My rank number (not used)
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MEM
  TMPDIR   Temporary directory to run scale-rm_pp
  DATADIR  Directory of SCALE data files
  SCPCALL  Called from which script? (fcst/cycle)
  COPYTOPO
  CATALOGUE

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
MEM="$1"; shift
TMPDIR="$1"; shift
DATADIR="$1"; shift   ###### no use
SCPCALL="${1:-cycle}"; shift
COPYTOPO="$1"; shift
CATALOGUE="$1"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

CONVERT_TOPO='.false.'
if [ "$TOPO_FORMAT" != 'prep' ]; then
  CONVERT_TOPO='.true.'
fi

CONVERT_LANDUSE='.false.'
if [ "$LANDUSE_FORMAT" != 'prep' ]; then
  CONVERT_LANDUSE='.true.'
fi

if ((DISK_MODE_TOPO_LANDUSE_DB == 2)); then
  DATADIR=$TMPDAT_S
else
  DATADIR=$TMPDAT_L
fi

USE_NESTING='.false.'
OFFLINE='.true.'
if ((BDY_FORMAT == 1)) && [ "$TOPO_FORMAT" != 'prep' ]; then
  USE_NESTING='.true.'
fi

if ((LANDUSE_UPDATE == 1)); then
  LANDUSE_OUT_BASENAME="$TMPOUT/${STIME}/landuse/landuse"
else
  LANDUSE_OUT_BASENAME="$TMPOUT/const/landuse/landuse"
fi

if [ "$SCPCALL" = 'cycle' ]; then
  IO_LOG_DIR='scale_pp'
else
  IO_LOG_DIR="${SCPCALL}_scale_pp"
fi

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_pp | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_LOG\"," \
        -e "/!--TOPO_OUT_BASENAME--/a TOPO_OUT_BASENAME = \"$TMPOUT/const/topo/topo\"," \
        -e "/!--LANDUSE_OUT_BASENAME--/a LANDUSE_OUT_BASENAME = \"${LANDUSE_OUT_BASENAME}\"," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--CONVERT_TOPO--/a CONVERT_TOPO = $CONVERT_TOPO," \
        -e "/!--CONVERT_LANDUSE--/a CONVERT_LANDUSE = $CONVERT_LANDUSE," \
        -e "/!--CNVTOPO_name--/a CNVTOPO_name = \"$TOPO_FORMAT\"," \
        -e "/!--GTOPO30_IN_DIR--/a GTOPO30_IN_DIR = \"$DATADIR/topo/GTOPO30/Products\"," \
        -e "/!--DEM50M_IN_DIR--/a DEM50M_IN_DIR = \"$DATADIR/topo/DEM50M/Products\"," \
        -e "/!--CNVLANDUSE_name--/a CNVLANDUSE_name = '$LANDUSE_FORMAT'," \
        -e "/!--GLCCv2_IN_DIR--/a GLCCv2_IN_DIR = \"$DATADIR/landuse/GLCCv2/Products\"," \
        -e "/!--LU100M_IN_DIR--/a LU100M_IN_DIR = \"$DATADIR/landuse/LU100M/Products\"," \
        -e "/!--COPYTOPO_IN_BASENAME--/a COPYTOPO_IN_BASENAME = \"${COPYTOPO}\"," \
        -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${CATALOGUE}\"," \
        -e "/!--USE_NESTING--/a USE_NESTING = $USE_NESTING," \
        -e "/!--OFFLINE--/a OFFLINE = $OFFLINE," \
    > $TMPDIR/pp.conf

#===============================================================================

exit 0
