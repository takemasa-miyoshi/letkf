#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE model run.
#  August  2014  created,  Guo-Yuan Lien
#  October 2014  modified, Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 12)); then
  cat >&2 << EOF

[pre_scale.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK MEM INIT OCEAN BDY TOPO LANDUSE STIME FCSTLEN FCSTINT HISTINT TMPDIR [SCPCALL BDY_STIME]

  MYRANK   My rank number (not used)
  MEM      Name of the ensemble member
  INIT     Basename of SCALE initial files
  OCEAN    Basename of SCALE initial ocean files
  BDY      Basename of SCALE boundary files
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  FCSTLEN  Forecast length (second)
  FCSTINT  Output interval of restart files (second)
  HISTINT  Output interval of history files (second)
  TMPDIR   Temporary directory to run the model
  SCPCALL
  BDY_STIME  (format: YYYYMMDDHHMMSS)

EOF
  exit 1
fi

MYRANK="$1"; shift
MEM="$1"; shift
INIT="$1"; shift
OCEAN="$1"; shift
BDY="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
STIME="$1"; shift
FCSTLEN="$1"; shift
FCSTINT="$1"; shift
HISTINT="$1"; shift
TMPDIR="$1"; shift
SCPCALL="${1:-cycle}"; shift
BDY_STIME="${1:-$STIME}"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

BS_YYYY=${BDY_STIME:0:4}
BS_MM=${BDY_STIME:4:2}
BS_DD=${BDY_STIME:6:2}
BS_HH=${BDY_STIME:8:2}
BS_II=${BDY_STIME:10:2}
BS_SS=${BDY_STIME:12:2}

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

if [ "${TOPO_TARGZ}" = 'T' ] ; then
  tmp_length=${#TOPO}
  tmp_length=$((tmp_length-4)) # cut "topo"
  UNCOMP_DIR="$(echo $TOPO | cut -c 1-${tmp_length} )"

  tar zxvf  ${TOPO}.tar.gz -C $UNCOMP_DIR > /dev/null
fi

if [ "${LANDUSE_TARGZ}" = 'T' ] ; then
  tmp_length=${#LANDUSE}
  tmp_length=$((tmp_length-7)) # cut "landuse"
  UNCOMP_DIR="$(echo $LANDUSE | cut -c 1-${tmp_length} )"

  tar zxvf  ${LANDUSE}.tar.gz -C $UNCOMP_DIR > /dev/null
fi

if [ "$OCEAN" = '-' ]; then
  OCEAN=$INIT
fi

if [ "$SCPCALL" = 'cycle' ]; then
  IO_LOG_DIR='scale'
else
  IO_LOG_DIR="${SCPCALL}_scale"
fi

#===============================================================================

TMPSUBDIR=$(basename "$(cd "$TMPDIR" && pwd)")

cat $TMPDAT/conf/config.nml.scale | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_LOG\"," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--TIME_DURATION--/a TIME_DURATION = ${FCSTLEN}.D0," \
        -e "/!--TIME_DT_ATMOS_RESTART--/a TIME_DT_ATMOS_RESTART = ${FCSTINT}.D0," \
        -e "/!--TIME_DT_OCEAN_RESTART--/a TIME_DT_OCEAN_RESTART = ${FCSTINT}.D0," \
        -e "/!--TIME_DT_LAND_RESTART--/a TIME_DT_LAND_RESTART = ${FCSTINT}.D0," \
        -e "/!--TIME_DT_URBAN_RESTART--/a TIME_DT_URBAN_RESTART = ${FCSTINT}.D0," \
        -e "/!--RESTART_IN_BASENAME--/a RESTART_IN_BASENAME = \"${INIT}\"," \
        -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${TMPSUBDIR}\/restart\"," \
        -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"${TOPO}\"," \
        -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${LANDUSE}\"," \
        -e "/!--ATMOS_BOUNDARY_IN_BASENAME--/a ATMOS_BOUNDARY_IN_BASENAME = \"${BDY}\"," \
        -e "/!--ATMOS_BOUNDARY_START_DATE--/a ATMOS_BOUNDARY_START_DATE = $BS_YYYY, $BS_MM, $BS_DD, $BS_HH, $BS_II, $BS_SS," \
        -e "/!--ATMOS_BOUNDARY_UPDATE_DT--/a ATMOS_BOUNDARY_UPDATE_DT = $BDYINT.D0," \
        -e "/!--OCEAN_RESTART_IN_BASENAME--/a OCEAN_RESTART_IN_BASENAME = \"${OCEAN}\"," \
        -e "/!--HISTORY_DEFAULT_BASENAME--/a HISTORY_DEFAULT_BASENAME = \"${TMPSUBDIR}\/history\"," \
        -e "/!--HISTORY_DEFAULT_TINTERVAL--/a HISTORY_DEFAULT_TINTERVAL = ${HISTINT}.D0," \
        -e "/!--MONITOR_OUT_BASENAME--/a MONITOR_OUT_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_monitor\"," \
        -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPDAT}/land/param.bucket.conf\"," \
        -e "/!--ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_GASPARA_IN_FILENAME = \"${TMPDAT}/rad/PARAG.29\"," \
        -e "/!--ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_AEROPARA_IN_FILENAME = \"${TMPDAT}/rad/PARAPC.29\"," \
        -e "/!--ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME--/a ATMOS_PHY_RD_MSTRN_HYGROPARA_IN_FILENAME = \"${TMPDAT}/rad/VARDATA.RM29\"," \
        -e "/!--ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME--/a ATMOS_PHY_RD_PROFILE_CIRA86_IN_FILENAME = \"${TMPDAT}/rad/cira.nc\"," \
        -e "/!--ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME--/a ATMOS_PHY_RD_PROFILE_MIPAS2001_IN_BASENAME = \"${TMPDAT}/rad/MIPAS\"," \
    > $TMPDIR/run.conf

#===============================================================================

exit 0
