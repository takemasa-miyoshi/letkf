#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of SCALE boundary creation.
#  October 2014  created,  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 12)); then
  cat >&2 << EOF

[pre_scale_init.sh] Prepare a temporary directory for SCALE model run.

Usage: $0 MYRANK TOPO LANDUSE BDYORG STIME MKINIT MEM MEM_BDY TMPDIR BDY_TIME_LIST NUMBER_OF_TSTEPS NUMBER_OF_SKIP_TSTEPS [SCPCALL]

  MYRANK   My rank number (not used)
  TOPO     Basename of SCALE topography files
  LANDUSE  Basename of SCALE land use files
  BDYORG   Path of the source boundary files
           SCALE history: XXX
           WRF: Basename of WRF files
  STIME    Start time (format: YYYYMMDDHHMMSS)
  MKINIT   Make initial condition as well?
            0: No
            1: Yes
  MEM      Name of the ensemble member
  MEM_BDY  Name of the ensemble member of the boundary data source
  TMPDIR   Temporary directory to run scale-rm_init
  BDY_TIME_LIST
  NUMBER_OF_TSTEPS
  NUMBER_OF_SKIP_TSTEPS
  SCPCALL

EOF
  exit 1
fi

MYRANK="$1"; shift
TOPO="$1"; shift
LANDUSE="$1"; shift
BDYORG="$1"; shift
STIME="$1"; shift
MKINIT="$1"; shift
MEM="$1"; shift
MEM_BDY="$1"; shift
TMPDIR="$1"; shift
BDY_TIME_LIST="$1"; shift
NUMBER_OF_TSTEPS="$1"; shift
NUMBER_OF_SKIP_TSTEPS="$1"; shift
SCPCALL="${1:-cycle}"

S_YYYY=${STIME:0:4}
S_MM=${STIME:4:2}
S_DD=${STIME:6:2}
S_HH=${STIME:8:2}
S_II=${STIME:10:2}
S_SS=${STIME:12:2}

historybaselen=7

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*

TMPSUBDIR=$(basename "$(cd "$TMPDIR" && pwd)")

RESTART_OUTPUT='.false.'
if ((MKINIT == 1 || (OCEAN_INPUT == 1 && OCEAN_FORMAT == 99))); then
  RESTART_OUTPUT='.true.'
fi

if ((BDY_FORMAT == 1)); then
  FILETYPE_ORG='SCALE-RM'
  USE_NESTING='.true.'
  LATLON_CATALOGUE_FNAME="$BDYORG/latlon_domain_catalogue.txt"
elif ((BDY_FORMAT == 2)); then
  FILETYPE_ORG='WRF-ARW'
  USE_NESTING='.false.'
  LATLON_CATALOGUE_FNAME=
else
  echo "[Error] $0: Unsupport boundary file types" >&2
  exit 1
fi

NUMBER_OF_FILES=0
for time_bdy in $BDY_TIME_LIST; do
  NUMBER_OF_FILES=$((NUMBER_OF_FILES+1))
done

i=0
for time_bdy in $BDY_TIME_LIST; do
  if ((NUMBER_OF_FILES <= 1)); then
    file_number=''
  else
    file_number="_$(printf %05d $i)"
  fi
  if ((BDY_ROTATING == 1)); then
    bdyorg_path="${BDYORG}/${STIME}/${MEM_BDY}"
  else
    bdyorg_path="${BDYORG}/const/${MEM_BDY}"
  fi
  if ((BDY_FORMAT == 1)); then
    if [ -s "${bdyorg_path}/${time_bdy}/history.pe000000.nc" ]; then
      for ifile in $(cd ${bdyorg_path}/${time_bdy} ; ls history*.nc 2> /dev/null); do
        ln -fs "${bdyorg_path}/${time_bdy}/${ifile}" $TMPDIR/bdydata${file_number}${ifile:$historybaselen}
      done
    else
      echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/${time_bdy}/history.*.nc'."
      exit 1
    fi
  elif ((BDY_FORMAT == 2)); then
    if [ -s "${bdyorg_path}/wrfout_${time_bdy}" ]; then
      ln -fs "${bdyorg_path}/wrfout_${time_bdy}" $TMPDIR/bdydata${file_number}
    else
      echo "[Error] $0: Cannot find source boundary file '${bdyorg_path}/wrfout_${time_bdy}'."
      exit 1
    fi
  fi
  i=$((i+1))
done

if [ "$SCPCALL" = 'cycle' ]; then
  IO_LOG_DIR='scale_init'
else
  IO_LOG_DIR="${SCPCALL}_scale_init"
fi

mkdir -p $TMPOUT/${STIME}/bdy/${MEM_BDY}

#===============================================================================

cat $TMPDAT/conf/config.nml.scale_init | \
    sed -e "/!--IO_LOG_BASENAME--/a IO_LOG_BASENAME = \"$TMPOUT/${STIME}/log/${IO_LOG_DIR}/${MEM}_LOG\"," \
        -e "/!--TIME_STARTDATE--/a TIME_STARTDATE = $S_YYYY, $S_MM, $S_DD, $S_HH, $S_II, $S_SS," \
        -e "/!--RESTART_OUTPUT--/a RESTART_OUTPUT = $RESTART_OUTPUT," \
        -e "/!--RESTART_OUT_BASENAME--/a RESTART_OUT_BASENAME = \"${TMPSUBDIR}\/init\"," \
        -e "/!--TOPO_IN_BASENAME--/a TOPO_IN_BASENAME = \"${TOPO}\"," \
        -e "/!--LANDUSE_IN_BASENAME--/a LANDUSE_IN_BASENAME = \"${LANDUSE}\"," \
        -e "/!--LAND_PROPERTY_IN_FILENAME--/a LAND_PROPERTY_IN_FILENAME = \"${TMPDAT}/land/param.bucket.conf\"," \
        -e "/!--BASENAME_BOUNDARY--/a BASENAME_BOUNDARY = \"$TMPOUT/${STIME}/bdy/${MEM_BDY}/boundary\"," \
        -e "/!--BASENAME_ORG--/a BASENAME_ORG = \"${TMPSUBDIR}\/bdydata\"," \
        -e "/!--FILETYPE_ORG--/a FILETYPE_ORG = \"${FILETYPE_ORG}\"," \
        -e "/!--NUMBER_OF_FILES--/a NUMBER_OF_FILES = ${NUMBER_OF_FILES}," \
        -e "/!--NUMBER_OF_TSTEPS--/a NUMBER_OF_TSTEPS = ${NUMBER_OF_TSTEPS}," \
        -e "/!--NUMBER_OF_SKIP_TSTEPS--/a NUMBER_OF_SKIP_TSTEPS = ${NUMBER_OF_SKIP_TSTEPS}," \
        -e "/!--BOUNDARY_UPDATE_DT--/a BOUNDARY_UPDATE_DT = $BDYINT.D0," \
        -e "/!--LATLON_CATALOGUE_FNAME--/a LATLON_CATALOGUE_FNAME = \"${LATLON_CATALOGUE_FNAME}\"," \
        -e "/!--USE_NESTING--/a USE_NESTING = $USE_NESTING," \
        -e "/!--OFFLINE--/a OFFLINE = .true.," \
    > $TMPDIR/init.conf

#===============================================================================

exit 0
