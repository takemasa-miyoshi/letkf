#!/bin/bash
#===============================================================================
#
#  Script to prepare the directory of obsope run; for each node.
#  December 2014  created  Guo-Yuan Lien
#
#===============================================================================

. config.main

if (($# < 10)); then
  cat >&2 << EOF

[pre_obsope_node.sh] 

Usage: $0 MYRANK STIME ATIME TMPDIR OBSDIR MEM_NODES MEM_NP SLOT_START SLOT_END SLOT_BASE [MEMBERSEQ]

  MYRANK      My rank number (not used)
  STIME       
  ATIME       Analysis time (format: YYYYMMDDHHMMSS)
  TMPDIR      Temporary directory to run the program
  OBSDIR      Directory of SCALE data files
  MEM_NODES   Number of nodes for a member
  MEM_NP      Number of processes for a member
  SLOT_START  Start observation timeslots
  SLOT_END    End observation timeslots
  SLOT_BASE   The base slot
  MEMBERSEQ

EOF
  exit 1
fi

MYRANK="$1"; shift
STIME="$1"; shift
ATIME="$1"; shift
TMPDIR="$1"; shift
OBSDIR="$1"; shift
MEM_NODES="$1"; shift
MEM_NP="$1"; shift
SLOT_START="$1"; shift
SLOT_END="$1"; shift
SLOT_BASE="$1"; shift
MEMBERSEQ=${1:-$MEMBER}

#===============================================================================

# Moved to init_all_node.sh
#-- H08 --
#if [ -e "$TMPDAT/rttov/rtcoef_himawari_8_ahi.dat" ]; then
#  ln -fs $TMPDAT/rttov/rtcoef_himawari_8_ahi.dat $TMPDIR
#fi
#if [ -e "$TMPDAT/rttov/sccldcoef_himawari_8_ahi.dat" ]; then
#  ln -fs $TMPDAT/rttov/sccldcoef_himawari_8_ahi.dat $TMPDIR
#fi

OBS_IN_NAME_LIST=
for iobs in $(seq $OBSNUM); do
  if [ "${OBSNAME[$iobs]}" != '' ]; then
#    ln -fs $OBSDIR/${OBSNAME[$iobs]}_${ATIME}.dat $TMPDIR/${OBSNAME[$iobs]}.dat
    OBS_IN_NAME_LIST="${OBS_IN_NAME_LIST}'$OBSDIR/${OBSNAME[$iobs]}_${ATIME}.dat', "
  fi
done

OBSDA_RUN_LIST=
for iobs in $(seq $OBSNUM); do
  if [ -n "${OBSOPE_SEPARATE[$iobs]}" ] && ((${OBSOPE_SEPARATE[$iobs]} == 1)); then
    OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.true., "
  else
    OBSDA_RUN_LIST="${OBSDA_RUN_LIST}.false., "
  fi
done

#===============================================================================

cat $TMPDAT/conf/config.nml.obsope | \
    sed -e "/!--MEMBER--/a MEMBER = $MEMBERSEQ," \
        -e "/!--OBS_IN_NUM--/a OBS_IN_NUM = $OBSNUM," \
        -e "/!--OBS_IN_NAME--/a OBS_IN_NAME = $OBS_IN_NAME_LIST" \
        -e "/!--OBSDA_RUN--/a OBSDA_RUN = $OBSDA_RUN_LIST" \
        -e "/!--OBSDA_OUT_BASENAME--/a OBSDA_OUT_BASENAME = '${TMPOUT}/${ATIME}/obsgues/@@@@/obsda.ext'," \
        -e "/!--HISTORY_IN_BASENAME--/a HISTORY_IN_BASENAME = '${TMPOUT}/${STIME}/hist/@@@@/history'," \
        -e "/!--SLOT_START--/a SLOT_START = $SLOT_START," \
        -e "/!--SLOT_END--/a SLOT_END = $SLOT_END," \
        -e "/!--SLOT_BASE--/a SLOT_BASE = $SLOT_BASE," \
        -e "/!--SLOT_TINTERVAL--/a SLOT_TINTERVAL = $LTIMESLOT.D0," \
        -e "/!--NNODES--/a NNODES = $NNODES," \
        -e "/!--PPN--/a PPN = $PPN," \
        -e "/!--MEM_NODES--/a MEM_NODES = $MEM_NODES," \
        -e "/!--MEM_NP--/a MEM_NP = $MEM_NP," \
    > $TMPDIR/obsope.conf

# These parameters are not important for obsope
cat $TMPDAT/conf/config.nml.scale >> $TMPDIR/obsope.conf

#===============================================================================

exit 0
