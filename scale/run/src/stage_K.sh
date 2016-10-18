#!/bin/bash
#===============================================================================
#
#  Create stage-in and stage-out scripts for the K-computer
#
#===============================================================================

. config.main

STAGING_DIR="$1"; shift
PROGNAME="$1"

#===============================================================================
# stage-in: Files in TMPDAT directory

if [ -s "$STAGING_DIR/stagein.dat" ]; then
  irank=0

  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    ftype="$(echo $line | cut -d '|' -s -f3)"
    if [ "$ftype" = 's' ]; then
      TMPDAT_STGtmp=$TMPDAT_S_STG
    else
      TMPDAT_STGtmp=$TMPDAT_STG
      if ((TMPDAT_MODE == 3)); then
        ftype='l'
      else
        ftype='s'
      fi
    fi
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          if [ "$ftype" = 'l' ]; then
            echo "#PJM --stgin-dir \"rank=* ${source} %r:${TMPDAT_STGtmp}/${destin} recursive=10\""
          else
#            echo "#PJM --stgin-dir \"rank=0 ${source} 0:${TMPDAT_STGtmp}/${destin} recursive=10\""
            echo "#PJM --stgin-dir \"rank=${irank} ${source} ${irank}:${TMPDAT_STGtmp}/${destin} recursive=10\""
            irank=$((irank+1))
            if ((irank >= $NNODES)); then
              irank=0
            fi
          fi
        else
          echo "#PJM --stgin-dir \"${source} ${TMPDAT_STGtmp}/${destin} recursive=10\""
        fi
      else
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          if [ "$ftype" = 'l' ]; then
            echo "#PJM --stgin \"rank=* ${source} %r:${TMPDAT_STGtmp}/${destin}\""
          else
#            echo "#PJM --stgin \"rank=0 ${source} 0:${TMPDAT_STGtmp}/${destin}\""
            echo "#PJM --stgin \"rank=${irank} ${source} ${irank}:${TMPDAT_STGtmp}/${destin}\""
            irank=$((irank+1))
            if ((irank >= $NNODES)); then
              irank=0
            fi
          fi
        else
          echo "#PJM --stgin \"${source} ${TMPDAT_STGtmp}/${destin}\""
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.dat" | sort | uniq
fi

#-------------------------------------------------------------------------------
# stage-in: Files in TMPOUT directory

if [ -s "$STAGING_DIR/stagein.out" ]; then
  irank=0

  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          if ((TMPOUT_MODE == 3)); then
            echo "#PJM --stgin-dir \"rank=* ${source} %r:${TMPOUT_STG}/${destin} recursive=10\""
          else
#            echo "#PJM --stgin-dir \"rank=0 ${source} 0:${TMPOUT_STG}/${destin} recursive=10\""
            echo "#PJM --stgin-dir \"rank=${irank} ${source} ${irank}:${TMPOUT_STG}/${destin} recursive=10\""
            irank=$((irank+1))
            if ((irank >= $NNODES)); then
              irank=0
            fi
          fi
        else
          echo "#PJM --stgin-dir \"${source} ${TMPOUT_STG}/${destin} recursive=10\""
        fi
      else
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          if ((TMPOUT_MODE == 3)); then
            echo "#PJM --stgin \"rank=* ${source} %r:${TMPOUT_STG}/${destin}\""
          else
#            echo "#PJM --stgin \"rank=0 ${source} 0:${TMPOUT_STG}/${destin}\""
            echo "#PJM --stgin \"rank=${irank} ${source} ${irank}:${TMPOUT_STG}/${destin}\""
            irank=$((irank+1))
            if ((irank >= $NNODES)); then
              irank=0
            fi
          fi
        else
          echo "#PJM --stgin \"${source} ${TMPOUT_STG}/${destin}\""
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.out" | sort | uniq
fi

#-------------------

i=0
while [ -s "$STAGING_DIR/stagein.out.$((i+1))" ]; do
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      if [ -d "$source" ]; then
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          echo "#PJM --stgin-dir \"rank=${i} ${source} ${i}:${TMPOUT_STG}/${destin} recursive=10\""
        else
          echo "#PJM --stgin-dir \"${source} ${TMPOUT_STG}/${destin} recursive=10\""
        fi
      else
        if [ "$STG_TYPE" = 'K_rankdir' ]; then
          echo "#PJM --stgin \"rank=${i} ${source} ${i}:${TMPOUT_STG}/${destin}\""
        else
          echo "#PJM --stgin \"${source} ${TMPOUT_STG}/${destin}\""
        fi
      fi
    fi
  done < "$STAGING_DIR/stagein.out.$((i+1))" | sort | uniq
  i=$((i+1))
done

#-------------------------------------------------------------------------------
# stage-in: nodefiles

if [ "$STG_TYPE" = 'K_rankdir' ]; then
  echo "#PJM --stgin-dir \"rank=* $TMPS/node %r:./node\""
#  echo "#PJM --stgin-dir \"rank=0 $TMPS/node 0:./node\""
else
  echo "#PJM --stgin-dir \"$TMPS/node ./node\""
fi

#-------------------------------------------------------------------------------
# stage-in: scripts

if [ "$STG_TYPE" = 'K_rankdir' ]; then
  echo "#PJM --stgin \"rank=* $TMPS/config.main %r:./config.main\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.rc %r:./config.rc\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/config.${PROGNAME} %r:./config.${PROGNAME}\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/${PROGNAME}.sh %r:./${PROGNAME}.sh\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/${PROGNAME}_step.sh %r:./${PROGNAME}_step.sh\""
  echo "#PJM --stgin \"rank=* $SCRP_DIR/src/* %r:./src/\""
else
  echo "#PJM --stgin \"$TMPS/config.main ./config.main\""
  echo "#PJM --stgin \"$SCRP_DIR/config.rc ./config.rc\""
  echo "#PJM --stgin \"$SCRP_DIR/config.${PROGNAME} ./config.${PROGNAME}\""
  echo "#PJM --stgin \"$SCRP_DIR/${PROGNAME}.sh ./${PROGNAME}.sh\""
  echo "#PJM --stgin \"$SCRP_DIR/${PROGNAME}_step.sh ./${PROGNAME}_step.sh\""
  echo "#PJM --stgin \"$SCRP_DIR/src/* ./src/\""
fi

#===============================================================================

#if ((SIMPLE_STGOUT <= 1)); then
##------

#-------------------------------------------------------------------------------
# stage-out: Files in TMPOUT directory

  if [ -s "$STAGING_DIR/stageout.out" ]; then
    irank=0

    while read line; do
      destin="$(echo $line | cut -d '|' -s -f1)"
      source="$(echo $line | cut -d '|' -s -f2)"
      ftype="$(echo $line | cut -d '|' -s -f3)"
      if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
        if [ "$ftype" = 'd' ] || [ "$ftype" = 'drm' ]; then
          if [ "$STG_TYPE" = 'K_rankdir' ]; then
            if ((TMPOUT_MODE == 3)); then
              echo "#PJM --stgout-dir \"rank=* %r:${TMPOUT_STG}/${source} ${destin} recursive=10\""
            else
#              echo "#PJM --stgout-dir \"rank=0 0:${TMPOUT_STG}/${source} ${destin} recursive=10\""
              echo "#PJM --stgout-dir \"rank=${irank} ${irank}:${TMPOUT_STG}/${source} ${destin} recursive=10\""
              irank=$((irank+1))
              if ((irank >= $NNODES)); then
                irank=0
              fi
            fi
          else
            echo "#PJM --stgout-dir \"${TMPOUT_STG}/${source} ${destin} recursive=10\""
          fi
        else
          if [ "$STG_TYPE" = 'K_rankdir' ]; then
            if ((TMPOUT_MODE == 3)); then
              echo "#PJM --stgout \"rank=* %r:${TMPOUT_STG}/${source} ${destin}\""
            else
#              echo "#PJM --stgout \"rank=0 0:${TMPOUT_STG}/${source} ${destin}\""
              echo "#PJM --stgout \"rank=${irank} ${irank}:${TMPOUT_STG}/${source} ${destin}\""
              irank=$((irank+1))
              if ((irank >= $NNODES)); then
                irank=0
              fi
            fi
          else
            echo "#PJM --stgout \"${TMPOUT_STG}/${source} ${destin}\""
          fi
        fi
      fi
    done < "$STAGING_DIR/stageout.out" # | sort | uniq
  fi

  #-------------------

  i=0
  while [ -s "$STAGING_DIR/stageout.out.$((i+1))" ]; do
    while read line; do
      destin="$(echo $line | cut -d '|' -s -f1)"
      source="$(echo $line | cut -d '|' -s -f2)"
      ftype="$(echo $line | cut -d '|' -s -f3)"
      if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
        if [ "$ftype" = 'd' ] || [ "$ftype" = 'drm' ]; then
          if [ "$STG_TYPE" = 'K_rankdir' ]; then
            echo "#PJM --stgout-dir \"rank=${i} ${i}:${TMPOUT_STG}/${source} ${destin} recursive=10\""
          else
            echo "#PJM --stgout-dir \"${TMPOUT_STG}/${source} ${destin} recursive=10\""
          fi
        else
          if [ "$STG_TYPE" = 'K_rankdir' ]; then
            echo "#PJM --stgout \"rank=${i} ${i}:${TMPOUT_STG}/${source} ${destin}\""
          else
            echo "#PJM --stgout \"${TMPOUT_STG}/${source} ${destin}\""
          fi
        fi
      fi
    done < "$STAGING_DIR/stageout.out.$((i+1))" # | sort | uniq
    i=$((i+1))
  done

#-------------------------------------------------------------------------------
# stage-out: standard log files

##------
#else # [ SIMPLE_STGOUT <= 1 ]
##------

#-------------------------------------------------------------------------------
# stage-out: everything

#  ALLOUTDIR="$OUTDIR/everything"

#  if [ "$STG_TYPE" = 'K_rankdir' ]; then
##    echo "#PJM --stgout-dir \"rank=* %r:./log $ALLOUTDIR/log recursive=10,stgout=all\""
#    echo "#PJM --stgout-dir \"rank=* %r:./run $ALLOUTDIR/run recursive=10,stgout=all\""
#    echo "#PJM --stgout-dir \"rank=* %r:./out $ALLOUTDIR/out recursive=10,stgout=all\""
#    echo "#PJM --stgout \"rank=* %r:./* $ALLOUTDIR/\""
#  else
##    echo "#PJM --stgout-dir \"./log $ALLOUTDIR/log recursive=10,stgout=all\""
#    echo "#PJM --stgout-dir \"./run $ALLOUTDIR/run recursive=10,stgout=all\""
#    echo "#PJM --stgout-dir \"./out $ALLOUTDIR/out recursive=10,stgout=all\""
#    echo "#PJM --stgout \"./* $ALLOUTDIR/\""
#  fi

##------
#fi

#===============================================================================

exit 0
