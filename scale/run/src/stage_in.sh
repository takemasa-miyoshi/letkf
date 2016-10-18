#!/bin/bash
#===============================================================================
#
#  Built-in stage-in script
#
#===============================================================================

. config.main

MYRANK="$1"   # a: run on the server node, stage in all files

if [ "$MYRANK" = '-' ]; then
  # If myrank is not passed using the first argument, determine myrank in another way.
  MYRANK=$MPI_DRANK
#  MYRANK=$(cat $NODEFILE_DIR/node | sed -ne "/$(hostname)/=") 
fi

#-------------------------------------------------------------------------------
# Files in TMPDAT directory

if [ "$MYRANK" = 'a' ] ||
   (((TMPDAT_MODE == 2 && MYRANK == 0) || TMPDAT_MODE == 3)); then

  if [ -s "$STAGING_DIR/stagein.dat" ]; then
    while read line; do
      source="$(echo $line | cut -d '|' -s -f1)"
      destin="$(echo $line | cut -d '|' -s -f2)"
      ftype="$(echo $line | cut -d '|' -s -f3)"
      if [ "$ftype" = 's' ]; then
        if ((MYRANK > 0)); then
          break
        fi
        TMPDATtmp=$TMPDAT_S
      else
        TMPDATtmp=$TMPDAT
      fi
      if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
        mkdir -p "$(dirname ${TMPDATtmp}/${destin})"
        if ((SCP_THREAD > 1)); then
          while (($(jobs -p | wc -l) >= SCP_THREAD)); do
            sleep 1s
          done
          $SCP -r "${SCP_HOSTPREFIX}${source}" "${TMPDATtmp}/${destin}" &
        else
          $SCP -r "${SCP_HOSTPREFIX}${source}" "${TMPDATtmp}/${destin}"
        fi
      fi
    done < "$STAGING_DIR/stagein.dat" | sort | uniq
    if ((SCP_THREAD > 1)); then
      wait
    fi
  fi

fi

#-------------------------------------------------------------------------------
# Files in TMPOUT directory

filelist=
if [ "$MYRANK" = 'a' ]; then
  filelist="$(ls $STAGING_DIR/stagein.out* 2> /dev/null)"
elif ((TMPOUT_MODE >= 2)); then
  if (((TMPOUT_MODE == 2 && MYRANK == 0) || TMPOUT_MODE == 3)); then
    filelist="$(ls $STAGING_DIR/stagein.out 2> /dev/null)"
  fi
  filelist="$filelist $(ls $STAGING_DIR/stagein.out.$((MYRANK+1)) 2> /dev/null)"
fi

for ifile in $filelist; do
  while read line; do
    source="$(echo $line | cut -d '|' -s -f1)"
    destin="$(echo $line | cut -d '|' -s -f2)"
    if [ ! -z "$source" ] && [ ! -z "$destin" ]; then
      mkdir -p "$(dirname ${TMPOUT}/${destin})"
      if ((SCP_THREAD > 1)); then
        while (($(jobs -p | wc -l) >= SCP_THREAD)); do
          sleep 1s
        done
        $SCP -r "${SCP_HOSTPREFIX}${source}" "${TMPOUT}/${destin}" &
      else
        $SCP -r "${SCP_HOSTPREFIX}${source}" "${TMPOUT}/${destin}"
      fi
    fi
  done < "$ifile" | sort | uniq
  if ((SCP_THREAD > 1)); then
    wait
  fi
done

#===============================================================================

exit 0
