#!/bin/bash
#===============================================================================
#
#  Run before the built-in stage-in script
#
#===============================================================================

. config.main
. src/func_util.sh

MYRANK="$1"

if [ "$MYRANK" = '-' ]; then
  # If myrank is not passed using the first argument, determine myrank in another way.
  MYRANK=$MPI_DRANK
#  MYRANK=$(cat $NODEFILE_DIR/node | sed -ne "/$(hostname)/=") 
fi

#-------------------------------------------------------------------------------
# Clear TMPDAT, TMPOUT directories

if (((TMPDAT_MODE == 2 && MYRANK == 0) || TMPDAT_MODE == 3)); then
  safe_init_tmpdir $TMPDAT || exit $?
  if [ "$TMPDAT_S" != "$TMPDAT" ] && ((MYRANK == 0)); then
    safe_init_tmpdir $TMPDAT_S || exit $?
  fi
fi

if (((TMPOUT_MODE == 2 && MYRANK == 0) || TMPOUT_MODE == 3)); then
  safe_init_tmpdir $TMPOUT || exit $?
fi

#===============================================================================

exit 0
