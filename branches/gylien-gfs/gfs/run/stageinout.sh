#!/bin/bash
#===============================================================================
#
#  Copy files between the head node and local disks on computing nodes
#  April 2013, Guo-Yuan Lien
#
#  *Require source 'configure.sh' and 'distribution.sh'
#   and run 'distribute_*' function first.
#
#===============================================================================

function stagein {
#-------------------------------------------------------------------------------
# $SHAREDISK = 0: Copy necessary files to local temporary directories on computing nodes.
# $SHAREDISK = 1: Link necessary files to temporary directories in the shared disk.
#
# Usage: stagein RUNDIR1 RUNDIR2 LOUTDIR [COPY_ALL]
#
#   RUNDIR1   Temperory run directory level 1
#   RUNDIR2   Temperory run directory level 2
#   LOUTDIR   Local output directory
#   COPY_ALL  0: Only copy files marked by 'ne' (do not clear entire directories)
#             1: Clear entire directories and then copy all files
#             (default: 1)
#
# Input files:
#   [$TMPMPI/stagein/run.1]    list of files to be copied/linked to $RUNDIR1
#   [$TMPMPI/stagein/run.2]    list of files to be copied/linked to $RUNDIR2
#   [$TMPMPI/stagein/out.XXX]  list of files with relative paths to $OUTDIR
#                               to be copied to local disk $LOUTDIR
#-------------------------------------------------------------------------------

if [ "$#" -lt 3 ]; then
  echo "[Error] function 'stagein': Insufficient arguments." 1>&2
  exit 1
fi

RUNDIR1="$1"
RUNDIR2="$2"
LOUTDIR="$3"
COPY_ALL=${4:-1}

tmpnode="$TMPMPI/node"
tmpstagein="$TMPMPI/stagein"

#-------------------------------------------------------------------------------

cd $TMPMPI

cat > stagein.sh << EOF
if [ -s "$tmpstagein/run.1" ]; then
  mkdir -p $RUNDIR1
  if [ "$COPY_ALL" = '1' ]; then
    rm -fr $RUNDIR1/*
  fi
  cd $RUNDIR1
  while read line; do
    flag="\$(echo \$line | cut -d '|' -s -f1)"
    source="\$(echo \$line | cut -d '|' -s -f2)"
    destin="\$(echo \$line | cut -d '|' -s -f3)"
    if [ '$COPY_ALL' = '1' ] || [ "\$flag" = 'ne' ]; then
      if [ ! -z "\$source" ] && [ ! -z "\$destin" ]; then
        mkdir -p "\$(dirname \$destin)"
EOF
if [ "$SHAREDISK" = '0' ]; then
  cat >> stagein.sh << EOF
        $CPCOMM -r "${HOSTPREFIX}\${source}" "\$destin"
EOF
else
  cat >> stagein.sh << EOF
        rm -f "\$destin"
        ln -s "\$source" "\$destin"
EOF
fi
cat >> stagein.sh << EOF
      fi
    fi
  done < "$tmpstagein/run.1"
fi

if [ -s "$tmpstagein/run.2" ]; then
  mkdir -p $RUNDIR2
  if [ "$COPY_ALL" = '1' ] && [ "$RUNDIR1" != "$RUNDIR2" ]; then
    rm -fr $RUNDIR2/*
  fi
  cd $RUNDIR2
  while read line; do
    flag="\$(echo \$line | cut -d '|' -s -f1)"
    source="\$(echo \$line | cut -d '|' -s -f2)"
    destin="\$(echo \$line | cut -d '|' -s -f3)"
    if [ '$COPY_ALL' = '1' ] || [ "\$flag" = 'ne' ]; then
      if [ ! -z "\$source" ] && [ ! -z "\$destin" ]; then
        mkdir -p "\$(dirname \$destin)"
EOF
if [ "$SHAREDISK" = '0' ]; then
  cat >> stagein.sh << EOF
        $CPCOMM -r "${HOSTPREFIX}\${source}" "\$destin"
EOF
else
  cat >> stagein.sh << EOF
        rm -f "\$destin"
        ln -s "\$source" "\$destin"
EOF
fi
cat >> stagein.sh << EOF
      fi
    fi
  done < "$tmpstagein/run.2"
fi
EOF

if [ "$SHAREDISK" = '0' ]; then
  cat >> stagein.sh << EOF

fullname="\$(ls $tmpstagein/out.\$(hostname)* 2> /dev/null)"
if [ -s "\$fullname" ]; then
  mkdir -p $LOUTDIR
  cd $LOUTDIR
  while read line; do
    flag="\$(echo \$line | cut -d '|' -s -f1)"
    file="\$(echo \$line | cut -d '|' -s -f2)"
    if [ '$COPY_ALL' = '1' ] || [ "\$flag" = 'ne' ]; then
      if [ ! -z "\$file" ]; then
        mkdir -p "\$(dirname \$file)"
        $CPCOMM "${HOSTPREFIX}${OUTDIR}/\${file}" "\$file"
      fi
    fi
  done < "\$fullname"
fi
EOF
fi

if [ "$SHAREDISK" = '0' ]; then
  nnodes=`cat $tmpnode/machinefile.node | wc -l`
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes bash stagein.sh &
  wait
else
  bash stagein.sh
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function stageout {
#-------------------------------------------------------------------------------
# Collect output files from local disks on computing nodes.
#
# Usage: stageout LOUTDIR [BACKGROUND]
#
#   LOUTDIR     Local output directory
#   BACKGROUND  Send this job to background
#               0: No,  do not send to background
#               1: Yes, send to background.
#                  Create a temporary file and delete it when finishing.
#               (default: No)
#
# Input files:
#   [$TMPMPI/stageout/out.XXX]  list of files with relative paths to $OUTDIR
#                                to be copied from local disk $LOUTDIR
#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  echo "[Error] function 'stageout': Insufficient arguments." 1>&2
  exit 1
fi

LOUTDIR="$1"
BACKGROUND="${2:-0}"

tmpnode="$TMPMPI/node"
tmpstageout="$TMPMPI/stageout"

#-------------------------------------------------------------------------------
if [ "$SHAREDISK" = '0' ]; then
#-------------------------------------------------------------------------------

mkdir -p $tmpstageout
cd $tmpstageout
for ifile in `ls out.*`; do
  while read line; do
    flag=`echo $line | cut -d '|' -s -f1`
    file=`echo $line | cut -d '|' -s -f2`
    if [ "$flag" = 'cp' ] || [ "$flag" = 'mv' ]; then
      mkdir -p "$OUTDIR/$(dirname $file)"
    fi
  done < "$ifile"
done

#-------------------------------------------------------------------------------

cd $TMPMPI

cat > stageout.sh << EOF
fullname="\$(ls $tmpstageout/out.\$(hostname)* 2> /dev/null)"
if [ -s "\$fullname" ]; then
  mkdir -p $LOUTDIR
  cd $LOUTDIR
  while read line; do
    flag="\$(echo \$line | cut -d '|' -s -f1)"
    file="\$(echo \$line | cut -d '|' -s -f2)"
    if [ -e "\$file" ]; then
      if [ "\$flag" = 'cp' ] || [ "\$flag" = 'mv' ]; then
        $CPCOMM "\$file" "${HOSTPREFIX}${OUTDIR}/\${file}"
      fi
      if [ "\$flag" = 'mv' ] || [ "\$flag" = 'rm' ]; then
        rm -f "\$file"
      fi
    fi
  done < "\$fullname"
fi
EOF

nnodes=`cat $tmpnode/machinefile.node | wc -l`
if [ "$BACKGROUND" = '0' ]; then
  $MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes bash stageout.sh &
  wait
else
  cat > stageout_bg.sh << EOF
cd $TMPMPI
touch 'stageout_bg.running'
$MPIBIN/mpiexec -machinefile $tmpnode/machinefile.node -n $nnodes bash stageout.sh
rm -f 'stageout_bg.running'
now=\$(date +'%Y-%m-%d %H:%M:%S')
echo "[\$now] Collect outputs (stageout) <finish on background>" 1>&2
EOF
  bash stageout_bg.sh &
fi

#-------------------------------------------------------------------------------
fi
#-------------------------------------------------------------------------------
}

#===============================================================================
