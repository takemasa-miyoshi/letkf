#!/bin/bash
#===============================================================================
#
#  Adaptively distribute members on nodes
#  April 2013, Guo-Yuan Lien
#
#  *Require source 'configure.sh' first.
#
#===============================================================================

function distribute_da_cycle {
#-------------------------------------------------------------------------------
# Distribute members on nodes for DA cycling run.
#
# Usage: distribute_da_cycle NODEFILE [WRITE_FILES] [DISPLAY]
#
#   NODEFILE     machinefile for mpiexec
#   WRITE_FILES  output machinefiles or not
#                0: No
#                1: Yes
#                (default: Yes)
#   DISPLAY      display distribution information or not
#                0: No
#                1: Yes
#                (default: No)
#
# Output files:
#   [$TMPMPI/node/machinefile]          machinefile for mpiexec
#   [$TMPMPI/node/machinefile.node]     machinefile for mpiexec (one core per node)
#   [$TMPMPI/node/machinefile.MMM.XXX]  machinefile of each member for mpiexec
#                                       (MMM = gfs, gsi)
#                                       (XXX = 001, 002, ...)
#
# Return variables:
#   $mmean (= $MEMBER+1)
#   $msprd (= $MEMBER+2)
#   $nnodes
#   $ppn
#   $totalnp
#   $mem_nodes_gfs
#   $mem_nodes_gsi
#   $mem_np_gfs
#   $mem_np_gsi
#   $node   [1, ..., $nnodes  ]
#   $node_m [1, ..., $MEMBER+2]
#   $name_m [1, ..., $MEMBER+2]
#   $nodes_gfs_m[1,...,$MEMBER]
#   $nodes_gsi_m[1,...,$MEMBER]
#-------------------------------------------------------------------------------

if [ "$#" -lt 1 ]; then
  echo "[Error] function 'distribute_da_cycle': Insufficient arguments." 1>&2
  exit 1
fi

local NODEFILE="$1"
local WRITE_FILES=${2:-1}
local DISPLAY=${3:-0}

if [ ! -s "$NODEFILE" ]; then
  echo "[Error] $0: Can not find \$NODEFILE '$NODEFILE'" 1>&2
  exit 1
fi

tmpnode="$TMPMPI/node"

mmean=$((MEMBER+1))
msprd=$((MEMBER+2))

#-------------------------------------------------------------------------------
# Read $NODEFILE and determine
#   $nnodes
#   $ppn
#   $totalnp
#   $node[1,...,$nnodes]

local nodelist=`cat $NODEFILE | sort | uniq`

local n=0
for inode in $nodelist; do
  n=$((n+1))
  node[$n]=$inode

  ippn=`cat $NODEFILE | grep $inode | wc -l`
  if [ "$n" -eq 1 ]; then
    ppn=$ippn
  elif [ "$ippn" -ne "$ppn" ]; then
    echo "[Error] $0: The number of cores per node must be a constant" 1>&2
    exit 1
  fi
done
nnodes=$n
totalnp=$((ppn*nnodes))

#-------------------------------------------------------------------------------
# If the minimum cores required to run a GFS or GSI member cannot be afforded, then exit.

if [ "$MIN_NP_GFS" -gt "$ppn" ] || [ "$MIN_NP_GSI" -gt "$ppn" ] && [ "$SHAREDISK" = '0' ]; then
  echo "[Error] When using local disks, minimum cores required to run a GFS or GSI cannot exceed cores per node." 1>&2
  exit 1
fi
if [ "$MIN_NP_GFS" -gt "$totalnp" ] || [ "$MIN_NP_GSI" -gt "$totalnp" ]; then
  echo "[Error] Minimum cores required to run a GFS or GSI exceed total avaiable cores." 1>&2
  exit 1
fi

#-------------------------------------------------------------------------------
# Find optimal number of cores ($mem_np) based on the ensemble size ($MEMBER)
# and the total nodes ($nnodes) and cores per node ($ppn).
# The 'optimal' here means minimum $mem_np but still can occupy all avaiable cores at once.

if [ "$nnodes" -ge "$MEMBER" ]; then
  if [ "$SHAREDISK" = '1' ]; then
    mem_nodes=$((nnodes/MEMBER))
  else
    mem_nodes=1
  fi
  mem_np=$((ppn*mem_nodes))
else
  mem_nodes=1
  mempn=$(((MEMBER-1)/nnodes+1))
  if [ "$mempn" -gt "$ppn" ]; then
    mem_np=1
  else
    mem_np=$((ppn/mempn))
  fi
fi

#-------------------------------------------------------------------------------
# Limited to the minimum cores required to run a GFS member ($MIN_NP_GFS)
#        and the maximum cores to run a GFS member ($MAX_NP_GFS),
# and require it occupy full nodes if multiple nodes are used,
# determine the numbers of cores ($mem_np_gfs) and nodes ($mem_nodes_gfs) used to run a GFS member
#        from the optimal number of cores per member ($mem_np),

if [ "$mem_np" -lt "$MIN_NP_GFS" ]; then
  mem_np_gfs=$MIN_NP_GFS
else
  mem_np_gfs=$mem_np
fi
if [ "$mem_np_gfs" -gt "$MAX_NP_GFS" ]; then
  mem_np_gfs=$MAX_NP_GFS
fi
if [ "$mem_np_gfs" -gt "$ppn" ]; then
  mem_nodes_gfs=$(((mem_np_gfs-1)/ppn+1))
  mem_np_gfs=$((ppn*mem_nodes_gfs))
else
  mem_nodes_gfs=1
fi

#-------------------------------------------------------------------------------
# Same as the above section, but for GSI usage ($mem_np_gsi, $mem_nodes_gsi).

if [ "$mem_np" -lt "$MIN_NP_GSI" ]; then
  mem_np_gsi=$MIN_NP_GSI
else
  mem_np_gsi=$mem_np
fi
if [ "$mem_np_gsi" -gt "$MAX_NP_GSI" ]; then
  mem_np_gsi=$MAX_NP_GSI
fi
if [ "$mem_np_gsi" -gt "$ppn" ]; then
  mem_nodes_gsi=$(((mem_np_gsi-1)/ppn+1))
  mem_np_gsi=$((ppn*mem_nodes_gsi))
else
  mem_nodes_gsi=1
fi

#-------------------------------------------------------------------------------
# Determine
#   $node_m[1,...,$MEMBER+2]
#   $name_m[1,...,$MEMBER+2]
# Create files
#   'machinefile'
#   'machinefile.gfs.XXX' (when $mem_nodes_gfs = 1)
#   'machinefile.gsi.XXX' (when $mem_nodes_gsi = 1)

mkdir -p $tmpnode
if [ "$WRITE_FILES" = '1' ]; then
  rm -f $tmpnode/machinefile*
fi

local m=0
local p=0
while [ "$m" -lt "$MEMBER" ]; do
  for i in `seq $ppn`; do
    for n in `seq $nnodes`; do
      m=$((m+1))
      p=$((p+1))
      if [ "$m" -le "$MEMBER" ]; then
        node_m[$m]=${node[$n]}
        name_m[$m]=`printf '%03d' $m`
        nodes_gfs_m[$m]=${node[$n]}
        nodes_gsi_m[$m]=${node[$n]}

        if [ "$mem_nodes_gfs" -eq 1 ] && [ "$WRITE_FILES" = '1' ]; then
          for ii in `seq $mem_np_gfs`; do
            echo ${node[$n]} >> $tmpnode/machinefile.gfs.${name_m[$m]}
          done
        fi
        if [ "$mem_nodes_gsi" -eq 1 ] && [ "$WRITE_FILES" = '1' ]; then
          for ii in `seq $mem_np_gsi`; do
            echo ${node[$n]} >> $tmpnode/machinefile.gsi.${name_m[$m]}
          done
        fi
      fi
      if [ "$p" -le "$totalnp" ] && [ "$WRITE_FILES" = '1' ]; then
        echo ${node[$n]} >> $tmpnode/machinefile
      fi
    done
  done
done

node_m[$mmean]=${node[1]}
name_m[$mmean]='mean'
node_m[$msprd]=${node[1]}
name_m[$msprd]='sprd'

#-------------------------------------------------------------------------------
# When $mem_nodes_gfs > 1, determine
#   $nodes_gfs_m[1,...,$MEMBER]
# Create file
#   'machinefile.gfs.XXX'

if [ "$mem_nodes_gfs" -gt 1 ]; then
  n=1
  for m in `seq $MEMBER`; do
    nodes_gfs_m[$m]=''
    n2=$((n+mem_nodes_gfs-1))
    if [ "$n2" -gt "$nnodes" ]; then
      n=1
      n2=$((n+mem_nodes_gfs-1))
    fi
    for nn in `seq $n $n2`; do
      nodes_gfs_m[$m]="${nodes_gfs_m[$m]}${node[$nn]} "
    done
    if [ "$WRITE_FILES" = '1' ]; then
      for i in `seq $ppn`; do
        for nn in `seq $n $n2`; do
          echo ${node[$nn]} >> $tmpnode/machinefile.gfs.${name_m[$m]}
        done
      done
    fi
    n=$((n+mem_nodes_gfs))
  done
fi

#-------------------------------------------------------------------------------
# When $mem_nodes_gsi > 1, determine
#   $nodes_gsi_m[1,...,$MEMBER]
# Create files
#   'machinefile.gsi.XXX'

if [ "$mem_nodes_gsi" -gt 1 ]; then
  n=1
  for m in `seq $MEMBER`; do
    nodes_gsi_m[$m]=''
    n2=$((n+mem_nodes_gfs-1))
    if [ "$n2" -gt "$nnodes" ]; then
      n=1
      n2=$((n+mem_nodes_gfs-1))
    fi
    for nn in `seq $n $n2`; do
      nodes_gsi_m[$m]="${nodes_gsi_m[$m]}${node[$nn]} "
    done
    if [ "$WRITE_FILES" = '1' ]; then
      for i in `seq $ppn`; do
        for nn in `seq $n $n2`; do
          echo ${node[$nn]} >> $tmpnode/machinefile.gsi.${name_m[$m]}
        done
      done
    fi
    n=$((n+mem_nodes_gsi))
  done
fi

#-------------------------------------------------------------------------------
# Determine
#   $nodes_gsi_m[$MEMBER+1]
# Create files
#   'machinefile.gsi.mean' (used to run GSI for ensemble mean)

if [ "$WRITE_FILES" = '1' ]; then
  if [ "$SHAREDISK" = '0' ]; then
    mean_np_gsi=$ppn
    if [ "$MAX_NP_GSI" -lt "$mean_np_gsi" ]; then
      mean_np_gsi=$MAX_NP_GSI
    fi
    nodes_gsi_m[$mmean]=${node_m[$mmean]}
    for i in `seq $mean_np_gsi`; do
      echo ${node_m[$mmean]} >> $tmpnode/machinefile.gsi.${name_m[$mmean]}
    done
  else
    mean_np_gsi=$totalnp
    if [ "$MAX_NP_GSI" -lt "$mean_np_gsi" ]; then
      mean_np_gsi=$MAX_NP_GSI
    fi
    if [ "$mean_np_gsi" -gt "$ppn" ]; then
      mean_nodes_gsi=$(((mean_np_gsi-1)/ppn+1))
      nodes_gsi_m[$mmean]=''
      for nn in `seq 1 $mean_nodes_gsi`; do
        nodes_gsi_m[$mmean]="${nodes_gsi_m[$mmean]}${node[$nn]} "
      done
      for i in `seq $ppn`; do
        for nn in `seq 1 $mean_nodes_gsi`; do
          echo ${node[$nn]} >> $tmpnode/machinefile.gsi.${name_m[$mmean]}
        done
      done
    else
      nodes_gsi_m[$mmean]=${node_m[$mmean]}
      for i in `seq $mean_np_gsi`; do
        echo ${node_m[$mmean]} >> $tmpnode/machinefile.gsi.${name_m[$mmean]}
      done
    fi
  fi
fi

#-------------------------------------------------------------------------------
# Create file
#   'machinefile.node'

if [ "$WRITE_FILES" = '1' ]; then
  for n in `seq $nnodes`; do
    echo ${node[$n]} >> $tmpnode/machinefile.node
  done
fi

#-------------------------------------------------------------------------------
# Display the current settings (when $DISPLAY = 1)

if [ "$DISPLAY" = '1' ]; then
  echo "  Nodes used:           $nnodes"
  for n in `seq $nnodes`; do
    echo "    ${node[$n]}"
  done
  echo
  echo "  Cores per node:       $ppn"
  echo "  Total cores:          $totalnp"
  echo
  if [ "$SHAREDISK" = '0' ]; then
    echo "  Local disks on nodes are used to store runtime files"
  else
    echo "  A share disk is used to store runtime files"
  fi
  echo
  if [ "$mem_nodes_gfs" -eq 1 ]; then
    echo "  Cores per GFS run:    $mem_np_gfs (1 node)"
  else
    echo "  Cores per GFS run:    $mem_np_gfs ($mem_nodes_gfs nodes)"
  fi
  if [ "$mem_nodes_gsi" -eq 1 ]; then
    echo "  Cores per GSI run:    $mem_np_gsi (1 node)"
  else
    echo "  Cores per GSI run:    $mem_np_gsi ($mem_nodes_gsi nodes)"
  fi
  echo
  echo "  Ensemble size:        $MEMBER"
  if [ "$SHAREDISK" = '0' ]; then
    for m in `seq $msprd`; do
      echo "    ${name_m[$m]} stored on ${node_m[$m]}"
    done
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================

function distribute_fcst {
#-------------------------------------------------------------------------------
# Distribute members on nodes for ensemble forecasts.
#
# Usage: distribute_fcst NODEFILE MEMBERS [WRITE_FILES] [DISPLAY] [CYCLES]
#
#   NODEFILE     machinefile for mpiexec
#   MEMBERS      List of forecast members
#   WRITE_FILES  output machinefiles or not
#                0: No
#                1: Yes
#                (default: Yes)
#   DISPLAY      display distribution information or not
#                0: No
#                1: Yes
#                (default: No)
#   CYCLES       Number of forecast cycles run in parallel
#                (default: 1)
#
# Output files:
#   [$TMPMPI/node/machinefile.node]          machinefile for mpiexec (one core per node)
#   [$TMPMPI/node/machinefile.gfs.YYYY.XXX]  machinefile of each member for mpiexec
#                                             (YYYY = 0001, 0002, 0003, ... : cycles)
#                                             (XXX = 001, 002, ..., mean : members)
#
# Return variables:
#   $nnodes
#   $ppn
#   $totalnp
#   $fmember
#   $fmembertot
#   $mem_nodes_gfs
#   $mem_np_gfs
#   $node   [1, ..., $nnodes    ]
#   $node_m [1, ..., $fmembertot]
#   $name_m [1, ..., $fmembertot]
#   $cycle_m[1, ..., $fmembertot]
#   $nodes_gfs_m[1, ..., $fmembertot]
#-------------------------------------------------------------------------------

if [ "$#" -lt 2 ]; then
  echo "[Error] function 'distribute_fcst': Insufficient arguments." 1>&2
  exit 1
fi

local NODEFILE="$1"
local MEMBERS="$2"
local WRITE_FILES=${3:-1}
local DISPLAY=${4:-0}
local CYCLES=${5:-1}

if [ ! -s "$NODEFILE" ]; then
  echo "[Error] $0: Can not find \$NODEFILE '$NODEFILE'" 1>&2
  exit 1
fi

tmpnode="$TMPMPI/node"

#-------------------------------------------------------------------------------
# Read $NODEFILE and determine
#   $nnodes
#   $ppn
#   $totalnp
#   $node[1,...,$nnodes]

local nodelist=`cat $NODEFILE | sort | uniq`

local n=0
for inode in $nodelist; do
  n=$((n+1))
  node[$n]=$inode

  ippn=`cat $NODEFILE | grep $inode | wc -l`
  if [ "$n" -eq 1 ]; then
    ppn=$ippn
  elif [ "$ippn" -ne "$ppn" ]; then
    echo "[Error] $0: The number of cores per node must be a constant" 1>&2
    exit 1
  fi
done
nnodes=$n
totalnp=$((ppn*nnodes))

#-------------------------------------------------------------------------------
# Determine
#   $fmember
#   $fmembertot

fmember=`echo $MEMBERS | wc -w`
fmembertot=$((fmember * CYCLES))

#-------------------------------------------------------------------------------
# If the minimum cores required to run a GFS member cannot be afforded, then exit.

if [ "$MIN_NP_GFS" -gt "$ppn" ] && [ "$SHAREDISK" = '0' ]; then
  echo "[Error] When using local disks, minimum cores required to run a GFS or GSI cannot exceed cores per node." 1>&2
  exit 1
fi
if [ "$MIN_NP_GFS" -gt "$totalnp" ]; then
  echo "[Error] Minimum cores required to run a GFS or GSI exceed total avaiable cores." 1>&2
  exit 1
fi

#-------------------------------------------------------------------------------
# Find optimal number of cores ($mem_np) based on the ensemble size ($MEMBER)
# and the total nodes ($nnodes) and cores per node ($ppn).
# The 'optimal' here means minimum $mem_np but still can occupy all avaiable cores at once.

if [ "$nnodes" -ge "$fmembertot" ]; then
  if [ "$SHAREDISK" = '1' ]; then
    mem_nodes=$((nnodes/fmembertot))
  else
    mem_nodes=1
  fi
  mem_np=$((ppn*mem_nodes))
else
  mem_nodes=1
  mempn=$(((fmembertot-1)/nnodes+1))
  if [ "$mempn" -gt "$ppn" ]; then
    mem_np=1
  else
    mem_np=$((ppn/mempn))
  fi
fi

#-------------------------------------------------------------------------------
# Limited to the minimum cores required to run a GFS member ($MIN_NP_GFS)
#        and the maximum cores to run a GFS member ($MAX_NP_GFS),
# and require it occupy full nodes if multiple nodes are used,
# determine the numbers of cores ($mem_np_gfs) and nodes ($mem_nodes_gfs) used to run a GFS member
#        from the optimal number of cores per member ($mem_np),

if [ "$mem_np" -lt "$MIN_NP_GFS" ]; then
  mem_np_gfs=$MIN_NP_GFS
else
  mem_np_gfs=$mem_np
fi
if [ "$mem_np_gfs" -gt "$MAX_NP_GFS" ]; then
  mem_np_gfs=$MAX_NP_GFS
fi
if [ "$mem_np_gfs" -gt "$ppn" ]; then
  mem_nodes_gfs=$(((mem_np_gfs-1)/ppn+1))
  mem_np_gfs=$((ppn*mem_nodes_gfs))
else
  mem_nodes_gfs=1
fi

#-------------------------------------------------------------------------------
# Determine
#   $node_m [1,...,$fmembertot]
#   $name_m [1,...,$fmembertot]
#   $cycle_m[1,...,$fmembertot]
# Create files
#   'machinefile'
#   'machinefile.gfs.YYYY.XXX' (when $mem_nodes_gfs = 1)

mkdir -p $tmpnode
if [ "$WRITE_FILES" = '1' ]; then
  rm -f $tmpnode/machinefile*
fi

local m=0
local p=0
local mm=0
local c=1
while [ "$m" -lt "$fmembertot" ]; do
  for i in `seq $ppn`; do
    for n in `seq $nnodes`; do
      m=$((m+1))
      p=$((p+1))
      mm=$((mm+1))
      if [ "$mm" -gt "$fmember" ]; then
        mm=1
        c=$((c+1))
      fi
      if [ "$m" -le "$fmembertot" ]; then
        node_m[$m]=${node[$n]}
        name_m[$m]=`echo $MEMBERS | cut -d ' ' -f$mm`
        if [ "${name_m[$m]}" != 'mean' ]; then
          name_m[$m]=`printf '%03d' $((10#${name_m[$m]}))`
        fi
        nodes_gfs_m[$m]=${node[$n]}
        cycle_m[$m]=`printf '%04d' $c`

        if [ "$mem_nodes_gfs" -eq 1 ] && [ "$WRITE_FILES" = '1' ]; then
          for ii in `seq $mem_np_gfs`; do
            echo ${node[$n]} >> $tmpnode/machinefile.gfs.${cycle_m[$m]}.${name_m[$m]}
          done
        fi
      fi
      if [ "$p" -le "$totalnp" ] && [ "$WRITE_FILES" = '1' ]; then
        echo ${node[$n]} >> $tmpnode/machinefile
      fi
    done
  done
done

#-------------------------------------------------------------------------------
# When $mem_nodes_gfs > 1, determine
#   $nodes_gfs_m[1,...,$MEMBER]
# Create file
#   'machinefile.gfs.XXX'

if [ "$mem_nodes_gfs" -gt 1 ]; then
  n=1
  for m in `seq $fmembertot`; do
    nodes_gfs_m[$m]=''
    n2=$((n+mem_nodes_gfs-1))
    if [ "$n2" -gt "$nnodes" ]; then
      n=1
      n2=$((n+mem_nodes_gfs-1))
    fi
    for nn in `seq $n $n2`; do
      nodes_gfs_m[$m]="${nodes_gfs_m[$m]}${node[$nn]} "
    done
    if [ "$WRITE_FILES" = '1' ]; then
      for i in `seq $ppn`; do
        for nn in `seq $n $n2`; do
          echo ${node[$nn]} >> $tmpnode/machinefile.gfs.${cycle_m[$m]}.${name_m[$m]}
        done
      done
    fi
    n=$((n+mem_nodes_gfs))
  done
fi

#-------------------------------------------------------------------------------
# Create file
#   'machinefile.node'

if [ "$WRITE_FILES" = '1' ]; then
  for n in `seq $nnodes`; do
    echo ${node[$n]} >> $tmpnode/machinefile.node
  done
fi

#-------------------------------------------------------------------------------

if [ "$DISPLAY" = '1' ]; then
  echo "  Nodes used:           $nnodes"
  for n in `seq $nnodes`; do
    echo "    ${node[$n]}"
  done
  echo
  echo "  Cores per node:       $ppn"
  echo "  Total cores:          $totalnp"
  echo
  if [ "$SHAREDISK" = '0' ]; then
    echo "  Local disks on nodes are used to store runtime files"
  else
    echo "  A share disk is used to store runtime files"
  fi
  echo
  if [ "$mem_nodes_gfs" -eq 1 ]; then
    echo "  Cores per GFS run:    $mem_np_gfs (1 node)"
  else
    echo "  Cores per GFS run:    $mem_np_gfs ($mem_nodes_gfs nodes)"
  fi
fi

#-------------------------------------------------------------------------------
}

#===============================================================================
