#!/bin/sh
#PJM -N cycle_d3_100m_100mem_225
#PJM -s
#PJM --rsc-list "node=72720"
#PJM --rsc-list "elapse=01:20:00"
#PJM --rsc-list "rscgrp=huge"
##PJM --rsc-list "node-quota=29G"
##PJM --mpi "shape=72720"
#PJM --mpi "proc=72720"
#PJM --mpi assign-online-node
#PJM --stg-transfiles all
#PJM --mpi "use-rankdir"

#PJM --stgin  "rank=* ./input_tar/%06r.tar %r:./input.tar"

#PJM --stgout-dir "rank=* %r:./out/20130713060030/gues /volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/set_3nest_d3_100m_100mem_720p/d3_100m_100mem_225/20130713060030/gues recursive=5"
#PJM --stgout-dir "rank=* %r:./out/20130713060030/anal /volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/set_3nest_d3_100m_100mem_720p/d3_100m_100mem_225/20130713060030/anal recursive=5"
#PJM --stgout-dir "rank=* %r:./out/20130713060000/log/scale /volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/set_3nest_d3_100m_100mem_720p/d3_100m_100mem_225/20130713060000/log/scale recursive=3"
#PJM --stgout-dir "rank=* %r:./out/20130713060030/log/obsope /volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/set_3nest_d3_100m_100mem_720p/d3_100m_100mem_225/20130713060030/log/obsope recursive=3"
#PJM --stgout-dir "rank=* %r:./out/20130713060030/log/letkf /volume63/data/hp150019/gylien/exp/scale-letkf/BDA_case130713/set_3nest_d3_100m_100mem_720p/d3_100m_100mem_225/20130713060030/log/letkf recursive=3"

. /work/system/Env_base_1.2.0-20-1
export OMP_NUM_THREADS=8
export PARALLEL=8

mpiexec -of-proc tmp /work/system/bin/msh "/bin/tar -xf ./input.tar"

./cycle.sh "20130713060000" "20130713060000" "all" "1" "5" || exit $?
