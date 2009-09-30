#!/bin/sh
LETKFDIR=/S/home00/G4015/y0266/enkf/afes/letkf
DATADIR=/S/data04/G4015/y0266/alera2
YMDH=1982010100
OBS=test
EXP=M10
rm -rf test.nqs2*
rm -rf Log
mkdir Log
cat << EOF > test.nqs2
#!/bin/sh
#PBS -q L
#PBS -b 1
#PBS -l elapstim_req=00:20:00
#PBS -l filecap_job=10gb
#PBS -I "${LETKFDIR}/letkf008.m01,ALL:./"
#PBS -I "${DATADIR}/data/fort.21,ALL:./"
#PBS -I "${DATADIR}/data/fort.41,ALL:./"
#PBS -I "${DATADIR}/${OBS}/obs/${YMDH}.dat,ALL:./obs01.dat"
#PBS -I "${DATADIR}/${OBS}/${EXP}/gues/%03{1-8}n/${YMDH}.grd,%{(0)@8}n:./gs01%03{1-8}n.grd"
#PBS -O "${DATADIR}/${OBS}/${EXP}/gues/mean/${YMDH}.grd,0:./gues_me.grd"
#PBS -O "${DATADIR}/${OBS}/${EXP}/gues/sprd/${YMDH}.grd,0:./gues_sp.grd"
#PBS -O "${DATADIR}/${OBS}/${EXP}/anal/%03{1-8}n/${YMDH}.grd,%{(0)@8}n:./anal%03{1-8}n.grd"
#PBS -O "${DATADIR}/${OBS}/${EXP}/anal/mean/${YMDH}.grd,0:./anal_me.grd"
#PBS -O "${DATADIR}/${OBS}/${EXP}/anal/sprd/${YMDH}.grd,0:./anal_sp.grd"
#PBS -O "${LETKFDIR}/run/Log/,%{(0)@8}n:./NOUT-%03{0-7}n"
#PBS -v F_SETBUF06=0
#PBS -v MPIPROGINF=ALL_DETAIL
#PBS -v F_PROGINF=DETAIL
#PBS -v F_FTRACE=YES

set -vx
ls -l
mpirun -nnp 8 letkf008.m01
ls -l
EOF
qsub test.nqs2
