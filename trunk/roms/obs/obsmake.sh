#!/bin/sh
set -e
CDIR=`pwd`
cd ..
ROMS=`pwd`
TRUEDIR=$ROMS/DATA/nature
OBSDIR=$ROMS/DATA/obs
PGM=obsmake.s01
# Initial time
IT=0001
# Final time
ET=0120
# start
cd $CDIR
mkdir -p $OBSDIR
ln -s $ROMS/model/bc/ee6_grd.nc grd.nc
# main loop
while test $IT -le $ET
do
ln -s $TRUEDIR/$IT.nc true.nc

./$PGM

mv obs.dat $OBSDIR/$IT.dat
rm true.nc

IT=`expr $IT + 1`
if test $IT -lt 1000
then
IT=0$IT
fi
if test $IT -lt 100
then
IT=0$IT
fi
if test $IT -lt 10
then
IT=0$IT
fi
done
rm grd.nc
echo "NORMAL END"
