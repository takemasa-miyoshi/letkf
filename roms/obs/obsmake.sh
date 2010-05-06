#!/bin/sh
set -e
CDIR=`pwd`
cd ..
ROMS=`pwd`
source $ROMS/../common/timeinc.sh
TRUEDIR=$ROMS/DATA/nature
EXP=5x20
ASCIIOBS=$ROMS/DATA/$EXP/insitu_obs
LETKFOBS=$ROMS/DATA/$EXP/letkf_obs
PGM=obsmake.s01
# Initial time
IY=2004
IM=01
ID=01
IH=00
# Final time
EY=2004
EM=03
ED=01
EH=00
# start
cd $CDIR
mkdir -p $ASCIIOBS
mkdir -p $LETKFOBS
for OBTYPE in hfradar_uv ship_sst profiles
do
mkdir -p $ASCIIOBS/$OBTYPE
mkdir -p $LETKFOBS/$OBTYPE
done
cp station.tbl $ASCIIOBS
cp obserr.tbl $ASCIIOBS
cp station.tbl $LETKFOBS
cp obserr.tbl $LETKFOBS
ln -fs $ROMS/model/bc/ee6_grd.nc grd.nc
# main loop
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
echo "PROCESSING $IY/$IM/$ID/$IH"
ln -fs $TRUEDIR/$IY$IM$ID${IH}_rst.nc true.nc

./$PGM

mv hrf.ascii $ASCIIOBS/hfradar_uv/$IY$IM$ID${IH}00_hfr.dat
mv hrf.letkf $LETKFOBS/hfradar_uv/$IY$IM$ID${IH}_obs.dat
mv sst.ascii $ASCIIOBS/ship_sst/$IY$IM$ID${IH}.asc
mv sst.letkf $LETKFOBS/ship_sst/$IY$IM$ID${IH}_obs.dat
mv prof.ascii $ASCIIOBS/profiles/fort.$IY$IM$ID${IH}
mv prof.letkf $LETKFOBS/profiles/$IY$IM$ID${IH}_obs.dat
rm true.nc
### 6hr increment
TY=`timeinc6hr $IY $IM $ID $IH | cut -c1-4`
TM=`timeinc6hr $IY $IM $ID $IH | cut -c5-6`
TD=`timeinc6hr $IY $IM $ID $IH | cut -c7-8`
TH=`timeinc6hr $IY $IM $ID $IH | cut -c9-10`
IY=$TY
IM=$TM
ID=$TD
IH=$TH

done
rm grd.nc
echo "NORMAL END"
