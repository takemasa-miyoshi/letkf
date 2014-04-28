#!/bin/sh
### Running on LINUX
#
# DATA SOURCE
#
DATASRC=2 #1:NOAA/NOMADS, 2:UCAR/DSS
#
# COMPILER
#
F90=pgf90
FFLAGS="-byteswapio" #big-endian IO
#
# TOOLS: grabbufr, BUFRLIB
#
GRABBUFR=$HOME/tools/grabbufr
BUFRLIB=$HOME/tools/BUFRLIB
# Initial time
IY=2008
IM=09
ID=01
IH=00
# Final time
EY=2008
EM=09
ED=01
EH=00
#
# DIR
#
CDIR=`pwd`
cd ..
WRF=`pwd`
OBSDIR=$WRF/DATA/obs
DLDIR=$OBSDIR/downloads
SAVEDIR=$OBSDIR/prepbufr
#
source $WRF/../common/timeinc.sh
#
# WKDIR
#
WKDIR=$OBSDIR/tmp
rm -rf $WKDIR
mkdir -p $WKDIR
cd $WKDIR
#
# Settings for DATASRC
#
if test $DATASRC -eq 1
then
DATAURL="http://nomads.ncdc.noaa.gov/data/gdas"
WGET="wget"
DLDIR=$OBSDIR/downloads/nomads
SAVEDIR=$OBSDIR/prepbufr/nomads
mkdir -p $DLDIR
mkdir -p $SAVEDIR
else
EMAIL=$1
PASSWD=$2
if test $# -ne 2
then
echo "USAGE: $0 EMAIL PASSWD"
exit
fi
DATAURL="https://dss.ucar.edu/datazone/dsszone/ds337.0"
WGET="wget --no-check-certificate"
DLDIR=$OBSDIR/downloads/ucar
SAVEDIR=$OBSDIR/prepbufr/ucar
mkdir -p $DLDIR
mkdir -p $SAVEDIR
cd $DLDIR
rm -f auth.dss_ucar_edu
OPT1='-O /dev/null --save-cookies auth.dss_ucar_edu --post-data'
OPT2="email=${EMAIL}&passwd=${PASSWD}&action=login"
$WGET $OPT1="$OPT2" https://dss.ucar.edu/cgi-bin/login
OPT="-N --load-cookies auth.dss_ucar_edu"
cd $WKDIR
fi
#
# build grabbufr
#
cp $GRABBUFR/grabbufr.f .
cp $GRABBUFR/spbufr.f .
$F90 -o grabbufr grabbufr.f spbufr.f
#
# build decoder
#
cp $WRF/../common/SFMT.f90 .
cp $WRF/../common/common.f90 .
cp $WRF/common/common_wrf.f90 .
cp $WRF/common/common_obs_wrf.f90 .
cp $CDIR/dec_prepbufr.f90 .
$F90 $FFLAGS -o decoder SFMT.f90 common.f90 common_wrf.f90 common_obs_wrf.f90 dec_prepbufr.f90 -L$BUFRLIB -lbufrlib
#
# Cycle run ### MAIN LOOP ###
#
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
echo " >>>"
echo " >>> BEGIN COMPUTATION OF $IY/$IM/$ID/$IH"
echo " >>>"
#
# GET NCEP PREPBUFR
#
cd $DLDIR
if test $DATASRC -eq 1
then
DATAF=$IY$IM$ID$IH.prepbufr.nr
else
if test $IY$IM -ge 200807
then
NR=nr
DATAF=prepbufr.gdas.$IY$IM$ID.t${IH}z.$NR
else
NR=wo40
DATAF=prepbufr.gdas.$IY$IM$ID$IH.$NR
fi
fi
if test ! -f $DATAF
then
if test $DATASRC -eq 1
then
$WGET $DATAURL/$IY$IM/$IY$IM$ID/gdas1.t${IH}z.prepbufr.nr
mv gdas1.t${IH}z.prepbufr.nr $IY$IM$ID$IH.prepbufr.nr
else
$WGET $OPT $DATAURL/prepbufr.$IY$IM$ID.$NR.tar.gz
tar zxvf prepbufr.$IY$IM$ID.$NR.tar.gz
rm prepbufr.$IY$IM$ID.$NR.tar.gz
mv $IY$IM$ID.$NR/* .
rmdir $IY$IM$ID.$NR
rm -rf $IY$IM$ID.$NR.tar.gz
fi
fi
#
# DECODE
#
cd $WKDIR
ln -s $DLDIR/$DATAF prepbufr.tmp
wc -c prepbufr.tmp | ./grabbufr prepbufr.tmp prepbufr.in
time ./decoder
#
# SAVE
#
touch fort.87
touch fort.88
touch fort.89
touch fort.90
touch fort.91
touch fort.92
touch fort.93
mkdir -p $SAVEDIR/obs$IY$IM$ID$IH
mv fort.87 $SAVEDIR/obs$IY$IM$ID$IH/t-3.dat
mv fort.88 $SAVEDIR/obs$IY$IM$ID$IH/t-2.dat
mv fort.89 $SAVEDIR/obs$IY$IM$ID$IH/t-1.dat
mv fort.90 $SAVEDIR/obs$IY$IM$ID$IH/t.dat
mv fort.91 $SAVEDIR/obs$IY$IM$ID$IH/t+1.dat
mv fort.92 $SAVEDIR/obs$IY$IM$ID$IH/t+2.dat
mv fort.93 $SAVEDIR/obs$IY$IM$ID$IH/t+3.dat
#
# Date change ### MAIN LOOP END ###
#
TY=`timeinc6hr $IY $IM $ID $IH | cut -c1-4`
TM=`timeinc6hr $IY $IM $ID $IH | cut -c5-6`
TD=`timeinc6hr $IY $IM $ID $IH | cut -c7-8`
TH=`timeinc6hr $IY $IM $ID $IH | cut -c9-10`
IY=$TY
IM=$TM
ID=$TD
IH=$TH
done
rm -f $DLDIR/auth.dss_ucar_edu
echo "NORMAL END"
