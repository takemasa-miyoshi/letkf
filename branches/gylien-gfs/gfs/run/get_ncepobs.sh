#!/bin/bash
#===============================================================================
#
#  Download NCEP conventional observation data from UCAR/DSS
#   -- adapted from Takemasa Miyoshi's LETKF google code,
#      March 2013, Guo-Yuan Lien
#
#===============================================================================

if [ -f configure.sh ]; then
  . configure.sh
else
  echo "[Error] $0: 'configure.sh' does not exist." 1>&2
  exit 1
fi
. datetime.sh

if [ "$#" -lt 3 ]; then
  cat 1>&2 << EOF

[get_ncepobs.sh] Download NCEP conventional observation data.
                 *use settings in 'configure.sh'

Usage: $0 EMAIL PASSWD STIME [ETIME] [IF_DECODE]

  EMAIL      UCAR/DSS account (register at http://rda.ucar.edu )
  PASSWD     UCAR/DSS account password
  STIME      Start time (format: YYYYMMDD)
  ETIME      End   time (format: YYYYMMDD)
             (default: same as STIME)
  IF_DECODE  Convert PREPBUFR to LETKF obs format or not
             0: No,  store only PREPBUFR format
             1: Yes, decode to LETKF obs format
             (default: Yes)

EOF
  exit 1
fi

EMAIL="$1"
PASSWD="$2"
STIME=$(datetime $3)
ETIME=$(datetime ${4:-$STIME})
IF_DECODE=${5:-1}

DATAURL="http://rda.ucar.edu/data/ds337.0/tarfiles"
LOGINURL="https://rda.ucar.edu/cgi-bin/login"
WGET="wget --no-check-certificate"
AUTH="auth.rda_ucar_edu"
OPTLOGIN="-O /dev/null --save-cookies $AUTH --post-data=\"email=${EMAIL}&passwd=${PASSWD}&action=login\""
OPT="-N --load-cookies $AUTH"

tmpsubdir="gfs-letkf_${USER}_${SYSNAME}_get_ncepobs"
tmprun="$TMP1/${tmpsubdir}"

#===============================================================================

mkdir -p $tmprun
rm -fr $tmprun/*
mkdir -p $tmprun/download
cd $tmprun
if [ "$IF_DECODE" = '1' ]; then
  cp $DIR/obs/dec_prepbufr .
fi

cd $tmprun/download
$WGET $OPTLOGIN $LOGINURL

time=$STIME
while [ "$time" -le "$ETIME" ]; do

  yyyy=${time:0:4}
  mm=${time:4:2}
  dd=${time:6:2}
  if [ "$yyyy$mm" -ge 200807 ]; then
    DATAF="prepbufr.$yyyy$mm$dd.nr.tar.gz"
  else
    DATAF="prepbufr.$yyyy$mm$dd.wo40.tar.gz"
  fi

  cd $tmprun/download
  $WGET $OPT "${DATAURL}/${DATAF}"
  tar xzf $DATAF
  rm -f $DATAF

  cd $tmprun
  for hh in '00' '06' '12' '18'; do
    timef="$yyyy$mm$dd$hh"
    echo
    echo "[${timef}]"
    echo

    if [ "$yyyy$mm" -ge 200807 ]; then
      mv -f download/$yyyy$mm$dd.nr/prepbufr.gdas.$yyyy$mm$dd.t${hh}z.nr \
            prepbufr.gdas.${timef}.nr
    else
      mv -f download/$yyyy$mm$dd.wo40/prepbufr.gdas.${timef}.wo40 \
            prepbufr.gdas.${timef}.nr
    fi
    wc -c prepbufr.gdas.${timef}.nr | $BUFRBIN/grabbufr prepbufr.gdas.${timef}.nr prepbufr.in

    if [ "$IF_DECODE" = '1' ]; then
      time ./dec_prepbufr
      touch fort.87
      touch fort.88
      touch fort.89
      touch fort.90
      touch fort.91
      touch fort.92
      touch fort.93
      mkdir -p $OBS/obs${timef}
      mv fort.87 $OBS/obs${timef}/t-3.dat
      mv fort.88 $OBS/obs${timef}/t-2.dat
      mv fort.89 $OBS/obs${timef}/t-1.dat
      mv fort.90 $OBS/obs${timef}/t.dat
      mv fort.91 $OBS/obs${timef}/t+1.dat
      mv fort.92 $OBS/obs${timef}/t+2.dat
      mv fort.93 $OBS/obs${timef}/t+3.dat
    fi

    mkdir -p $OBSNCEP/obs${timef}
    mv -f prepbufr.gdas.${timef}.nr $OBSNCEP/obs${timef}
#    mv -f prepbufr.gdas.${timef}.nr $OBSNCEP/obs${timef}/gdas1.t${hh}z.prepbufr.nr
  done

time=$(datetime $time 1 d)
done

rm -rf $tmprun

#===============================================================================

exit 0
