#!/bin/bash
#===============================================================================
#
#  Download NCEP conventional observation data from UCAR/DSS
#   -- adapted from Takemasa Miyoshi's LETKF google code,
#      March 2013,             Guo-Yuan Lien
#      January 2015, modified, Guo-Yuan Lien
#
#===============================================================================

cd "$(dirname "$0")"
myname=$(basename "$0")
myname1=${myname%.*}

#===============================================================================
# Configuration

. config.main
res=$? && ((res != 0)) && exit $res

. src/func_datetime.sh
. src/func_util.sh

#-------------------------------------------------------------------------------

USAGE="
[$myname] Download NCEP conventional observation data.

Usage: $myname EMAIL PASSWD STIME [ETIME] [IF_DECODE]

  EMAIL      UCAR/DSS account (register at http://rda.ucar.edu )
  PASSWD     UCAR/DSS account password
  STIME      Start time (format: YYYYMMDD)
  ETIME      End   time (format: YYYYMMDD)
             (default: same as STIME)
  IF_DECODE  Convert PREPBUFR to LETKF obs format or not
             0: No,  store only PREPBUFR format
             1: Yes, decode to LETKF obs format
             (default: Yes)
"

#-------------------------------------------------------------------------------

if [ "$1" == '-h' ] || [ "$1" == '--help' ]; then
  echo "$USAGE"
  exit 0
fi
if (($# < 3)); then
  echo "$USAGE" >&2
  exit 1
fi

EMAIL="$1"
PASSWD="$2"
STIME=$(datetime $3)
ETIME=$(datetime ${4:-$STIME})
IF_DECODE=${5:-1}

#-------------------------------------------------------------------------------

DATAURL="http://rda.ucar.edu/data/ds337.0/tarfiles"
LOGINURL="https://rda.ucar.edu/cgi-bin/login"
WGET="wget --no-check-certificate"
AUTH="auth.rda_ucar_edu"
OPTLOGIN="-O /dev/null --save-cookies $AUTH --post-data=\"email=${EMAIL}&passwd=${PASSWD}&action=login\""
OPT="-N --load-cookies $AUTH"

#===============================================================================

safe_init_tmpdir $TMPS
cd $TMPS

if ((IF_DECODE == 1)); then
  cp $OBSUTIL_DIR/dec_prepbufr .
fi

mkdir -p download
cd download
$WGET $OPTLOGIN $LOGINURL

time=$STIME
while ((time <= ETIME)); do

  yyyy=${time:0:4}
  mm=${time:4:2}
  dd=${time:6:2}
  if (("$yyyy$mm" >= 200807)); then
    DATAF="prepbufr.$yyyy$mm$dd.nr.tar.gz"
  else
    DATAF="prepbufr.$yyyy$mm$dd.wo40.tar.gz"
  fi

  cd $TMPS/download
  #$WGET $OPT "${DATAURL}/${DATAF}"
  $WGET $OPT "${DATAURL}/${yyyy}/${DATAF}"
  tar xzf $DATAF
  rm -f $DATAF

  cd $TMPS
  for hh in '00' '06' '12' '18'; do
    timef="$yyyy$mm$dd${hh}0000"
    echo
    echo "[${timef}]"
    echo

    if (("$yyyy$mm" >= 200807)); then
      mv -f download/$yyyy$mm$dd.nr/prepbufr.gdas.$yyyy$mm$dd.t${hh}z.nr \
            prepbufr.gdas.${timef}.nr
    else
      mv -f download/$yyyy$mm$dd.wo40/prepbufr.gdas.$yyyy$mm$dd$hh.wo40 \
            prepbufr.gdas.${timef}.nr
    fi
    wc -c prepbufr.gdas.${timef}.nr | $BUFRBIN/grabbufr prepbufr.gdas.${timef}.nr prepbufr.in

    if ((IF_DECODE == 1)); then
      time ./dec_prepbufr
      touch fort.90
      mkdir -p $OBS
      mv fort.90 $OBS/obs_${timef}.dat
    fi

    mv -f prepbufr.gdas.${timef}.nr $OBSNCEP/obs_${timef}
  done

time=$(datetime $time 1 d)
done

safe_rm_tmpdir $TMPS

#===============================================================================

echo

exit 0
