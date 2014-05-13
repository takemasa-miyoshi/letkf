#!/bin/bash
#===============================================================================
#
#  Script to run GSI in observer mode.
#  August    2013  created, Runhua Yang
#  September 2013  modified, Guo-Yuan Lien
#
#===============================================================================

if [ "$#" -lt 6 ]; then
  cat 1>&2 << EOF

[run_gsi.sh] Prepare a temporary directory for GSI run,
             and may run the model.

Usage: $0 TMPDIR EXECDIR FIXDIR FIXCRTMDIR GFSRUNDIR OBSDIR
          [MODE] [WINDOW] [THIN_OPT] [MPIEXEC] [NP] [MACHINE]

  TMPDIR      Temporary directory to run GSI
  EXECDIR     Directory of GSI executable files
  FIXDIR      Directory of GSI fix files
  FIXCRTMDIR  Directory of CRTM fix files
  GFSRUNDIR   Directory of GFS sig/sfc outputs (used as the background)
  OBSDIR      Directory of observation files
              (MODE=1: NCEP BUFR; MODE=2: obs_input* from a GSI run)
  MODE        GSI observer mode
              1: Run for ensemble mean
              2: Run for each members
              (default: 1)
  WINDOW      Assimilation window (+/- hour)
              (default: 3.0)
  THIN_OPT    Thinning option (refer to the definition in configure.sh)
              (default: 3)
  MPIEXEC     Directory of the mpiexec program
              'no': do not run GFS model (just prepare the temporary directory)
              'se': run GFS model serially
              (default: 'no')
  NP          Number of cores
  MACHINE     machinefile for mpiexec

EOF
  exit 1
fi

TMPDIR="$1"
EXECDIR="$2"
FIXDIR="$3"
FIXCRTMDIR="$4"
GFSRUNDIR="$5"
OBSDIR="$6"
MODE="${7:-1}"
WINDOW="${8:-3.0}"
THIN_OPT="${9:-3}"
MPIEXEC="${10:-no}"
NP="${11:-1}"
MACHINE="$12"

GSIEXEC="${EXECDIR}/global_gsi"
SIGHDR="${EXECDIR}/global_sighdr"

#-------------------------------------------------------------------------------

nslots=7
for is in `seq $nslots`; do
  fh=$((3+is-1))
  fhh_sl[$is]=`printf '%02d' $fh`
  if [ "$fh" -eq 6 ]; then
    baseslot=$is
  fi
done

#===============================================================================

sigbase="${GFSRUNDIR}/SIG.F${fhh_sl[$baseslot]}"

JCAP=`echo jcap | $SIGHDR $sigbase`
LEVS=`echo levs | $SIGHDR $sigbase`
LONR=`echo lonr | $SIGHDR $sigbase`
LATR=`echo latr | $SIGHDR $sigbase`
JCAP_B=$JCAP
NLON_A=$LONR
NLAT_A=$(($LATR+2))
IGEN=0
DELTIM=$((3600/(JCAP/20)))
ENDIANNESS=Big_Endian

#----------------------------------------------------------------
# Following part does not change with ensemble members
#----------------------------------------------------------------
# Set fixed files
#   berror   = forecast model background error statistics
#   specoef  = CRTM spectral coefficients
#   trncoef  = CRTM transmittance coefficients
#   emiscoef = CRTM coefficients for IR sea surface emissivity model
#   aerocoef = CRTM coefficients for aerosol effects
#   cldcoef  = CRTM coefficients for cloud effects
#   satinfo  = text file with information about assimilation of brightness temperatures
#   satangl  = angle dependent bias correction file (fixed in time)
#   pcpinfo  = text file with information about assimilation of prepcipitation rates
#   ozinfo   = text file with information about assimilation of ozone data
#   errtable = text file with obs error for conventional data (optional)
#   convinfo = text file with information about assimilation of conventional data
#   bufrtable= text file ONLY needed for single obs test (oneobstest=.true.)
#   bftab_sst= bufr table for sst ONLY needed for sst retrieval (retrieval=.true.)

anavinfo=${FIXDIR}/global_anavinfo.l64.txt
berror=${FIXDIR}/${ENDIANNESS}/global_berror.l${LEVS}y${NLAT_A}.f77
emiscoef=${FIXCRTMDIR}/EmisCoeff/${ENDIANNESS}/EmisCoeff.bin
aercoef=${FIXCRTMDIR}/AerosolCoeff/${ENDIANNESS}/AerosolCoeff.bin
cldcoef=${FIXCRTMDIR}/CloudCoeff/${ENDIANNESS}/CloudCoeff.bin
satinfo=${FIXDIR}/global_satinfo.txt
scaninfo=${FIXDIR}/global_scaninfo.txt
satangl=${FIXDIR}/global_satangbias.txt
pcpinfo=${FIXDIR}/global_pcpinfo.txt
ozinfo=${FIXDIR}/global_ozinfo.txt
if [ "$THIN_OPT" = '4' ]; then
  convinfo=${FIXDIR}/global_convinfo.thin.t${JCAP}.txt
else
  convinfo=${FIXDIR}/global_convinfo.txt
fi
atmsbeamdat=${FIXDIR}/atms_beamwidth.txt
errtable=${FIXDIR}/prepobs_errtable.global
bufrtable=${FIXDIR}/prepobs_prep.bufrtable
bftab_sst=${FIXDIR}/bufrtab.012


#expid=${expnm}.enmn.clnscript
#datges=$TOPDIR/guess
#tmpdir=$WORKDIR/tmp_sigmap/${expid}
#savdir=$WORKDIR/sigmap/${expid}
#SELECTOBS=$tmpdir

#rm -rf $tmpdir
#mkdir -p $tmpdir
#cd $tmpdir
#rm -rf core*

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*
cd $TMPDIR

# Copy executable and fixed files to $TMPDIR
ln -fs $GSIEXEC .

ln -fs $anavinfo anavinfo
ln -fs $berror   berror_stats
ln -fs $emiscoef EmisCoeff.bin
ln -fs $aercoef  AerosolCoeff.bin
ln -fs $cldcoef  CloudCoeff.bin
ln -fs $satangl  satbias_angle
ln -fs $satinfo  satinfo
ln -fs $scaninfo scaninfo
ln -fs $pcpinfo  pcpinfo
ln -fs $ozinfo   ozinfo
ln -fs $convinfo convinfo
ln -fs $atmsbeamdat atms_beamwidth.txt
ln -fs $errtable errtable
ln -fs $bufrtable prepobs_prep.bufrtable
ln -fs $bftab_sst bftab_sstphr

# Copy CRTM coefficient files based on entries in satinfo file
#nsatsen=`cat satinfo | $wc -l`
nsatsen=0                        #do not use satellite 
isatsen=1
while [ "$isatsen" -le "$nsatsen" ]; do
  flag=`head -n $isatsen satinfo | tail -1 | cut -c1-1`
  if [ "$flag" != "!" ]; then
    satsen=`head -n $isatsen satinfo | tail -1 | cut -f 2 -d ' '`
    if [ ! -s "${satsen}.SpcCoeff.bin" ]; then
      ln -fs ${FIXCRTMDIR}/SpcCoeff/${ENDIANNESS}/${satsen}.SpcCoeff.bin .
      ln -fs ${FIXCRTMDIR}/TauCoeff/${ENDIANNESS}/${satsen}.TauCoeff.bin .
    fi
  fi
  isatsen=$((isatsen+1))
done

PREPB_SATWND=".true."

# use prepbufr only
if [ "$MODE" -eq 1 ]; then
  ln -fs ${OBSDIR}/prepbufr.gdas.nr prepbufr
else
  ln -fs ${OBSDIR}/obs_input.* .
fi

# Copy bias correction, atmospheric and surface files
#ln -fs $datges/${prefix_tbc}.abias              ./satbias_in
#ln -fs $datges/${prefix_tbc}.satang             ./satbias_angle

for is in `seq $nslots`; do
  if [ -s "${GFSRUNDIR}/SIG.F${fhh_sl[$is]}" ] && [ -s "${GFSRUNDIR}/SFC.F${fhh_sl[$is]}" ]; then
    ln -fs ${GFSRUNDIR}/SIG.F${fhh_sl[$is]} sigf${fhh_sl[$is]}
    ln -fs ${GFSRUNDIR}/SFC.F${fhh_sl[$is]} sfcf${fhh_sl[$is]}
  fi
done

#===============================================================================

cat << EOF > gsiparm.anl
 &SETUP
   miter=0,niter(1)=0,niter(2)=0,
EOF
if [ "$MODE" -eq 1 ]; then
  echo "   lread_obs_save=.true.,lread_obs_skip=.false.," >> gsiparm.anl
else
  echo "   lread_obs_save=.false.,lread_obs_skip=.true.," >> gsiparm.anl
fi
cat << EOF >> gsiparm.anl
   lwrite_predterms=.true.,lwrite_peakwt=.true.,reduce_diag=.true.,passive_bc=.false.,
   niter_no_qc(1)=50,niter_no_qc(2)=0,
   write_diag(1)=.true.,write_diag(2)=.false.,write_diag(3)=.false.,
   qoption=2,
   gencode=$IGEN,factqmin=5.0,factqmax=5.0,deltim=$DELTIM,
   ndat=64,iguess=-1,
   oneobtest=.false.,retrieval=.false.,l_foto=.false.,
   use_pbl=.false.,use_compress=.true.,nsig_ext=12,gpstop=50.,
   use_gfs_nemsio=.false.,
   use_prepb_satwnd=$PREPB_SATWND,
 /
 &GRIDOPTS
   JCAP_B=$JCAP_B,JCAP=$JCAP,NLAT=$NLAT_A,NLON=$NLON_A,nsig=$LEVS,hybrid=.true.,
   regional=.false.,nlayers(63)=3,nlayers(64)=6,
 /
 &BKGERR
   vs=0.7,
   hzscl=1.7,0.8,0.5,
   hswgt=0.45,0.3,0.25,
   bw=0.0,norsp=4,
   bkgv_flowdep=.false.,bkgv_rewgtfct=1.5,
 /
 &ANBKGERR
   anisotropic=.false.,
 /
 &JCOPTS
   ljcdfi=.false.,alphajc=0.0,ljcpdry=.true.,bamp_jcpdry=5.0e7,
 /
 &STRONGOPTS
   jcstrong=.false.,nstrong=0,nvmodes_keep=0,period_max=0.,period_width=0,
   jcstrong_option=2,baldiag_full=.false.,baldiag_inc=.false.,
 /
 &OBSQC
   dfact=0.75,dfact1=3.0,noiqc=.true.,oberrflg=.false.,c_varqc=0.02,
   use_poq7=.true.,qc_noirjaco3_pole=.true.,
   tcp_width=60.0,tcp_ermin=2.0,tcp_ermax=12.0,
 /
 &OBS_INPUT
   dmesh(1)=145.0,dmesh(2)=150.0,time_window_max=${WINDOW},
   dfile(01)='prepbufr',  dtype(01)='ps',        dplat(01)=' ',       dsis(01)='ps',                 dval(01)=0.0, dthin(01)=0, dsfcalc(01)=0,
   dfile(02)='prepbufr'   dtype(02)='t',         dplat(02)=' ',       dsis(02)='t',                  dval(02)=0.0, dthin(02)=0, dsfcalc(02)=0,
   dfile(03)='prepbufr',  dtype(03)='q',         dplat(03)=' ',       dsis(03)='q',                  dval(03)=0.0, dthin(03)=0, dsfcalc(03)=0,
   dfile(04)='prepbufr',  dtype(04)='pw',        dplat(04)=' ',       dsis(04)='pw',                 dval(04)=0.0, dthin(04)=0, dsfcalc(04)=0,
   dfile(05)='satwndbufr',dtype(05)='uv',        dplat(05)=' ',       dsis(05)='uv',                 dval(05)=0.0, dthin(05)=0, dsfcalc(05)=0,
   dfile(06)='prepbufr',  dtype(06)='uv',        dplat(06)=' ',       dsis(06)='uv',                 dval(06)=0.0, dthin(06)=0, dsfcalc(06)=0,
   dfile(07)='prepbufr',  dtype(07)='spd',       dplat(07)=' ',       dsis(07)='spd',                dval(07)=0.0, dthin(07)=0, dsfcalc(07)=0,
   dfile(08)='prepbufr',  dtype(08)='dw',        dplat(08)=' ',       dsis(08)='dw',                 dval(08)=0.0, dthin(08)=0, dsfcalc(08)=0,
   dfile(09)='radarbufr', dtype(09)='rw',        dplat(09)=' ',       dsis(09)='rw',                 dval(09)=0.0, dthin(09)=0, dsfcalc(09)=0,
   dfile(10)='prepbufr',  dtype(10)='sst',       dplat(10)=' ',       dsis(10)='sst',                dval(10)=0.0, dthin(10)=0, dsfcalc(10)=0,
   dfile(11)='gpsrobufr', dtype(11)='gps_bnd',   dplat(11)=' ',       dsis(11)='gps',                dval(11)=0.0, dthin(11)=0, dsfcalc(11)=0,
   dfile(12)='ssmirrbufr',dtype(12)='pcp_ssmi',  dplat(12)='dmsp',    dsis(12)='pcp_ssmi',           dval(12)=0.0, dthin(12)=-1,dsfcalc(12)=0,
   dfile(13)='tmirrbufr', dtype(13)='pcp_tmi',   dplat(13)='trmm',    dsis(13)='pcp_tmi',            dval(13)=0.0, dthin(13)=-1,dsfcalc(13)=0,
   dfile(14)='sbuvbufr',  dtype(14)='sbuv2',     dplat(14)='n16',     dsis(14)='sbuv8_n16',          dval(14)=0.0, dthin(14)=0, dsfcalc(14)=0,
   dfile(15)='sbuvbufr',  dtype(15)='sbuv2',     dplat(15)='n17',     dsis(15)='sbuv8_n17',          dval(15)=0.0, dthin(15)=0, dsfcalc(15)=0,
   dfile(16)='sbuvbufr',  dtype(16)='sbuv2',     dplat(16)='n18',     dsis(16)='sbuv8_n18',          dval(16)=0.0, dthin(16)=0, dsfcalc(16)=0,
   dfile(17)='hirs3bufr', dtype(17)='hirs3',     dplat(17)='n17',     dsis(17)='hirs3_n17',          dval(17)=0.0, dthin(17)=1, dsfcalc(17)=0,
   dfile(18)='hirs4bufr', dtype(18)='hirs4',     dplat(18)='metop-a', dsis(18)='hirs4_metop-a',      dval(18)=0.0, dthin(18)=1, dsfcalc(18)=1,
   dfile(19)='gimgrbufr', dtype(19)='goes_img',  dplat(19)='g11',     dsis(19)='imgr_g11',           dval(19)=0.0, dthin(19)=1, dsfcalc(19)=0,
   dfile(20)='gimgrbufr', dtype(20)='goes_img',  dplat(20)='g12',     dsis(20)='imgr_g12',           dval(20)=0.0, dthin(20)=1, dsfcalc(20)=0,
   dfile(21)='airsbufr',  dtype(21)='airs',      dplat(21)='aqua',    dsis(21)='airs281SUBSET_aqua', dval(21)=0.0, dthin(21)=1, dsfcalc(21)=1,
   dfile(22)='amsuabufr', dtype(22)='amsua',     dplat(22)='n15',     dsis(22)='amsua_n15',          dval(22)=0.0, dthin(22)=1, dsfcalc(22)=1,
   dfile(23)='amsuabufr', dtype(23)='amsua',     dplat(23)='n18',     dsis(23)='amsua_n18',          dval(23)=0.0, dthin(23)=1, dsfcalc(23)=1,
   dfile(24)='amsuabufr', dtype(24)='amsua',     dplat(24)='metop-a', dsis(24)='amsua_metop-a',      dval(24)=0.0, dthin(24)=1, dsfcalc(24)=1,
   dfile(25)='airsbufr',  dtype(25)='amsua',     dplat(25)='aqua',    dsis(25)='amsua_aqua',         dval(25)=0.0, dthin(25)=1, dsfcalc(25)=1,
   dfile(26)='amsubbufr', dtype(26)='amsub',     dplat(26)='n17',     dsis(26)='amsub_n17',          dval(26)=0.0, dthin(26)=1, dsfcalc(26)=1,
   dfile(27)='mhsbufr',   dtype(27)='mhs',       dplat(27)='n18',     dsis(27)='mhs_n18',            dval(27)=0.0, dthin(27)=1, dsfcalc(27)=1,
   dfile(28)='mhsbufr',   dtype(28)='mhs',       dplat(28)='metop-a', dsis(28)='mhs_metop-a',        dval(28)=0.0, dthin(28)=1, dsfcalc(28)=1,
   dfile(29)='ssmitbufr', dtype(29)='ssmi',      dplat(29)='f14',     dsis(29)='ssmi_f14',           dval(29)=0.0, dthin(29)=1, dsfcalc(29)=0,
   dfile(30)='ssmitbufr', dtype(30)='ssmi',      dplat(30)='f15',     dsis(30)='ssmi_f15',           dval(30)=0.0, dthin(30)=1, dsfcalc(30)=0,
   dfile(31)='amsrebufr', dtype(31)='amsre_low', dplat(31)='aqua',    dsis(31)='amsre_aqua',         dval(31)=0.0, dthin(31)=1, dsfcalc(31)=0,
   dfile(32)='amsrebufr', dtype(32)='amsre_mid', dplat(32)='aqua',    dsis(32)='amsre_aqua',         dval(32)=0.0, dthin(32)=1, dsfcalc(32)=0,
   dfile(33)='amsrebufr', dtype(33)='amsre_hig', dplat(33)='aqua',    dsis(33)='amsre_aqua',         dval(33)=0.0, dthin(33)=1, dsfcalc(33)=0,
   dfile(34)='ssmisbufr', dtype(34)='ssmis',     dplat(34)='f16',     dsis(34)='ssmis_f16',          dval(34)=0.0, dthin(34)=1, dsfcalc(34)=0,
   dfile(35)='gsnd1bufr', dtype(35)='sndrd1',    dplat(35)='g12',     dsis(35)='sndrD1_g12',         dval(35)=0.0, dthin(35)=1, dsfcalc(35)=0,
   dfile(36)='gsnd1bufr', dtype(36)='sndrd2',    dplat(36)='g12',     dsis(36)='sndrD2_g12',         dval(36)=0.0, dthin(36)=1, dsfcalc(36)=0,
   dfile(37)='gsnd1bufr', dtype(37)='sndrd3',    dplat(37)='g12',     dsis(37)='sndrD3_g12',         dval(37)=0.0, dthin(37)=1, dsfcalc(37)=0,
   dfile(38)='gsnd1bufr', dtype(38)='sndrd4',    dplat(38)='g12',     dsis(38)='sndrD4_g12',         dval(38)=0.0, dthin(38)=1, dsfcalc(38)=0,
   dfile(39)='gsnd1bufr', dtype(39)='sndrd1',    dplat(39)='g11',     dsis(39)='sndrD1_g11',         dval(39)=0.0, dthin(39)=1, dsfcalc(39)=0,
   dfile(40)='gsnd1bufr', dtype(40)='sndrd2',    dplat(40)='g11',     dsis(40)='sndrD2_g11',         dval(40)=0.0, dthin(40)=1, dsfcalc(40)=0,
   dfile(41)='gsnd1bufr', dtype(41)='sndrd3',    dplat(41)='g11',     dsis(41)='sndrD3_g11',         dval(41)=0.0, dthin(41)=1, dsfcalc(41)=0,
   dfile(42)='gsnd1bufr', dtype(42)='sndrd4',    dplat(42)='g11',     dsis(42)='sndrD4_g11',         dval(42)=0.0, dthin(42)=1, dsfcalc(42)=0,
   dfile(43)='gsnd1bufr', dtype(43)='sndrd1',    dplat(43)='g13',     dsis(43)='sndrD1_g13',         dval(43)=0.0, dthin(43)=1, dsfcalc(43)=0,
   dfile(44)='gsnd1bufr', dtype(44)='sndrd2',    dplat(44)='g13',     dsis(44)='sndrD2_g13',         dval(44)=0.0, dthin(44)=1, dsfcalc(44)=0,
   dfile(45)='gsnd1bufr', dtype(45)='sndrd3',    dplat(45)='g13',     dsis(45)='sndrD3_g13',         dval(45)=0.0, dthin(45)=1, dsfcalc(45)=0,
   dfile(46)='gsnd1bufr', dtype(46)='sndrd4',    dplat(46)='g13',     dsis(46)='sndrD4_g13',         dval(46)=0.0, dthin(46)=1, dsfcalc(46)=0,
   dfile(47)='iasibufr',  dtype(47)='iasi',      dplat(47)='metop-a', dsis(47)='iasi616_metop-a',    dval(47)=0.0, dthin(47)=1, dsfcalc(47)=1,
   dfile(48)='gomebufr',  dtype(48)='gome',      dplat(48)='metop-a', dsis(48)='gome_metop-a',       dval(48)=0.0, dthin(48)=2, dsfcalc(48)=0,
   dfile(49)='omibufr',   dtype(49)='omi',       dplat(49)='aura',    dsis(49)='omi_aura',           dval(49)=0.0, dthin(49)=2, dsfcalc(49)=0,
   dfile(50)='sbuvbufr',  dtype(50)='sbuv2',     dplat(50)='n19',     dsis(50)='sbuv8_n19',          dval(50)=0.0, dthin(50)=0, dsfcalc(50)=0,
   dfile(51)='hirs4bufr', dtype(51)='hirs4',     dplat(51)='n19',     dsis(51)='hirs4_n19',          dval(51)=0.0, dthin(51)=1, dsfcalc(51)=1,
   dfile(52)='amsuabufr', dtype(52)='amsua',     dplat(52)='n19',     dsis(52)='amsua_n19',          dval(52)=0.0, dthin(52)=1, dsfcalc(52)=1,
   dfile(53)='mhsbufr',   dtype(53)='mhs',       dplat(53)='n19',     dsis(53)='mhs_n19',            dval(53)=0.0, dthin(53)=1, dsfcalc(53)=1,
   dfile(54)='tcvitl'     dtype(54)='tcp',       dplat(54)=' ',       dsis(54)='tcp',                dval(54)=0.0, dthin(54)=0, dsfcalc(54)=0,
   dfile(55)='seviribufr',dtype(55)='seviri',    dplat(55)='m08',     dsis(55)='seviri_m08',         dval(55)=0.0, dthin(55)=1, dsfcalc(55)=0,
   dfile(56)='seviribufr',dtype(56)='seviri',    dplat(56)='m09',     dsis(56)='seviri_m09',         dval(56)=0.0, dthin(56)=1, dsfcalc(56)=0,
   dfile(57)='seviribufr',dtype(57)='seviri',    dplat(57)='m10',     dsis(57)='seviri_m10',         dval(57)=0.0, dthin(57)=1, dsfcalc(57)=0,
   dfile(58)='hirs4bufr', dtype(58)='hirs4',     dplat(58)='metop-b', dsis(58)='hirs4_metop-b',      dval(58)=0.0, dthin(58)=1, dsfcalc(58)=0,
   dfile(59)='amsuabufr', dtype(59)='amsua',     dplat(59)='metop-b', dsis(59)='amsua_metop-b',      dval(59)=0.0, dthin(59)=1, dsfcalc(59)=0,
   dfile(60)='mhsbufr',   dtype(60)='mhs',       dplat(60)='metop-b', dsis(60)='mhs_metop-b',        dval(60)=0.0, dthin(60)=1, dsfcalc(60)=0,
   dfile(61)='iasibufr',  dtype(61)='iasi',      dplat(61)='metop-b', dsis(61)='iasi616_metop-b',    dval(61)=0.0, dthin(61)=1, dsfcalc(61)=0,
   dfile(62)='gomebufr',  dtype(62)='gome',      dplat(62)='metop-b', dsis(62)='gome_metop-b',       dval(62)=0.0, dthin(62)=2, dsfcalc(62)=0,
   dfile(63)='atmsbufr',  dtype(63)='atms',      dplat(63)='npp',     dsis(63)='atms_npp',           dval(63)=0.0, dthin(63)=1, dsfcalc(63)=0,
   dfile(64)='crisbufr',  dtype(64)='cris',      dplat(64)='npp',     dsis(64)='cris_npp',           dval(64)=0.0, dthin(64)=1, dsfcalc(64)=0,
 /
 &SUPEROB_RADAR
 /
 &LAG_DATA
 /
 &HYBRID_ENSEMBLE
 /
 &RAPIDREFRESH_CLDSURF
   dfi_radar_latent_heat_time_period=30.0,
 /
 &CHEM
 /
 &SINGLEOB_TEST
 /
EOF

#===============================================================================

if [ "$MPIEXEC" != 'no' ]; then
  if [ "$MPIEXEC" = 'sr' ]; then
    ./global_gsi > gsi.log 2>&1
  elif [ -s "$MACHINE" ]; then
    $MPIEXEC -n $NP -machinefile $MACHINE ./global_gsi > gsi.log 2>&1
  else
    $MPIEXEC -n $NP ./global_gsi > gsi.log 2>&1
  fi
fi

#===============================================================================

exit 0
