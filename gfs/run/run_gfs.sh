#!/bin/bash
#===============================================================================
#
#  Script to run GFS model.
#  October 2012  simplified from 'exglobal_fcst.sh.sms.sh', Guo-Yuan Lien
#  April   2013  modified, Guo-Yuan Lien
#
#===============================================================================

if [ "$#" -lt 7 ]; then
  cat 1>&2 << EOF

[run_gfs.sh] Prepare a temporary directory for GFS model run,
             and may run the model.

Usage: $0 SIGI SFCI FHMAX FHOUT TMPDIR EXECDIR FIXDIR [MPIEXEC] [NP] [MACHINE]

  SIGI     Initial condition: GFS sigma file
  SFCI     Initial condition: GFS surface file
  FHMAX    Forecast length (hour)
  FHOUT    Output interval (hour)
  TMPDIR   Temporary directory to run the model
  EXECDIR  Directory of GFS executable files
  FIXDIR   Directory of GFS fix files
  MPIEXEC  Directory of the mpiexec program
           'no': do not run GFS model (just prepare the temporary directory)
           'se': run GFS model serially
           (default: 'no')
  NP       Number of cores
  MACHINE  machinefile for mpiexec

EOF
  exit 1
fi

SIGI="$1"
SFCI="$2"
FHMAX="$3"
FHOUT="$4"
TMPDIR="$5"
EXECGLOBAL="$6"
FIXGLOBAL="$7"
MPIEXEC="${8:-no}"
NP="${9:-1}"
MACHINE="$10"

FCSTEXEC="${EXECGLOBAL}/global_fcst"
SIGHDR="${EXECGLOBAL}/global_sighdr"

#===============================================================================

JCAP=`echo jcap | $SIGHDR $SIGI`
LEVS=`echo levs | $SIGHDR $SIGI`
LONR=`echo lonr | $SIGHDR $SIGI`
LATR=`echo latr | $SIGHDR $SIGI`
LONF=`echo lonf | $SIGHDR $SIGI`
LATG=`echo latf | $SIGHDR $SIGI`
NTRAC=`echo ntrac | $SIGHDR $SIGI`
NMTVR=14
LSOIL=4
NTOZ=2
NTCW=3
NCLD=1
NGPTC=30
######
if [ "$JCAP" -eq 62 ]; then
  DELTIM=1200
else
  DELTIM=$((3600/(JCAP/20)))
fi
######

IDVC=`echo idvc | $SIGHDR $SIGI`
if [ $IDVC = 1 ] ; then
 HYBRID=.false.
 GEN_COORD_HYBRID=.false.
elif [ $IDVC = 2 ] ; then
 HYBRID=.true.
 GEN_COORD_HYBRID=.false.
elif [ $IDVC = 3 ] ; then
 HYBRID=.false.
 GEN_COORD_HYBRID=.true.
fi

#===============================================================================

mkdir -p $TMPDIR
rm -fr $TMPDIR/*
cd $TMPDIR

ln -fs $FCSTEXEC .

ln -fs ${FIXGLOBAL}/global_co2con.l${LEVS}.f77 fort.15
ln -fs ${FIXGLOBAL}/global_mtnvar.t${JCAP}.f77 fort.24
ln -fs ${FIXGLOBAL}/global_tbthe.f77 fort.27
ln -fs ${FIXGLOBAL}/global_o3prdlos.f77 fort.28
ln -fs ${FIXGLOBAL}/global_cldtune.f77 fort.43
ln -fs ${FIXGLOBAL}/global_o3clim.txt fort.48
ln -fs ${FIXGLOBAL}/global_sfc_emissivity_idx.txt sfc_emissivity_idx.txt
ln -fs ${FIXGLOBAL}/global_orography.t$JCAP.grb orography
ln -fs ${FIXGLOBAL}/global_climaeropac_global.txt aerosol.dat
for m in 01 02 03 04 05 06 07 08 09 10 11 12; do
  if [ -s "${FIXGLOBAL}/global_aeropac3a.m$m.txt" ]; then
    ln -fs ${FIXGLOBAL}/global_aeropac3a.m$m.txt aeropac3a.m$m
  fi
done
for y in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012; do
  if [ -s "${FIXGLOBAL}/global_co2historicaldata_$y.txt" ]; then
    ln -fs ${FIXGLOBAL}/global_co2historicaldata_$y.txt co2historicaldata_$y.txt
  fi
done

ln -fs $SIGI sig_ini
ln -fs $SFCI sfc_ini
  
#===============================================================================

cat << EOF > gfs_namelist.rc

#nam_gfs +++++++++++++++++++++++++++
NLUNIT:                  35
DELTIM:                  ${DELTIM}.0
FHROT:                   0
NAMELIST:                gfs_namelist
TOTAL_MEMBER:            1
GRIB_INPUT:              0
PE_MEMBER01:             0

# For stachastic purturbed runs -  added by dhou and wyang
HH_INCREASE:             $FHMAX
HH_FINAL:                $FHMAX
HH_START:                00
ADVANCECOUNT_SETUP:      0

SPS_PARM1:               0.005 10.0 0.005 10.0 0.0 0.0 0.0 0.0 0.0 0.0
SPS_PARM2:               0.105 0.03 0.12 42.0 0.0 0.0 0.0 0.0 0.0 0.0
SPS_PARM3:               0.2 0.34 -0.34 3.0 0.0 0.0 0.0 0.0 0.0 0.0


#ESMF_State_Namelist +++++++++++++++

ENS_SPS:                          .false.
HOUTASPS:                         10000

IDATE1_IMPORT:                    .false.
Z_IMPORT:                         .false.
PS_IMPORT:                        .false.
VOR_IMPORT:                       .false.
DIV_IMPORT:                       .false.
TEMP_IMPORT:                      .false.
Q_IMPORT:                         .false.
OZ_IMPORT:                        .false.
SCLD_IMPORT:                      .false.
TRIEO_IMPORT:                     .false.

IDATE1_EXPORT:                    .false.
Z_EXPORT:                         .false.
PS_EXPORT:                        .false.
VOR_EXPORT:                       .false.
DIV_EXPORT:                       .false.
TEMP_EXPORT:                      .false.
Q_EXPORT:                         .false.
OZ_EXPORT:                        .false.
SCLD_EXPORT:                      .false.
TRIEO_EXPORT:                     .true.

# Surface state.
#---------------
OROGRAPHY_IMPORT:                 .false.
T_SKIN_IMPORT:                    .false.
SOIL_MOIS_IMPORT:                 .false.
SNOW_DEPTH_IMPORT:                .false.
SOIL_T_IMPORT:                    .false.
DEEP_SOIL_T_IMPORT:               .false.
ROUGHNESS_IMPORT:                 .false.
CONV_CLOUD_COVER_IMPORT:          .false.
CONV_CLOUD_BASE_IMPORT:           .false.
CONV_CLOUD_TOP_IMPORT:            .false.
ALBEDO_VISIBLE_SCATTERED_IMPORT:  .false.
ALBEDO_VISIBLE_BEAM_IMPORT:       .false.
ALBEDO_NEARIR_SCATTERED_IMPORT:   .false.
ALBEDO_NEARIR_BEAM_IMPORT:        .false.
SEA_LEVEL_ICE_MASK_IMPORT:        .false.
VEGETATION_COVER_IMPORT:          .false.
CANOPY_WATER_IMPORT:              .false.
M10_WIND_FRACTION_IMPORT:         .false.
VEGETATION_TYPE_IMPORT:           .false.
SOIL_TYPE_IMPORT:                 .false.
ZENEITH_ANGLE_FACSF_IMPORT:       .false.
ZENEITH_ANGLE_FACWF_IMPORT:       .false.
UUSTAR_IMPORT:                    .false.
FFMM_IMPORT:                      .false.
FFHH_IMPORT:                      .false.
SEA_ICE_THICKNESS_IMPORT:         .false.
SEA_ICE_CONCENTRATION_IMPORT:     .false.
TPRCP_IMPORT:                     .false.
SRFLAG_IMPORT:                    .false.
ACTUAL_SNOW_DEPTH_IMPORT:         .false.
LIQUID_SOIL_MOISTURE_IMPORT:      .false.
VEGETATION_COVER_MIN_IMPORT:      .false.
VEGETATION_COVER_MAX_IMPORT:      .false.
SLOPE_TYPE_IMPORT:                .false.
SNOW_ALBEDO_MAX_IMPORT:           .false.

OROGRAPHY_EXPORT:                 .false.
T_SKIN_EXPORT:                    .false.
SOIL_MOIS_EXPORT:                 .false.
SNOW_DEPTH_EXPORT:                .false.
SOIL_T_EXPORT:                    .false.
DEEP_SOIL_T_EXPORT:               .false.
ROUGHNESS_EXPORT:                 .false.
CONV_CLOUD_COVER_EXPORT:          .false.
CONV_CLOUD_BASE_EXPORT:           .false.
CONV_CLOUD_TOP_EXPORT:            .false.
ALBEDO_VISIBLE_SCATTERED_EXPORT:  .false.
ALBEDO_VISIBLE_BEAM_EXPORT:       .false.
ALBEDO_NEARIR_SCATTERED_EXPORT:   .false.
ALBEDO_NEARIR_BEAM_EXPORT:        .false.
SEA_LEVEL_ICE_MASK_EXPORT:        .false.
VEGETATION_COVER_EXPORT:          .false.
CANOPY_WATER_EXPORT:              .false.
M10_WIND_FRACTION_EXPORT:         .false.
VEGETATION_TYPE_EXPORT:           .false.
SOIL_TYPE_EXPORT:                 .false.
ZENEITH_ANGLE_FACSF_EXPORT:       .false.
ZENEITH_ANGLE_FACWF_EXPORT:       .false.
UUSTAR_EXPORT:                    .false.
FFMM_EXPORT:                      .false.
FFHH_EXPORT:                      .false.
SEA_ICE_THICKNESS_EXPORT:         .false.
SEA_ICE_CONCENTRATION_EXPORT:     .false.
TPRCP_EXPORT:                     .false.
SRFLAG_EXPORT:                    .false.
ACTUAL_SNOW_DEPTH_EXPORT:         .false.
LIQUID_SOIL_MOISTURE_EXPORT:      .false.
VEGETATION_COVER_MIN_EXPORT:      .false.
VEGETATION_COVER_MAX_EXPORT:      .false.
SLOPE_TYPE_EXPORT:                .false.
SNOW_ALBEDO_MAX_EXPORT:           .false.

EOF

cat  > gfs_namelist <<EOF
&nam_mrf
  FHOUT=$FHOUT,
  FHMAX=$FHMAX,
  IGEN=0,
  DELTIM=$DELTIM,
  FHRES=24,
  FHZER=6,
  FHLWR=3,
  FHSWR=1,
  FHROT=0,
  FHDFI=0,
  FHCYC=0,
  ntrac=$NTRAC,
  jcap=$JCAP,
  levs=$LEVS,
  lonf=$LONF,
  lonr=$LONR,
  latg=$LATG,
  latr=$LATR,
  ntoz=$NTOZ,
  ntcw=$NTCW,
  ncld=$NCLD,
  lsoil=$LSOIL,
  nmtvr=14,
  ngptc=30,
  hybrid=$HYBRID,
  tfiltc=0.85,
  gen_coord_hybrid=$GEN_COORD_HYBRID,
  FHOUT_HF=1,
  FHMAX_HF=0,
  IEMS=0,
  ISOL=0,
  IAER=111,
  ICO2=1,
  gfsio_in=.false.,
  gfsio_out=.false.,
!  ncw=20,120,
!  flgmin=0.180,0.220,
!  phigs_d=70.0,
!  random_clds=.false.,
!  redgg_a=.false.,
!  redgg_b=.false.,
!  newsas=.false.,
!  sashal=.false.,
!  zhao_mic=.false.,
/
&TRACER_CONSTANT
/
&SOIL_VEG
  LPARAM = .false.
/
&NAMSFC
  FNGLAC="${FIXGLOBAL}/global_glacier.2x2.grb",
  FNMXIC="${FIXGLOBAL}/global_maxice.2x2.grb",
  FNTSFC="${FIXGLOBAL}/cfs_oi2sst1x1monclim19822001.grb",
  FNSNOC="${FIXGLOBAL}/global_snoclim.1.875.grb",
  FNZORC="${FIXGLOBAL}/global_zorclim.1x1.grb",
  FNALBC="${FIXGLOBAL}/global_albedo4.1x1.grb",
  FNAISC="${FIXGLOBAL}/cfs_ice1x1monclim19822001.grb",
  FNTG3C="${FIXGLOBAL}/global_tg3clim.2.6x1.5.grb",
  FNVEGC="${FIXGLOBAL}/global_vegfrac.0.144.decpercent.grb",
  FNVETC="${FIXGLOBAL}/global_vegtype.1x1.grb",
  FNSOTC="${FIXGLOBAL}/global_soiltype.1x1.grb",
  FNSMCC="${FIXGLOBAL}/global_soilmcpc.1x1.grb",
  FNMSKH="${FIXGLOBAL}/seaice_newland.grb",
  FNTSFA="",
  FNACNA="",
  FNSNOA="",
  FNVMNC="${FIXGLOBAL}/global_shdmin.0.144x0.144.grb",
  FNVMXC="${FIXGLOBAL}/global_shdmax.0.144x0.144.grb",
  FNSLPC="${FIXGLOBAL}/global_slope.1x1.grb",
  FNABSC="${FIXGLOBAL}/global_snoalb.1x1.grb",
  LDEBUG=.false.,
  FSMCL(2)=99999,
  FSMCL(3)=99999,
  FSMCL(4)=99999,
  FTSFS=90,
  FAISS=99999,
  FSNOL=99999,
  FSICL=99999,
  FTSFL=99999,
  FAISL=99999,
  FVETL=99999,
  FSOTL=99999,
  FvmnL=99999,
  FvmxL=99999,
  FSLPL=99999,
  FABSL=99999,
  FSNOS=99999,
  FSICS=99999,
/
EOF

#===============================================================================

if [ "$MPIEXEC" != 'no' ]; then
# -- Only for IBM machines??
#  export MP_BINDPROC="yes"
#  export MEMORY_AFFINIY="MCM"
#  export MP_SYNC_QP="yes"
#  export MP_SHARED_MEMORY="NO"
#  export MP_COREFILE_FORMAT="lite"
#  export XLSMPOPTS="parthds=1:stack=128000000"
#  export BIND_TASKS="no"

  if [ "$MPIEXEC" = 'sr' ]; then
    ./global_fcst > gfs.log 2>&1
  elif [ -s "$MACHINE" ]; then
    $MPIEXEC -n $NP -machinefile $MACHINE ./global_fcst > gfs.log 2>&1
  else
    $MPIEXEC -n $NP ./global_fcst > gfs.log 2>&1
  fi
fi

#===============================================================================

exit 0
