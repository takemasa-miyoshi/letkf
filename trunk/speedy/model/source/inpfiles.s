# Script to link fortran units to input files; $1 = resolution (t21, t30, ..)

SB=../data/bc/$1/clim
SC=../data/bc/$1/anom

 ln -s $SB/orog_lsm_alb.${1}.grd         fort.20
 ln -s $SB/sst_8190clim.${1}.sea.grd     fort.21
 ln -s $SB/seaice_8190clim.${1}.sea.grd  fort.22
 ln -s $SB/skt_8190clim.${1}.land.grd    fort.23
 ln -s $SB/sndep_8190clim.${1}.land.grd  fort.24
 ln -s $SB/veget.${1}.land.grd           fort.25
 ln -s $SB/soilw_8190clim.${1}.land.grd  fort.26

#cp    $SC/sst_anom_1950_1999_mean8190.${1}.grd  fort.30
 cp    $SC/sst_anom_7990.${1}.grd	         fort.30
#cp    $SC/sst_anom_8190.${1}.sea.grd	         fort.30
#cp    $SC/elnino_anom_p1.${1}.grd	         fort.30
	
#cp    $HOME/speedy/hflux/clim_hflux_582.grd	 fort.31

