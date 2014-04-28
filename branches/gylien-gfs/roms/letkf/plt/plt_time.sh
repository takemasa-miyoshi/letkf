#!/bin/sh
#  plot bred vectors
set -e
OBS=5x20_new
EXP=LETKF20H50V100
VARLIST="u v temp salt"
ZLIST="40 30 20 10"
OUTPUT=figs/$OBS/$EXP
mkdir -p $OUTPUT
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
cd $CDIR
for VAR in $VARLIST
do
case "$VAR" in
  "u"    ) UNIT=mps;XMAX=110;YMAX=114;;
  "v"    ) UNIT=mps;XMAX=111;YMAX=113;;
  "temp" ) UNIT=K  ;XMAX=111;YMAX=114;;
  "salt" ) UNIT=PSU;XMAX=111;YMAX=114;;
  *      ) echo "ABEND";exit;;
esac
XMIN=1
YMIN=1
for Z in $ZLIST
do
END=`ls $NATURE/200401*_rst.nc | wc -l`
END1=`expr $END - 1`
X0=`expr $XMIN - 1` || test 1 -eq 1
X1=`expr $XMAX - 1`
Y0=`expr $YMIN - 1` || test 1 -eq 1
Y1=`expr $YMAX - 1`
Z1=`expr $Z - 1`
cat << EOF > tmp.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
;wks=gsn_open_wks("x11","gsun02n")
wks=gsn_open_wks("eps","out")
;gsn_define_colormap(wks,"testcmap")
;gsn_define_colormap(wks,"gui_default")
gsn_define_colormap(wks,"rainbow")
;gsn_define_colormap(wks,"gsltod")
files1=systemfunc("ls $NATURE/200401*_rst.nc")
files2=systemfunc("ls $ROMS/DATA/$OBS/$EXP/gues/mean/200401*_rst.nc")
files3=systemfunc("ls $ROMS/DATA/$OBS/$EXP/anal/mean/200401*_rst.nc")
files4=systemfunc("ls $ROMS/DATA/$OBS/$EXP/gues/sprd/200401*_rst.nc")
files5=systemfunc("ls $ROMS/DATA/$OBS/$EXP/anal/sprd/200401*_rst.nc")
xlab=fspan(-770,21,114)
ylab=fspan(0,0.25*$END1,$END)

data=addfiles(files1,"r")
z=addfiles_GetVar(data,files1,"$VAR")
z1=z(0:$END1,$Z1,$Y0:$Y1,$X0:$X1)
z1&time=ylab
;printVarSummary(z1)
delete(z)
delete(data)
;
data=addfiles(files2,"r")
z=addfiles_GetVar(data,files2,"$VAR")
z2=z(0:$END1,$Z1,$Y0:$Y1,$X0:$X1)
z2&time=ylab
;printVarSummary(z2)
delete(z)
delete(data)
;
data=addfiles(files3,"r")
z=addfiles_GetVar(data,files3,"$VAR")
z3=z(0:$END1,$Z1,$Y0:$Y1,$X0:$X1)
z3&time=ylab
delete(z)
delete(data)
;
data=addfiles(files4,"r")
z=addfiles_GetVar(data,files4,"$VAR")
z4=z(0:$END1,$Z1,$Y0:$Y1,$X0:$X1)
z4&time=ylab
delete(z)
delete(data)
;
data=addfiles(files5,"r")
z=addfiles_GetVar(data,files5,"$VAR")
z5=z(0:$END1,$Z1,$Y0:$Y1,$X0:$X1)
z5&time=ylab
delete(z)
delete(data)
;
rms=new((/4,$END/),typeof(z1))
do i=0,$END1
rms(0,i)=sqrt(avg((z1(i,:,:)-z2(i,:,:))^2))
rms(1,i)=sqrt(avg((z1(i,:,:)-z3(i,:,:))^2))
rms(2,i)=sqrt(avg((z4(i,:,:))^2))
rms(3,i)=sqrt(avg((z5(i,:,:))^2))
end do
res=True
res@tiXAxisString="TIME [day]"
res@tiYAxisString="$VAR RMSE [$UNIT]"
res@xyLineColors=(/"blue","red","blue","red"/)
res@xyDashPatterns=(/0,0,1,1/)
res@pmLegendDisplayMode="Always"
res@pmLegendSide="Top"
res@pmLegendParallelPosF=0.8
res@pmLegendOrthogonalPosF=-0.35
res@pmLegendWidthF=0.2
res@pmLegendHeight=0.1
res@lgLegendFontHeightF=0.7
res@lgPerimOn=False
res@tmXBMode="Manual"
res@tmXBTickSpacingF=10
res@xyExplicitLegendLabels=(/"RMSE_G","RMSE_A","SPRD_G","SPRD_A"/)
plot=gsn_csm_xy(wks,ylab,rms,res)
end
EOF
ncl tmp.ncl
rm tmp.ncl
mv out.eps $OUTPUT/time_${VAR}_$Z.eps
done
done
echo "NORMAL END"
