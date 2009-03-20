#!/bin/sh
#  plot bred vectors
set -e
END=120
VAR=temp
UNIT=K
XMIN=1
XMAX=111
YMIN=1
YMAX=114
Z=40
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
EXP=letkf001
END1=`expr $END - 1`
X0=`expr $XMIN - 1` || test 1 -eq 1
X1=`expr $XMAX - 1`
Y0=`expr $YMIN - 1` || test 1 -eq 1
Y1=`expr $YMAX - 1`
Z1=`expr $Z - 1`
cd $CDIR
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
files1=systemfunc("ls $NATURE/*.nc")
files2=systemfunc("ls $ROMS/DATA/$EXP/gues/mean/*.nc")
files3=systemfunc("ls $ROMS/DATA/$EXP/anal/mean/*.nc")
files4=systemfunc("ls $ROMS/DATA/$EXP/gues/sprd/*.nc")
files5=systemfunc("ls $ROMS/DATA/$EXP/anal/sprd/*.nc")
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
display out.eps
mv out.eps figs/time_${EXP}_$VAR.eps
echo "NORMAL END"
