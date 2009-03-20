#!/bin/sh
#  plot bred vectors
set -e
TIME=0237
MEM=001
XMIN=-770
XMAX=0
DEPTH=-250.0
Y=57
CDIR=`pwd`
cd ../..
ROMS=`pwd`
NATURE=$ROMS/DATA/nature
GRD=$ROMS/model/bc/ee6_grd.nc
BVEXP=bv_ubar_5E-3
cd $CDIR
cat << EOF > tmp.ncl
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;wks=gsn_open_wks("x11","gsun02n")
wks=gsn_open_wks("eps","out")
;gsn_define_colormap(wks,"testcmap")
;gsn_define_colormap(wks,"gui_default")
gsn_define_colormap(wks,"rainbow")
;gsn_define_colormap(wks,"gsltod")
data1=addfile("$NATURE/$TIME.nc","r")
data2=addfile("$ROMS/DATA/$BVEXP/bv/$MEM/$TIME.nc","r")
grddata=addfile("$GRD","r")
xlab=fspan(-770,21,114)
ylab=fspan(0,791,114)
h1=grddata->h($Y,:)
h=dble2flt(h1)
cs=(/ -0.921962804129276, -0.783682417800398, -0.666141723246888, -0.566230075492245, -0.481303371858989, -0.409114077354825, -0.347751744988103, -0.295592456917874, -0.251255848428496, -0.213568577393686, -0.181533272470966, -0.154302138258862, -0.131154518892926, -0.111477826314426, -0.0947513284896102, -0.0805323685443561, -0.0684446501117523, -0.0581682788712608, -0.0494312967346868, -0.0420024846361897, -0.0356852434564318, -0.0303123911431393, -0.0257417383369477, -0.0218523254140864, -0.0185412213614628, -0.0157207997682656, -0.0133164198448193, -0.0112644510982221, -0.00951058938897752, -0.00800841980405428, -0.0067181883136592, -0.00560574970434628, -0.00464166394613383, -0.00380041707951966, -0.00305974600362987, -0.00240004929641906, -0.00180386847462664, -0.00125542596534104, -0.000740207561862293, -0.000244578313804319 /)

z1=data1->temp(0,:,$Y,:)
z1&xi_rho=xlab
;printVarSummary(z1)
zeta1=data1->zeta(0,$Y,:)
zdim=dimsizes(z1)
x=new(zdim,typeof(z1))
z=new(zdim,typeof(z1))
nlev=40
do k=0,nlev-1
sc=(k+1.0-nlev-0.5)/nlev
do i=0,113
  z(k,i)=zeta1(i)+(zeta1(i)+h(i))*(10.0*sc+cs(k)*h(i))/(h(i)+10.0)
end do
x(k,:)=xlab
end do
res1=True
res1@tiMainString="NATURE(shade), BV(contour)"
res1@tiYAxisString="DEPTH [m]"
res1@tiXAxisString="DISTANCE FROM COAST [km]"
res1@gsnDraw=False
res1@gsnFrame=False
res1@gsnSpreadColors=True
res1@cnFillOn=True
res1@cnLinesOn=False
;res1@cnCellFillEdgeColor=2
;res1@cnFillMode="CellFill"
res1@cnLevelSelectionMode="ManualLevels"
res1@cnMaxLevelValF=20
res1@cnMinLevelValF=7
res1@cnLevelSpacingF=0.1
res1@cnLineLabelsOn=False
res1@lbLabelAutoStride=True
res1@lbOrientation="Vertical"
res1@trXMinF=$XMIN
res1@trXMaxF=$XMAX
res1@trYMinF=$DEPTH
res1@tmYLMode="Manual"
res1@tmYLTickSpacingF=50.0
res1@trYMaxF=max(z)
res1@trGridType="TriangularMesh"
res1@sfXArray=x
res1@sfYArray=z
plot1=gsn_csm_contour(wks,z1,res1)

z2=data2->temp(0,:,$Y,:)
z2&xi_rho=xlab
z2=z2/max(abs(z2))
;printVarSummary(z2)
res2=True
res2@tiMainString=""
res2@gsnDraw=False
res2@gsnFrame=False
res2@gsnLeftString=""
res2@gsnRightString=""
res2@tiXAxisString=""
res2@tiYAxisString=""
res2@cnLineLabelsOn=False
res2@cnLevelSpacingF=.05
res2@gsnContourZeroLineThicknessF=0
res2@gsnContourNegLineDashPattern=1
res2@trXMinF=$XMIN
res2@trXMaxF=$XMAX
res2@trYMinF=$DEPTH
res1@tmYLMode="Manual"
res1@tmYLTickSpacingF=50.0
res2@trYMaxF=max(z)
res2@trGridType="TriangularMesh"
res2@sfXArray=x
res2@sfYArray=z
plot2=gsn_csm_contour(wks,z2,res2)
;plot2=ColorShadeLeGeContour(plot2,-0.01,"cyan",0.01,"salmon")

overlay(plot1,plot2)

draw(plot1)
frame(wks)
end
EOF
ncl tmp.ncl
rm tmp.ncl
mv out.eps figs/bvvert$BVEXP.eps
display figs/bvvert$BVEXP.eps
echo "NORMAL END"
