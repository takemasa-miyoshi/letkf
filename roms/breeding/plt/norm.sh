#!/bin/sh
#  plot fcst norm
set -e
COLOR=color
CDIR=`pwd`
cd ../..
ROMS=`pwd`
BVEXP=bv_sst_5E-1
#BVEXP=breeding
RESCALE=0.5
#
cd $CDIR
awk '{print NR/4, $1}' $ROMS/DATA/$BVEXP/rescale001.dat > 001.dat
awk '{print NR/4, $1}' $ROMS/DATA/$BVEXP/rescale002.dat > 002.dat
awk '{print NR/4, $1}' $ROMS/DATA/$BVEXP/rescale003.dat > 003.dat
cat << EOF > tmp.plt
set key bottom right
set size 0.6
set xrange [0.5:60]
set xlabel "day"
set ylabel "norm [m/s]"
set title "fcst bv norm = sqrt(ubar^2 + vbar^2)"
plot $RESCALE notitle, \
"001.dat" ti "001" w l, "002.dat" ti "002" w l, "003.dat" ti "003" w l
set terminal postscript eps enhanced $COLOR
set output "output.eps"
plot $RESCALE notitle, \
"001.dat" ti "001" w l, "002.dat" ti "002" w l, "003.dat" ti "003" w l
set terminal png
set output "output.png"
plot $RESCALE notitle, \
"001.dat" ti "001" w l, "002.dat" ti "002" w l, "003.dat" ti "003" w l
EOF
gnuplot -persist tmp.plt
rm tmp.plt
rm 001.dat
rm 002.dat
rm 003.dat
mv output.png $CDIR/figs/norm$BVEXP.png
mv output.eps $CDIR/figs/norm$BVEXP.eps
echo "NORMAL END"
