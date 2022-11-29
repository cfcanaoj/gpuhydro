reset

pngflag=1
if(pngflag==1)set terminal push
if(pngflag==1)set terminal pngcairo dashed enhanced font "Helvetica,20"
if(pngflag==1)set encoding utf8

set style line 1 lt 1 lw 6 lc rgb "#ff2800" # universal design red 
set style line 2 lt 1 lw 6 lc rgb "#0041ff" # universal design blue
set style line 3 lt 1 lw 6 lc rgb "#35a16B" # universal design green
set style line 4 lt 1 lw 6 lc rgb "#faf500" # universal design yellow
set style line 5 lt 1 lw 6 lc rgb "#66ccff" # universal design sky-blue,azure
set style line 6 lt 1 lw 6 lc rgb "#ff99a0" # universal design pink
set style line 7 lt 1 lw 6 lc rgb "#ff9900" # universal design orange
set style line 8 lt 1 lw 6 lc rgb "#9a0079" # universal design purple
set style line 9 lt 1 lw 6 lc rgb "#663300" # universal design brown
 
set style line 91 lt 1 lw 2 lc rgb "black" # 
set style line 92 lt 2 lw 6 lc rgb "black" #

# input file


input= sprintf("output/spc%05d.dat",ifnum)

##########################################
# Evolusiton of shock and Gain Radius
##########################################

set log 
set format y "10^{%L}"

set xlabel "Wave number"
set xrange [1:100]

####################
# kinetic energy
####################
outputfile= sprintf("figures/ksp%05d.png",ifnum)
if(pngflag==1)set output outputfile

set ylabel "Kinetic energy"
set yrange [*:*]

set key right top

plot NaN notitle \
, input  u 1:2  notitle w l ls 1  \
#, 6*x**(-3.0/2.0) title "-3/2" w l ls 91


####################
# Kinetic helicity
####################

outputfile= sprintf("figures/hks%05d.png",ifnum)
if(pngflag==1)set output outputfile

set log 
set format y "10^{%L}"

set ylabel "Enstrophy"

plot NaN notitle \
, input  u 1:3  notitle w l ls 1  \
#, 0.3*(x/10)**(-5.0/3.0) title "-5/3" w l ls 91

##########
# 
##########
outputfile= sprintf("figures/hmm%05d.png",ifnum)
if(pngflag==1)set output outputfile

set log 

set ylabel "Magnetic helicity mimic"

set log 
set format y "10^{%L}"

plot NaN notitle \
, input  u 1:4  title "2D MHD" w l ls 1  \
#, 5.0e-7*(x/10)**(-3.0/2.0) title "-3/2" w l ls 91


##########
# 
##########
outputfile= sprintf("figures/hcr%05d.png",ifnum)
if(pngflag==1)set output outputfile

set log 

set ylabel "Cross helicity"

set log 
set format y "10^{%L}"

plot NaN notitle \
, input  u 1:5  title "2D MHD" w l ls 1  \
#, 1.0e-3*(x/10)**(-3.0/2.0) title "-3/2" w l ls 91

##########
# 
##########

outputfile= sprintf("figures/msp%05d.png",ifnum)
if(pngflag==1)set output outputfile

set log 

set ylabel "Magnetic energy"

set log 
set format y "10^{%L}"

plot NaN notitle \
, input  u 1:6  title "2D MHD" w l ls 1  \
#, 1.0e-3*(x/10)**(-3.0/2.0) title "-3/2" w l ls 91


reset
set term pop
