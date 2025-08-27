#####################
# initialize
######################

set tmargin 0
set bmargin 0
set lmargin 0
set rmargin 0

set size ratio -1
set view map
unset key

#set palette rgbformulae 21,22,23 # black red yellow white

#set palette functions sqrt(gray), gray**3, sin(gray*2*pi) # black blue red yellow
#set palette rgb 33,13,10 # blue-green-yellow-red

#set palette define (-1.0 "blue", 0.0 "white", 1.0 "red")
#set palette define (1.0 "black", 1.0 "black")

set palette rgb 33,13,10 # rainbow
#set palette define (0.9 "dark-magenta", 1.0 "dark-magenta", 1.14 "navy", 1.28 "blue", 1.42 "dark-green", 1.56 "light-green", 1.7 "yellow", 1.84 "orange" ,2.0 "dark-red", 2.1 "dark-red")

set style line 10 lt  1 lw 4 lc rgb "white" #

set style line 1 lt 1 lw 1 lc rgb "black" #

##########################################
# parameters
##########################################

srange=0.5

##########################################
# mainloop
##########################################

# PNG
if (exist("ifnum")==0 ) set term push
set term pngcairo enhanced font "Helvetica, 12" size 600,600
# crop 

if (exist("ifnum")==0 ) ifnum=100

ifnames = sprintf("output/vor%05d.dat",ifnum)
ifnamev = sprintf("meddata/Sc%05d.xsv",ifnum)

command = sprintf("ls %s 1> /dev/null 2> /dev/null ; echo $? ",ifnames)
flag=0
flag=system(command)
print  ifnames." found"

command = sprintf(" head -n 1 %s | sed 's/#//' ",ifnames)
time   = system(command)
print "time= ".time

# Position of color bar
set colorbox horizontal user origin 0.175, 0.87 size 0.71, 0.04
#set cbtics 1.0
#set cbtics offset 0,3.2

maxcb=20.0
set cbrange [-maxcb:maxcb]

ofname = sprintf("figures/vor%05d.png",ifnum)
set output ofname

set size 1.0

set origin 0.05,0.0
set xlabel "X" offset 0,0
set xtics 0.25
set ylabel "Y" offset 0,0
set ytics 0.25

set xrange [-0.5:0.5]
set yrange [-0.5:0.5]

vr=srange/10

set label 1 time at screen 0.65, screen 0.845

set label 2 "PLM" at screen 0.3, screen 0.845

splot  \
  ifnames u ($1):($2):3 w pm3d  \

ofname = sprintf("figures/jcd%05d.png",ifnum)
set output ofname

maxcb=20.0
set cbrange [-maxcb:maxcb]

splot  \
  ifnames u ($1):($2):4 w pm3d  \

unset label 1

reset
set term pop
