##########################################
# Makefile for Turbulence Studies
##########################################

##########################################
# Programs
##########################################

vis=Visualization.x

######################	
# complile options
######################
fc=ifort -extend-source
#fopt=-g -traceback -check all -fpe0

HDF5 = /home/takiwkkz/hdf5
FLAGHDF =  -I${HDF5}/include -L${HDF5}/lib -lhdf5 -lhdf5_fortran 

##########################################
# Timeseqencial data and figures
##########################################

######################
# directory
######################

dira := output
dirf := figures
dirm := movies

countfile=control.dat

######################
# files
######################

filename = ${dird}/unf00001.dat
fileexists = $(shell ls | grep ${filename})
ifeq (${fileexists}, ${filename})
    BIN= $(shell ls ${dird}/unf*.dat)
else
    BIN= ${filename}
endif

VOR   = $(patsubst ${dird}/unf%.dat,${dira}/vor%.dat,$(BIN))
VORPNG= $(patsubst ${dira}/vor%.dat,${dirf}/vor%.png,$(VOR))
KINPNG= $(patsubst ${dira}/vor%.dat,${dirf}/kin%.png,$(VOR))
JCDPNG= $(patsubst ${dira}/vor%.dat,${dirf}/jcd%.png,$(VOR))
MAGPNG= $(patsubst ${dira}/vor%.dat,${dirf}/mag%.png,$(VOR))

filenameS = ${dira}/spc00001.dat
fileexistsS = $(shell ls | grep ${filenameS})
ifeq (${fileexistsS}, ${filenameS})
    SPC= $(shell ls ${dira}/spc*.dat)
else
    SPC= ${filename}
endif

KSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/ksp%.png,$(SPC))
HKSPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hks%.png,$(SPC))
KSPPNG= $(patsubst ${dira}/spc%.dat,${dirf}/ksp%.png,$(SPC))
HCRPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hcr%.png,$(SPC))
HMMPNG= $(patsubst ${dira}/spc%.dat,${dirf}/hmm%.png,$(SPC))

##########################################
# Movies
##########################################
##

# 2D snaps evolution
vormovie=${dirm}/anivor.mp4
kinmovie=${dirm}/anikin.mp4
jcdmovie=${dirm}/anijcd.mp4
magmovie=${dirm}/animag.mp4

2Dsnapmovies = ${vormovie} ${kinmovie} ${jcdmovie} ${magmovie}

# spectrum evolution
kspmovie=${dirm}/aniksp.mp4
hksmovie=${dirm}/anihks.mp4
hcrmovie=${dirm}/anihcr.mp4
hmmmovie=${dirm}/anihmm.mp4
mspmovie=${dirm}/animsp.mp4

specmovies = ${kspmovie} ${hksmovie} ${hcrmovie} ${hmmmovie} ${mspmovie}

moviefile = ${2Dsnapmovies} ${specmovies}

##########################################
# Final spectrum
##########################################

specfile=k-E_k.png

##########################################
# Time evolution
##########################################

timefile=t-E.png

##########################################
# Proceadures
##########################################

all: ${moviefile} ${specfile} ${timefile}

.PHONY: all clean allclean

##########################################
# Time evolution
##########################################
${timefile}: TimeEvolution.plt t-prof.dat
	gnuplot TimeEvolution.plt

t-prof.dat: MakeTimeseq.sh ${TOT}
	./MakeTimeseq.sh
#################
# ffmpeg
#################

movie: ${moviefile}

${vormovie}: MakeMovie.sh ${VORPNG}
	./MakeMovie.sh vor

${kinmovie}: MakeMovie.sh ${KINPNG}
	./MakeMovie.sh kin

${jcdmovie}: MakeMovie.sh ${JCDPNG}
	./MakeMovie.sh jcd

${magmovie}: MakeMovie.sh ${MAGPNG}
	./MakeMovie.sh mag

${kspmovie}: MakeMovie.sh ${KSPPNG}
	./MakeMovie.sh ksp

${hksmovie}: MakeMovie.sh ${HKSPNG}

#################
# gnuplot
#################

########################
# spectrum
########################

2Dsnaps: ${VORPNG}

${VORPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor vor.plt

${JCDPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor vor.plt

${KINPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor vor.plt

${MAGPNG}: vor.plt MakeSnap.sh ${VOR}
	./MakeSnap.sh vor vor.plt

spectrum: ${specfile}

${specfile}: spectrum-fin.plt
	gnuplot spectrum-fin.plt

specs: ${KSPPNG} ${HSKPNG} 

${KSPPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spc spectrum.plt

${HSKPNG}: spectrum.plt MakeSnap.sh ${SPC}
	./MakeSnap.sh spc spectrum.plt

#################
# visualization
#################

${vis}: Visualization.f90
	${fc} ${fopt} ${FLAGHDF} $< -o ${vis}

${VOR}: ${vis} ${BIN} ${countfile}
	./${vis}

#################
# analysis
#################

count-number: ${countfile}

${countfile}: CountBindata.sh
	./CountBindata.sh

#################
# clean up
#################
allclean:
	rm -fr ${dirm} ${dirf} ${dira}
