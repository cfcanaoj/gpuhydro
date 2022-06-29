dird := meddata
dirf := figures
dirm := movies

MOVIES := ${dirm}/anidnx.mp4

XSS= ${dird}/Sc00001.xss $(shell ls ${dird}/Sc*.xss)
DNPNG=$(patsubst ${dird}/Sc%.xss,${dirf}/dnx%.png,$(XSS))

exe=simplehyd.x

all: ${exe} ${MOVIES}

.PHONY: all clean allclean

#################
# ffmpeg
#################
${dirm}/anidnx.mp4: MakeMovie.sh ${DNPNG} Init.sh
	./Init.sh
	./MakeMovie.sh dnx

#################
# gnuplot
#################

${DNPNG}: dn2dx.plt MakeSnap.sh ${XSS} Init.sh
	./Init.sh
	rm -fr ${dirf}/dnx*.png
	./MakeSnap.sh

${XSS}: ${exe} Init.sh
	./Init.sh
	rm -fr ${dird}/Sc*.xss
	./${exe}

fc=ifort -extend-source
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exe}: main.F
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
