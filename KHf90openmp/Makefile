
exe=kh.x

all: ${exe}

.PHONY: all clean allclean

fc= ftn
foptopenmp = -qopenmp
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exe}: main.f90
	${fc} ${fopt} ${foptopenmp} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

