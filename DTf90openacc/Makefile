dird := bindata

exe=Simulation.x
ana=Analysis.x

all: ${exe}

.PHONY: all clean allclean

#################
# Fortran
#################

fc= nvfortran
foptopenacc = -Minfo=accel -acc -mcmodel=medium
fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

${exe}: Simulation.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${exe}

${ana}: Analysis.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${ana}

#################
# clean up
#################
clean:
	rm -f ${exe} ${ana} *.o *.mod *~

allclean:
	rm -fr ${dirm} ${dirf} ${dird} ${dirf}
