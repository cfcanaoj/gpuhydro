dirbin := bindata

exe=Simulation.x
ana=Analysis.x

all: ${exe}

.PHONY: all clean allclean

#################
# Fortran
#################

fc= mpif90
foptopenacc = -Minfo=accel -acc -mcmodel=medium
fopt = -g -traceback -O2
#fopt = -g -traceback -check all -fpe0

.SUFFIXES: .f90 .o

.f90.o: 
	${fc} ${fopt} ${foptopenacc} -c $< 

OBJS=main.o mpiroutines.o

mpiroutines.o:
main.o: main.f90  mpiroutines.o 

${exe}: ${OBJS}
	${fc} ${fopt} ${foptopenacc} ${OBJS} -o ${exe}

${ana}: Analysis.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${ana}

#################
# clean up
#################
clean:
	rm -f ${exe} ${ana} *.o *.mod *~

allclean:
	rm -fr ${dirbin}
