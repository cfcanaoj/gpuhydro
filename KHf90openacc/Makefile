dird := snapshots

exe := kh.x

all: ${exe}

.PHONY: all clean allclean

##################################
# setup for CfCA g00
##################################
fc= mpif90 -mcmodel=medium
foptopenacc = -Minfo=accel -acc
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0
##################################

.SUFFIXES: .f90 .o

.f90.o: 
	${fc} ${fopt} ${foptopenacc} -c $< 

OBJS=main.o mpiroutines.o

mpiroutines.o:
main.o: main.f90 mpiroutines.o 

${exe}: ${OBJS}
	${fc} ${fopt} ${foptopenacc} ${OBJS} -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dird}/*.xss
