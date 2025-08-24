dird := snapshots

exe  := kh.x

all: ${exe}

.PHONY: all clean allclean

fc= nvfortran
foptopenacc = -Minfo=accel -acc
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exe}: main.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
