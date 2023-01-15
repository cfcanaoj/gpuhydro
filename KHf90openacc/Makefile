dird := snapshots

exeori=kh_ori.x
exepra=kh_pra.x

all: ${exeori} ${exepra}

.PHONY: all clean allclean

fc= nvfortran
foptopenacc = -Minfo=accel -acc
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exeori}: main_ori.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${exeori}

${exepra}: main_pra.f90
	${fc} ${fopt} ${foptopenacc} $< -o ${exepra}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
