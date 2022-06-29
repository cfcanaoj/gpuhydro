dird := snapshots

exe=kh.x

all: ${exe}

.PHONY: all clean allclean

fc= nvfortran
fopt=-g -traceback -O2
#fopt=-g -traceback -check all -fpe0

${exe}: main.f90
	${fc} ${fopt} $< -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

allclean:
	rm -f ${dirm}/*.wmv ${dirf}/*.png ${dird}/*.xss
