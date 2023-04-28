
exe=kh.x

all: ${exe}

.PHONY: all clean allclean

##################################
# setup for CfCA XC50 begin
##################################
fc= ifort -mcmodel=medium
foptopenmp = -qopenmp -qopt-report=5
#foptopenmp = -parallel
fopt = -g -traceback -O2
#fopt=-g -traceback -check all -fpe0

# setup for CfCA XC50 end
##################################

##################################
# setup for Fugaku
##################################
#fc= mpifrtpx -Free
#foptopenmp = -Kopenmp
#fopt = -Kfast
# setup for Fugaku end
##################################

mpiroutines.o: mpiroutines.f90 
	${fc} ${fopt} ${foptopenmp} -c $< -o mpirompiroutines.o

main.o: main.f90 
	${fc} ${fopt} ${foptopenmp} -c $< -o main.o

${exe}: main.o mpiroutines.o
	${fc} ${fopt} ${foptopenmp} main.o mpirompiroutines.o -o ${exe}

#################
# clean up
#################
clean:
	rm -f ${exe} *.o *.mod *~

